// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "laplace.h"

#include <strings.h>

#include "gamma.h"
#include "beta.h"
#include "distmodelfactory.h"
#include "matrix.h"
#include "utils.h"
#include "compoundmodel.h"
#include "distgamma.h"
#include "distlognormal.h"
#include "distbeta.h"
#include "distbetaprime.h"
#include "distmultinormal.h"


#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

// Number of points of the log-spaced profile used by the boundary fits.
#define LAPLACE_GRID 10
// Smallest point of that profile, so that log(x) stays finite for a MAP at 0.
#define LAPLACE_GRID_FLOOR 1.0e-8

// Map a flattened element index to its owning Parameter and local index, so a
// single vector Parameter (e.g. all branch lengths / all ratios) is treated as
// one entry per scalar element rather than collapsing onto element 0.
static Parameter* laplace_element(Parameters* ps, size_t e, size_t* local){
	for(size_t i = 0; i < Parameters_count(ps); i++){
		Parameter* p = Parameters_at(ps, i);
		size_t sz = Parameter_size(p);
		if(e < sz){ *local = e; return p; }
		e -= sz;
	}
	*local = 0;
	return NULL;
}

// First derivative of model->logP wrt one element. The step has an absolute floor
// so that it does not vanish at 0, and a forward difference is used when a central
// one would step below the parameter's lower bound (e.g. a branch length at 0).
static double laplace_first_derivative(Model* model, Parameter* p, size_t local){
	double v = Parameter_value_at(p, local);
	double h = 1.0e-5 * fmax(fabs(v), 1.0e-3);
	double lower = p->cnstr != NULL ? Parameter_lower(p) : -INFINITY;

	Parameter_set_value_at(p, v + h, local);
	double fp = model->logP(model);
	double d;
	if(v - h > lower){
		Parameter_set_value_at(p, v - h, local);
		double fm = model->logP(model);
		d = (fp - fm)/(2.0*h);
	}
	else{
		Parameter_set_value_at(p, v, local);
		d = (fp - model->logP(model))/h;
	}
	Parameter_set_value_at(p, v, local);
	model->logP(model);
	return d;
}

// Flattened first derivatives of `model` over every scalar element of `ps`.
static void laplace_first_derivatives(Model* model, Parameters* ps, double* grad){
	size_t dim = Parameters_size(ps);
	for(size_t k = 0; k < dim; k++){
		size_t local;
		Parameter* p = laplace_element(ps, k, &local);
		grad[k] = laplace_first_derivative(model, p, local);
	}
}

// Gamma log density in the shape/rate parameterization. Unlike log(gsl pdf) it
// does not underflow, and it is finite at x = 0 for an exponential (shape 1).
static double laplace_log_gamma(double x, double shape, double rate){
	double lp = shape*log(rate) - gsl_sf_lngamma(shape) - rate*x;
	if(shape != 1.0) lp += (shape - 1.0)*log(x);
	return lp;
}

// Profile model->logP along a log-spaced grid x[0..N) from max(map, floor) to hi,
// for the element `local` of `p`, whose value is x + shift. On return y[j] holds
// logP(x[j]) - logP(x[0]) and the parameter is restored to map + shift.
static void laplace_profile(Model* model, Parameter* p, size_t local, double map,
                            double shift, double hi, double* x, double* y){
	double lo = fmax(map, LAPLACE_GRID_FLOOR);
	hi = fmax(hi, 2.0*lo);
	log_spaced_spaced_vector2(x, lo, hi, LAPLACE_GRID);
	for(size_t j = 0; j < LAPLACE_GRID; j++){
		Parameter_set_value_at(p, x[j] + shift, local);
		y[j] = model->logP(model);
	}
	Parameter_set_value_at(p, map + shift, local);
	model->logP(model);
	for(size_t j = LAPLACE_GRID; j-- > 0;){
		y[j] -= y[0];
	}
}

// Gamma fit for an element whose Laplace (mode and curvature) fit is unusable.
// Two cases are handled; anything else returns false and leaves shape/rate alone.
//
//  - The MAP sits on or near the zero boundary, or logP is not concave there: an
//    exponential matching the slope, shape 1 and rate |f'(m)|. For a likelihood
//    that is smooth at 0 this is accurate; fitting a shape < 1 instead is badly
//    conditioned, since log q(m) carries (shape - 1) log(m) with m close to 0.
//  - A small MAP with a spuriously large variance: the mode is fixed at m (so
//    rate = (shape - 1)/m) and the shape is fitted in (1, 100] by least squares
//    of the gamma log density on the log posterior profile, both relative to the
//    first grid point. With u = log x_j - log x_0 and v = x_j - x_0 the gamma log
//    density difference is (shape - 1)(u - v/m), linear in shape - 1, so the fit
//    has a closed form. If it fails the exponential above is used.
static bool laplace_boundary_gamma(Model* model, Parameter* p, size_t local,
                                   double map, double shift, double dlogP,
                                   double d2logP, double* shape, double* rate){
	double rate0 = map * -d2logP;
	double shape0 = rate0*map + 1.0;
	bool boundary = map < 1.e-6 || d2logP >= 0;
	bool flat = !boundary && shape0/(rate0*rate0) > 0.1 && map < 0.0001;
	if(!boundary && !flat) return false;

	*shape = 1.0;
	*rate = fmax(fabs(dlogP), 1.0e-8);
	if(boundary) return true;

	double x[LAPLACE_GRID];
	double y[LAPLACE_GRID];
	laplace_profile(model, p, local, map, shift, 0.5, x, y);
	double num = 0;
	double den = 0;
	for(size_t j = 1; j < LAPLACE_GRID; j++){
		double w = log(x[j]) - log(x[0]) - (x[j] - x[0])/map;
		num += y[j]*w;
		den += w*w;
	}
	double c = num/den;
	if(c > 0){
		*shape = 1.0 + (c > 99.0 ? 99.0 : c);
		*rate = (*shape - 1.0)/map;
	}
	return true;
}

// Gamma Laplace approximation of each element of `ps` (element value = x + shift).
// Returns -sum log q(map) and writes the fitted parameters into `dm` when non-NULL,
// following its shape/rate or shape/scale parameterization.
static double laplace_gamma(Laplace* laplace, Parameters* ps, double shift,
                            DistributionModel* dm){
	Model* posterior = laplace->model;
	size_t dim = Parameters_size(ps);
	double* alphas = dvector(dim);
	double* betas = dvector(dim);
	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);

	posterior->logP(posterior);
	posterior->hessian(posterior, ps, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, ps, dlogP_all);

	double logQ = 0;
	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(ps, i, &local);
		double map = Parameter_value_at(param, local) - shift;
		double d2logP = d2logP_all[i];
		// rate = -f''(m) * m, shape = rate * m + 1
		double rate = map * -d2logP;
		double shape = rate*map + 1;
		laplace_boundary_gamma(posterior, param, local, map, shift, dlogP_all[i],
		                       d2logP, &shape, &rate);

		logQ += laplace_log_gamma(map, shape, rate);
		alphas[i] = shape;
		betas[i] = rate;
		if (dm != NULL && dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
			betas[i] = 1.0/rate;
		}
	}

	if(dm != NULL){
		Parameter_set_values(Parameters_at(dm->parameters, 0), alphas);
		Parameter_set_values(Parameters_at(dm->parameters, 1), betas);
	}

	free(alphas);
	free(betas);
	free(d2logP_all);
	free(dlogP_all);
	return -logQ;
}

// Lognormal Laplace approximation of each element of `ps`:
//	sigma = sqrt(-1/(f''(m)*m^2)), mu = log(m) + sigma^2
// Elements for which it is unusable (MAP near 0, logP not concave, or mu > 5) fall
// back to a gamma. Since `dm` is a lognormal, those elements are stored as the
// lognormal with the same mean and variance as the fitted gamma.
static double laplace_lognormal(Laplace* laplace, Parameters* ps, DistributionModel* dm){
	Model* posterior = laplace->model;
	size_t dim = Parameters_size(ps);
	double* mus = dvector(dim);
	double* sigmas = dvector(dim);
	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);

	posterior->logP(posterior);
	posterior->hessian(posterior, ps, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, ps, dlogP_all);

	double logQ = 0;
	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(ps, i, &local);
		double map = Parameter_value_at(param, local);
		double d2logP = d2logP_all[i];

		double sigma = sqrt(-1.0/(d2logP*map*map));
		double mu = log(map) + sigma*sigma;
		if (map < 1.e-6 || d2logP >= 0 || mu > 5) {
			double rate = map * -d2logP;
			double shape = rate*map + 1;
			laplace_boundary_gamma(posterior, param, local, map, 0.0, dlogP_all[i],
			                       d2logP, &shape, &rate);
			logQ += laplace_log_gamma(map, shape, rate);

			double mean = shape/rate;
			double s2 = log1p(1.0/shape); // log(1 + var/mean^2)
			sigma = sqrt(s2);
			mu = log(mean) - s2/2.0;
		}
		else{
			logQ += log(gsl_ran_lognormal_pdf(map, mu, sigma));
		}
		mus[i] = mu;
		sigmas[i] = sigma;
	}

	if(dm != NULL){
		Parameter_set_values(Parameters_at(dm->parameters, 0), mus);
		Parameter_set_values(Parameters_at(dm->parameters, 1), sigmas);
	}
	free(mus);
	free(sigmas);
	free(d2logP_all);
	free(dlogP_all);
	return -logQ;
}

// Beta Laplace approximation of each element of `ps` (bounds (0,1)). Matching the
// mode m = (alpha-1)/(alpha+beta-2) and the curvature
//	f''(m) = -(alpha-1)/m^2 - (beta-1)/(1-m)^2 = -s/(m(1-m)),  s = alpha+beta-2
// gives s = -f''(m) m (1-m), alpha = 1 + m s, beta = 1 + (1-m) s.
static double laplace_beta(Laplace* laplace, Parameters* ps, DistributionModel* dm){
	Model* posterior = laplace->model;
	size_t dim = Parameters_size(ps);
	double* alphas = dvector(dim);
	double* betas = dvector(dim);
	double* d2logP_all = dvector(dim);

	posterior->logP(posterior);
	posterior->hessian(posterior, ps, HESSIAN_DIAGONAL, d2logP_all);

	double logQ = 0;
	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(ps, i, &local);
		double map = Parameter_value_at(param, local);
		double s = -d2logP_all[i]*map*(1.0 - map);
		if (s <= 0) {
			fprintf(stderr, "Beta Laplace: logP is not concave at %s[%zu] = %g\n",
			        Parameter_name(param), local, map);
		}
		double alpha = 1.0 + map*s;
		double beta = 1.0 + (1.0 - map)*s;

		logQ += log(gsl_ran_beta_pdf(map, alpha, beta));
		alphas[i] = alpha;
		betas[i] = beta;
	}

	if(dm != NULL){
		Parameter_set_values(Parameters_at(dm->parameters, 0), alphas);
		Parameter_set_values(Parameters_at(dm->parameters, 1), betas);
	}
	free(alphas);
	free(betas);
	free(d2logP_all);
	return -logQ;
}

double calculate_laplace_beta(Laplace* laplace){
	DistributionModel* dm = laplace->refdist != NULL ? laplace->refdist->obj : NULL;
	double logP = laplace->model->logP(laplace->model);
	logP += laplace_beta(laplace, laplace->parameters, dm);
	printf("Beta Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma(Laplace* laplace){
	DistributionModel* dm = laplace->refdist != NULL ? laplace->refdist->obj : NULL;
	double logP = laplace->model->logP(laplace->model);
	logP += laplace_gamma(laplace, laplace->parameters, 0.0, dm);
	printf("Gamma Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_lognormal(Laplace* laplace){
	DistributionModel* dm = laplace->refdist != NULL ? laplace->refdist->obj : NULL;
	double logP = laplace->model->logP(laplace->model);
	logP += laplace_lognormal(laplace, laplace->parameters, dm);
	printf("Lognormal Laplace: %f\n", logP);
	return logP;
}

// Laplace approximation in log space with a full covariance: y = log(x) is
// approximated by a multivariate normal. The second-order expansion of
//	f(y) = logP(exp(y)) + sum y_i          (sum y_i is the log-Jacobian of x = e^y)
// around y* = log(map) is integrated exactly:
//	log Z = f(y*) + d/2 log(2 pi) - 1/2 log|A| + 1/2 g' A^-1 g,  A = -H
// with g and H the gradient and Hessian of f at y*. The gradient term matters
// because the posterior mode in x is not the mode in y; with it the estimate is
// exact whenever f is quadratic in y (a lognormal posterior). The matching normal
// has mean y* + A^-1 g and covariance A^-1. By the chain rule
//	g_i  = x_i f'_i + 1
//	H_ii = x_i f'_i + x_i^2 f''_ii,   H_ij = x_i x_j f''_ij  (i != j)
// Elements below 1e-6 are moved to 1e-4 first, since log(0) is undefined.
double calculate_laplace_multivariate_normal(Laplace* laplace){
	Model* model = laplace->model;
	Parameters* ps = laplace->parameters;
	DistributionModel* dm = laplace->refdist != NULL ? laplace->refdist->obj : NULL;
	size_t n = Parameters_size(ps);

	double* orig = dvector(n);
	double* map = dvector(n);
	for (size_t k = 0; k < n; k++) {
		size_t local;
		Parameter* p = laplace_element(ps, k, &local);
		orig[k] = map[k] = Parameter_value_at(p, local);
		if (map[k] < 1.0e-6) {
			map[k] += 1.0e-4;
			Parameter_set_value_at(p, map[k], local);
		}
	}

	double logP = model->logP(model);
	double* grad = dvector(n);
	double* Hnat = dvector(n*n);
	laplace_first_derivatives(model, ps, grad);
	model->hessian(model, ps, HESSIAN_FULL, Hnat);

	gsl_matrix* A = gsl_matrix_alloc(n, n);
	gsl_vector* g = gsl_vector_alloc(n);
	gsl_vector* delta = gsl_vector_alloc(n);
	double fy = logP;
	for (size_t i = 0; i < n; i++) {
		double xi = map[i];
		fy += log(xi);
		gsl_vector_set(g, i, xi*grad[i] + 1.0);
		gsl_matrix_set(A, i, i, -(xi*grad[i] + Hnat[i*n + i]*xi*xi));
		for (size_t j = i+1; j < n; j++) {
			double Aij = -Hnat[i*n + j]*xi*map[j];
			gsl_matrix_set(A, i, j, Aij);
			gsl_matrix_set(A, j, i, Aij);
		}
	}
	free(grad);
	free(Hnat);

	double logZ = NAN;
	gsl_error_handler_t* handler = gsl_set_error_handler_off();
	int status = gsl_linalg_cholesky_decomp1(A);
	gsl_set_error_handler(handler);

	if (status != GSL_SUCCESS) {
		fprintf(stderr, "Multivariatenormal Laplace: the negative Hessian in log space"
		        " is not positive definite\n");
	}
	else {
		double logdetA = 0;
		for (size_t i = 0; i < n; i++) {
			logdetA += 2.0*log(gsl_matrix_get(A, i, i));
		}
		gsl_linalg_cholesky_solve(A, g, delta);
		double gAg = 0;
		gsl_blas_ddot(g, delta, &gAg);
		logZ = fy + 0.5*n*log(2.0*M_PI) - 0.5*logdetA + 0.5*gAg;

		if (dm != NULL) {
			Parameter* mu = Parameters_at(dm->parameters, 0);
			Parameter* cov = Parameters_at(dm->parameters, 1);
			double* values = dvector(n*n);
			for (size_t i = 0; i < n; i++) {
				values[i] = log(map[i]) + gsl_vector_get(delta, i);
			}
			Parameter_set_values(mu, values);

			// covariance = A^-1, then either stored whole or as its packed Cholesky factor
			gsl_linalg_cholesky_invert(A);
			if (Parameter_size(cov) == n*n) {
				for (size_t i = 0; i < n; i++) {
					for (size_t j = 0; j < n; j++) {
						values[i*n + j] = gsl_matrix_get(A, i, j);
					}
				}
			}
			else {
				gsl_linalg_cholesky_decomp1(A);
				size_t row = 0;
				for (size_t i = 0; i < n; i++) {
					for (size_t j = 0; j <= i; j++) {
						values[row++] = gsl_matrix_get(A, i, j);
					}
				}
			}
			Parameter_set_values(cov, values);
			free(values);
		}
	}

	for (size_t k = 0; k < n; k++) {
		size_t local;
		Parameter* p = laplace_element(ps, k, &local);
		Parameter_set_value_at(p, orig[k], local);
	}
	model->logP(model);

	printf("Multivariatenormal Laplace: %f\n", logZ);

	free(orig);
	free(map);
	gsl_vector_free(g);
	gsl_vector_free(delta);
	gsl_matrix_free(A);
	return logZ;
}

// Beta prime Laplace approximation:
//	alpha = 1 - f''(m) * (m^2) * (m + 1)
//	beta = -f''(m) * m * (m + 1) - 1
// When beta < 0 an exponential-like shape is used (alpha = 1, beta = |f'(m)| - 1).
// If that beta is below 2, the mode is fixed at m (alpha = m (beta + 1) + 1) and
// beta is fitted to the log posterior profile. With u = log x_j - log x_0 and
// w = log(1 + x_j) - log(1 + x_0) the log density difference is
// (beta + 1)(m u - (m + 1) w), linear in beta + 1, so the fit has a closed form.
double calculate_laplace_betaprime(Laplace* laplace){
	Model* posterior = laplace->model;
	DistributionModel* dm = laplace->refdist != NULL ? laplace->refdist->obj : NULL;
	Parameters* ps = laplace->parameters;
	size_t dim = Parameters_size(ps);
	double* alphas = dvector(dim);
	double* betas = dvector(dim);
	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);
	double logP = posterior->logP(posterior);
	posterior->hessian(posterior, ps, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, ps, dlogP_all);

	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(ps, i, &local);
		double map = Parameter_value_at(param, local);
		double d2logP = d2logP_all[i];
		double alpha = 1.0 - d2logP*(map*map)*(map + 1.0);
		double beta = -d2logP*map*(map + 1.0) - 1.0;

		if (beta < 0) {
			beta = fabs(dlogP_all[i]) - 1;
			alpha = 1;

			if (beta < 2) {
				double x[LAPLACE_GRID];
				double y[LAPLACE_GRID];
				laplace_profile(posterior, param, local, map, 0.0, 0.5, x, y);
				double num = 0;
				double den = 0;
				for (size_t j = 1; j < LAPLACE_GRID; j++) {
					double z = map*(log(x[j]) - log(x[0]))
					           - (map + 1.0)*(log1p(x[j]) - log1p(x[0]));
					num += y[j]*z;
					den += z*z;
				}
				beta = num/den - 1.0;
				beta = beta < 2 ? 2 : (beta > 100 ? 100 : beta);
				alpha = map*(beta + 1) + 1;
			}
		}
		logP -= dlogbetaprime(map, alpha, beta);
		alphas[i] = alpha;
		betas[i] = beta;
	}

	if(dm != NULL){
		Parameter_set_values(Parameters_at(dm->parameters, 0), alphas);
		Parameter_set_values(Parameters_at(dm->parameters, 1), betas);
	}
	free(alphas);
	free(betas);
	free(d2logP_all);
	free(dlogP_all);

	printf("Beta' Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma_from_mcmc(Laplace* laplace){
	Model* posterior = laplace->model;
	Model* empirical = laplace->empirical;
	DistributionModel* dm = empirical->obj;
	const double *alphas = Parameter_values(Parameters_at(dm->parameters, 0));
	const double *betas = Parameter_values(Parameters_at(dm->parameters, 1));
	bool scale = dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE;

	for (size_t i = 0; i < Parameters_size(laplace->parameters); i++) {
		size_t local;
		Parameter* p = laplace_element(laplace->parameters, i, &local);
		double rate = scale ? 1.0/betas[i] : betas[i];
		double mode = (alphas[i]-1.0)/rate;
		if (mode < 0) {
			mode = 1.e-6;
		}
		Parameter_set_value_at(p, mode, local);
	}
	double logP = posterior->logP(posterior);
	double logQ = empirical->logP(empirical);
	double logLaplace = logP - logQ;

	printf("Gamma MCMC Laplace: %f\n", logLaplace);
	return logLaplace;
}

double calculate_laplace_lognormal_from_mcmc(Laplace* laplace){
	Model* posterior = laplace->model;
	Model* empirical = laplace->empirical;
	DistributionModel* dm = empirical->obj;
	const double *mus = Parameter_values(Parameters_at(dm->parameters, 0));
	const double *sigmas = Parameter_values(Parameters_at(dm->parameters, 1));

	for (size_t i = 0; i < Parameters_size(laplace->parameters); i++) {
		size_t local;
		Parameter* p = laplace_element(laplace->parameters, i, &local);
		double mode = exp(mus[i]-sigmas[i]*sigmas[i]);
		Parameter_set_value_at(p, mode, local);
	}
	double logP = posterior->logP(posterior);
	double logQ = empirical->logP(empirical);
	double logLaplace = logP - logQ;

	printf("Lognormal MCMC Laplace: %f\n", logLaplace);
	return logLaplace;
}

// The distribution of each block of parameters is given by a reference compound
// model (possibly a mixture of families); its distributions receive the fits.
double calculate_laplace(Laplace* laplace){
	CompoundModel* cm = laplace->refdist->obj;
	double logP = laplace->model->logP(laplace->model);
	for (int i = 0; i < cm->count; i++) {
		DistributionModel* dm = cm->models[i]->obj;
		double logQ;
		if (dm->type == DISTRIBUTION_GAMMA) {
			logQ = laplace_gamma(laplace, dm->x, dm->shift, dm);
			printf("Gamma Laplace: %f\n", logQ);
		}
		else if (dm->type == DISTRIBUTION_BETA) {
			logQ = laplace_beta(laplace, dm->x, dm);
			printf("Beta Laplace: %f\n", logQ);
		}
		else if (dm->type == DISTRIBUTION_LOGNORMAL) {
			logQ = laplace_lognormal(laplace, dm->x, dm);
			printf("Lognormal Laplace: %f\n", logQ);
		}
		else{
			fprintf(stderr, "Distribution not supported by laplace\n");
			exit(2);
		}
		logP += logQ;
	}
	printf("Laplace: %f\n", logP);
	return logP;
}

void _free_Laplace(Laplace* laplace){
	free_Parameters(laplace->parameters);
	laplace->model->free(laplace->model);
	if(laplace->refdist != NULL)laplace->refdist->free(laplace->refdist);
    if(laplace->empirical != NULL) laplace->empirical->free(laplace->empirical);
	free(laplace);
}

// The model whose log density is approximated, given by "target".
static Model* laplace_get_model(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	char* ref = get_json_node_value_string(node, "target");
	if (ref == NULL) {
		fprintf(stderr, "Laplace %s requires a \"target\"\n", id);
		exit(1);
	}
	Model* model = Hashtable_get(hash, ref+1);
	if (model == NULL) {
		fprintf(stderr, "Laplace %s: model %s not found\n", id, ref);
		exit(1);
	}
	model->ref_count++;
	return model;
}

Laplace* new_Laplace_from_json2(json_node* node, Hashtable* hash){
    static const json_field schema[] = {
        {"distribution", JSON_REQUIRED, JSON_OBJECT},
        {"target", JSON_REQUIRED, JSON_STRING},
        {"x", JSON_REQUIRED, JSON_ANY},
    };
    json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

    char* id = get_json_node_value_string(node, "id");
    json_node* dist_node = get_json_node(node, "distribution");
    json_node* x_node = get_json_node(node, "x");
    Laplace* laplace = malloc(sizeof(Laplace));
    laplace->model = laplace_get_model(node, hash);
    laplace->parameters = NULL;
    laplace->refdist = NULL;
    laplace->empirical = NULL;
    laplace->free = _free_Laplace;

    char* d_string = get_json_node_value_string(dist_node, "distribution");
    if(strcasecmp(d_string, "gamma") == 0){
        laplace->refdist = new_GammaDistributionModel_from_json(dist_node, hash);
    }
    else if(strcasecmp(d_string, "lognormal") == 0){
        laplace->refdist = new_LogNormalDistributionModel_from_json(dist_node, hash);
    }
    else if(strcasecmp(d_string, "beta") == 0){
        laplace->refdist = new_BetaDistributionModel_from_json(dist_node, hash);
    }
    else if(strcasecmp(d_string, "betaprime") == 0){
        laplace->refdist = new_BetaPrimeDistributionModel_from_json(dist_node, hash);
    }
    else if(strcasecmp(d_string, "multivariatenormal") == 0){
        laplace->refdist = new_MultivariateNormalDistributionModel_from_json(dist_node, hash);
    }
    else{
        fprintf(stderr, "Distribution unknown: %s in Laplace with ID: %s\n", d_string, id);
        exit(1);
    }

    DistributionModel* dm = laplace->refdist->obj;
    if (dm->type == DISTRIBUTION_GAMMA) {
        laplace->calculate = calculate_laplace_gamma;
    }
    else if (dm->type == DISTRIBUTION_LOGNORMAL) {
        laplace->calculate = calculate_laplace_lognormal;
    }
    else if (dm->type == DISTRIBUTION_BETA) {
        laplace->calculate = calculate_laplace_beta;
    }
    else if (dm->type == DISTRIBUTION_BETA_PRIME) {
        laplace->calculate = calculate_laplace_betaprime;
    }
    else if (dm->type == DISTRIBUTION_NORMAL_MULTIVARIATE) {
        laplace->calculate = calculate_laplace_multivariate_normal;
    }
    else{
        fprintf(stderr, "Distribution unknown in Laplace with ID: %s\n", id);
        exit(1);
    }

    laplace->parameters = distmodel_get_x(id, x_node, hash);

    return laplace;
}

static Model* laplace_get_empirical(json_node* empirical, Hashtable* hash){
	Model* empiricalDist = NULL;
	if (empirical->node_type == MJSON_OBJECT) {
		empiricalDist = new_DistributionModel_from_json(empirical, hash);
		char* id = get_json_node_value_string(empirical, "id");
		Hashtable_add(hash, id, empiricalDist);
	}
	else if(empirical->node_type == MJSON_STRING){
		char* ref = (char*)empirical->value;
		empiricalDist = Hashtable_get(hash, ref+1);
		empiricalDist->ref_count++;
	}
	else{
		fprintf(stderr, "Laplace: \"empirical\" should be a distribution or a reference\n");
		exit(10);
	}
	return empiricalDist;
}

Laplace* new_Laplace_from_json(json_node* node, Hashtable* hash){
	json_node* dist_node = get_json_node(node, "distribution");
	if(dist_node != NULL && dist_node->node_type != MJSON_STRING){
		return new_Laplace_from_json2(node, hash);
	}

	static const json_field schema[] = {
	    {"distribution", JSON_OPTIONAL, JSON_STRING},
	    {"empirical", JSON_OPTIONAL, JSON_ANY},
	    {"target", JSON_REQUIRED, JSON_STRING},
	    {"parameters", JSON_OPTIONAL, JSON_ANY},
	    {"ref", JSON_OPTIONAL, JSON_STRING},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

	char* id = get_json_node_value_string(node, "id");
	char* dist_string = get_json_node_value_string(node, "distribution");
	char* ref = get_json_node_value_string(node, "ref");
	Laplace* laplace = malloc(sizeof(Laplace));
	laplace->model = laplace_get_model(node, hash);
	laplace->parameters = new_Parameters(1);
	laplace->refdist = NULL;
	laplace->empirical = NULL;
	laplace->free = _free_Laplace;
	json_node* empirical = get_json_node(node, "empirical");
	if (ref != NULL) {
		laplace->refdist = Hashtable_get(hash, ref+1);
		if (laplace->refdist == NULL) {
			fprintf(stderr, "Laplace %s: reference %s not found\n", id, ref);
			exit(1);
		}
		laplace->refdist->ref_count++;
	}
	// distribution is inferred from the reference distribution
	// it can be a mixture of distributions
	if (dist_string == NULL) {
		if (laplace->refdist == NULL || laplace->refdist->type != MODEL_COMPOUND) {
			fprintf(stderr, "Laplace %s requires a \"distribution\" or a compound \"ref\"\n", id);
			exit(1);
		}
		laplace->calculate = calculate_laplace;
		return laplace;
	}

	get_parameters_references(node, hash, laplace->parameters);
	if (Parameters_count(laplace->parameters) == 0) {
		fprintf(stderr, "Laplace %s requires \"parameters\"\n", id);
		exit(1);
	}

	if(strcasecmp(dist_string, "gamma") == 0){
		if(empirical != NULL){
			laplace->empirical = laplace_get_empirical(empirical, hash);
			laplace->calculate = calculate_laplace_gamma_from_mcmc;
		}
		else{
			laplace->calculate = calculate_laplace_gamma;
		}
	}
	else if(strcasecmp(dist_string, "lognormal") == 0){
		if(empirical != NULL){
			laplace->empirical = laplace_get_empirical(empirical, hash);
			laplace->calculate = calculate_laplace_lognormal_from_mcmc;
		}
		else{
			laplace->calculate = calculate_laplace_lognormal;
		}
	}
	else if(strcasecmp(dist_string, "betaprime") == 0){
		laplace->calculate = calculate_laplace_betaprime;
	}
	else if(strcasecmp(dist_string, "multivariate") == 0){
		laplace->calculate = calculate_laplace_multivariate_normal;
	}
	else if(strcasecmp(dist_string, "beta") == 0){
		laplace->calculate = calculate_laplace_beta;
	}
	else{
		fprintf(stderr, "Laplace distribution not available %s\n", dist_string);
		exit(13);
	}

	return laplace;
}
