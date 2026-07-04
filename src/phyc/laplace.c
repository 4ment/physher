// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "laplace.h"

#include <strings.h>

#include "gamma.h"
#include "beta.h"
#include "distmodelfactory.h"
#include "matrix.h"
#include "utils.h"
#include "optimizer.h"
#include "compoundmodel.h"
#include "distgamma.h"
#include "distlognormal.h"
#include "distbeta.h"
#include "distbetaprime.h"
#include "distmultinormal.h"


#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

struct laplace_data_t{
	double map;
	double* x;
	double* y;
	double*yy;
	size_t N;
};

// shape is fixed and scale is optimized
static double _func_gamma_fixed_shape( Parameters *params, double *grad, void *data ){
	struct laplace_data_t* d = (struct laplace_data_t*)data;
	double alpha = Parameters_value(params, 0);
	double beta = d->map;
	double sumYY = 0;
	for (size_t i = 0; i < d->N; i++) {
		d->yy[i] = gsl_ran_gamma_pdf(d->x[i], alpha, 1.0/beta);
		sumYY += d->yy[i];
	}
	double maxY = dmax_vector(d->yy, d->N);
	double sum = 0;
	for (size_t i = 0; i < d->N; i++) {
		sum += pow((d->yy[i]-maxY) - d->y[i], 2);
	}
	return sum;
}

// scale is optimized and shape is constrained
static double _func_gamma_fixed_mode( Parameters *params, double *grad, void *data ){
	struct laplace_data_t* d = (struct laplace_data_t*)data;
	double alpha = Parameters_value(params, 0);
	// mode = (alpha-1)/beta for alpha >= 1
	double beta = (alpha - 1.0)/d->map;
	double sumYY = 0;
	for (size_t i = 0; i < d->N; i++) {
		d->yy[i] = gsl_ran_gamma_pdf(d->x[i], alpha, 1.0/beta);
		sumYY += d->yy[i];
	}
	double maxY = dmax_vector(d->yy, d->N);
	double sum = 0;
	for (size_t i = 0; i < d->N; i++) {
		//sum += pow(d->yy[i]/sumYY/d->N - d->y[i], 2);
		sum += pow((d->yy[i]-maxY) - d->y[i], 2);
	}
	return sum;
}

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

// Flattened first derivatives of `model` over every scalar element of `ps`.
static void laplace_first_derivatives(Model* model, Parameters* ps, double eps, double* grad){
	size_t k = 0;
	for(size_t i = 0; i < Parameters_count(ps); i++){
		Parameter* p = Parameters_at(ps, i);
		Model_first_derivatives(model, p, eps, grad + k);
		k += Parameter_size(p);
	}
}

// Parameters are modeled with a beta distribution B(alpha,beta)
// Bounds of the parameters are (0,1)
// Shoud check out the Kumaraswamy distribution
double calculate_laplace_beta(Laplace* laplace){
	/*
	a = m*m/(1-m)/(1-m)
	b = 1 - m*m*f''(m)
	c = 1 - m
	d = 2*m - 1
	x = (b*c + d)/m
	beta = (a*c - x*m)/(a*c - m)
	alpha = (beta - 1)*a + b
	*/
	
	Model* posterior = laplace->model;
	Model* refdist = laplace->refdist;
	DistributionModel* dm = NULL;
	double* alphas = NULL;
	double* betas = NULL;

	size_t dim = Parameters_size(laplace->parameters);
	if(refdist != NULL){
		dm = refdist->obj;
		alphas = malloc(sizeof(double)*dim);
		betas = malloc(sizeof(double)*dim);
	}
	double* d2logP_all = dvector(dim);
	double logP = posterior->logP(posterior);
	posterior->hessian(posterior, laplace->parameters, HESSIAN_DIAGONAL, d2logP_all);

	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(laplace->parameters, i, &local);
		double map = Parameter_value_at(param, local);
		double d2logP = d2logP_all[i];
		double a = map*map/(1.0 - map)/(1.0 - map);
		double b = 1.0 - map*map*d2logP;
		double c = 1.0 - map;
		double d = 2.0*map - 1.0;
		double x = (b*c + d)/map;
		double beta = (a*c - x*map)/(a*c - map);
		double alpha = (beta - 1.0)*a + b;

		// should handle small values B(alpha, 1) (exponential shape)
		// should handle large values B(1, Beta) (mirror-image exponential shape)
		// check that we don't get alpha <= 1 or beta <= 1

		if (map < 1.e-6) {

		}
		else if (1.0-map < 1.e-6) {

		}
		if (alpha <= 1 || beta <=1 ) {

		}

		logP -= log(gsl_ran_beta_pdf(map, alpha, beta));

		if(dm != NULL){
			alphas[i] = alpha;
			betas[i] = beta;
		}
	}

	if(dm != NULL){
		Parameter_set_values(Parameters_at(dm->parameters, 0), alphas);
		Parameter_set_values(Parameters_at(dm->parameters, 1), betas);
	}
	free(alphas);
	free(betas);
	free(d2logP_all);

	return logP;
}

double calculate_laplace_beta2(Laplace* laplace, DistributionModel* dm){
	/*
	 a = m*m/(1-m)/(1-m)
	 b = 1 - m*m*f''(m)
	 c = 1 - m
	 d = 2*m - 1
	 x = (b*c + d)/m
	 beta = (a*c - x*m)/(a*c - m)
	 alpha = (beta - 1)*a + b
	 */
	
	Model* posterior = laplace->model;
	double logP = 0;//posterior->logP(posterior);
	Parameters* parameters = dm->x;
	size_t dim = Parameters_size(parameters);
	double* alphas = malloc(sizeof(double)*dim);
	double* betas = malloc(sizeof(double)*dim);
	double* d2logP_all = dvector(dim);
	posterior->logP(posterior); // initialize
	posterior->hessian(posterior, parameters, HESSIAN_DIAGONAL, d2logP_all);

	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* parameter = laplace_element(parameters, i, &local);
		double map = Parameter_value_at(parameter, local);
		double d2logP = d2logP_all[i];
		double a = map*map/(1.0 - map)/(1.0 - map);
		double b = 1.0 - map*map*d2logP;
		double c = 1.0 - map;
		double d = 2.0*map - 1.0;
		double x = (b*c + d)/map;
		double beta = (a*c - x*map)/(a*c - map);
		double alpha = (beta - 1.0)*a + b;

		// should handle small values B(alpha, 1) (exponential shape)
		// should handle large values B(1, Beta) (mirror-image exponential shape)
		// check that we don't get alpha <= 1 or beta <= 1

		if (map < 1.e-6) {

		}
		else if (1.0-map < 1.e-6) {

		}
		if (alpha <= 1 || beta <=1 ) {

		}

		logP -= log(gsl_ran_beta_pdf(map, alpha, beta));

		alphas[i] = alpha;
		betas[i] = beta;
	}

	Parameter_set_values(Parameters_at(dm->parameters, 0), alphas);
	Parameter_set_values(Parameters_at(dm->parameters, 1), betas);

	free(d2logP_all);
	free(alphas);
	free(betas);

	printf("Beta Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma(Laplace* laplace){
	// beta = rate = -f''(m) * m
	// alpha = shape = rate * m + 1
	Model* posterior = laplace->model;
	Model* refdist = laplace->refdist;
	DistributionModel* dm = NULL;
	double* alphas = NULL;
	double* betas = NULL;

	size_t dim = Parameters_size(laplace->parameters);
	if(refdist != NULL){
		dm = refdist->obj;
		alphas = malloc(sizeof(double)*dim);
		betas = malloc(sizeof(double)*dim);
	}

	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);

	double logP = posterior->logP(posterior);
	posterior->hessian(posterior, laplace->parameters, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, laplace->parameters, 1e-5, dlogP_all);
	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(laplace->parameters, i, &local);
		double map = Parameter_value_at(param, local);
		double d2logP = d2logP_all[i];
		double rate = map * -d2logP;
		double shape = rate*map + 1;

		// Very small branch -> exponential shape
		if (map < 1.e-6 || d2logP >= 0) {
			double dlogP = dlogP_all[i];
			shape = 1;
			rate = fabs(dlogP);

			log_spaced_spaced_vector2(x, map, 0.5, N);

			for (size_t j = 1; j < N; j++) {
				Parameter_set_value_at(param, x[j], local);
				y[j] = laplace->model->logP(laplace->model);
			}
			Parameter_set_value_at(param, map, local);
			y[0] = laplace->model->logP(laplace->model);

			double maxY = y[0];
			for (int j = 0; j < N; j++) {
				y[j] -= maxY;
			}
            double lower = 0.001;
            double upper = 1;
			double guess = 1.0 - 0.001;
			Parameters* ps = new_Parameters(1);
			Parameters_move(ps, new_Parameter("", guess, new_Constraint(lower, upper)));

			struct laplace_data_t data = {rate, x, y, yy, N};
			double fx = _func_gamma_fixed_shape(ps, NULL, &data);

			Parameters_set_value(ps, 0, lower);
			double fa = _func_gamma_fixed_shape(ps, NULL, &data);
			Parameters_set_value(ps, 0, upper);
			double fb = _func_gamma_fixed_shape(ps, NULL, &data);
			Parameters_set_value(ps, 0, guess);

			if(fa > fx && fx < fb){
				Optimizer* opt = new_Optimizer(OPT_BRENT);
				opt_set_data(opt, &data);
				opt_set_objective_function(opt, _func_gamma_fixed_shape);
				opt_set_parameters(opt, ps);
				double min;
				opt_optimize(opt, ps, &min);

				shape = Parameters_value(ps, 0);

				free_Optimizer(opt);
			}
			free_Parameters(ps);
		}
		// Small branch with a maximum and spurious large variance
		else if(shape/(rate*rate) > 0.1 && map < 0.0001){
			double dlogP = dlogP_all[i];
			shape = 1;
			rate = fabs(dlogP);

			log_spaced_spaced_vector2(x, map, 0.5, N);

			for (size_t j = 1; j < N; j++) {
				Parameter_set_value_at(param, x[j], local);
				y[j] = laplace->model->logP(laplace->model);
			}
			Parameter_set_value_at(param, map, local);
			y[0] = laplace->model->logP(laplace->model);

			double maxY = y[0];
			for (int j = 0; j < N; j++) {
				y[j] -= maxY;
			}

			double guess = 1.0 + 0.001;
			Parameters* ps = new_Parameters(1);
			Parameters_move(ps, new_Parameter("", guess, new_Constraint(1, 100)));

			struct laplace_data_t data = {map, x, y, yy, N};
			double fx = _func_gamma_fixed_mode(ps, NULL, &data);

			Parameters_set_value(ps, 0, 1);
			double fa = _func_gamma_fixed_mode(ps, NULL, &data);
			Parameters_set_value(ps, 0, 100);
			double fb = _func_gamma_fixed_mode(ps, NULL, &data);
			Parameters_set_value(ps, 0, guess);

			if(fa > fx && fx < fb){
				Optimizer* opt = new_Optimizer(OPT_BRENT);
				opt_set_data(opt, &data);
				opt_set_objective_function(opt, _func_gamma_fixed_mode);
				opt_set_parameters(opt, ps);
				double min;
				opt_optimize(opt, ps, &min);

				shape = Parameters_value(ps, 0);
				rate = (shape - 1)/map;

				free_Optimizer(opt);
			}
			free_Parameters(ps);
		}

		logP -= log(gsl_ran_gamma_pdf(map, shape, 1.0/rate));

		if(dm != NULL){
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
				rate = 1.0/rate;
			}
			alphas[i] = shape;
			betas[i] = rate;
		}
	}

	if(dm != NULL){
		Parameter_set_values(Parameters_at(dm->parameters, 0), alphas);
		Parameter_set_values(Parameters_at(dm->parameters, 1), betas);
	}

	free(x);
	free(y);
	free(yy);
	free(d2logP_all);
	free(dlogP_all);
	free(alphas);
	free(betas);

	printf("Gamma Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma2(Laplace* laplace, DistributionModel* dm){
	// beta = rate = -f''(m) * m
	// alpha = shape = rate * m + 1
	Model* posterior = laplace->model;
	Parameters* parameters = dm->x;
	size_t dim = Parameters_size(parameters);
	double* alphas = malloc(sizeof(double)*dim);
	double* betas = malloc(sizeof(double)*dim);

	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);

	double logP = 0;
	posterior->logP(posterior);// initialize
	posterior->hessian(posterior, parameters, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, parameters, 1e-5, dlogP_all);

	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* parameter = laplace_element(parameters, i, &local);
		double map_orig = Parameter_value_at(parameter, local);
		double map = map_orig - dm->shift;
		double dlogP = dlogP_all[i];
		double d2logP = d2logP_all[i];
		double rate = map * -d2logP;
		double shape = rate*map + 1;
		
		// Very small branch -> exponential shape
		if (map < 1.e-6 || d2logP >= 0) {
			dlogP = dlogP_all[i];
			shape = 1;
			rate = fabs(dlogP);

			log_spaced_spaced_vector2(x, map, 0.5, N);

			for (size_t j = 1; j < N; j++) {
				Parameter_set_value_at(parameter, x[j] + dm->shift, local);
				y[j] = posterior->logP(posterior);
			}
			Parameter_set_value_at(parameter, map_orig, local);
			y[0] = posterior->logP(posterior);
			
			double maxY = y[0];
			for (int j = 0; j < N; j++) {
				y[j] -= maxY;
			}
			double lower = 0.001;
			double upper = 1;
			double guess = 1.0 - 0.001;
			Parameters* ps = new_Parameters(1);
			Parameters_move(ps, new_Parameter("", guess, new_Constraint(lower, upper)));
			
			struct laplace_data_t data = {rate, x, y, yy, N};
			double fx = _func_gamma_fixed_shape(ps, NULL, &data);
			
			Parameters_set_value(ps, 0, lower);
			double fa = _func_gamma_fixed_shape(ps, NULL, &data);
			Parameters_set_value(ps, 0, upper);
			double fb = _func_gamma_fixed_shape(ps, NULL, &data);
			Parameters_set_value(ps, 0, guess);
			
			if(fa > fx && fx < fb){
				Optimizer* opt = new_Optimizer(OPT_BRENT);
				opt_set_data(opt, &data);
				opt_set_objective_function(opt, _func_gamma_fixed_shape);
				opt_set_parameters(opt, ps);
				double min;
				opt_optimize(opt, ps, &min);
				
				shape = Parameters_value(ps, 0);
				
				free_Optimizer(opt);
			}
			free_Parameters(ps);
		}
		// Small branch with a maximum and spurious large variance
		else if(shape/(rate*rate) > 0.1 && map < 0.0001){
			dlogP = dlogP_all[i];
			shape = 1;
			rate = fabs(dlogP);

			log_spaced_spaced_vector2(x, map, 0.5, N);

			for (size_t j = 1; j < N; j++) {
				Parameter_set_value_at(parameter, x[j] + dm->shift, local);
				y[j] = posterior->logP(posterior);
			}
			Parameter_set_value_at(parameter, map_orig, local);
			y[0] = posterior->logP(posterior);
			
			double maxY = y[0];
			for (int j = 0; j < N; j++) {
				y[j] -= maxY;
			}
			
			double guess = 1.0 + 0.001;
			Parameters* ps = new_Parameters(1);
			Parameters_move(ps, new_Parameter("", guess, new_Constraint(1, 100)));
			
			struct laplace_data_t data = {map, x, y, yy, N};
			double fx = _func_gamma_fixed_mode(ps, NULL, &data);
			
			Parameters_set_value(ps, 0, 1);
			double fa = _func_gamma_fixed_mode(ps, NULL, &data);
			Parameters_set_value(ps, 0, 100);
			double fb = _func_gamma_fixed_mode(ps, NULL, &data);
			Parameters_set_value(ps, 0, guess);
			
			if(fa > fx && fx < fb){
				Optimizer* opt = new_Optimizer(OPT_BRENT);
				opt_set_data(opt, &data);
				opt_set_objective_function(opt, _func_gamma_fixed_mode);
				opt_set_parameters(opt, ps);
				double min;
				opt_optimize(opt, ps, &min);
				
				shape = Parameters_value(ps, 0);
				rate = (shape - 1)/map;
				
				free_Optimizer(opt);
			}
			free_Parameters(ps);
		}
		printf("map: %f dlogP: %e d2logP: %f [%f,%f]\n",map, dlogP, d2logP, shape, rate);

		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
			rate = 1.0/rate;
		}
		alphas[i] = shape;
		betas[i] = rate;
		logP -= log(gsl_ran_gamma_pdf(map, shape, rate));
	}

	Parameter_set_values(Parameters_at(dm->parameters, 0), alphas);
	Parameter_set_values(Parameters_at(dm->parameters, 1), betas);

	free(alphas);
	free(betas);
	free(d2logP_all);
	free(dlogP_all);

	free(x);
	free(y);
	free(yy);

	printf("Gamma Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_multivariate_normal(Laplace* laplace){
    Model* refdist = laplace->refdist;
    DistributionModel* dm = NULL;
    if(refdist != NULL){
        dm = refdist->obj;
    }

	// dim is the number of flattened scalar elements: a single vector Parameter
	// (e.g. all branch lengths) contributes Parameter_size entries, not one.
	size_t paramCount = Parameters_size(laplace->parameters);
	gsl_matrix* H = gsl_matrix_alloc(paramCount, paramCount);
	gsl_vector* mu = gsl_vector_alloc(paramCount);
	gsl_matrix * L = gsl_matrix_alloc (paramCount, paramCount);
	gsl_permutation* perm = gsl_permutation_alloc (paramCount);
	gsl_vector* work = gsl_vector_alloc(paramCount);
	gsl_set_error_handler_off();
	double logP = laplace->model->logP(laplace->model);
	double epsilon = 0.0001;

	// Gather the flattened MAP values, first derivatives, and the natural-space
	// Hessian. The Laplace approximation is built in log space, so the chain rule
	// is applied below: with y = log(x),
	//   d2/dy_i^2     = x_i f'_i        + x_i^2 f''_ii
	//   d2/dy_i dy_j  =                   x_i x_j f''_ij   (i != j)
	double* map = dvector(paramCount);
	double* grad = dvector(paramCount);
	double* Hnat = dvector(paramCount * paramCount);

	size_t k = 0;
	for (size_t p = 0; p < Parameters_count(laplace->parameters); p++) {
		Parameter* param = Parameters_at(laplace->parameters, p);
		Model_first_derivatives(laplace->model, param, 1e-5, grad + k);
		for (size_t e = 0; e < Parameter_size(param); e++) {
			double v = Parameter_value_at(param, e);
			if (v < 1.0e-6) v += epsilon;
			map[k++] = v;
		}
	}

	laplace->model->hessian(laplace->model, laplace->parameters, HESSIAN_FULL, Hnat);

	for (size_t i = 0; i < paramCount; i++) {
		double mapi = map[i];
		double Hii = grad[i]*mapi + Hnat[i*paramCount + i]*mapi*mapi;
		gsl_matrix_set(H, i, i, Hii);
		gsl_vector_set(mu, i, log(mapi));

		for (size_t j = i+1; j < paramCount; j++) {
			double Hij = Hnat[i*paramCount + j] * mapi * map[j];
			gsl_matrix_set(H, i, j, Hij);
			gsl_matrix_set(H, j, i, Hij);
		}
	}
	free(grad);
	free(Hnat);
	int signum;
	gsl_linalg_LU_decomp (H, perm, &signum);

	gsl_linalg_LU_invert (H, perm, L);
	for (int i = 0; i < paramCount; i++) {
		for (int j = 0; j < paramCount; j++) {
			gsl_matrix_set(L, i, j, -gsl_matrix_get(L, i, j));
		}
	}
	gsl_linalg_cholesky_decomp1(L);
	
	double logQ = 0;
	gsl_ran_multivariate_gaussian_log_pdf(mu, mu, L, &logQ, work);
	
	for (size_t i = 0; i < paramCount; i++) {
		// log-density of the Jacobian of x = exp(y): -sum log(x_i)
		logQ -= log(map[i]);
	}

    if(refdist != NULL){
		Parameter* mu = Parameters_at(dm->parameters, 0);
		Parameter* sigma = Parameters_at(dm->parameters, 1);
		size_t row = 0;
		for (size_t i = 0; i < paramCount; i++) {
            dm->tempx[i] = map[i];
		}
        if(paramCount*paramCount == Parameter_size(sigma)){
            for (size_t i = 0; i < paramCount; i++) {
                for (size_t j = 0; j < paramCount; j++) {
                    dm->tempx[row++] = gsl_matrix_get(H, i, j);
                }
            }
        }
        else{
            for(size_t i = 0; i < paramCount; i++){
                for(size_t j = 0; j <= i; j++){
                    dm->tempx[row++] = gsl_matrix_get(L, i, j);
                }
            }
			Parameter_set_values(sigma, dm->tempx);
        }
    }
	printf("Multivariatenormal Laplace: %f logQ: %f %d\n", logP - logQ, logQ, gsl_ran_multivariate_gaussian_log_pdf(mu, mu, L, &logQ, work));

	free(map);
	gsl_vector_free(work);
	gsl_matrix_free(H);
	gsl_matrix_free(L);
	gsl_vector_free(mu);
	gsl_permutation_free(perm);
	return logP - logQ;
}

double calculate_laplace_lognormal(Laplace* laplace){
	//	sigma = sqrt(-1/(f''(m)*m^2))
	//	mu    = log(m)+sigma^2
	Model* posterior = laplace->model;
    Model* refdist = laplace->refdist;
    DistributionModel* dm = NULL;
	double* mus = NULL;
	double* sigmas = NULL;

	size_t dim = Parameters_size(laplace->parameters);
    if(refdist != NULL){
        dm = refdist->obj;
		mus = malloc(sizeof(double)*dim);
		sigmas = malloc(sizeof(double)*dim);
    }

	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);
	double logP = posterior->logP(posterior);
	posterior->hessian(posterior, laplace->parameters, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, laplace->parameters, 1e-5, dlogP_all);
	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));

	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(laplace->parameters, i, &local);
		double map = Parameter_value_at(param, local);
		double d2logP = d2logP_all[i];

		double sigma = sqrt(-1.0/(d2logP*map*map));
		double mu = log(map) + sigma*sigma;
		if (map < 1.e-6 || d2logP >= 0 || mu > 5) {
			double rate = map * -d2logP;
			double shape = rate*map + 1;
			// Very small branch -> exponential shape
			if (map < 1.e-6 || d2logP >= 0) {
				double dlogP = dlogP_all[i];
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameter_set_value_at(param, x[j], local);
					y[j] = laplace->model->logP(laplace->model);
				}
				Parameter_set_value_at(param, map, local);
				y[0] = laplace->model->logP(laplace->model);
				
				double maxY = y[0];
				for (int j = 0; j < N; j++) {
					y[j] -= maxY;
				}
				double lower = 0.001;
				double upper = 1;
				double guess = 1.0 - 0.001;
				Parameters* ps = new_Parameters(1);
				Parameters_move(ps, new_Parameter("", guess, new_Constraint(lower, upper)));
				
				struct laplace_data_t data = {rate, x, y, yy, N};
				double fx = _func_gamma_fixed_shape(ps, NULL, &data);
				
				Parameters_set_value(ps, 0, lower);
				double fa = _func_gamma_fixed_shape(ps, NULL, &data);
				Parameters_set_value(ps, 0, upper);
				double fb = _func_gamma_fixed_shape(ps, NULL, &data);
				Parameters_set_value(ps, 0, guess);
				
				if(fa > fx && fx < fb){
					Optimizer* opt = new_Optimizer(OPT_BRENT);
					opt_set_data(opt, &data);
					opt_set_objective_function(opt, _func_gamma_fixed_shape);
					opt_set_parameters(opt, ps);
					double min;
					opt_optimize(opt, ps, &min);
					
					shape = Parameters_value(ps, 0);
					
					free_Optimizer(opt);
				}
				free_Parameters(ps);
			}
			// Small branch with a maximum and spurious large variance
			else if(shape/(rate*rate) > 0.1 && map < 0.0001){
				double dlogP = dlogP_all[i];
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameter_set_value_at(param, x[j], local);
					y[j] = laplace->model->logP(laplace->model);
				}
				Parameter_set_value_at(param, map, local);
				y[0] = laplace->model->logP(laplace->model);
				
				double maxY = y[0];
				for (int j = 0; j < N; j++) {
					y[j] -= maxY;
				}
				
				double guess = 1.0 + 0.001;
				Parameters* ps = new_Parameters(1);
				Parameters_move(ps, new_Parameter("", guess, new_Constraint(1, 100)));
				
				struct laplace_data_t data = {map, x, y, yy, N};
				double fx = _func_gamma_fixed_mode(ps, NULL, &data);
				
				Parameters_set_value(ps, 0, 1);
				double fa = _func_gamma_fixed_mode(ps, NULL, &data);
				Parameters_set_value(ps, 0, 100);
				double fb = _func_gamma_fixed_mode(ps, NULL, &data);
				Parameters_set_value(ps, 0, guess);
				
				if(fa > fx && fx < fb){
					Optimizer* opt = new_Optimizer(OPT_BRENT);
					opt_set_data(opt, &data);
					opt_set_objective_function(opt, _func_gamma_fixed_mode);
					opt_set_parameters(opt, ps);
					double min;
					opt_optimize(opt, ps, &min);
					
					shape = Parameters_value(ps, 0);
					rate = (shape - 1)/map;
					
					free_Optimizer(opt);
				}
				free_Parameters(ps);
			}
			logP -= log(gsl_ran_gamma_pdf(map, shape, 1.0/rate));
			continue;
		}
		logP -= log(gsl_ran_lognormal_pdf(map, mu, sigma));
        
        if(dm != NULL){
            mus[i] = mu;
			sigmas[i] = sigma;
        }
	}

	if(dm != NULL){
		Parameter_set_values(Parameters_at(dm->parameters, 0), mus);
		Parameter_set_values(Parameters_at(dm->parameters, 1), sigmas);
	}
	free(mus);
	free(sigmas);
	free(d2logP_all);
	free(dlogP_all);

	free(x);
	free(y);
	free(yy);

	printf("Lognormal Laplace: %f\n", logP);
	return logP;
}


double calculate_laplace_lognormal2(Laplace* laplace, DistributionModel* dm){
	//	sigma = sqrt(-1/(f''(m)*m^2))
	//	mu    = log(m)+sigma^2
	Model* posterior = laplace->model;
	Parameters* parameters = dm->x;
	size_t dim = Parameters_size(parameters);
	double logP = 0;
	posterior->logP(posterior);
	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	double* mus = malloc(sizeof(double)*dim);
	double* sigmas = malloc(sizeof(double)*dim);
	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);
	posterior->hessian(posterior, parameters, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, parameters, 1e-5, dlogP_all);

	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* parameter = laplace_element(parameters, i, &local);
		double map = Parameter_value_at(parameter, local);
		double d2logP = d2logP_all[i];

		double sigma = sqrt(-1.0/(d2logP*map*map));
		double mu = log(map) + sigma*sigma;
		if (map < 1.e-6 || d2logP >= 0 || mu > 5) {
			double rate = map * -d2logP;
			double shape = rate*map + 1;
			// Very small branch -> exponential shape
			if (map < 1.e-6 || d2logP >= 0) {
				double dlogP = dlogP_all[i];
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameter_set_value_at(parameter, x[j], local);
					y[j] = posterior->logP(posterior);
				}
				Parameter_set_value_at(parameter, map, local);
				y[0] = posterior->logP(posterior);
				
				double maxY = y[0];
				for (int j = 0; j < N; j++) {
					y[j] -= maxY;
				}
				double lower = 0.001;
				double upper = 1;
				double guess = 1.0 - 0.001;
				Parameters* ps = new_Parameters(1);
				Parameters_move(ps, new_Parameter("", guess, new_Constraint(lower, upper)));
				
				struct laplace_data_t data = {rate, x, y, yy, N};
				double fx = _func_gamma_fixed_shape(ps, NULL, &data);
				
				Parameters_set_value(ps, 0, lower);
				double fa = _func_gamma_fixed_shape(ps, NULL, &data);
				Parameters_set_value(ps, 0, upper);
				double fb = _func_gamma_fixed_shape(ps, NULL, &data);
				Parameters_set_value(ps, 0, guess);
				
				if(fa > fx && fx < fb){
					Optimizer* opt = new_Optimizer(OPT_BRENT);
					opt_set_data(opt, &data);
					opt_set_objective_function(opt, _func_gamma_fixed_shape);
					opt_set_parameters(opt, ps);
					double min;
					opt_optimize(opt, ps, &min);
					
					shape = Parameters_value(ps, 0);
					
					free_Optimizer(opt);
				}
				free_Parameters(ps);
			}
			// Small branch with a maximum and spurious large variance
			else if(shape/(rate*rate) > 0.1 && map < 0.0001){
				double dlogP = dlogP_all[i];
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameter_set_value_at(parameter, x[j], local);
					y[j] = posterior->logP(posterior);
				}
				Parameter_set_value_at(parameter, map, local);
				y[0] = posterior->logP(posterior);
				
				double maxY = y[0];
				for (int j = 0; j < N; j++) {
					y[j] -= maxY;
				}
				
				double guess = 1.0 + 0.001;
				Parameters* ps = new_Parameters(1);
				Parameters_move(ps, new_Parameter("", guess, new_Constraint(1, 100)));
				
				struct laplace_data_t data = {map, x, y, yy, N};
				double fx = _func_gamma_fixed_mode(ps, NULL, &data);
				
				Parameters_set_value(ps, 0, 1);
				double fa = _func_gamma_fixed_mode(ps, NULL, &data);
				Parameters_set_value(ps, 0, 100);
				double fb = _func_gamma_fixed_mode(ps, NULL, &data);
				Parameters_set_value(ps, 0, guess);
				
				if(fa > fx && fx < fb){
					Optimizer* opt = new_Optimizer(OPT_BRENT);
					opt_set_data(opt, &data);
					opt_set_objective_function(opt, _func_gamma_fixed_mode);
					opt_set_parameters(opt, ps);
					double min;
					opt_optimize(opt, ps, &min);
					
					shape = Parameters_value(ps, 0);
					rate = (shape - 1)/map;
					
					free_Optimizer(opt);
				}
				free_Parameters(ps);
			}
			logP -= log(gsl_ran_gamma_pdf(map, shape, 1.0/rate));
			continue;
		}
		logP -= log(gsl_ran_lognormal_pdf(map, mu, sigma));
		
		mus[i] = mu;
		sigmas[i] = sigma;
	}
	

	Parameter_set_values(Parameters_at(dm->parameters, 0), mus);
	Parameter_set_values(Parameters_at(dm->parameters, 1), sigmas);

	free(mus);
	free(sigmas);
	free(d2logP_all);
	free(dlogP_all);
	free(x);
	free(y);
	free(yy);

	printf("Lognormal Laplace: %f\n", logP);
	return logP;
}

double _func_betaprime( Parameters *params, double *grad, void *data ){
	struct laplace_data_t* d = (struct laplace_data_t*)data;
	double beta = Parameters_value(params, 0);
	double alpha = d->map*(beta+1.0) + 1.0;
	double sumYY = 0;
	for (size_t i = 0; i < d->N; i++) {
		d->yy[i] = dbetaprime(d->x[i], alpha, beta);
		sumYY += d->yy[i];
	}
	double sum = 0;
	for (size_t i = 0; i < d->N; i++) {
		sum += pow(d->yy[i]/sumYY/d->N - d->y[i], 2);
	}
	return sum;
}

double calculate_laplace_betaprime(Laplace* laplace){
//	alpha = 1 - f''(m) * (m^2) * (m + 1)
//	beta = -f''(m) * m * (m + 1) - 1
	Model* posterior = laplace->model;
    Model* refdist = laplace->refdist;
    DistributionModel* dm = NULL;double* alphas = NULL;
	double* betas = NULL;

	size_t dim = Parameters_size(laplace->parameters);
    if(refdist != NULL){
        dm = refdist->obj;
		alphas = malloc(sizeof(double)*dim);
		betas = malloc(sizeof(double)*dim);
    }
	double* d2logP_all = dvector(dim);
	double* dlogP_all = dvector(dim);
	double logP = posterior->logP(posterior);
	posterior->hessian(posterior, laplace->parameters, HESSIAN_DIAGONAL, d2logP_all);
	laplace_first_derivatives(posterior, laplace->parameters, 1e-5, dlogP_all);
	for (size_t i = 0; i < dim; i++) {
		size_t local;
		Parameter* param = laplace_element(laplace->parameters, i, &local);
		double map = Parameter_value_at(param, local);
		double d2logP = d2logP_all[i];
		double alpha = 1.0 - d2logP*(map*map)*(map + 1.0);
		double beta = -d2logP*map*(map + 1.0) - 1.0;

        if (beta < 0) {
			double dlogP = dlogP_all[i];
			beta = fabs(dlogP) - 1;
			alpha = 1;

			if (beta < 2) {
				double* x = log_spaced_spaced_vector(map, 0.5, 10);
				double y[10];
				double sumY = 0;
				for (int j = 0; j < 10; j++) {
					Parameter_set_value_at(param, x[j], local);
					y[j] = laplace->model->logP(laplace->model);
					sumY += y[j];
				}
				double minY = dmin_vector(y, 10);
				for (int j = 0; j < 10; j++) {
					y[j] -= minY;
					y[j] /= sumY/9;
				}
				Parameters* ps = new_Parameters(1);
				Parameters_move(ps, new_Parameter("", 2, new_Constraint(2, 100)));
				double yy[9];
				struct laplace_data_t data = {map, x, y, yy, 9};
				Optimizer* opt = new_Optimizer(OPT_BRENT);
				opt_set_data(opt, &data);
				opt_set_objective_function(opt, _func_betaprime);
				opt_set_parameters(opt, ps);
				double min;
				opt_optimize(opt, ps, &min);

				beta = Parameters_value(ps, 0);
				alpha = map*(beta+1) + 1;
				Parameter_set_value_at(param, map, local);

				free_Optimizer(opt);
				free_Parameters(ps);
				free(x);
			}
		}
		logP -= dlogbetaprime(map, alpha, beta);

        if(dm != NULL){
            alphas[i] = alpha;
			betas[i] = beta;
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

	printf("Beta' Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma_from_mcmc(Laplace* laplace){
    Model* posterior = laplace->model;
    Model* empirical = laplace->empirical;
    DistributionModel* dm = empirical->obj;
	const double *alphas = Parameter_values(Parameters_at(dm->parameters, 0));
	const double *betas = Parameter_values(Parameters_at(dm->parameters, 1));

	for (int i = 0; i< Parameters_count(laplace->parameters); i++) {
//		printf("mean: %f MLE: %f mode: %f diff: %f [%f %f] var: %f skew: %f\n", alpha/beta, Parameters_value(laplace->parameters, i), (alpha-1.0)/beta,
//			   Parameters_value(laplace->parameters, i)- (alpha-1.0)/beta, alpha, beta, alpha/beta/beta, 2.0/(sqrt(alpha)));
		double mode = (alphas[i]-1.0)/betas[i];
		if (mode < 0) {
			mode = 1.e-6;
		}
		Parameters_set_value(laplace->parameters, i, mode);
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
	
	for (int i = 0; i< Parameters_count(laplace->parameters); i++) {
		//		printf("mean: %f MLE: %f mode: %f diff: %f [%f %f] var: %f skew: %f\n", alpha/beta, Parameters_value(laplace->parameters, i), exp(mu-sigma*sigma),
		//			   Parameters_value(laplace->parameters, i)- exp(mu-sigma*sigma), mu, sigma, (exp(sigma*sigma)-1.0)*exp(2.0*mu+sigma*sigma), exp(sigma*sigma + 2.0)*sqrt(exp(sigma*sigma) - 1.0));
		double mode = exp(mus[i]-sigmas[i]*sigmas[i]);
		Parameters_set_value(laplace->parameters, i, mode);
	}
	double logP = posterior->logP(posterior);
	double logQ = empirical->logP(empirical);
	double logLaplace = logP - logQ;
	
	printf("Lognormal MCMC Laplace: %f\n", logLaplace);
	return logLaplace;
}

double calculate_laplace(Laplace* laplace){
	CompoundModel* cm = laplace->refdist->obj;
	double logP = 0;
	for (int i = 0; i < cm->count; i++) {
		DistributionModel* dm = cm->models[i]->obj;
		if (dm->type == DISTRIBUTION_GAMMA) {
			logP += calculate_laplace_gamma2(laplace, dm);
		}
		else if (dm->type == DISTRIBUTION_BETA) {
			logP += calculate_laplace_beta2(laplace, dm);
		}
		else if (dm->type == DISTRIBUTION_LOGNORMAL) {
			logP += calculate_laplace_lognormal2(laplace, dm);
		}
		else{
			fprintf(stderr, "Distribution not supported by laplace\n");
			exit(2);
		}
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

Laplace* new_Laplace_from_json2(json_node* node, Hashtable* hash){
    static const json_field schema[] = {
        {"distribution", JSON_REQUIRED, JSON_OBJECT | JSON_STRING},
        {"model", JSON_REQUIRED, JSON_STRING},
        {"x", JSON_REQUIRED, JSON_ANY},
    };
    json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
    
    char* id = get_json_node_value_string(node, "id");
    char* ref_model = get_json_node_value_string(node, "model");
    json_node* dist_node = get_json_node(node, "distribution");
    json_node* x_node = get_json_node(node, "x");
    Laplace* laplace = malloc(sizeof(Laplace));
    laplace->model = Hashtable_get(hash, ref_model+1);
    laplace->model->ref_count++;
    laplace->parameters = NULL;
    laplace->refdist = NULL;
    laplace->empirical = NULL;
    laplace->free = _free_Laplace;
    
    
    if(dist_node->node_type == MJSON_STRING){
        char* r = get_json_node_value_string(node, "distribution");
        laplace->refdist = Hashtable_get(hash, r+1);
        laplace->refdist->ref_count++;
    }
    else{
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
            laplace->refdist = new_BetaDistributionModel_from_json(dist_node, hash);
        }
        else if(strcasecmp(d_string, "multivariatenormal") == 0){
            laplace->refdist = new_MultivariateNormalDistributionModel_from_json(dist_node, hash);
        }
        else{
            fprintf(stderr, "Distribution unknown: %s in Laplace with ID: %s", d_string, id);
            exit(1);
        }
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
        fprintf(stderr, "Distribution unknown in Laplace with ID: %s", id);
        exit(1);
    }
    
    laplace->parameters = distmodel_get_x(id, x_node, hash);
    
    return laplace;
}


Laplace* new_Laplace_from_json(json_node* node, Hashtable* hash){
    json_node* dist_node = get_json_node(node, "distribution");
    if(dist_node->node_type != MJSON_STRING){
        return new_Laplace_from_json2(node, hash);
    }

	static const json_field schema[] = {
	    {"distribution", JSON_OPTIONAL, JSON_ANY},
	    {"empirical", JSON_OPTIONAL, JSON_ANY},
	    {"model", JSON_OPTIONAL, JSON_ANY},
	    {"parameters", JSON_OPTIONAL, JSON_ANY},
	    {"ref", JSON_OPTIONAL, JSON_ANY},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	
	char* ref_model = get_json_node_value_string(node, "model");
	char* dist_string = get_json_node_value_string(node, "distribution");
	json_node* ref = get_json_node(node, "ref");
	Laplace* laplace = malloc(sizeof(Laplace));
	laplace->model = Hashtable_get(hash, ref_model+1);
	laplace->model->ref_count++;
	laplace->parameters = new_Parameters(1);
	laplace->refdist = NULL;
	laplace->empirical = NULL;
	laplace->free = _free_Laplace;
	json_node* empirical = get_json_node(node, "empirical");
	if (ref != NULL) {
		char* r = get_json_node_value_string(node, "ref");
		laplace->refdist = Hashtable_get(hash, r+1);
		laplace->refdist->ref_count++;
	}
	// distribution is inferred from the reference distribution
	// it can be a mixture of distributions
	if (dist_string == NULL) {
		laplace->calculate = calculate_laplace;
		return laplace;
	}
	
	get_parameters_references(node, hash, laplace->parameters);
	
	if(strcasecmp(dist_string, "gamma") == 0){
		if(empirical != NULL){
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
				exit(10);
			}
			laplace->empirical = empiricalDist;
			laplace->calculate = calculate_laplace_gamma_from_mcmc;
		}
		else{
			laplace->calculate = calculate_laplace_gamma;
		}
	}
	else if(strcasecmp(dist_string, "lognormal") == 0){
		if(empirical != NULL){
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
				exit(10);
			}
			laplace->empirical = empiricalDist;
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
