//
//  laplace.c
//  physher
//
//  Created by Mathieu Fourment on 17/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "laplace.h"

#include <strings.h>

#include "gamma.h"
#include "beta.h"
#include "distmodel.h"
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
	if(refdist != NULL){
		dm = refdist->obj;
	}
	double logP = posterior->logP(posterior);
	
	for (int i = 0; i < Parameters_count(laplace->parameters); i++) {
		double map = Parameters_value(laplace->parameters, i);
		double d2logP = posterior->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
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
			Parameters_set_value(dm->parameters[0], i, alpha);
			Parameters_set_value(dm->parameters[1], i, beta);
		}
	}
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
	
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameter* parameter = Parameters_at(parameters, i);
		printf("beta %s\n",parameter->name);
		double map = Parameter_value(parameter);
		double dlogP;// = posterior->dlogP(laplace->model, parameter);
		double d2logP = Model_second_derivative(posterior, parameter, &dlogP, 1e-5);// posterior->d2logP(laplace->model, parameter);
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
		
		Parameters_set_value(dm->parameters[0], i, alpha);
		Parameters_set_value(dm->parameters[1], i, beta);
//		printf("map: %f dlogP: %f d2logP: %f [%f,%f]\n",map, dlogP, d2logP, alpha, beta);
	}
	
	printf("Beta Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma(Laplace* laplace){
	// beta = rate = -f''(m) * m
	// alpha = shape = rate * m + 1
	Model* posterior = laplace->model;
	Model* refdist = laplace->refdist;
	DistributionModel* dm = NULL;
	if(refdist != NULL){
		dm = refdist->obj;
	}
	
	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	
	double logP = posterior->logP(posterior);
	for (int i = 0; i < Parameters_count(laplace->parameters); i++) {
		double map = Parameters_value(laplace->parameters, i);
		double d2logP = laplace->model->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
		double rate = map * -d2logP;
		double shape = rate*map + 1;
		
		// Very small branch -> exponential shape
		if (map < 1.e-6 || d2logP >= 0) {
			double dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
			shape = 1;
			rate = fabs(dlogP);
			
			log_spaced_spaced_vector2(x, map, 0.5, N);
			
			for (size_t j = 1; j < N; j++) {
				Parameters_set_value(laplace->parameters, i, x[j]);
				y[j] = laplace->model->logP(laplace->model);
			}
			Parameters_set_value(laplace->parameters, i, map);
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
			double dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
			shape = 1;
			rate = fabs(dlogP);
			
			log_spaced_spaced_vector2(x, map, 0.5, N);
			
			for (size_t j = 1; j < N; j++) {
				Parameters_set_value(laplace->parameters, i, x[j]);
				y[j] = laplace->model->logP(laplace->model);
			}
			Parameters_set_value(laplace->parameters, i, map);
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
			Parameters_set_value(dm->parameters[0], i, shape);
			Parameters_set_value(dm->parameters[1], i, rate);
		}
	}
	
	free(x);
	free(y);
	free(yy);
	
	printf("Gamma Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma2(Laplace* laplace, DistributionModel* dm){
	// beta = rate = -f''(m) * m
	// alpha = shape = rate * m + 1
	Model* posterior = laplace->model;
	Parameters* parameters = dm->x;
	
	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	
	double logP = 0;
	posterior->logP(posterior);// initialize
	
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameter* parameter = Parameters_at(parameters, i);
		printf("gamma %s\n",parameter->name);
		double map_orig = Parameter_value(parameter);
		double map = map_orig - dm->shift;
		double dlogP;
		double d2logP = posterior->d2logP(posterior, parameter);
//		double d2logP = Model_second_derivative(posterior, parameter, &dlogP, 1e-5);
		double rate = map * -d2logP;
		double shape = rate*map + 1;
		
		// Very small branch -> exponential shape
		if (map < 1.e-6 || d2logP >= 0) {
			double dlogP = posterior->dlogP(posterior, parameter);
			shape = 1;
			rate = fabs(dlogP);
			
			log_spaced_spaced_vector2(x, map, 0.5, N);
			
			for (size_t j = 1; j < N; j++) {
				Parameter_set_value(parameter, x[j] + dm->shift);
				y[j] = posterior->logP(posterior);
			}
			Parameter_set_value(parameter, map_orig);
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
			dlogP = posterior->dlogP(posterior, parameter);
			shape = 1;
			rate = fabs(dlogP);
			
			log_spaced_spaced_vector2(x, map, 0.5, N);
			
			for (size_t j = 1; j < N; j++) {
				Parameter_set_value(parameter, x[j] + dm->shift);
				y[j] = posterior->logP(posterior);
			}
			Parameter_set_value(parameter, map_orig);
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

		Parameters_set_value(dm->parameters[0], i, shape);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
			rate = 1.0/rate;
		}
		Parameters_set_value(dm->parameters[1], i, rate);
		logP -= log(gsl_ran_gamma_pdf(map, shape, rate));
	}
	
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
    
	size_t paramCount = Parameters_count(laplace->parameters);
	gsl_matrix* H = gsl_matrix_alloc(paramCount, paramCount);
	gsl_vector* mu = gsl_vector_alloc(paramCount);
	gsl_matrix * L = gsl_matrix_alloc (paramCount, paramCount);
	gsl_permutation* perm = gsl_permutation_alloc (paramCount);
	gsl_vector* work = gsl_vector_alloc(paramCount);
	gsl_set_error_handler_off();
	double logP = laplace->model->logP(laplace->model);
	double epsilon = 0.0001;

	for (int i = 0; i < paramCount; i++) {
		double mapi = Parameters_value(laplace->parameters, i);
		if (mapi < 1.0e-6) {
			mapi += epsilon;
		}
		double dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
		double d2logP = laplace->model->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
		double Hii = dlogP*mapi + d2logP*mapi*mapi;
		//printf("%f %f %f %f\n",mapi,dlogP,d2logP, Hii);
		gsl_matrix_set(H, i, i, Hii);
		gsl_vector_set(mu, i, log(mapi));
		
		for (int j = i+1; j < paramCount; j++) {
			double mapj = Parameters_value(laplace->parameters, j);
			if (mapj < 1.0e-6) {
				mapj += epsilon;
			}
//			double didj = Model_mixed_derivative(laplace->model, Parameters_at(laplace->parameters, i), Parameters_at(laplace->parameters, j));
			double didj = laplace->model->ddlogP(laplace->model, Parameters_at(laplace->parameters, i), Parameters_at(laplace->parameters, j));
//			printf("%f %f %f\n", didj, Parameters_value(laplace->parameters, i), Parameters_value(laplace->parameters, j));
			double Hij = didj * mapi * mapj;
			gsl_matrix_set(H, i, j, Hij);
			gsl_matrix_set(H, j, i, Hij);
		}
	}
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
		double val = Parameters_value(laplace->parameters, i);
		if (val < 1.0e-6) {
			val += epsilon;
		}
		logQ -= log(val);
//		printf("%f %e\n", log(Parameters_value(laplace->parameters, i)), Parameters_value(laplace->parameters, i));
	}
	
    if(refdist != NULL){
        if(paramCount*paramCount == Parameters_count(dm->parameters[1])){
            for (int i = 0; i < paramCount; i++) {
                Parameters_set_value(dm->parameters[0], i, Parameters_value(laplace->parameters, i));
                for (int j = 0; j < paramCount; j++) {
                    double sigma = gsl_matrix_get(H, i, j);
                    Parameters_set_value(dm->parameters[1], i*paramCount+j, sigma);
                }
            }
        }
        else{
            size_t row = 0;
            for(int i = 0; i < paramCount; i++){
                for(int j = 0; j <= i; j++){
                    Parameters_set_value(dm->parameters[1], row, gsl_matrix_get(L, i, j));
                    row++;
                }
            }
        }
    }
	printf("Multivariatenormal Laplace: %f logQ: %f %d\n", logP - logQ, logQ, gsl_ran_multivariate_gaussian_log_pdf(mu, mu, L, &logQ, work));

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
    if(refdist != NULL){
        dm = refdist->obj;
    }
    
	double logP = posterior->logP(posterior);
	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	
	for (int i = 0; i < Parameters_count(laplace->parameters); i++) {
		double map = Parameters_value(laplace->parameters, i);
		double d2logP = laplace->model->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
		
		double sigma = sqrt(-1.0/(d2logP*map*map));
		double mu = log(map) + sigma*sigma;
		if (map < 1.e-6 || d2logP >= 0 || mu > 5) {
			double rate = map * -d2logP;
			double shape = rate*map + 1;
			// Very small branch -> exponential shape
			if (map < 1.e-6 || d2logP >= 0) {
				double dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameters_set_value(laplace->parameters, i, x[j]);
					y[j] = laplace->model->logP(laplace->model);
				}
				Parameters_set_value(laplace->parameters, i, map);
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
				double dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameters_set_value(laplace->parameters, i, x[j]);
					y[j] = laplace->model->logP(laplace->model);
				}
				Parameters_set_value(laplace->parameters, i, map);
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
            Parameters_set_value(dm->parameters[0], i, mu);
            Parameters_set_value(dm->parameters[1], i, sigma);
        }
	}
	
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
	double logP = 0;
	posterior->logP(posterior);
	int N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameter* parameter = Parameters_at(parameters, i);
		double map = Parameter_value(parameter);
		double d2logP = posterior->d2logP(posterior, parameter);
		
		double sigma = sqrt(-1.0/(d2logP*map*map));
		double mu = log(map) + sigma*sigma;
		if (map < 1.e-6 || d2logP >= 0 || mu > 5) {
			double rate = map * -d2logP;
			double shape = rate*map + 1;
			// Very small branch -> exponential shape
			if (map < 1.e-6 || d2logP >= 0) {
				double dlogP = posterior->dlogP(posterior, parameter);
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameter_set_value(parameter, x[j]);
					y[j] = posterior->logP(posterior);
				}
				Parameter_set_value(parameter, map);
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
				double dlogP = posterior->dlogP(posterior, parameter);
				shape = 1;
				rate = fabs(dlogP);
				
				log_spaced_spaced_vector2(x, map, 0.5, N);
				
				for (size_t j = 1; j < N; j++) {
					Parameter_set_value(parameter, x[j]);
					y[j] = posterior->logP(posterior);
				}
				Parameter_set_value(parameter, map);
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
		
		Parameters_set_value(dm->parameters[0], i, mu);
		Parameters_set_value(dm->parameters[1], i, sigma);
	}
	
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
    DistributionModel* dm = NULL;
    if(refdist != NULL){
        dm = refdist->obj;
    }
	double logP = posterior->logP(posterior);
	for (int i = 0; i < Parameters_count(laplace->parameters); i++) {
		double map = Parameters_value(laplace->parameters, i);
		double d2logP = laplace->model->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
		double alpha = 1.0 - d2logP*(map*map)*(map + 1.0);
		double beta = -d2logP*map*(map + 1.0) - 1.0;

        if (beta < 0) {
			double dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
			beta = fabs(dlogP) - 1;
			alpha = 1;

			if (beta < 2) {
				double* x = log_spaced_spaced_vector(map, 0.5, 10);
				double y[10];
				double sumY = 0;
				for (int j = 0; j < 10; j++) {
					Parameters_set_value(laplace->parameters, i, x[j]);
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
				Parameters_set_value(laplace->parameters, i, map);

				free_Optimizer(opt);
				free_Parameters(ps);
				free(x);
			}
		}
		logP -= dlogbetaprime(map, alpha, beta);
        
        if(dm != NULL){
            Parameters_set_value(dm->parameters[0], i, alpha);
            Parameters_set_value(dm->parameters[1], i, beta);
        }
	}
	
	printf("Beta' Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_gamma_from_mcmc(Laplace* laplace){
    Model* posterior = laplace->model;
    Model* empirical = laplace->empirical;
    DistributionModel* dm = empirical->obj;
	
	for (int i = 0; i< Parameters_count(laplace->parameters); i++) {
		double alpha = Parameters_value(dm->parameters[0], i);
		double beta = Parameters_value(dm->parameters[1], i);
//		printf("mean: %f MLE: %f mode: %f diff: %f [%f %f] var: %f skew: %f\n", alpha/beta, Parameters_value(laplace->parameters, i), (alpha-1.0)/beta,
//			   Parameters_value(laplace->parameters, i)- (alpha-1.0)/beta, alpha, beta, alpha/beta/beta, 2.0/(sqrt(alpha)));
		double mode = (alpha-1.0)/beta;
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
	
	for (int i = 0; i< Parameters_count(laplace->parameters); i++) {
		double mu = Parameters_value(dm->parameters[0], i);
		double sigma = Parameters_value(dm->parameters[1], i);
		//		printf("mean: %f MLE: %f mode: %f diff: %f [%f %f] var: %f skew: %f\n", alpha/beta, Parameters_value(laplace->parameters, i), exp(mu-sigma*sigma),
		//			   Parameters_value(laplace->parameters, i)- exp(mu-sigma*sigma), mu, sigma, (exp(sigma*sigma)-1.0)*exp(2.0*mu+sigma*sigma), exp(sigma*sigma + 2.0)*sqrt(exp(sigma*sigma) - 1.0));
		double mode = exp(mu-sigma*sigma);
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
    char* allowed[] = {
        "distribution",
        "model",
        "x"
    };
    json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
    
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

	char* allowed[] = {
		"distribution",
		"empirical",
		"model",
		"parameters",
		"ref"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
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
