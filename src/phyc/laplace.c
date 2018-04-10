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
#include "optimize.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

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
			Parameters_set_value(dm->parameters, i*2, shape);
			Parameters_set_value(dm->parameters, i*2+1, rate);
		}
	}
	
	free(x);
	free(y);
	free(yy);
	
	printf("Gamma Laplace: %f\n", logP);
	return logP;
}

double calculate_laplace_lognormal(Laplace* laplace){
	//	sigma = sqrt(-1/(f''(m)*m^2))
	//	mu    = log(m)+sigma^2
	Model* posterior = laplace->model;
	double logP = posterior->logP(posterior);
	for (int i = 0; i < Parameters_count(laplace->parameters); i++) {
		double map = Parameters_value(laplace->parameters, i);
		double d2logP = laplace->model->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
		
		double sigma = sqrt(-1.0/(d2logP*map*map));
		double mu = log(map) + sigma*sigma;
		if (map < 1.e-4 || d2logP >= 0 || mu > 5) {
			double dlogP = posterior->dlogP(posterior, Parameters_at(laplace->parameters, i));
			logP -= log(gsl_ran_gamma_pdf(map, 1.0, 1.0/fabs(dlogP)));
		}
		else{
			logP -= log(gsl_ran_lognormal_pdf(map, mu, sigma));
		}
	}
	
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
	double logP = posterior->logP(posterior);
	for (int i = 0; i < Parameters_count(laplace->parameters); i++) {
		double map = Parameters_value(laplace->parameters, i);
		double d2logP = laplace->model->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
//		printf("%d] %f %f %f\n", i, d2logP, Model_second_derivative(posterior, Parameters_at(laplace->parameters, i), NULL, 0.0001), map);
		double alpha = 1.0 - d2logP*(map*map)*(map + 1.0);
		double beta = -d2logP*map*(map + 1.0) - 1.0;
//		printf("%f %f %f %s %f %f\n", map, laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i)), d2logP, Parameters_name(laplace->parameters, i), alpha, beta);
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
	}
	
	printf("Beta' Laplace: %f\n", logP);
	return logP;
}

void _free_Laplace(Laplace* laplace){
	free_Parameters(laplace->parameters);
	laplace->model->free(laplace->model);
	if(laplace->refdist != NULL)laplace->refdist->free(laplace->refdist);
	free(laplace);
}

Laplace* new_Laplace_from_json(json_node* node, Hashtable* hash){
	char* ref_model = get_json_node_value_string(node, "model");
	char* dist_string = get_json_node_value_string(node, "distribution");
	json_node* ref = get_json_node(node, "ref");
	Laplace* laplace = malloc(sizeof(Laplace));
	laplace->model = Hashtable_get(hash, ref_model+1);
	laplace->model->ref_count++;
	laplace->parameters = new_Parameters(1);
	laplace->refdist = NULL;
	get_parameters_references(node, hash, laplace->parameters);
	if (ref != NULL) {
		char* r = get_json_node_value_string(node, "ref");
		laplace->refdist = Hashtable_get(hash, r+1);
		laplace->refdist->ref_count++;
	}
	if(strcasecmp(dist_string, "gamma") == 0){
		laplace->calculate = calculate_laplace_gamma;
	}
	else if(strcasecmp(dist_string, "lognormal") == 0){
		laplace->calculate = calculate_laplace_lognormal;
	}
	else if(strcasecmp(dist_string, "betaprime") == 0){
		laplace->calculate = calculate_laplace_betaprime;
	}
	laplace->free = _free_Laplace;
	return laplace;
}
