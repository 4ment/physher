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
#include "distmodel.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

double calculate_laplace_gamma(Laplace* laplace){
	// beta = rate = -f''(m) * m
	// alpha = shape = rate * m + 1
	Model* posterior = laplace->model;
	Model* refdist = laplace->refdist;
	DistributionModel* dm = NULL;
	if(refdist != NULL){
		dm = refdist->obj;
	}
	
	double logP = posterior->logP(posterior);
	for (int i = 0; i < Parameters_count(laplace->parameters); i++) {
		double map = Parameters_value(laplace->parameters, i);
//		printf("d2logP: %e %f (%f)\n", d2logP, map, Model_second_derivative(posterior, Parameters_at(laplace->parameters, i), NULL, 0.001));
		double d2logP = laplace->model->d2logP(laplace->model, Parameters_at(laplace->parameters, i));
		
		if (map < 1.e-6 || d2logP >= 0) {
			double dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
			logP -= log(gsl_ran_gamma_pdf(map, 1.0, -1.0/dlogP));
//			printf("dlogP: %f %f\n", dlogP, Model_first_derivative(posterior, Parameters_at(laplace->parameters, i), 0.001));
			if(dm != NULL){
				Parameters_set_value(dm->parameters, i*2, 1.0);
				Parameters_set_value(dm->parameters, i*2+1, -dlogP);
			}
		}
		else{
			double rate = map * -d2logP;
			//logP += dloggamma(map, map*rate+1.0, rate); // does not always work :(
			logP -= log(gsl_ran_gamma_pdf(map, map*rate+1.0, 1.0/rate));
//			printf("%f map: %f mu: %f sigma: %f d2:%f %f\n", log(gsl_ran_gamma_pdf(map, map*rate+1.0, 1.0/rate)), map, map*rate+1.0, rate, d2logP, Model_second_derivative(posterior, Parameters_at(laplace->parameters, i), NULL, 0.00001));
			if(dm != NULL){
				Parameters_set_value(dm->parameters, i*2, map*rate+1.0);
				Parameters_set_value(dm->parameters, i*2+1, rate);
			}
		}
	}
	
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
//			printf("%e %f %f %s\n", map, dlogP, d2logP, Parameters_name(laplace->parameters, i));
//			if (beta < 0) {
//				double x = map*1000;//0.2*map + map;
//				Parameters_set_value(laplace->parameters, i, x);
//				dlogP = laplace->model->dlogP(laplace->model, Parameters_at(laplace->parameters, i));
//				beta = (dlogP*x*x + (dlogP+1.0)*x - map)/(map - x);
//				//alpha = map*(beta + 1.0) + 1.0;
//				alpha = ((dlogP*x*x + dlogP*x + 1.0)*map - x)/(map - x);
//				Parameters_set_value(laplace->parameters, i, map);
//				printf("%e %f dlogP: %f %f\n", alpha, beta, dlogP, (alpha-1.0)*log(map) - (alpha+beta)*log(1.0+map) - gsl_sf_lnbeta(alpha, beta));
//				double temp[20] = {0.00000001,0.0000001,0.0000002,0.0000005,0.000001,0.000002, 0.000005, 0.000010, 0.00001, 0.0001, 0.001, 0.003, 0.005, 0.01, 0.05, 0.06, 0.1, 0.2, 0.3, 0.5};
//				for (int j = 0; j < 20; j++) {
//					Parameters_set_value(laplace->parameters, i, temp[j]);
//					printf("%e,", laplace->model->logP(laplace->model));
//				}
//				printf("\n");
//				exit(1);
//			}
//			else{
//				double temp[20] = {0.00000001,0.0000001,0.0000002,0.0000005,0.000001,0.000002, 0.000005, 0.000010, 0.00001, 0.0001, 0.001, 0.003, 0.005, 0.01, 0.05, 0.06, 0.1, 0.2, 0.3, 0.5};
//				for (int j = 0; j < 20; j++) {
//					Parameters_set_value(laplace->parameters, i, temp[j]);
//					printf("%e,", laplace->model->logP(laplace->model));
//				}
//				printf("\n");
//				exit(1);
//			}
		}
		logP -= (alpha-1.0)*log(map) - (alpha+beta)*log(1.0+map) - gsl_sf_lnbeta(alpha, beta);
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
