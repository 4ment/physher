//
//  klpq.c
//  physher
//
//  Created by Mathieu Fourment on 8/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "klpq.h"

#include <string.h>

#include "matrix.h"
#include "gaussian.h"
#include "transforms.h"
#include "klqp.h"


#include "bound.h"
#include "compoundmodel.h"
#include "distmodelfactory.h"
#include "treelikelihood.h"
#include "matrix.h"


double _klpq_snis_logP(Bound* klpq){
	size_t infCount = 0;
    Model** distributions = &klpq->variational;
	size_t distCount = 1;
	if(klpq->variational->type == MODEL_COMPOUND){
		CompoundModel* cm = klpq->variational->obj;
		distributions = cm->models;
		distCount = cm->count;
	}

	double* logw = dvector(klpq->samples);
	double sum = -DBL_MAX;
	double estimate = 0;

    for(size_t i = 0; i < klpq->samples; i++){
		klpq->variational->sample(klpq->variational);
		double logQ = klpq->variational->logP(klpq->variational);
		double logP = klpq->joint->logP(klpq->joint);
		if(!isinf(logP) && !isnan(logP)){
			logw[i] = logP - logQ;
		}
		else{
			infCount++;
			i--;
		}
		if (infCount == 10) {
			free(logw);
			return NAN;
		}
		sum = logaddexp(logw[i], sum);
    }
	for (size_t i = 0; i < klpq->samples; i++) {
		estimate += exp(logw[i] - sum) * logw[i];
	}
	free(logw);
    return estimate;
}

double _klpq_snis_gradient(Bound* klpq, Parameters* parameters){
	size_t infCount = 0;
    Model** distributions = &klpq->variational;
	size_t distCount = 1;
	if(klpq->variational->type == MODEL_COMPOUND){
		CompoundModel* cm = klpq->variational->obj;
		distributions = cm->models;
		distCount = cm->count;
	}
	size_t totalSize = 0;
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		totalSize += Parameter_size(Parameters_at(parameters, i));
	}
	double* grads = dvector(totalSize*klpq->samples); 

	double* logw = dvector(klpq->samples);
	double sum = -DBL_MAX;
	double estimate = 0;
	size_t index = 0;

    for(size_t i = 0; i < klpq->samples; i++){
		klpq->variational->sample(klpq->variational);
		double logQ = klpq->variational->logP(klpq->variational);
		double logP = klpq->joint->logP(klpq->joint);

		Parameters_zero_grad(parameters);
		klpq->variational->gradient(klpq->variational, parameters);

		for(size_t k = 0; k < Parameters_count(parameters); k++){
			Parameter* p = Parameters_at(parameters, k);
			for(size_t j = 0; j < Parameter_size(p); j++){
				grads[index++] = p->grad[j];
			}
		}

		if(!isinf(logP) && !isnan(logP)){
			logw[i] = logP - logQ;
		}
		else{
			infCount++;
			i--;
		}
		if (infCount == 10) {
			free(logw);
			return NAN;
		}
		sum = logaddexp(logw[i], sum);
    }

	Parameters_zero_grad(parameters);
	index = 0;
	for (size_t i = 0; i < klpq->samples; i++) {
		double w = exp(logw[i] - sum);
		for(size_t k = 0; k < Parameters_count(parameters); k++){
			Parameter* p = Parameters_at(parameters, k);
			for(size_t j = 0; j < Parameter_size(p); j++){
				p->grad[j] -= grads[index++] * w;
			}
		}
		estimate += w * logw[i];
	}
	free(logw);
	free(grads);
    return -estimate;
}

Model* new_KLpqBound_from_json(json_node* node, Hashtable* hash) {
    Model* model = new_AbstractBoundModel_from_json(node, hash);
    Bound* bound = model->obj;
    bound->logP = _klpq_snis_logP;
    bound->gradient = _klpq_snis_gradient;
    return model;
}

double variational_klpq_elbo(variational_t* var){
    double elbo = 0;
    Model* posterior = var->posterior;
    int inf_count = 0;
	double sum = -DBL_MAX;
	double* logw = dvector(var->elbo_samples); // log weights
    
    for (int i = 0; i < var->elbo_samples; i++) {
        double jacobian = 0.0;
        double logQ = 0.0;
        for(int j = 0; j < var->block_count; j++){
            variational_block_t* block = var->blocks[j];
            block->sample1(block, &jacobian);
            logQ += block->logQ(block, Vector_data(block->etas));
            
        }
        double logP = posterior->logP(posterior);
		if(!isinf(logP)){
			logw[i] = logP -logQ + jacobian;
			sum = logaddexp(logw[i], sum);
		}
        else inf_count++;
        
        if(isinf(elbo) || isnan(elbo)){
			free(logw);
            return elbo;
        }
    }
	
	for (int i = 0; i < var->elbo_samples; i++) {
		elbo += exp(logw[i]-sum) * logw[i];
	}
	free(logw);
    return -elbo;
}

double klpq_normal_meanfield(variational_t* var){
	if (var->initialized == false) {
		klqp_meanfield_normal_init(var);
		var->initialized = true;
	}
	double elbo = 0;
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;

	int inf_count = 0;
	double* logw = dvector(var->elbo_samples); // log weights
	double* samples = dvector(dim);
	double sum = -DBL_MAX;
	
	for (int i = 0; i < var->elbo_samples; i++) {
		double jacobian = 0.0;
		
		for (int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
			
			if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
			
			double mu;
			double sigma = exp(Parameter_value(var_p_sigma));
			
			if(Parameter_estimate(var_p_mu)){
				mu = Parameter_value(var_p_mu);
			}
			else{
				//				if(Parameter_value(p) < 1.0e-6){
				mu = -10 + sigma * sigma;
				//				}
				//				else{
				//					mu = log(Parameter_value(p)) + sigma*sigma;
				//				}
			}
			double zeta = rnorm() * sigma + mu;
			double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), &jacobian);
			Parameter_set_value(p, theta);
			samples[j] = theta;
		}
		
		double logP = posterior->logP(posterior) + jacobian;
		double logQ = var->logP(var, samples);
//		printf("%f %f\n", logP, logQ);
		logw[i] = logP-logQ;
		sum = logaddexp(logw[i], sum);
			
		if(isinf(elbo) || isnan(elbo)){
			free(logw);
			free(samples);
			return elbo;
		}
	}
	for (int i = 0; i < var->elbo_samples; i++) {
		elbo += exp(logw[i]-sum) * logw[i];
	}
//	print_dvector(logw, var->elbo_samples);
//	exit(11);
	free(logw);
	free(samples);
	return -elbo; // we want to minimize with a maximizer
}

void grad_klpq_normal_meanfield(variational_t* var, const Parameters* parameters, double* grads){
	if (var->initialized == false) {
		klqp_meanfield_normal_init(var);
		// save for later
		//		for (int i = 0; i < Parameters_count(var->parameters); i++) {
		//			Parameter_store(Parameters_at(var->parameters, i));
		//		}
		var->initialized = true;
	}
	
	Model* posterior = var->posterior;
	size_t dim = Parameters_count(var->parameters);
	
	double* samples = dvector(dim);
	double* zeta = dvector(dim*var->grad_samples);
	memset(grads, 0, sizeof(double)*dim*2);
	
	double* logw = dvector(var->grad_samples); // log weights
	double sum = -DBL_MAX;
	
	for (int i = 0; i < var->grad_samples; i++) {
		double jacobian = 0.0;
		for ( int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
			
			if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
			
			double mu;
			double sigma = exp(Parameter_value(var_p_sigma));
			
			if(Parameter_estimate(var_p_mu)){
				mu = Parameter_value(var_p_mu);
			}
			else{
				//				if(Parameter_value(p) < 1.0e-6){
				mu = -10 + sigma * sigma;
				//				}
				//				else{
				//					mu = log(Parameter_value(p)) + sigma*sigma;
				//				}
			}
			
			zeta[i*dim+j] = rnorm() * sigma + mu;
			
			double theta = inverse_transform(zeta[i*dim+j], Parameter_lower(p), Parameter_upper(p), &jacobian);
			//						printf("%f %f %f %d\n", mu, sigma, theta, Parameters_at(var->var_parameters, j)->estimate);
			Parameter_set_value(p, theta);
			samples[j] = theta;
		}
		
		double logP = posterior->logP(posterior) + jacobian;
		double logQ = var->logP(var, samples);
		
		logw[i] = logP-logQ;
		sum = logaddexp(logw[i], sum);
	}
	
	for (int i = 0; i < var->grad_samples; i++) {
		
		double w = exp(logw[i]-sum);
		
		for ( int j = 0; j < dim; j++) {
			Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
		
			if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
			
			double mu = Parameter_value(var_p_mu);
			double sigma = exp(Parameter_value(var_p_sigma));
			
			double sigma2 = sigma*sigma;
			double x = zeta[i*dim+j];
			grads[j] += w * (x - mu)/sigma2;
			grads[dim+j] += w * (mu*mu - sigma2 -2.0*mu*x + x*x)/(sigma2*sigma);
		}
	}
//	print_dvector(samples, dim);
//	print_dvector(logw, dim);
//	print_dvector(grads, 2*dim);
//	exit(11);
	free(samples);
	free(zeta);
	free(logw);
	
	if (var->file != NULL) {
		fprintf(var->file, "%zu", var->iter++);
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%e", Parameters_value(var->var_parameters, i));
		}
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%e", grads[i]);
		}
		fprintf(var->file, "\n");
	}
}

