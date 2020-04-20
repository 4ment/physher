//
//  distexp.c
//  physher
//
//  Created by Mathieu Fourment on 15/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distexp.h"

#include <string.h>
#include <strings.h>

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"

double DistributionModel_log_exp(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	dm->lp = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = Parameters_value(dm->parameters[0], i);
			dm->lp += log(lambda) - lambda * Parameters_value(dm->x, i);
		}
	}
	else{
		double lambda = Parameters_value(dm->parameters[0], 0);
		dm->lp = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			dm->lp -= lambda * Parameters_value(dm->x, i);
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_exp_mean(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	dm->lp = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double mean = Parameters_value(dm->parameters[0], i);
			dm->lp += -log(mean) - Parameters_value(dm->x, i)/mean;
		}
	}
	else{
		double mean = Parameters_value(dm->parameters[0], 0);
		dm->lp = -log(mean) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			dm->lp -=  Parameters_value(dm->x, i)/mean;
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_exp_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = Parameters_value(dm->parameters[0], i);
			logP += log(lambda) - lambda * values[i];
		}
	}
	else{
		double lambda = Parameters_value(dm->parameters[0], 0);
		logP = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP -= lambda * values[i];
		}
	}
	return logP;
}

double DistributionModel_log_exp_with_values_mean(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double mean = Parameters_value(dm->parameters[0], i);
			logP += -log(mean) - values[i]/mean;
		}
	}
	else{
		double mean = Parameters_value(dm->parameters[0], 0);
		logP = -log(mean) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP -= values[i]/mean;
		}
	}
	return logP;
}

double DistributionModel_dlog_exp(DistributionModel* dm, const Parameter* p){
    if (p == Parameters_at(dm->parameters[0], 0)) {
        double dlogf = Parameters_count(dm->x)/Parameters_value(dm->parameters[0], 0);
        for(int i = 0; i < Parameters_count(dm->x); i++){
            dlogf -= Parameters_value(dm->x, i);
        }
        return dlogf;
    }
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (p == Parameters_at(dm->x, i)) {
            if(Parameters_count(dm->parameters[0]) == 1){
                return -Parameters_value(dm->parameters[0], 0);
            }
            return -Parameters_value(dm->parameters[0], i);
		}
	}
	return 0;
}

double DistributionModel_dlog_exp_mean(DistributionModel* dm, const Parameter* p){
    if (p == Parameters_at(dm->parameters[0], 0)) {
        double mean = Parameters_value(dm->parameters[0], 0);
        double dlogf = -Parameters_count(dm->x)/mean;
        double mean2 = mean;
        for(int i = 0; i < Parameters_count(dm->x); i++){
            dlogf += Parameters_value(dm->x, i)/mean2;
        }
        return dlogf;
    }
    
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (p == Parameters_at(dm->x, i)) {
            if(Parameters_count(dm->parameters[0]) == 1){
                return -1.0/Parameters_value(dm->parameters[0], 0);
            }
            return -1.0/Parameters_value(dm->parameters[0], i);
		}
	}
	return 0;
}

static void DistributionModel_exp_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, 1.0/Parameters_value(dm->parameters[0], i));
		}
	}
	else{
		double mean = 1.0/Parameters_value(dm->parameters[0], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, mean);
		}
	}
}

static void DistributionModel_exp_sample_mean(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, Parameters_value(dm->parameters[0], i));
		}
	}
	else{
		double mean = Parameters_value(dm->parameters[0], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, mean);
		}
	}
}


static double DistributionModel_exp_sample_evaluate(DistributionModel* dm){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, 1.0/Parameters_value(dm->parameters[0], i));
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		double mean = 1.0/Parameters_value(dm->parameters[0], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, mean);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_log_exp(dm);
}

static double DistributionModel_exp_sample_evaluate_mean(DistributionModel* dm){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, Parameters_value(dm->parameters[0], i));
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		double mean = Parameters_value(dm->parameters[0], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, mean);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_log_exp_mean(dm);
}

DistributionModel* new_ExponentialDistributionModel_with_parameters(Parameters** parameters, Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, 1, x);
	dm->type = DISTRIBUTION_EXPONENTIAL;
	dm->parameterization = parameterization;
	if(parameterization == DISTRIBUTION_EXPONENTIAL_RATE){
		dm->logP = DistributionModel_log_exp;
		dm->logP_with_values = DistributionModel_log_exp_with_values;
		dm->dlogP = DistributionModel_dlog_exp;
		dm->sample = DistributionModel_exp_sample;
		dm->sample_evaluate = DistributionModel_exp_sample_evaluate;
	}
	else{
		dm->logP = DistributionModel_log_exp_mean;
		dm->logP_with_values = DistributionModel_log_exp_with_values_mean;
		dm->dlogP = DistributionModel_dlog_exp_mean;
		dm->sample = DistributionModel_exp_sample_mean;
		dm->sample_evaluate = DistributionModel_exp_sample_evaluate_mean;
	}
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
    dm->shift = 0;
	return dm;
}

Model* new_ExponentialDistributionModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters** parameters = NULL;
    size_t parameters_dim = 1;
    distribution_parameterization parameterization = DISTRIBUTION_EXPONENTIAL_MEAN;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        size_t paramCount = Parameters_count(x);
        parameters = malloc(sizeof(Parameters*));
        parameters[0] = new_Parameters(paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            Parameters_move(parameters[0], new_Parameter("mean", m, new_Constraint(0, INFINITY)));
            free_Vector(samples[i]);
        }
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        parameters = malloc(sizeof(Parameters*));
        parameters[0] = new_Parameters(Parameters_count(x));
        for (int i = 0; i < Parameters_count(x); i++) {
            Parameters_move(parameters[0], new_Parameter("mean", 1, new_Constraint(0, INFINITY)));
        }
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        
        if (strcasecmp(parameters_node->children[0]->key, "lambda") != 0 && strcasecmp(parameters_node->children[0]->key, "mean") != 0) {
            fprintf(stderr, "Normal distribution should be parametrized with mean or lambda\n");
            exit(13);
        }
        
        if (strcasecmp(parameters_node->children[0]->key, "lambda") == 0){
            parameterization = DISTRIBUTION_EXPONENTIAL_RATE;
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
    }
    
    DistributionModel* dm = new_ExponentialDistributionModel_with_parameters(parameters, x, parameterization);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free(parameters);
	
	return model;
}
