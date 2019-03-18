//
//  distexp.c
//  physher
//
//  Created by Mathieu Fourment on 15/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distexp.h"

#include <string.h>

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"

double DistributionModel_log_exp(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = Parameters_value(dm->parameters, i);
			logP += log(lambda) - lambda * Parameters_value(dm->x, i);
		}
	}
	else{
		double lambda = Parameters_value(dm->parameters, 0);
		logP = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP -= lambda * Parameters_value(dm->x, i);
		}
	}
	return logP;
}

double DistributionModel_log_exp_mean(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = 1.0/Parameters_value(dm->parameters, i);
			logP += log(lambda) - lambda * Parameters_value(dm->x, i);
		}
	}
	else{
		double lambda = 1.0/Parameters_value(dm->parameters, 0);
		logP = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP -= lambda * Parameters_value(dm->x, i);
		}
	}
	return logP;
}

double DistributionModel_log_exp_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = Parameters_value(dm->parameters, i);
			logP += log(lambda) - lambda * values[i];
		}
	}
	else{
		double lambda = Parameters_value(dm->parameters, 0);
		logP = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP -= log(lambda) * values[i];
		}
	}
	return logP;
}

double DistributionModel_log_exp_with_values_mean(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = 1.0/Parameters_value(dm->parameters, i);
			logP += log(lambda) - lambda * values[i];
		}
	}
	else{
		double lambda = 1.0/Parameters_value(dm->parameters, 0);
		logP = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP -= log(lambda) * values[i];
		}
	}
	return logP;
}

double DistributionModel_dlog_exp(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			return -Parameters_value(dm->parameters, 0);
		}
	}
	return 0;
}

double DistributionModel_dlog_exp_mean(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			return -1.0/Parameters_value(dm->parameters, 0);
		}
	}
	return 0;
}

double DistributionModel_ddlog_exp(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	return 0.0;
}

double DistributionModel_d2log_exp(DistributionModel* dm, const Parameter* p){
	return 0;
}

static void DistributionModel_exp_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, 1.0/Parameters_value(dm->parameters, i));
		}
	}
	else{
		double mean = 1.0/Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, mean);
		}
	}
}

static void DistributionModel_exp_sample_mean(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, Parameters_value(dm->parameters, i));
		}
	}
	else{
		double mean = Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_exponential(dm->rng, mean);
		}
	}
}


static double DistributionModel_exp_sample_evaluate(DistributionModel* dm){
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, 1.0/Parameters_value(dm->parameters, i));
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		double mean = 1.0/Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, mean);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_log_exp(dm);
}

static double DistributionModel_exp_sample_evaluate_mean(DistributionModel* dm){
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, Parameters_value(dm->parameters, i));
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		double mean = Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_exponential(dm->rng, mean);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_log_exp_mean(dm);
}

DistributionModel* new_ExponentialDistributionModel(const double lambda, const Parameters* x){
	Parameters* ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter("exp.lambda", lambda, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->type = DISTRIBUTION_EXPONENTIAL;
	dm->parameterization = DISTRIBUTION_EXPONENTIAL_RATE;
	dm->logP = DistributionModel_log_exp;
	dm->logP_with_values = DistributionModel_log_exp_with_values;
	dm->dlogP = DistributionModel_dlog_exp;
	dm->d2logP = DistributionModel_d2log_exp;
	dm->ddlogP = DistributionModel_ddlog_exp;
	dm->sample = DistributionModel_exp_sample;
	dm->sample_evaluate = DistributionModel_exp_sample_evaluate;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_ExponentialDistributionModel_with_parameters(Parameters* parameters, const Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
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
	dm->d2logP = DistributionModel_d2log_exp;
	dm->ddlogP = DistributionModel_ddlog_exp;
	return dm;
}

Model* new_ExponentialDistributionModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	
	Parameters* x = new_Parameters(1);
	get_parameters_references2(node, hash, x, "x");
	
	char* file = get_json_node_value_string(node, "file");
	Vector** samples = NULL;
	if (file != NULL){
		get_parameters_references2(node, hash, x, "x");
		size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
		samples = read_log_for_parameters_t(file, burnin, x);
	}
	
	Parameters* parameters = new_Parameters(1);
	get_parameters_references(node, hash, parameters);

	double mean = false;
	json_node* p_node = get_json_node(node, "parameters");
	if (Parameters_count(parameters) == 1 && strcasecmp(p_node->children[0]->key, "mean") == 0) {
		Parameters_set_value(parameters, 0, 1.0/Parameters_value(parameters, 0));
		mean = true;
	}
	char* parameterization_string = get_json_node_value_string(node, "parameterization");
	if(parameterization_string!= NULL && strcasecmp(parameterization_string, "mean") == 0){
		mean = true;
	}
	
	// empirical
	if (samples != NULL) {
		size_t paramCount = Parameters_count(x);
		
		for (int i = 0; i < paramCount; i++) {
			double* vec = Vector_data(samples[i]);
			double m = dmean(vec, Vector_length(samples[i]));
			if(mean){
				Parameters_move(parameters, new_Parameter("lambda", m, NULL));
			}
			else{
				Parameters_move(parameters, new_Parameter("lambda", 1.0/m, NULL));
			}
			free_Vector(samples[i]);
		}
		free(samples);
	}
	else{
		for (int i = 0; i < Parameters_count(parameters); i++) {
			Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
		}
	}
	distribution_parameterization param = mean ? DISTRIBUTION_EXPONENTIAL_MEAN : DISTRIBUTION_EXPONENTIAL_RATE;
	DistributionModel* dm = new_ExponentialDistributionModel_with_parameters(parameters, x, param);
	
	dm->shift = get_json_node_value_double(node, "shift", 0);
	Model* model = new_DistributionModel2(id, dm);
	
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	free_Parameters(parameters);
	free_Parameters(x);
	
	return model;
}