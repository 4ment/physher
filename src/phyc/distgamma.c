//
//  distgamma.c
//  physher
//
//  Created by Mathieu Fourment on 30/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distgamma.h"


#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"

double DistributionModel_log_gamma(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	dm->lp = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(gsl_ran_gamma_pdf(x, alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(gsl_ran_gamma_pdf(x, alpha, beta));
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_gamma_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			logP += log(gsl_ran_gamma_pdf(values[i] - dm->shift, alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(gsl_ran_gamma_pdf(values[i] - dm->shift, alpha, beta));
		}
	}
	return logP;
}

double DistributionModel_dlog_gamma(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters, 0);
			double beta = Parameters_value(dm->parameters, 1);
			if(Parameters_count(dm->parameters) > 2){
				alpha = Parameters_value(dm->parameters, i*2);
				beta = Parameters_value(dm->parameters, i*2+1);
			}
			double x = Parameter_value(p) - dm->shift;
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
				beta = 1.0/beta;
			}
			return (alpha-1.0)/x - beta;
			
		}
	}
	return 0;
}

double DistributionModel_d2log_gamma(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters, 0);
			if(Parameters_count(dm->parameters) > 2){
				alpha = Parameters_value(dm->parameters, i*2);
			}
			double x = Parameter_value(p) - dm->shift;
			return -(alpha-1.0)/x/x;
		}
	}
	return 0;
}

static void DistributionModel_gamma_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			samples[i] = gsl_ran_gamma(dm->rng, alpha, beta);
			samples[i] += dm->shift;
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_gamma(dm->rng, alpha, beta);
			samples[i] += dm->shift;
		}
	}
}

static double DistributionModel_gamma_sample_evaluate(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			double sample = gsl_ran_gamma(dm->rng, alpha, beta);
			Parameters_set_value(dm->x, i, sample + dm->shift);
			logP += log(gsl_ran_gamma_pdf(sample, alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_gamma(dm->rng, alpha, beta);
			Parameters_set_value(dm->x, i, sample + dm->shift);
			logP += log(gsl_ran_gamma_pdf(sample, alpha, beta));
		}
	}
	return logP;
}

DistributionModel* new_GammaDistributionModel(const double shape, const double rate, const Parameters* x){
	Parameters* ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter("gamma.shape", shape, NULL));
	Parameters_move(ps, new_Parameter("gamma.rate", rate, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->type = DISTRIBUTION_GAMMA;
	dm->parameterization = DISTRIBUTION_GAMMA_SHAPE_RATE;
	dm->logP = DistributionModel_log_gamma;
	dm->logP_with_values = DistributionModel_log_gamma_with_values;
	dm->dlogP = DistributionModel_dlog_gamma;
	dm->d2logP = DistributionModel_d2log_gamma;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_gamma_sample;
	dm->sample_evaluate = DistributionModel_gamma_sample_evaluate;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_GammaDistributionModel_with_parameters(Parameters* parameters, const Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_GAMMA;
	dm->logP = DistributionModel_log_gamma;
	dm->logP_with_values = DistributionModel_log_gamma_with_values;
	dm->dlogP = DistributionModel_dlog_gamma;
	dm->d2logP = DistributionModel_d2log_gamma;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_gamma_sample;
	dm->sample_evaluate = DistributionModel_gamma_sample_evaluate;
	dm->parameterization = parameterization;
	return dm;
}

Model* new_GammaDistributionModel_from_json(json_node* node, Hashtable* hash){
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
	
	bool scale = false;
	// empirical
	if (samples != NULL) {
		size_t paramCount = Parameters_count(x);
		
		for (int i = 0; i < paramCount; i++) {
			double* vec = Vector_data(samples[i]);
			double m = dmean(vec, Vector_length(samples[i]));
			double v = variance(vec, Vector_length(samples[i]), m);
			Parameters_move(parameters, new_Parameter("shape", m*m/v, NULL));
			Parameters_move(parameters, new_Parameter("rate", m/v, NULL));
			free_Vector(samples[i]);
		}
		free(samples);
	}
	else if(get_json_node(node, "parameters") == NULL){
		get_parameters_references2(node, hash, x, "x");
		for (int i = 0; i < Parameters_count(x); i++) {
			Parameters_move(parameters, new_Parameter("shape", 1, new_Constraint(0, INFINITY)));
			Parameters_move(parameters, new_Parameter("rate", 1, new_Constraint(0, INFINITY)));
		}
	}
	else{
		get_parameters_references(node, hash, parameters);
		json_node* x_node = get_json_node(node, "parameters");
		for (int i = 0; i < x_node->child_count; i++) {
			if (strcasecmp(x_node->children[i]->key, "scale") == 0) {
				scale = true;
			}
			else if (strcasecmp(x_node->children[i]->key, "rate") != 0 && strcasecmp(x_node->children[i]->key, "shape") != 0) {
				fprintf(stderr, "Gamma distribution should be parametrized with rate and (shape or scale)\n");
				exit(13);
			}
		}
		if (strcasecmp(x_node->children[0]->key, "shape") != 0) {
			Parameters_swap_index(parameters, 0, 1);
		}
		for (int i = 0; i < Parameters_count(parameters); i++) {
			Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
		}
	}
	
	distribution_parameterization param = scale ? DISTRIBUTION_GAMMA_SHAPE_SCALE: DISTRIBUTION_GAMMA_SHAPE_RATE;
	DistributionModel* dm = new_GammaDistributionModel_with_parameters(parameters, x, param);
	
	Model* model = new_DistributionModel2(id, dm);
	
	model->samplable = true;
	dm->shift = get_json_node_value_double(node, "shift", 0);
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	free_Parameters(parameters);
	free_Parameters(x);
	
	return model;
}
