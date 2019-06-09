//
//  distcauchy.c
//  physher
//
//  Created by Mathieu Fourment on 8/6/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distcauchy.h"

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"
#include "mathconstant.h"

double DistributionModel_log_cauchy(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	dm->lp = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(gsl_ran_cauchy_pdf(x, alpha));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(gsl_ran_cauchy_pdf(x, alpha));
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_cauchy_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			logP += log(gsl_ran_cauchy_pdf(values[i] - dm->shift, alpha));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(gsl_ran_cauchy_pdf(values[i] - dm->shift, alpha));
		}
	}
	return logP;
}

double DistributionModel_dlog_cauchy(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters, 0);
			double x = Parameter_value(p) - dm->shift;
			return -2.0*x/(alpha*alpha*(x*x/alpha/alpha + 1.0));
		}
	}
	return 0;
}

double DistributionModel_d2log_cauchy(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters, 0);
			if(Parameters_count(dm->parameters) > 2){
				alpha = Parameters_value(dm->parameters, i*2);
			}
			double x = Parameter_value(p) - dm->shift;
			double alpha2 = alpha*alpha;
			return -2.0*(alpha2 - x*x)/(alpha2*alpha2 + 2.0*alpha2*x*x + x*x*x*x);
		}
	}
	return 0;
}

static void DistributionModel_cauchy_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			samples[i] = gsl_ran_cauchy(dm->rng, alpha);
			samples[i] += dm->shift;
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_cauchy(dm->rng, alpha);
			samples[i] += dm->shift;
		}
	}
}

static double DistributionModel_cauchy_sample_evaluate(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double sample = gsl_ran_cauchy(dm->rng, alpha);
			Parameters_set_value(dm->x, i, sample + dm->shift);
			logP += log(gsl_ran_cauchy_pdf(sample, alpha));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_cauchy(dm->rng, alpha);
			Parameters_set_value(dm->x, i, sample + dm->shift);
			logP += log(gsl_ran_cauchy_pdf(sample, alpha));
		}
	}
	return logP;
}

DistributionModel* new_CauchyDistributionModel(const double location, const double scale, const Parameters* x){
	Parameters* ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter("cauchy.scale", scale, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->type = DISTRIBUTION_CAUCHY;
	dm->logP = DistributionModel_log_cauchy;
	dm->logP_with_values = DistributionModel_log_cauchy_with_values;
	dm->dlogP = DistributionModel_dlog_cauchy;
	dm->d2logP = DistributionModel_d2log_cauchy;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_cauchy_sample;
	dm->sample_evaluate = DistributionModel_cauchy_sample_evaluate;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_CauchyDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_CAUCHY;
	dm->logP = DistributionModel_log_cauchy;
	dm->logP_with_values = DistributionModel_log_cauchy_with_values;
	dm->dlogP = DistributionModel_dlog_cauchy;
	dm->d2logP = DistributionModel_d2log_cauchy;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_cauchy_sample;
	dm->sample_evaluate = DistributionModel_cauchy_sample_evaluate;
	return dm;
}

Model* new_CauchyDistributionModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	
	Parameters* x = new_Parameters(1);
	get_parameters_references2(node, hash, x, "x");
	
	char* file = get_json_node_value_string(node, "file");
	Parameters* parameters = new_Parameters(1);
	
	if(get_json_node(node, "parameters") == NULL){
		for (int i = 0; i < Parameters_count(x); i++) {
//			Parameters_move(parameters, new_Parameter("location", 1, new_Constraint(0, INFINITY)));
			Parameters_move(parameters, new_Parameter("scale", 1, new_Constraint(0, INFINITY)));
		}
	}
	else{
		get_parameters_references(node, hash, parameters);
		for (int i = 0; i < Parameters_count(parameters); i++) {
			Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
		}
	}
	
	DistributionModel* dm = new_CauchyDistributionModel_with_parameters(parameters, x);
	
	Model* model = new_DistributionModel2(id, dm);
	
	model->samplable = true;
	dm->shift = get_json_node_value_double(node, "shift", 0);
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	free_Parameters(parameters);
	free_Parameters(x);
	
	return model;
}
