//
//  distkumaraswamy.c
//  physher
//
//  Created by Mathieu Fourment on 19/8/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distkumaraswamy.h"

#include <strings.h>

#include <gsl/gsl_randist.h>

#include "parametersio.h"

double DistributionModel_log_kumaraswamy(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	dm->lp = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters, i*2);
			double b = Parameters_value(dm->parameters, i*2+1);
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(a*b*pow(x, a - 1.0) * pow( 1.0 - pow(x, a), b - 1.0));
		}
	}
	else{
		double a = Parameters_value(dm->parameters, 0);
		double b = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(a*b*pow(x, a - 1.0) * pow( 1.0 - pow(x, a), b - 1.0));
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_kumaraswamy_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters, i*2);
			double b = Parameters_value(dm->parameters, i*2+1);
			logP += log(a*b*pow(values[i], a - 1.0) * pow( 1.0 - pow(values[i], a), b - 1.0));
		}
	}
	else{
		double a = Parameters_value(dm->parameters, 0);
		double b = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(a*b*pow(values[i], a - 1.0) * pow( 1.0 - pow(values[i], a), b - 1.0));
		}
	}
	return logP;
}


double DistributionModel_dlog_kumaraswamy(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double a;
			double b;
			if(Parameters_count(dm->parameters) > 1){
				a = Parameters_value(dm->parameters, i*2);
				b = Parameters_value(dm->parameters, i*2+1);
			}
			else{
				a = Parameters_value(dm->parameters, 0);
				b = Parameters_value(dm->parameters, 1);
			}
			double x = Parameters_value(dm->x, i);
			return -(a - 1.0)/x + (b - 1.0)*a*pow(x, a - 1.0)/(pow(x, a) - 1.0);
		}
	}
	return 0;
}

// F_X(x) = p(X <= x)
double DistributionModel_kumaraswamy_inverse_CDF(double p, double a, double b){
	return pow(1.0 - pow(1.0 - p, 1.0/b), 1.0/a);
}

static void DistributionModel_kumaraswamy_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters, i*2);
			double b = Parameters_value(dm->parameters, i*2+1);
			samples[i] = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
		}
	}
	else{
		double a = Parameters_value(dm->parameters, 0);
		double b = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
		}
	}
}

static double DistributionModel_kumaraswamy_sample_evaluate(DistributionModel* dm){
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters, i*2);
			double b = Parameters_value(dm->parameters, i*2+1);
			double sample = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		double a = Parameters_value(dm->parameters, 0);
		double b = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_log_kumaraswamy(dm);
}

DistributionModel* new_KumaraswamyDistributionModel(const double lambda, const Parameters* x){
	Parameters* ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter("kumaraswamy.a", lambda, NULL));
	Parameters_move(ps, new_Parameter("kumaraswamy.b", lambda, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->type = DISTRIBUTION_KUMARASWAMY;
	dm->parameterization = 0;
	dm->logP = DistributionModel_log_kumaraswamy;
	dm->logP_with_values = DistributionModel_log_kumaraswamy_with_values;
	dm->dlogP = DistributionModel_dlog_kumaraswamy;
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_kumaraswamy_sample;
	dm->sample_evaluate = DistributionModel_kumaraswamy_sample_evaluate;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_KumaraswamyDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_KUMARASWAMY;
	dm->parameterization = 0;
	
	dm->logP = DistributionModel_log_kumaraswamy;
	dm->logP_with_values = DistributionModel_log_kumaraswamy_with_values;
	dm->dlogP = DistributionModel_dlog_kumaraswamy;
	dm->sample = DistributionModel_kumaraswamy_sample;
	dm->sample_evaluate = DistributionModel_kumaraswamy_sample_evaluate;
	
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
	return dm;
}

Model* new_KumaraswamyDistributionModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	
	Parameters* x = new_Parameters(1);
	get_parameters_references2(node, hash, x, "x");
	
	Parameters* parameters = new_Parameters(1);
	get_parameters_references(node, hash, parameters);
	
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
	}
	
	DistributionModel* dm = new_KumaraswamyDistributionModel_with_parameters(parameters, x);
	
	dm->shift = get_json_node_value_double(node, "shift", 0);
	Model* model = new_DistributionModel2(id, dm);
	
	model->samplable = true;
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	free_Parameters(parameters);
	free_Parameters(x);
	
	return model;
}
