//
//  gmrf.c
//  physher
//
//  Created by Mathieu Fourment on 18/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "gmrf.h"

#include <tgmath.h>

#include "distmodel.h"
#include "mathconstant.h"

double DistributionModel_log_gmrf(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double sum = 0;
	double precision = Parameters_value(dm->parameters, 0);
	
	for(int i = 1; i < fieldDimension; i++){
		double popSize = Parameters_value(x, i);
		double popSize1 = Parameters_value(x, i - 1);
		sum += pow(popSize1 - popSize, 2.0);
	}
	
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	
	return dm->lp;
}


double DistributionModel_log_gmrf_with_values(DistributionModel* dm, const double* values){
	size_t fieldDimension = Parameters_count(dm->x);
	double sum = 0;
	double precision = Parameters_value(dm->parameters, 0);
	
	for(size_t i = 1; i < fieldDimension; i++){
		sum += pow(values[i - 1] - values[i], 2.0);
	}
	
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_dlog_gmrf(DistributionModel* dm, const Parameter* p){
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double precision = Parameters_value(dm->parameters, 0);
	
	// precision
	if(p == Parameters_at(dm->parameters, 0)){
		double sum = 0;
		for(int i = 1; i < fieldDimension; i++){
			double popSize = Parameters_value(x, i);
			double popSize1 = Parameters_value(x, i - 1);
			sum += pow(popSize1 - popSize, 2.0);
		}
		return (fieldDimension - 1)/2.0/precision - sum/2.0;
	}
	
	// domain
	for (int i = 0; i < fieldDimension; i++) {
		if(p == Parameters_at(x, i)){
			// f(x) = (x_i - x_{i-1})^2
			// dfdx_i = 2x_i -2x_{i-1}
			if (i == 0) {
				return -(Parameters_value(x, i) - Parameters_value(x, i+1))*precision;
			}
			else if(i == fieldDimension-1){
				return -(Parameters_value(x, i) - Parameters_value(x, i-1))*precision;
			}
			
			// f(x) = (x_{i-1} - x_{i})^2 + (x_i - x_{i+1})^2
			// dfdx_i = 2x_i -2x_{i-1} + 2x_i -2x_{i+1}
			return -(2.0*Parameters_value(x, i) - Parameters_value(x, i-1) - Parameters_value(x, i+1))*precision;
		}
	}
	return 0;
}


double DistributionModel_d2log_gmrf(DistributionModel* dm, const Parameter* p){
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double precision = Parameters_value(dm->parameters, 0);
	
	// precision
	if(p == Parameters_at(dm->parameters, 0)){
		return -(fieldDimension - 1)/2.0/(precision*precision);
	}
	
	// domain
	for (int i = 0; i < fieldDimension; i++) {
		if(p == Parameters_at(x, i)){
			if (i == 0 || i == fieldDimension-1){
				return -precision;
			}
			return -2.0*precision;
		}
	}
	return 0;
}

double DistributionModel_ddlog_gmrf(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	
	// precision
	if(p1 == Parameters_at(dm->parameters, 0) && p2 == Parameters_at(dm->parameters, 0)){
		return DistributionModel_d2log_gmrf(dm, p1);
	}
	// precision and one of the x
	if(p1 == Parameters_at(dm->parameters, 0) || p2 == Parameters_at(dm->parameters, 0)){
		for (int i = 0; i < fieldDimension; i++) {
			if((p1 == Parameters_at(x, i) || p2 == Parameters_at(x, i)) ){
				if (i == 0) {
					return -(Parameters_value(x, i) - Parameters_value(x, i+1));
				}
				else if(i == fieldDimension-1){
					return -(Parameters_value(x, i) - Parameters_value(x, i-1));
				}
				return -(2.0*Parameters_value(x, i) - Parameters_value(x, i-1) - Parameters_value(x, i+1));
			}
		}
	}
	
	int first = -1;
	int second = -1;
	
	for (int i = 0; i < fieldDimension; i++) {
		if(p1 == Parameters_at(dm->x, i)) first = i;
		else if(p2 == Parameters_at(dm->x, i)) second = i;
	}
	
	if(first != -1){
		if (first == second) {
			return DistributionModel_d2log_gmrf(dm, p1);
		}
		// consecutive
		else if(abs(first - second) == 1){
			return Parameters_value(dm->parameters, 0);
		}
	}
	
	return 0;
}

static void DistributionModel_gmrf_sample(DistributionModel* dm, double* samples){
	fprintf(stderr, "DistributionModel_gmrf_sample not implemented\n");
	exit(2);
}

static double DistributionModel_gmrf_sample_evaluate(DistributionModel* dm){
	fprintf(stderr, "DistributionModel_gmrf_sample_evaluate not implemented\n");
	exit(2);
	return DistributionModel_log_gmrf(dm);
}

DistributionModel* new_GMRF_with_parameters(Parameters* parameters, const Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_GMRF;
	dm->parameterization = parameterization;
	dm->logP = DistributionModel_log_gmrf;
	dm->logP_with_values = DistributionModel_log_gmrf_with_values;
	dm->dlogP = DistributionModel_dlog_gmrf;
	dm->sample = DistributionModel_gmrf_sample;
	dm->sample_evaluate = DistributionModel_gmrf_sample_evaluate;
	dm->d2logP = DistributionModel_d2log_gmrf;
	dm->ddlogP = DistributionModel_ddlog_gmrf;
	return dm;
}

Model* new_GMRFModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	
	Parameters* x = new_Parameters(1);
	get_parameters_references2(node, hash, x, "x");
	
	Parameters* parameters = new_Parameters(1);
	get_parameters_references(node, hash, parameters);
	
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
	}
	
	DistributionModel* dm = new_GMRF_with_parameters(parameters, x, 0);
	
	Model* model = new_DistributionModel2(id, dm);
	
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	free_Parameters(parameters);
	free_Parameters(x);
	
	return model;
}