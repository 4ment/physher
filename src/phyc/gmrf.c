//
//  gmrf.c
//  physher
//
//  Created by Mathieu Fourment on 18/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "gmrf.h"

#include <math.h>
#include <strings.h>

#include "distmodel.h"
#include "mathconstant.h"
#include "matrix.h"
#include "demographicmodels.h"

#include <gsl/gsl_randist.h>

double DistributionModel_log_gmrf(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double sum = 0;
	double precision = Parameters_value(dm->parameters[0], 0);
	
	for(int i = 1; i < fieldDimension; i++){
		double popSize = Parameters_value(x, i);
		double popSize1 = Parameters_value(x, i - 1);
		sum += pow(popSize1 - popSize, 2.0);
	}
	
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	
	return dm->lp;
}


double DistributionModel_log_gmrf_time_aware(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	Coalescent* coal = dm->data;
	if(coal->need_update_intervals) coal->update_intervals(coal);
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double* intervals = dvector(fieldDimension);
	int j = 0;
	double interval = 0;
	for( int i = 0; i< coal->n; i++  ){
		interval += coal->times[i];
		if(coal->iscoalescent[i]){
			intervals[j++] = interval;
			interval = 0;
		}
	}
	double sum = 0;
	double precision = Parameters_value(dm->parameters[0], 0);
	
	for(int i = 1; i < fieldDimension; i++){
		double popSize = Parameters_value(x, i);
		double popSize1 = Parameters_value(x, i - 1);
		sum += pow(popSize1 - popSize, 2.0)*2.0/(intervals[i]+intervals[i-1]);
	}
	free(intervals);
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	
	return dm->lp;
}

double DistributionModel_log_gammarf(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double precision = Parameters_value(dm->parameters[0], 0);
	dm->lp = 0;
	for(int i = 1; i < fieldDimension; i++){
		dm->lp += log(gsl_ran_gamma_pdf(Parameters_value(x, i), precision, Parameters_value(x, i-1)/precision));
	}
	dm->need_update = false;
	
	return dm->lp;
}


double DistributionModel_log_gmrf_with_values(DistributionModel* dm, const double* values){
	size_t fieldDimension = Parameters_count(dm->x);
	double sum = 0;
	double precision = Parameters_value(dm->parameters[0], 0);
	
	for(size_t i = 1; i < fieldDimension; i++){
		sum += pow(values[i - 1] - values[i], 2.0);
	}
	
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	return dm->lp;
}


double DistributionModel_log_gmrf_with_values_time_aware(DistributionModel* dm, const double* values){
	Coalescent* coal = dm->data;
	if(coal->need_update_intervals) coal->update_intervals(coal);
	size_t fieldDimension = Parameters_count(dm->x);
	double sum = 0;
	double precision = Parameters_value(dm->parameters[0], 0);
	double* intervals = dvector(fieldDimension);
	int j = 0;
	double interval = 0;
	for( int i = 0; i< coal->n; i++  ){
		interval += coal->times[i];
		if(coal->iscoalescent[i]){
			intervals[j++] = interval;
			interval = 0;
		}
	}
	
	for(size_t i = 1; i < fieldDimension; i++){
		sum += pow(values[i - 1] - values[i], 2.0)*2.0/(intervals[i]+intervals[i-1]);
	}
	
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	free(intervals);
	return dm->lp;
}

double DistributionModel_dlog_gmrf(DistributionModel* dm, const Parameter* p){
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double precision = Parameters_value(dm->parameters[0], 0);
	
	// precision
	if(p == Parameters_at(dm->parameters[0], 0)){
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

double DistributionModel_dlog_gmrf_time_aware(DistributionModel* dm, const Parameter* p){
	Coalescent* coal = dm->data;
	if(coal->need_update_intervals) coal->update_intervals(coal);
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double* intervals = dvector(fieldDimension);
	int j = 0;
	double interval = 0;
	for( int i = 0; i< coal->n; i++  ){
		interval += coal->times[i];
		if(coal->iscoalescent[i]){
			intervals[j++] = interval;
			interval = 0;
		}
	}
	double precision = Parameters_value(dm->parameters[0], 0);
	double dlogP = 0;
	
	// precision
	if(p == Parameters_at(dm->parameters[0], 0)){
		double sum = 0;
		for(int i = 1; i < fieldDimension; i++){
			double popSize = Parameters_value(x, i);
			double popSize1 = Parameters_value(x, i - 1);
			sum += pow(popSize1 - popSize, 2.0)*2.0/(intervals[i]+intervals[i-1]);
		}
		dlogP = (fieldDimension - 1)/2.0/precision - sum/2.0;
	}
	else{
		// domain
		for (int i = 0; i < fieldDimension; i++) {
			if(p == Parameters_at(x, i)){
				// f(x) = (x_i - x_{i-1})^2
				// dfdx_i = 2x_i -2x_{i-1}
				if (i == 0) {
					dlogP = -(Parameters_value(x, i) - Parameters_value(x, i+1))*precision*2.0/(intervals[i]+intervals[i-1]);
				}
				else if(i == fieldDimension-1){
					dlogP = -(Parameters_value(x, i) - Parameters_value(x, i-1))*precision*2.0/(intervals[i]+intervals[i-1]);
				}
				else{
					// f(x) = (x_{i-1} - x_{i})^2 + (x_i - x_{i+1})^2
					// dfdx_i = 2x_i -2x_{i-1} + 2x_i -2x_{i+1}
					dlogP = -((Parameters_value(x, i) - Parameters_value(x, i-1))*2.0/(intervals[i]+intervals[i-1]) +
							  (Parameters_value(x, i) - Parameters_value(x, i+1))*2.0/(intervals[i]+intervals[i+1]))*precision;
				}
				break;
			}
		}
	}
	free(intervals);
	return dlogP;
}

double DistributionModel_d2log_gmrf(DistributionModel* dm, const Parameter* p){
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double precision = Parameters_value(dm->parameters[0], 0);
	
	// precision
	if(p == Parameters_at(dm->parameters[0], 0)){
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

double DistributionModel_d2log_gmrf_time_aware(DistributionModel* dm, const Parameter* p){
	Coalescent* coal = dm->data;
	if(coal->need_update_intervals) coal->update_intervals(coal);
	
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	double precision = Parameters_value(dm->parameters[0], 0);
	double d2logP = 0;
	double* intervals = dvector(fieldDimension);
	int j = 0;
	double interval = 0;
	for( int i = 0; i< coal->n; i++  ){
		interval += coal->times[i];
		if(coal->iscoalescent[i]){
			intervals[j++] = interval;
			interval = 0;
		}
	}
	
	// precision
	if(p == Parameters_at(dm->parameters[0], 0)){
		d2logP = -(fieldDimension - 1)/2.0/(precision*precision);
	}
	else{
		// domain
		for (int i = 0; i < fieldDimension; i++) {
			if(p == Parameters_at(x, i)){
				if (i == 0){
					d2logP =  -precision*2.0/(intervals[i]+intervals[i-1]);
				}
				else if(i == fieldDimension-1){
					d2logP =  -precision*2.0/(intervals[i]+intervals[i+1]);
				}
				else{
					d2logP = -(2.0/(intervals[i]+intervals[i-1]) + 2.0/(intervals[i]+intervals[i+1]))*precision;
				}
				break;
			}
		}
	}
	free(intervals);
	return d2logP;
}

double DistributionModel_ddlog_gmrf(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	
	// precision
	if(p1 == Parameters_at(dm->parameters[0], 0) && p2 == p1){
		return DistributionModel_d2log_gmrf(dm, p1);
	}
	// precision and one of the x
	if(p1 == Parameters_at(dm->parameters[0], 0) || p2 == p1){
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
			return Parameters_value(dm->parameters[0], 0);
		}
	}
	
	return 0;
}

double DistributionModel_ddlog_gmrf_time_aware(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	Coalescent* coal = dm->data;
	if(coal->need_update_intervals) coal->update_intervals(coal);
	
	Parameters* x = dm->x;
	size_t fieldDimension = Parameters_count(x);
	
	double* intervals = dvector(fieldDimension);
	int j = 0;
	double interval = 0;
	for( int i = 0; i< coal->n; i++  ){
		interval += coal->times[i];
		if(coal->iscoalescent[i]){
			intervals[j++] = interval;
			interval = 0;
		}
	}
	double ddlogP = 0;
	
	// precision
	if(p1 == Parameters_at(dm->parameters[0], 0) && p2 == p1){
		ddlogP = DistributionModel_d2log_gmrf(dm, p1);
	}
	// precision and one of the x
	else if(p1 == Parameters_at(dm->parameters[0], 0) || p2 == p1){
		for (int i = 0; i < fieldDimension; i++) {
			if((p1 == Parameters_at(x, i) || p2 == Parameters_at(x, i)) ){
				if (i == 0) {
					ddlogP = -(Parameters_value(x, i) - Parameters_value(x, i+1))*2.0/(intervals[i]+intervals[i+1]);
				}
				else if(i == fieldDimension-1){
					ddlogP = -(Parameters_value(x, i) - Parameters_value(x, i-1))*2.0/(intervals[i]+intervals[i-1]);
				}
				else{
					ddlogP = -((Parameters_value(x, i) - Parameters_value(x, i-1))*2.0/(intervals[i]+intervals[i-1]) +
							   (Parameters_value(x, i) - Parameters_value(x, i+1))*2.0/(intervals[i]+intervals[i+1]));
				}
				break;
			}
		}
	}
	else{
		int first = -1;
		int second = -1;
		
		for (int i = 0; i < fieldDimension; i++) {
			if(p1 == Parameters_at(dm->x, i)) first = i;
			else if(p2 == Parameters_at(dm->x, i)) second = i;
		}
		
		if(first != -1){
			if (first == second) {
				ddlogP = DistributionModel_d2log_gmrf(dm, p1);
			}
			// consecutive
			else if(abs(first - second) == 1){
				ddlogP = Parameters_value(dm->parameters[0], 0);
			}
		}
	}
	free(intervals);
	return ddlogP;
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

DistributionModel* new_GMRF_with_parameters(Parameters** parameters, Parameters* x, Coalescent* coalescent, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, 1, x);
	dm->type = DISTRIBUTION_GMRF;
	dm->parameterization = parameterization;
	dm->logP = DistributionModel_log_gmrf;
	dm->logP_with_values = DistributionModel_log_gmrf_with_values;
	dm->dlogP = DistributionModel_dlog_gmrf;
	dm->sample = DistributionModel_gmrf_sample;
	dm->sample_evaluate = DistributionModel_gmrf_sample_evaluate;
	dm->d2logP = DistributionModel_d2log_gmrf;
	dm->ddlogP = DistributionModel_ddlog_gmrf;
	dm->data = coalescent;
    dm->shift = INFINITY;
	return dm;
}

Model* new_GMRFModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	char* model_key = get_json_node_value_string(node, "tree");
    json_node* parameters_node = get_json_node(node, "parameters");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
	
    if (parameters_node->child_count != 1 && strcasecmp(parameters_node->children[0]->key, "precision") != 0) {
        fprintf(stderr, "GMRF should be parametrized with precision parameters\n");
        exit(13);
    }
    
    size_t parameters_dim = 1;
    Parameters** parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
	
	DistributionModel* dm = NULL;
	Model* model = NULL;
	
	if (model_key != NULL) {
		Model *m = Hashtable_get(hash, model_key+1);
		dm = new_GMRF_with_parameters(parameters, x, m->obj, 0);
		model = new_DistributionModel3(id, dm, m);
		dm->logP = DistributionModel_log_gmrf_time_aware;
		dm->dlogP = DistributionModel_dlog_gmrf_time_aware;
		dm->d2logP = DistributionModel_d2log_gmrf_time_aware;
		dm->ddlogP = NULL; //TODO
		m->listeners->add(m->listeners, model);
	}
	else{
		dm = new_GMRF_with_parameters(parameters, x, NULL, 0);
		model = new_DistributionModel2(id, dm);
	}
	
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	free_Parameters(parameters[0]);
    free(parameters);
	
	return model;
}
