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
	Parameter* xParameter = Parameters_at(dm->x, 0);
	const double* x = Parameter_values(xParameter);
	size_t fieldDimension = Parameter_size(xParameter);
	double sum = 0;
	double precision = Parameter_value(Parameters_at(dm->parameters, 0));
	
	for(size_t i = 1; i < fieldDimension; i++){
		sum += pow(x[i-1] - x[i], 2.0);
	}
	
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	
	return dm->lp;
}


double DistributionModel_log_gmrf_time_aware(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	Coalescent* coal = dm->data;
	if(coal->need_update_intervals) coal->update_intervals(coal);

	Parameter* xParameter = Parameters_at(dm->x, 0);
	const double* x = Parameter_values(xParameter);
	size_t fieldDimension = Parameter_size(xParameter);
	double* intervals = dvector(fieldDimension);
	size_t j = 0;
	double interval = 0;
	for( size_t i = 0; i < coal->n; i++ ){
		interval += coal->times[i];
		if(coal->iscoalescent[i]){
			intervals[j++] = interval;
			interval = 0;
		}
	}
	double sum = 0;
	double precision = Parameter_value(Parameters_at(dm->parameters, 0));
	
	for(size_t i = 1; i < fieldDimension; i++){
		sum += pow(x[i-1] - x[i], 2.0)*2.0/(intervals[i]+intervals[i-1]);
	}
	free(intervals);
	dm->lp = log(precision)*(fieldDimension - 1)/2.0 - sum*precision/2.0 - (fieldDimension - 1)/2.0 * LOG_2PI;
	dm->need_update = false;
	
	return dm->lp;
}

double DistributionModel_log_gammarf(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	
	Parameter* xParameter = Parameters_at(dm->x, 0);
	const double* x = Parameter_values(xParameter);
	size_t fieldDimension = Parameter_size(xParameter);
	double precision = Parameter_value(Parameters_at(dm->parameters, 0));
	dm->lp = 0;
	for(size_t i = 1; i < fieldDimension; i++){
		dm->lp += log(gsl_ran_gamma_pdf(x[i], precision, x[i-1]/precision));
	}
	dm->need_update = false;
	
	return dm->lp;
}

double DistributionModel_gmrf_gradient(DistributionModel* dm, const Parameters* parameters){
	Parameter* xParameter = Parameters_at(dm->x, 0);
	const double* x = Parameter_values(xParameter);
	size_t fieldDimension = Parameter_size(xParameter);
	Parameter* precisionParameter = Parameters_at(dm->parameters, 0);
	double precision = Parameter_value(precisionParameter);
	
	// precision
	Parameter *precisionx = Parameters_depends(parameters, precisionParameter);
	if(precisionx != NULL){
		double sum = 0;
		for(size_t i = 1; i < fieldDimension; i++){
			sum += pow(x[i-1] - x[i], 2.0);
		}
		double dlogP = (fieldDimension - 1)/2.0/precision - sum/2.0;
		precisionParameter->grad[0] += dlogP;
		if(precisionParameter != precisionx){
			precisionParameter->transform->backward(precisionParameter->transform, &dlogP);
		}
	}
	
	// domain
	Parameter *xx = Parameters_depends(parameters, xParameter);
	if(xx != NULL){
		double* grad = dvector(fieldDimension);
		// f(x) = (x_i - x_{i-1})^2
		// dfdx_i = 2x_i -2x_{i-1}
		grad[0] = -(x[0] - x[1])*precision;
		xParameter->grad[0] += grad[0];
		size_t i = 1;
		for (; i < fieldDimension-1; i++) {
			// f(x) = (x_{i-1} - x_{i})^2 + (x_i - x_{i+1})^2
			// dfdx_i = 2x_i -2x_{i-1} + 2x_i -2x_{i+1}
			grad[i] = -(2.0*x[i] - x[i-1] - x[i+1])*precision;
			xParameter->grad[i] += grad[i];
		}
		grad[i] += -(x[i] - x[i-1])*precision;
		xParameter->grad[i] += grad[i];
		if(xParameter != xx){
			xParameter->transform->backward(xParameter->transform, grad);
		}
		free(grad);
	}
	return 0;
}

double DistributionModel_gmrf_time_aware_gradient(DistributionModel* dm, const Parameters* parameters){
	Coalescent* coal = dm->data;
	if(coal->need_update_intervals) coal->update_intervals(coal);
	
	Parameter* xParameter = Parameters_at(dm->x, 0);
	const double* x = Parameter_values(xParameter);
	size_t fieldDimension = Parameter_size(xParameter);
	Parameter* precisionParameter = Parameters_at(dm->parameters, 0);
	double precision = Parameter_value(precisionParameter);

	double* intervals = dvector(fieldDimension);
	size_t j = 0;
	double interval = 0;
	for( size_t i = 0; i< coal->n; i++  ){
		interval += coal->times[i];
		if(coal->iscoalescent[i]){
			intervals[j++] = interval;
			interval = 0;
		}
	}

	Parameters* treeModelParameters = new_Parameters(1);
	Parameters* reparam = get_reparams(coal->tree);
	Node** nodes = Tree_nodes(coal->tree);
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		Parameter* parameter = Parameters_at(parameters, i);
		// ratios and root_height transformed
		if(parameter->model == MODEL_TREE_TRANSFORM){
			for(size_t j = 0; j < Parameters_count(reparam); j++){
				Parameter* xx = Parameters_depends(parameters, Parameters_at(reparam, j));
				if(xx != NULL) {
					Parameters_add(treeModelParameters, xx);
				}
			}
		}
		// heights
		else if(parameter->model == MODEL_TREE){
			Parameter* xx = Parameters_depends(parameters, nodes[parameter->id]->height);
			if(xx != NULL){
				Parameters_add(treeModelParameters, xx);
			}
		}
	}
	if(Parameters_count(treeModelParameters) > 0){
		size_t tipCount = Tree_tip_count(coal->tree);
		size_t nodeCount = Tree_node_count(coal->tree);
		double* heightGradient = dvector(nodeCount);
		j = 0;
		size_t previousInternalIndex = 0;
		size_t previousInternalIndex2 = 0;
		for( size_t i = 0; i< coal->n; i++  ){
			// f(x) = (x_i - x_{i-1})^2*2/(intervals[i]+intervals[i-1]) = (x_i - x_{i-1})^2*2/(time_i - time_{i-2})
			if(coal->iscoalescent[i]){
				if(j >= 1){
					double temp = -pow(x[j] - x[j-1], 2)*2/pow(intervals[j]+intervals[j-1], 2)*precision/2.0;
					heightGradient[coal->nodes[i]->index] += -temp;
					if( j > 1){
						heightGradient[previousInternalIndex2] += temp;
					}
				}
				previousInternalIndex2 = previousInternalIndex;
				previousInternalIndex = coal->nodes[i]->index;
			}
		}
		Tree_height_backward(coal->tree, treeModelParameters, heightGradient);
		free(heightGradient);
	}
	free_Parameters(treeModelParameters);

	// precision
	Parameter *precisionx = Parameters_depends(parameters, precisionParameter);
	if(precisionx != NULL){
		double sum = 0;
		for(size_t i = 1; i < fieldDimension; i++){
			sum += pow(x[i-1] - x[i], 2.0)*2.0/(intervals[i]+intervals[i-1]);
		}
		double dlogP = (fieldDimension - 1)/2.0/precision - sum/2.0;
		precisionParameter->grad[0] += dlogP;
		if(precisionParameter != precisionx){
			precisionParameter->transform->backward(precisionParameter->transform, &dlogP);
		}
	}
	
	// domain
	Parameter *xx = Parameters_depends(parameters, xParameter);
	if(xx != NULL){
		double* grad = dvector(fieldDimension);
		grad[0] = -(x[0] - x[1])*precision*2.0/(intervals[0]+intervals[1]);
		xParameter->grad[0] += grad[0];
		size_t i = 1;
		for (; i < fieldDimension-1; i++) {
			grad[i] = -((x[i] - x[i-1])*2.0/(intervals[i]+intervals[i-1]) + 
						(x[i] - x[i+1])*2.0/(intervals[i]+intervals[i+1]))*precision;
			xParameter->grad[i] += grad[i];
		}
		grad[i] += -(x[i] - x[i-1])*precision*2.0/(intervals[i]+intervals[i-1]);;
		xParameter->grad[i] += grad[i];
		if(xParameter != xx){
			xParameter->transform->backward(xParameter->transform, grad);
		}
		free(grad);
	}
	free(intervals);
	return 0;
}

static void DistributionModel_gmrf_sample(DistributionModel* dm){
	fprintf(stderr, "DistributionModel_gmrf_sample not implemented\n");
	exit(2);
}

DistributionModel* new_GMRF_with_parameters(Parameters* parameters, Parameter* x, Coalescent* coalescent, distribution_parameterization parameterization){
	Parameters* xx = new_Parameters(1);
	Parameters_add(xx, x);
	DistributionModel* dm = new_DistributionModel(parameters, xx);
	free_Parameters(xx);
	dm->type = DISTRIBUTION_GMRF;
	dm->parameterization = parameterization;
	dm->logP = DistributionModel_log_gmrf;
	dm->gradient2 = DistributionModel_gmrf_gradient;
	dm->sample = DistributionModel_gmrf_sample;
	dm->data = coalescent;
    dm->shift = INFINITY;
	return dm;
}

Model* new_GMRFModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	char* model_key = get_json_node_value_string(node, "tree");
    json_node* parameters_node = get_json_node(node, "parameters");
    
    json_node* x_node = get_json_node(node, "x");
    Parameter* x = new_Parameter_from_json(x_node, hash);
	
    if (parameters_node->child_count != 1 && strcasecmp(parameters_node->children[0]->key, "precision") != 0) {
        fprintf(stderr, "GMRF should be parametrized with precision parameters\n");
        exit(13);
    }
    
	json_node* precisionNode = get_json_node(parameters_node, "precision");
    Parameter* precision = new_Parameter_from_json(precisionNode, hash);
	Parameters* parameters = new_Parameters(1);
	Parameters_add(parameters, precision);
	Hashtable_add(hash, Parameter_name(precision), precision);
	
	DistributionModel* dm = NULL;
	Model* model = NULL;
	
	if (model_key != NULL) {
		Model *m = Hashtable_get(hash, model_key+1);
		dm = new_GMRF_with_parameters(parameters, x, m->obj, 0);
		model = new_DistributionModel3(id, dm, m);
		dm->logP = DistributionModel_log_gmrf_time_aware;
		dm->gradient2 = DistributionModel_gmrf_time_aware_gradient;
		m->listeners->add(m->listeners, model);
	}
	else{
		dm = new_GMRF_with_parameters(parameters, x, NULL, 0);
		model = new_DistributionModel2(id, dm);
	}
	
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
    free_Parameters(parameters);
	
	return model;
}
