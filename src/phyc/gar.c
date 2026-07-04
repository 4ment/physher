// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "gar.h"

#include <math.h>
#include <strings.h>

#include "distmodel.h"
#include "mathconstant.h"
#include "matrix.h"

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>


double DistributionModel_log_gar(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	Parameter* xParameter = Parameters_at(dm->x, 0);
	const double* x = Parameter_values(xParameter);
	size_t fieldDimension = Parameter_size(xParameter);
	double sum = 0;
	double shape = Parameter_value(Parameters_at(dm->parameters, 0));
	dm->lp = 0;
	for(size_t i = 1; i < fieldDimension; i++){
		dm->lp += log(gsl_ran_gamma_pdf(x[i], shape, x[i-1]/shape));
	}
	dm->need_update = false;
	
	return dm->lp;
}

void DistributionModel_gar_gradient(DistributionModel* dm, Parameters* parameters){
	Parameter* xParameter = Parameters_at(dm->x, 0);
	const double* x = Parameter_values(xParameter);
	size_t fieldDimension = Parameter_size(xParameter);
	Parameter* shapeParameter = Parameters_at(dm->parameters, 0);
	double shape = Parameter_value(shapeParameter);

	// shape
	Parameter *shapex = Parameters_depends(parameters, shapeParameter);
	if(shapex != NULL){
		double psiShape = gsl_sf_psi(shape);
		double dlogP = 0;
		for (size_t i = 1; i < fieldDimension; i++) {
			double xi   = x[i];
			double xim1 = x[i-1];

			// d log p_i / d shape
			dlogP += log((shape * xi) / xim1) + 1.0 - psiShape - xi / xim1;
	    }
		shapeParameter->grad[0] += dlogP;
		if(shapeParameter != shapex){
			shapeParameter->transform->backward(shapeParameter->transform, &dlogP);
		}
	}
	
	// domain
	Parameter *xx = Parameters_depends(parameters, xParameter);
	if(xx != NULL){
		double* grad = dvector(fieldDimension);
		memset(grad, 0, fieldDimension * sizeof(double));
		for (size_t i = 1; i < fieldDimension; i++) {
			double xi   = x[i];
			double xim1 = x[i-1];
			// d log p_i / d x_i
			grad[i] += (shape - 1.0) / xi - shape / xim1;
			// d log p_i / d x_{i-1}
			grad[i-1] += shape * (xi - xim1) / (xim1 * xim1);
	    }
		
		for (size_t i = 0; i < fieldDimension; i++) {
			xParameter->grad[i] += grad[i];
		}
		if(xParameter != xx){
			xParameter->transform->backward(xParameter->transform, grad);
		}
		free(grad);
	}
}


DistributionModel* new_GAR_with_parameters(Parameters* parameters, Parameter* x, distribution_parameterization parameterization){
	Parameters* xx = new_Parameters(1);
	Parameters_add(xx, x);
	DistributionModel* dm = new_DistributionModel(parameters, xx);
	free_Parameters(xx);
	dm->type = DISTRIBUTION_GAR;
	dm->parameterization = parameterization;
	dm->logP = DistributionModel_log_gar;
	dm->gradient = DistributionModel_gar_gradient;
	dm->sample = NULL;
    dm->shift = INFINITY;
	return dm;
}

Model* new_GARModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
    json_node* parameters_node = get_json_node(node, "parameters");
    
    json_node* x_node = get_json_node(node, "x");
    Parameter* x = new_Parameter_from_json(x_node, hash);
	
    if (parameters_node != NULL && parameters_node->child_count == 1 && strcasecmp(parameters_node->children[0]->key, "shape") != 0) {
        fprintf(stderr, "GAR should be parametrized with shape parameters\n");
        exit(13);
    }
	Parameter* shape = NULL;
	if(parameters_node == NULL){
		shape = new_Parameter("shape", 1.0, new_Constraint(0, INFINITY));
	}
	else{
		json_node* shapeNode = get_json_node(parameters_node, "shape");
		if(shapeNode == NULL){
			shape = new_Parameter("shape", 1.0, new_Constraint(0, INFINITY));
		}
		else{
			shape = new_Parameter_from_json(shapeNode, hash);
			Hashtable_add(hash, Parameter_name(shape), shape);
		}
	}
	Parameters* parameters = new_Parameters(1);
	Parameters_add(parameters, shape);
	
	DistributionModel* dm = new_GAR_with_parameters(parameters, x, DISTRIBUTION_GAMMA_SHAPE_RATE);
	Model* model = new_DistributionModel2(id, dm);
	
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
    free_Parameters(parameters);
	
	return model;
}
