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

#ifndef GSL_DISABLED
#include <gsl/gsl_randist.h>
#endif

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"


// parameter is lambda or mean
// out is always lambda
static void _expand_parameter_to_lambda(distribution_parameterization parameterization, Parameters* x, Parameter* parameter, double* out) {
    if (Parameter_size(parameter) == 1) {
        size_t dim = 0;
        for (size_t k = 0; k < Parameters_count(x); k++) {
            dim += Parameter_size(Parameters_at(x, k));
        }
        if (dim > 1) {
            double pValue = Parameter_value(parameter);
            if(parameterization == DISTRIBUTION_EXPONENTIAL_MEAN){
                for (size_t k = 0; k < dim; k++) {
                    out[k] = 1.0 / pValue;
                }
            }
            else{
                for (size_t k = 0; k < dim; k++) {
                    out[k] = pValue;
                }
            }
        }
    }
    else{
        if(parameterization == DISTRIBUTION_EXPONENTIAL_MEAN){
            for(size_t k = 0; k < Parameter_size(parameter); k++){
                out[k] = 1.0/Parameter_value_at(parameter, k);
            }
        }
        else{
            memcpy(out, Parameter_values(parameter), sizeof(double)*Parameter_size(parameter));
        }
    }
}

double DistributionModel_exp_log_prob(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    size_t dimX = Parameters_count(dm->x);
    _expand_parameter_to_lambda(dm->parameterization, dm->x, Parameters_at(dm->parameters, 0), dm->tempp);
    size_t index = 0;
    for(size_t j = 0; j < dimX; j++){
        Parameter* x = Parameters_at(dm->x, j);
        size_t dim = Parameter_size(x);
        const double* values = Parameter_values(x);
        for (size_t i = 0; i < dim; i++) {
            dm->lp += log(dm->tempp[index]) - dm->tempp[index] * values[i];
            index++;
        }
    }
    
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_exp_log_prob_grad(DistributionModel* dm, const Parameters* parameters){
    size_t dimX = Parameters_count(dm->x);
    Parameter* parameter = Parameters_at(dm->parameters, 0); // lambda or mean
    _expand_parameter_to_lambda(DISTRIBUTION_EXPONENTIAL_RATE, dm->x, parameter, dm->tempp); // dm->tempp contains lambda or mean
    size_t index = 0;

	for (size_t i = 0; i < dimX; i++) {
        Parameter* x = Parameters_at(dm->x, i);
        Parameter* xx = Parameters_depends(parameters, x); // check if it is x or x->transform->parameter
        // x is constrained and matches support of distribution (e.g. branch lengths)
        // xx is unconstrained (e.g. transformed branch lengths)
		if (xx != NULL) {
            size_t sizeX = Parameter_size(x);
            for (size_t j = 0; j < sizeX; j++) {
                dm->tempx[j] = -dm->tempp[index];
                x->grad[j] += dm->tempx[j];
                index++;
            }
            if(xx != x){
                x->transform->backward(x->transform, dm->tempx);
            }
		}
	}

    Parameter* parameterx = Parameters_depends(parameters, parameter);
    if(parameterx != NULL){
        index = 0;
        for (size_t i = 0; i < dimX; i++) {
            Parameter* x = Parameters_at(dm->x, i);
            const double* xValues = Parameter_values(x);
            size_t sizeX = Parameter_size(x);
            if(dm->parameterization == DISTRIBUTION_EXPONENTIAL_RATE){
                for (size_t j = 0; j < sizeX; j++) {
                    dm->tempx[index] = 1.0/dm->tempp[index] - xValues[j];
                    parameter->grad[index] += dm->tempx[index];
                    index++;
                }
            }
            else {
                for (size_t j = 0; j < sizeX; j++) {
                    dm->tempx[index] = (xValues[j] - dm->tempp[index])/(dm->tempp[index]*dm->tempp[index]);
                    parameter->grad[index] += dm->tempx[index];
                    index++;
                }
            }
        }
        if(parameterx != parameter){
            parameter->transform->backward(parameter->transform, dm->tempx);
        }
    }
	return 0;
}

void DistributionModel_exp_log_prob_hessian_diag(DistributionModel* dm, const Parameters* parameters){
    Parameter* parameter = Parameters_at(dm->parameters, 0);
    Parameter* parameterx = Parameters_depends(parameters, parameter);
    _expand_parameter_to_lambda(DISTRIBUTION_EXPONENTIAL_RATE, dm->x, parameter, dm->tempx); // dm->tempx contains lambda or mean
    
    if(parameterx != NULL){
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);

            if(dm->parameterization == DISTRIBUTION_EXPONENTIAL_RATE){
                for(size_t j = 0; j < dim; j++){
                    dm->tempp[index] = -1.0/(dm->tempx[index]*dm->tempx[index]);
                    parameter->grad[index] += dm->tempp[index];
                    index++;
                }
            }
            else{
                for(size_t j = 0; j < dim; j++){
                    dm->tempp[index] = (dm->tempx[index] - 2.0*xValues[j])/(dm->tempx[index]*dm->tempx[index]*dm->tempx[index]);
                    parameter->grad[index] += dm->tempp[index];
                    index++;
                }
            }
        }
        if(parameter != parameterx){
            parameter->transform->backward(parameter->transform, dm->tempp);
        }
    }

    // second derivative wrt x is zero
}

static void DistributionModel_exp_sample(DistributionModel* dm){
    Parameter* parameter = Parameters_at(dm->parameters, 0);
    size_t dimX = Parameters_count(dm->x);
    _expand_parameter_to_lambda(dm->parameterization, dm->x, parameter, dm->tempp);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(parameter) == 1){
        double lambdaValue = Parameter_value(parameter);
        if(dm->parameterization == DISTRIBUTION_EXPONENTIAL_MEAN){
            lambdaValue = 1.0 / lambdaValue; // convert mean to rate
        }
        
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
#ifndef GSL_DISABLED
                dm->tempx[i] = gsl_ran_exponential(dm->rng, 1.0/ lambdaValue);
#else
                dm->tempx[i] = rexp(lambdaValue);
#endif
            }
			Parameter_set_values(x, dm->tempx);
        }
    }
    // multiple parameter and len(x) >= 1 and \sum_i len(x_i) == len(parameter)
    else{
        size_t index = 0;
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
#ifndef GSL_DISABLED
                dm->tempx[i] = gsl_ran_exponential(dm->rng, 1.0/dm->tempp[index]);
#else
                dm->tempx[i] = rexp(dm->tempp[index]);
#endif
                index++;
            }
			Parameter_set_values(x, dm->tempx);
        }
    }
}

DistributionModel* new_ExponentialDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_EXPONENTIAL;
	dm->parameterization = parameterization;
    dm->log_prob = DistributionModel_exp_log_prob;
    dm->log_prob_grad = DistributionModel_exp_log_prob_grad;
    dm->log_prob_hessian_diag = DistributionModel_exp_log_prob_hessian_diag;
    dm->sample = DistributionModel_exp_sample;
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = INFINITY;
	return dm;
}

Model* new_ExponentialDistributionModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
	size_t paramCount = Parameter_size(Parameters_at(x, 0));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(1);
    Parameters_set_name2(parameters, "parameters.exp");
	Parameter* parameter = NULL;
    distribution_parameterization parameterization = DISTRIBUTION_EXPONENTIAL_RATE;

    // empirical
    if (file != NULL) {
        parameterization = DISTRIBUTION_EXPONENTIAL_MEAN;
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* pValues = malloc(sizeof(double)*paramCount);
        
        for (size_t i = 0; i < paramCount; i++) {
            const double* vec = Vector_data(samples[i]);
            pValues[i] = mean(vec, Vector_length(samples[i]));
            // Parameters_move(parameters[0], new_Parameter("mean", m, new_Constraint(0, INFINITY)));
            free_Vector(samples[i]);
        }

		parameter = new_Parameter2("mean", pValues, paramCount, new_Constraint(0, INFINITY));

		free(pValues);
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        parameterization = DISTRIBUTION_EXPONENTIAL_MEAN;
		double* pValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            pValues[i] = 1;
    	}
		parameter = new_Parameter2("mean", pValues, paramCount, new_Constraint(0, INFINITY));
		free(pValues);
	}
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        json_node* parameter_node = get_json_node(parameters_node, "lambda");
        if(parameter_node == NULL){
            parameter_node = get_json_node(parameters_node, "mean");
            if(parameter_node != NULL){
                parameterization = DISTRIBUTION_EXPONENTIAL_MEAN;
            }
        }

        if(parameter_node == NULL){
            fprintf(stderr, "Exponential distribution should be parametrized with mean or lambda\n");
            exit(13);
        }

        parameter = distmodel_parse_parameter(parameter_node, hash, "", 0, INFINITY);
    }
    
	Parameters_move(parameters, parameter);

    DistributionModel* dm = new_ExponentialDistributionModel_with_parameters(parameters, x, parameterization);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
	
	return model;
}
