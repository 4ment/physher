// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "distexp.h"

#include <string.h>
#include <strings.h>

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"


// The single distribution parameter is the rate lambda or the mean 1/lambda,
// depending on dm->parameterization. It is a shared scalar (mask == 0, mode 1) or
// one value per element of x (mask all-ones, mode 2); see distmodel.h.
static inline double _exp_rate(const DistributionModel* dm, double value){
    return dm->parameterization == DISTRIBUTION_EXPONENTIAL_MEAN ? 1.0/value : value;
}

double DistributionModel_log_exp(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    Parameter* parameter = Parameters_at(dm->parameters, 0);
    const double* values = Parameter_values(parameter);
    size_t mask = -(Parameter_size(parameter) != 1);
    size_t index = 0;
    for(size_t k = 0; k < Parameters_count(dm->x); k++){
        Parameter* x = Parameters_at(dm->x, k);
        const double* xValues = Parameter_values(x);
        for (size_t j = 0; j < Parameter_size(x); j++) {
            double lambda = _exp_rate(dm, values[index & mask]);
            dm->lp += log(lambda) - lambda * xValues[j];
            index++;
        }
    }
    dm->need_update = false;
    return dm->lp;
}

void DistributionModel_gradient_exp(DistributionModel* dm, Parameters* parameters){
    Parameter* parameter = Parameters_at(dm->parameters, 0);
    const double* values = Parameter_values(parameter);
    size_t mask = -(Parameter_size(parameter) != 1);

    // d/dx log p = -lambda
    size_t index = 0;
    for (size_t k = 0; k < Parameters_count(dm->x); k++) {
        Parameter* x = Parameters_at(dm->x, k);
        size_t sizeX = Parameter_size(x);
        // x is constrained and matches the support (e.g. branch lengths); xx may be
        // its unconstrained leaf (e.g. transformed branch lengths)
        Parameter* xx = Parameters_depends(parameters, x);
        if (xx != NULL) {
            for (size_t j = 0; j < sizeX; j++) {
                dm->tempx[j] = -_exp_rate(dm, values[(index + j) & mask]);
                x->grad[j] += dm->tempx[j];
            }
            if(xx != x){
                x->transform->backward(x->transform, dm->tempx);
            }
        }
        index += sizeX;
    }

    // d/dlambda log p = 1/lambda - x;  d/dmean log p = x/mean^2 - 1/mean
    // A shared scalar folds every element onto tempp[0] so that backward sees the
    // full sum.
    Parameter* px = Parameters_depends(parameters, parameter);
    if(px != NULL){
        memset(dm->tempp, 0, Parameter_size(parameter) * sizeof(double));
        index = 0;
        for (size_t k = 0; k < Parameters_count(dm->x); k++) {
            Parameter* x = Parameters_at(dm->x, k);
            const double* xValues = Parameter_values(x);
            for (size_t j = 0; j < Parameter_size(x); j++) {
                size_t idx = index & mask;
                double v = values[idx];
                double g;
                if(dm->parameterization == DISTRIBUTION_EXPONENTIAL_MEAN){
                    g = xValues[j]/(v*v) - 1.0/v;
                }
                else{
                    g = 1.0/v - xValues[j];
                }
                dm->tempp[idx] += g;
                parameter->grad[idx] += g;
                index++;
            }
        }
        if(px != parameter){
            parameter->transform->backward(parameter->transform, dm->tempp);
        }
    }
}

static void DistributionModel_exp_sample(DistributionModel* dm){
    const double* lambda = Parameter_values(Parameters_at(dm->parameters, 0));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
				dm->tempx[i] = gsl_ran_exponential(dm->rng, 1.0/ lambda[0]);
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
                dm->tempx[i] = gsl_ran_exponential(dm->rng, 1.0/ lambda[index]);
                index++;
            }
			Parameter_set_values(x, dm->tempx);
        }
    }
}

static void DistributionModel_exp_sample_mean(DistributionModel* dm){
 	const double* mean = Parameter_values(Parameters_at(dm->parameters, 0));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
				dm->tempx[i] = gsl_ran_exponential(dm->rng, mean[0]);
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
                dm->tempx[i] = gsl_ran_exponential(dm->rng, mean[index]);
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
    dm->logP = DistributionModel_log_exp;
    dm->gradient = DistributionModel_gradient_exp;
    dm->sample = DistributionModel_exp_sample;
	if(parameterization == DISTRIBUTION_EXPONENTIAL_RATE){
	}
	else{
		dm->sample = DistributionModel_exp_sample_mean;
	}
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
