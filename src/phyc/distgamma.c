//
//  distgamma.c
//  physher
//
//  Created by Mathieu Fourment on 30/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distgamma.h"

#include <strings.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"

// gls use gsl_ran_gamma_pdf(x, a, b) where x is the value, a is the shape and b is the scale

double DistributionModel_log_gamma(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* beta = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double betaValue = beta[0];
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			betaValue = 1.0/betaValue;
		}

        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(gsl_ran_gamma_pdf(values[i] - dm->shift, *alpha, betaValue));
            }
        }
    }
    // multiple parameter and len(x) >= 1 and \sum_i len(x_i) == len(parameter)
    else{
        size_t index = 0;
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                double betaValue = beta[index];
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    betaValue = 1.0/betaValue;
                }
                dm->lp = log(gsl_ran_gamma_pdf(values[i] - dm->shift, alpha[i], betaValue));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_normal_log_prob_gradient(DistributionModel* dm, const Parameters* parameters){
    size_t dimX = Parameters_count(dm->x);
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* beta = Parameters_at(dm->parameters, 1);
    double* alphaValues = NULL;
    double* betaValues = NULL;

    distmodel_expand_2parameters(dm->x, dm->parameters, &alphaValues, &betaValues);

    Parameter* alphax = Parameters_depends(parameters, alpha);
    Parameter* betax = Parameters_depends(parameters, beta);

    if(alphax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                double betaValue = betaValues[index];
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			        betaValue = 1.0/betaValue;
		        }
                dm->tempp[index] = log(xValues[j]) - log(betaValues[index]) - gsl_sf_psi(alphaValues[index]);
                alpha->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(alpha != alphax){
            alpha->transform->backward(alpha->transform, dm->tempp);
        }
    }

    if(betax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    dm->tempp[j] = alphaValues[index]/betaValues[index] - xValues[j];
                }
                else{
                    dm->tempp[index] = xValues[j]/(betaValues[index]*betaValues[index]) - alphaValues[index]/betaValues[index];
                }
                beta->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(beta != betax){
            beta->transform->backward(beta->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for(size_t k = 0; k < dimX; k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        if (xx != NULL) {
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            
            if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
                for(size_t j = 0; j < dim; j++){
                    dm->tempp[j] = (alphaValues[index] - 1.0)/xValues[j] - 1.0/betaValues[index];
                    x->grad[j] += dm->tempp[j];
                    index++;
                }
            }
            else{
                for(size_t j = 0; j < dim; j++){
                    dm->tempp[j] = (alphaValues[index] - 1.0)/xValues[j] - betaValues[index];
                    x->grad[j] += dm->tempp[j];
                    index++;
                }
            }
            if(x != xx){
                x->transform->backward(x->transform, dm->tempp);
            }
        }
        else{
            index += Parameter_size(x);
        }
    }

    if(Parameter_values(alpha) != alphaValues){
        free(alphaValues);
        free(betaValues);
    }
    return 0;
}

void DistributionModel_normal_log_prob_hessian_diag(DistributionModel* dm, const Parameters* parameters){
    size_t dimX = Parameters_count(dm->x);
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* beta = Parameters_at(dm->parameters, 1);
    double* alphaValues = NULL;
    double* betaValues = NULL;

    distmodel_expand_2parameters(dm->x, dm->parameters, &alphaValues, &betaValues);

    Parameter* alphax = Parameters_depends(parameters, alpha);
    Parameter* betax = Parameters_depends(parameters, beta);

    if(alphax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                double betaValue = betaValues[index];
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			        betaValue = 1.0/betaValue;
		        }
                dm->tempp[index] = - gsl_sf_psi_1(alphaValues[index]);
                alpha->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(alpha != alphax){
            alpha->transform->backward(alpha->transform, dm->tempp);
        }
    }

    if(betax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    dm->tempp[j] = -alphaValues[index]/(betaValues[index]*betaValues[index]);
                }
                else{
                    dm->tempp[index] = -2.0*xValues[j]/(betaValues[index]*betaValues[index]*betaValues[index]) + alphaValues[index]/(betaValues[index]*betaValues[index]);
                }
                beta->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(beta != betax){
            beta->transform->backward(beta->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for(size_t k = 0; k < dimX; k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        if (xx != NULL) {
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                dm->tempp[j] = -(alphaValues[index] - 1.0)/(xValues[j]*xValues[j]);
                x->grad[j] += dm->tempp[j];
                index++;
            }
            if(x != xx){
                x->transform->backward(x->transform, dm->tempp);
            }
        }
        else{
            index += Parameter_size(x);
        }
    }

    if(Parameter_values(alpha) != alphaValues){
        free(alphaValues);
        free(betaValues);
    }
}

static void DistributionModel_gamma_sample(DistributionModel* dm){
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* beta = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double betaValue = beta[0];
        if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
            betaValue = 1.0/betaValue;
        }
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_gamma(dm->rng, *alpha, betaValue);
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
            const double* values = Parameter_values(x);
            double betaValue = beta[index];
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				betaValue = 1.0/betaValue;
			}
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_gamma(dm->rng, alpha[index], betaValue);
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

DistributionModel* new_GammaDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_GAMMA;
	dm->log_prob = DistributionModel_log_gamma;
	dm->log_prob_grad = DistributionModel_normal_log_prob_gradient;
    dm->log_prob_hessian_diag = DistributionModel_normal_log_prob_hessian_diag;
	dm->sample = DistributionModel_gamma_sample;
	dm->parameterization = parameterization;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = INFINITY;
	return dm;
}

Model* new_GammaDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    Parameter* alpha = NULL;
    Parameter* beta = NULL;
    
    distribution_parameterization parameterization = DISTRIBUTION_GAMMA_SHAPE_RATE;
    bool scale = false;
    
    char* parameterization_string = get_json_node_value_string(node, "parameterization");
    if(parameterization_string!= NULL && strcasecmp(parameterization_string, "scale") == 0){
        scale = true;
    }

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* alphaValues = malloc(sizeof(double)*paramCount);
        double* betaValues = malloc(sizeof(double)*paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            const double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            alphaValues[i] = m*m/v;
            if (scale) {
                betaValues[i] = v/m;
            }
            else{
				betaValues[i] = m/v;
            }
            free_Vector(samples[i]);
        }
		
		alpha = new_Parameter2("alpha", alphaValues, paramCount, new_Constraint(0, INFINITY));
        beta = new_Parameter2("beta", betaValues, paramCount, new_Constraint(0, INFINITY));

        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* alphaValues = malloc(sizeof(double)*paramCount);
        double* betaValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            alphaValues[i] = 0;
            betaValues[i] = 1;
        }
        alpha = new_Parameter2("alpha", alphaValues, paramCount, new_Constraint(0, INFINITY));
        beta = new_Parameter2("beta", betaValues, paramCount, new_Constraint(0, INFINITY));
        
        free(alphaValues);
        free(betaValues);
    }
    else{
        json_node* parametersNode = get_json_node(node, "parameters");
        json_node* alphaNode = get_json_node(parametersNode, "shape");
        json_node* betaNode = get_json_node(parametersNode, "rate");
        if(betaNode == NULL){
            betaNode = get_json_node(parametersNode, "scale");
            parameterization = DISTRIBUTION_GAMMA_SHAPE_SCALE;
        }

        if(betaNode == NULL){
            fprintf(stderr, "Gamma distribution should be parametrized with shape and (scale or rate)\n");
            exit(13);
        }
        
        alpha = distmodel_parse_parameter(alphaNode, hash, "", 0, INFINITY);
        beta = distmodel_parse_parameter(betaNode, hash, "", 0, INFINITY);
    }

    Parameters_move(parameters, alpha);
    Parameters_move(parameters, beta);

    DistributionModel* dm = new_GammaDistributionModel_with_parameters(parameters, x, parameterization);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");

    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
