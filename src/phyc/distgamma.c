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

// NOTE: GSL uses shape and scale parameterization

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
                dm->lp += log(gsl_ran_gamma_pdf(values[i] - dm->shift, alpha[i], betaValue));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_gamma_gradient2(DistributionModel* dm, const Parameters* parameters){
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* beta = Parameters_at(dm->parameters, 1);
    
    const double* alphaValues = Parameter_values(alpha);
    const double* betaValues = Parameter_values(beta);
    size_t mask = -(Parameter_size(alpha) != 1);

    Parameter* alphax = Parameters_depends(parameters, alpha);
    Parameter* betax = Parameters_depends(parameters, beta);
    if(alphax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            memset(dm->tempp, 0, dim * sizeof(double));
            for(size_t j = 0; j < dim; j++){
                double alphaValue = alphaValues[index & mask];
                double betaValue = betaValues[index & mask];
                size_t idx = j & mask;
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    dm->tempp[idx] = log(xValues[j]) + log(betaValue) - gsl_sf_psi(alphaValue);
                }
                else{
                    dm->tempp[idx] = log(xValues[j]) - log(betaValue) - gsl_sf_psi(alphaValue);
                }
                alpha->grad[idx] += dm->tempp[idx];
                index++;
            }
        }
        if(alpha != alphax){
            alpha->transform->backward(alpha->transform, dm->tempp);
        }
    }

    if(betax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            memset(dm->tempp, 0, dim * sizeof(double));
            for(size_t j = 0; j < dim; j++){
                double alphaValue = alphaValues[index & mask];
                double betaValue = betaValues[index & mask];
                size_t idx = j & mask;
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    dm->tempp[idx] = alphaValue/betaValue - xValues[j];
                }
                else{
                    dm->tempp[idx] =  (xValues[j] -alphaValue*betaValue)/(betaValue*betaValue);
                }
                beta->grad[idx] += dm->tempp[idx];
                index++;
            }
        }

        if(beta != betax){
            beta->transform->backward(beta->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for(size_t k = 0; k < Parameters_count(dm->x); k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        size_t sizeX = Parameter_size(x);
        if (xx != NULL) {
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < sizeX; j++){
                double alphaValue = alphaValues[index & mask];
                double betaValue = betaValues[index & mask];
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
                    betaValue = 1.0/betaValue;
                }
                dm->tempp[index] = (alphaValue-1.0)/xValues[j] - betaValue;
                x->grad[index] += dm->tempp[index];
                index++;
            }
            if(x != xx){
                x->transform->backward(x->transform, dm->tempp);
            }
        }
        else{
            index += sizeX;
        }
    }
    return 0;
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

// static double DistributionModel_gamma_sample_evaluate(DistributionModel* dm){
// 	double logP = 0;
// 	if(Parameters_count(dm->parameters[0]) > 1){
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double alpha = Parameters_value(dm->parameters[0], i);
// 			double beta = Parameters_value(dm->parameters[1], i);
// 			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
// 				beta = 1.0/beta;
// 			}
// 			double sample = gsl_ran_gamma(dm->rng, alpha, beta);
// 			Parameters_set_value(dm->x, i, sample + dm->shift);
// 			logP += log(gsl_ran_gamma_pdf(sample, alpha, beta));
// 		}
// 	}
// 	else{
// 		double alpha = Parameters_value(dm->parameters[0], 0);
// 		double beta = Parameters_value(dm->parameters[1], 1);
// 		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
// 			beta = 1.0/beta;
// 		}
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double sample = gsl_ran_gamma(dm->rng, alpha, beta);
// 			Parameters_set_value(dm->x, i, sample + dm->shift);
// 			logP += log(gsl_ran_gamma_pdf(sample, alpha, beta));
// 		}
// 	}
// 	return logP;
// }

DistributionModel* new_GammaDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_GAMMA;
	dm->logP = DistributionModel_log_gamma;
    dm->gradient2 = DistributionModel_gamma_gradient2;
	// dm->logP_with_values = DistributionModel_log_gamma_with_values;
    // if(Parameters_count(parameters[0]) > 1){
    // }
    // else{
	dm->sample = DistributionModel_gamma_sample;
	// dm->sample_evaluate = DistributionModel_gamma_sample_evaluate;
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
