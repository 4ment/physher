//
//  distbeta.c
//  physher
//
//  Created by mathieu on 17/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distbeta.h"

#include <strings.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"


double DistributionModel_beta_log_prob(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* beta = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(gsl_ran_beta_pdf(values[i], *alpha, *beta));
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
                dm->lp += log(gsl_ran_beta_pdf(values[i], alpha[index], beta[index]));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_beta_log_prob_gradient(DistributionModel* dm, const Parameters* parameters){
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* beta = Parameters_at(dm->parameters, 1);
    double* alphaValues = NULL;
    double* betaValues = NULL;

    distmodel_expand_2parameters(dm->x, dm->parameters, &alphaValues, &betaValues);

    Parameter* alphax = Parameters_depends(parameters, alpha);
    Parameter* betax = Parameters_depends(parameters, beta);

    if(alphax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                dm->tempp[index] = log(xValues[j])  - gsl_sf_psi(alphaValues[index]) + gsl_sf_psi(alphaValues[index] + betaValues[index]);
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
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                dm->tempp[index] =  log(xValues[j])  - gsl_sf_psi(betaValues[index]) + gsl_sf_psi(alphaValues[index] + betaValues[index]);
                beta->grad[index] += dm->tempp[index];
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
        if (xx != NULL) {
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
        
            for(size_t j = 0; j < dim; j++){
                dm->tempp[j] = (alphaValues[index] - 1)/xValues[j] - (betaValues[index] - 1)/(1-xValues[j]);
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
    return 0;
}

static void DistributionModel_beta_sample(DistributionModel* dm){
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* beta = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_beta(dm->rng, *alpha, *beta);
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
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_beta(dm->rng, alpha[index], beta[index]);
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

DistributionModel* new_BetaDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_BETA;
    dm->log_prob = DistributionModel_beta_log_prob;
    dm->log_prob_grad = DistributionModel_beta_log_prob_gradient;
    dm->sample = DistributionModel_beta_sample;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = 1;
    return dm;
}

Model* new_BetaDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = new_Parameters(1);
    distmodel_get_parameters(x_node, hash, x);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));
    Parameter* alpha = NULL;
    Parameter* beta = NULL;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* alphaValues = malloc(sizeof(double)*paramCount);
        double* betaValues = malloc(sizeof(double)*paramCount);
        
        for (size_t i = 0; i < paramCount; i++) {
            const double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            alphaValues[i] = m*(m*(1.0 - m)/v - 1.0);
            betaValues[i] = (1.0 - m)*(m*(1.0 - m)/v - 1.0);
            free_Vector(samples[i]);
        }

        alpha = new_Parameter2("alpha", alphaValues, paramCount, new_Constraint(0, INFINITY));
        beta = new_Parameter2("beta", betaValues, paramCount, new_Constraint(0, INFINITY));

        free(alphaValues);
        free(betaValues);
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
    else {
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "alpha") != 0 && strcasecmp(parameters_node->children[i]->key, "beta") != 0) {
                fprintf(stderr, "Beta distribution should be parametrized with alpha and beta\n");
                exit(13);
            }
        }

        json_node* alpha_node = get_json_node(parameters_node, "alpha");
        json_node* beta_node = get_json_node(parameters_node, "beta");
        alpha = new_Parameter_from_json(alpha_node, hash);
        beta = new_Parameter_from_json(beta_node, hash);
    }
    
    Parameters_add(parameters, alpha);
    Parameters_add(parameters, beta);
    
    DistributionModel* dm = new_BetaDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");

    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
