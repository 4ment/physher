//
//  distlognormal.c
//  physher
//
//  Created by mathieu on 17/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distlognormal.h"

#include <strings.h>

#ifndef GSL_DISABLED
#include <gsl/gsl_randist.h>
#endif
#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"
#include "gaussian.h"
#include "lognormal.h"


double DistributionModel_lognormal_log_prob(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* mu = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigma = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
#ifdef GSL_DISABLED
            dm->lp += log(dlnorm(values[i], *mu, *sigma));
#else
            dm->lp += log(gsl_ran_lognormal_pdf(values[i], *mu, *sigma));
                
#endif
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
#ifdef GSL_DISABLED
                dm->lp += log(dlnorm(values[i], mu[index], sigma[index]));
#else
                dm->lp += log(gsl_ran_lognormal_pdf(values[i], mu[index], sigma[index]));
#endif
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_lognormal_log_prob_gradient(DistributionModel* dm,
                                             const Parameters* parameters) {
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* sigma = Parameters_at(dm->parameters, 1);
    double* muValues = NULL;
    double* sigmaValues = NULL;

    distmodel_expand_2parameters(dm->x, dm->parameters, &muValues, &sigmaValues);

    Parameter* mux = Parameters_depends(parameters, mu);
    Parameter* sigmax = Parameters_depends(parameters, sigma);

    if (mux != NULL) {
        size_t index = 0;
        for (size_t k = 0; k < Parameters_count(dm->x); k++) {
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for (size_t j = 0; j < dim; j++) {
                dm->tempp[index] = (log(xValues[j]) - muValues[index]) / pow(sigmaValues[index], 2.0);
                mu->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if (mux != mu) {
            mu->transform->backward(mu->transform, dm->tempp);
        }
    }

    if (sigmax != NULL) {
        size_t index = 0;
        for (size_t k = 0; k < Parameters_count(dm->x); k++) {
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            double z = 0.0;
            for (size_t j = 0; j < dim; j++) {
                z = (log(xValues[j]) - muValues[index]) / sigmaValues[index];
                dm->tempp[index] = -1.0 / sigmaValues[index] + z*z / sigmaValues[index];
                sigma->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if (sigmax != sigma) {
            sigma->transform->backward(sigma->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for (size_t k = 0; k < Parameters_count(dm->x); k++) {
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        size_t sizeX = Parameter_size(x);
        if (xx != NULL) {
            const double* xValues = Parameter_values(x);
            for (size_t j = 0; j < sizeX; j++) {
                dm->tempp[j] = -1.0 / xValues[j] - (log(xValues[j]) - muValues[index]) / (sigmaValues[index]*sigmaValues[index] * xValues[j]);
                x->grad[j] += dm->tempp[j];
                index++;
            }
            if (x != xx) {
                x->transform->backward(x->transform, dm->tempp);
            }
        } else {
            index += sizeX;
        }
    }

    if (Parameter_values(mu) != muValues) {
        free(muValues);
        free(sigmaValues);
    }
    return 0;
}

void DistributionModel_lognormal_log_prob_hessian_diag(DistributionModel* dm, const Parameters* parameters){
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* sigma = Parameters_at(dm->parameters, 1);
    double* muValues = NULL;
    double* sigmaValues = NULL;

    distmodel_expand_2parameters(dm->x, dm->parameters, &muValues, &sigmaValues);

    Parameter* mux = Parameters_depends(parameters, mu);
    Parameter* sigmax = Parameters_depends(parameters, sigma);

    if(mux != NULL){
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                dm->tempp[index] = -1.0/(sigmaValues[index]*sigmaValues[index]);
                mu->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(mu != mux){
            mu->transform->backward(mu->transform, dm->tempp);
        }
    }

    if(sigmax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            double z = 0.0;
            for(size_t j = 0; j < dim; j++){ 
                z = (log(xValues[j]) - muValues[index] )/sigmaValues[index];
                dm->tempp[index] = (1.0-3.0*z*z)/(sigmaValues[index]*sigmaValues[index]);
                sigma->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(sigma != sigmax){
            sigma->transform->backward(sigma->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for(size_t k = 0; k < Parameters_count(dm->x); k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        if (xx != NULL) {
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            double z = 0.0;
            for(size_t j = 0; j < dim; j++){
                z = (log(xValues[j]) - muValues[index] )/sigmaValues[index];
                dm->tempp[j] = z/(xValues[j]*xValues[j]*sigmaValues[index]);
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

    if(Parameter_values(mu) != muValues){
        free(muValues);
        free(sigmaValues);
    }
}

static void DistributionModel_lognormal_sample(DistributionModel* dm){
    const double* mu = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigma = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
#ifdef GSL_DISABLED
                dm->tempx[i] = rlnorm(*mu, *sigma);
#else
                dm->tempx[i] = gsl_ran_lognormal(dm->rng, *mu, *sigma);
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
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
#ifdef GSL_DISABLED
                dm->tempx[i] = rlnorm(mu[index], sigma[index]);
#else
                dm->tempx[i] = gsl_ran_lognormal(dm->rng, mu[index], sigma[index]);
#endif
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

static void DistributionModel_lognormal_rsample(DistributionModel* dm){
    const double* mu = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigma = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t sizeX = Parameter_size(x);
            const double* values = Parameter_values(x);
            double* temp = dvector(sizeX);
            for (size_t i = 0; i < sizeX; i++) {
#ifdef GSL_DISABLED
                dm->tempx[i] = rnorm();
#else
                dm->tempx[i] = gsl_ran_gaussian(dm->rng, 1.0);
#endif
                temp[i] = exp(dm->tempx[i] * *sigma + *mu);
            }
            Parameter_set_values(x, temp);
            free(temp);
        }
    }
    // multiple parameter and len(x) >= 1 and \sum_i len(x_i) == len(parameter)
    else{
        size_t index = 0;
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t sizeX = Parameter_size(x);
            const double* values = Parameter_values(x);
            double* temp = dvector(sizeX);
            for (size_t i = 0; i < sizeX; i++) {
#ifdef GSL_DISABLED
                dm->tempx[index] = rnorm();
#else
                dm->tempx[i] = gsl_ran_gaussian(dm->rng, 1.0);
#endif
                temp[i] = exp(dm->tempx[index] + mu[index]);
                index++;
            }
            Parameter_set_values(x, temp);
            free(temp);
        }
    }
}

// Calculate the gradient of the log PDF at x wrt its parameters where x was sampled using the reparmetrization trick
// eta ~ Normal(0, 1)
// x = exp(μ + σ*eta)
// dL(x)/dμ = dL(x)/dx * dx/dμ = dL(x)/dx * x
// dL/dσ = dL/dx * dx/dσ = dL/dx * eta*x
static void DistributionModel_lognormal_reparam_backprop(DistributionModel* dm) {
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* sigma = Parameters_at(dm->parameters, 1);
    const double* muValues = Parameter_values(mu);
    const double* sigmaValues = Parameter_values(sigma);
    size_t dimX = Parameters_count(dm->x);
    size_t index = 0;
    for (size_t j = 0; j < dimX; j++) {
        Parameter* x = Parameters_at(dm->x, j);
        const double* xValues = Parameter_values(x);
        size_t dim = Parameter_size(x);
        for (size_t i = 0; i < dim; i++) {
            mu->grad[index] += x->grad[i] * xValues[j];
            dm->tempp[i] = x->grad[i] * xValues[j] * dm->tempx[index];
            sigma->grad[index] += dm->tempp[i];
            index++;
        }
        if (sigma->transform != NULL) {
            sigma->transform->backward(sigma->transform, dm->tempp);
        }
    }
}

static double DistributionModel_lognormal_entropy(DistributionModel* dm) {
    // Entropy: log_{2}(sqrt(2\pi) \sigma e^{\mu + 1/2})
    const double* muValues = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigmaValues = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimP = Parameter_size(Parameters_at(dm->parameters, 1));
    double entropy = 0;
    for (size_t i = 0; i < dimP; i++) {
        entropy += log2(sigmaValues[i]) + muValues[i];
    }
    return entropy + log2(sqrt(2 * M_PI)) * dimP + dimP / 2;
}

static void DistributionModel_lognormal_entropy_gradient(DistributionModel* dm,
                                                         const Parameters* parameters) {
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* sigma = Parameters_at(dm->parameters, 1);
    Parameter* mux = Parameters_depends(parameters, mu);
    Parameter* sigmax = Parameters_depends(parameters, sigma);
    const double* sigmaValues = Parameter_values(sigma);
    size_t dimP = Parameter_size(mu);
    if (mux != NULL) {
        for (size_t i = 0; i < dimP; i++) {
            dm->tempp[i] = 1.0;
            mu->grad[i] += 1.0;
        }
        if (mux != mu) {
            mu->transform->backward(mu->transform, dm->tempp);
        }
    }
    if (sigmax != NULL) {
        for (size_t i = 0; i < dimP; i++) {
            dm->tempp[i] = 1.0 / sigmaValues[i];
            sigma->grad[i] += dm->tempp[i];
        }
        if (sigmax != sigma) {
            sigma->transform->backward(sigma->transform, dm->tempp);
        }
    }
}

DistributionModel* new_LogNormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_LOGNORMAL;
    dm->log_prob = DistributionModel_lognormal_log_prob;
    dm->log_prob_grad = DistributionModel_lognormal_log_prob_gradient;
    dm->reparam_backprop = DistributionModel_lognormal_reparam_backprop;
    dm->log_prob_hessian_diag = DistributionModel_lognormal_log_prob_hessian_diag;
    dm->sample = DistributionModel_lognormal_sample;
    dm->rsample = DistributionModel_lognormal_rsample;
    dm->entropy = DistributionModel_lognormal_entropy;
    dm->entropy_grad = DistributionModel_lognormal_entropy_gradient;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = INFINITY;
    return dm;
}

Model* new_LogNormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    Parameter* mu = NULL;
    Parameter* sigma = NULL;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* muValues = malloc(sizeof(double)*paramCount);
        double* sigmaValues = malloc(sizeof(double)*paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            //TODO: that coming from normal
            const double* vec = Vector_data(samples[i]);
            muValues[i] = mean(vec, Vector_length(samples[i]));
            sigmaValues[i] = sqrt(variance(vec, Vector_length(samples[i]), muValues[i]));
            free_Vector(samples[i]);
        }

        mu = new_Parameter2("mu", muValues, paramCount, new_Constraint(-INFINITY, INFINITY));
        sigma = new_Parameter2("sigma", sigmaValues, paramCount, new_Constraint(0, INFINITY));

        free(muValues);
        free(sigmaValues);
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* muValues = malloc(sizeof(double)*paramCount);
        double* sigmaValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            muValues[i] = 0;
            sigmaValues[i] = 1;
        }
        mu = new_Parameter2("mu", muValues, paramCount, new_Constraint(-INFINITY, INFINITY));
        sigma = new_Parameter2("sigma", sigmaValues, paramCount, new_Constraint(0, INFINITY));
        
        free(muValues);
        free(sigmaValues);
    }
    else{
        json_node* parametersNode = get_json_node(node, "parameters");
        json_node* muNode = get_json_node(parametersNode, "mu");
        json_node* sigmaNode = get_json_node(parametersNode, "sigma");

        if(muNode == NULL || sigmaNode == NULL){
            fprintf(stderr, "LogNormal distribution should be parametrized with mean and sigma\n");
            exit(13);
        }
    
        mu = distmodel_parse_parameter(muNode, hash, "", -INFINITY, INFINITY);
        sigma = distmodel_parse_parameter(sigmaNode, hash, "", 0.0, INFINITY);
    }

    Parameters_move(parameters, mu);
    Parameters_move(parameters, sigma);
    
    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
