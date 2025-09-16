//
//  distnormal.c
//  physher
//
//  Created by mathieu on 24/1/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distnormal.h"

#include <gsl/gsl_randist.h>
#include <math.h>
#include <strings.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"
#include "transforms.h"
#include "gaussian.h"

#define LOG_TWO_PI (log(2.0)+log(M_PI))

double dnormal(double x, double mu, double sigma){
    return gsl_ran_gaussian_pdf(x - mu, sigma);
}

double rnormal(const gsl_rng *rng, double mu, double sigma){
    return gsl_ran_gaussian(rng, sigma) + mu;
}

double DistributionModel_normal_logP(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* mu = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigma = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double sigmaValue = sigma[0];
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigmaValue = sqrt(1.0/sigmaValue);
        }
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(gsl_ran_gaussian_pdf(values[i] - *mu, sigmaValue));
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
                double sigmaValue = sigma[index];
                if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                    sigmaValue = sqrt(1.0/sigmaValue);
                }
                dm->lp += log(gsl_ran_gaussian_pdf(values[i] - mu[index], sigmaValue));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

static void _expand_parameters(Parameters* x, Parameters* parameters, double** muValues,
                               double** sigmaValues) {
    Parameter* mu = Parameters_at(parameters, 0);
    Parameter* sigma = Parameters_at(parameters, 1);
    *muValues = Parameter_values(mu);
    *sigmaValues = Parameter_values(sigma);

    if (Parameter_size(mu) == 1) {
        size_t dim = 0;
        for (size_t k = 0; k < Parameters_count(x); k++) {
            dim += Parameter_size(Parameters_at(x, k));
        }
        if (dim > 1) {
            *muValues = dvector(dim);
            *sigmaValues = dvector(dim);
            double muValue = Parameter_value(mu);
            double sigmaValue = Parameter_value(sigma);
            for (size_t k = 0; k < dim; k++) {
                *muValues[k] = muValue;
                *sigmaValues[k] = sigmaValue;
            }
        }
    }
}

double DistributionModel_normal_gradient2(DistributionModel* dm, const Parameters* parameters){
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* sigma = Parameters_at(dm->parameters, 1);
    double* muValues = NULL;
    double* sigmaValues = NULL;

    _expand_parameters(dm->x, dm->parameters, &muValues, &sigmaValues);

    Parameter* mux = Parameters_depends(parameters, mu);
    Parameter* sigmax = Parameters_depends(parameters, sigma);

    if(mux != NULL){
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                double sigmaValue = sigmaValues[index];
                if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                    sigmaValue = sqrt(1.0/sigmaValue);
                }
                dm->tempp[j] = (xValues[j] - muValues[index])/(sigmaValue*sigmaValue);
                mu->grad[j] += dm->tempp[j];
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
            for(size_t j = 0; j < dim; j++){
                double sigmaValue = sigmaValues[index];
                if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                    dm->tempp[j] = 1.0/(2.0*sigmaValue) - pow(xValues[j] - muValues[index], 2.0)/2.0;
                }
                else{
                    dm->tempp[j] =  (muValues[index]*muValues[index] - sigmaValue*sigmaValue - 2.0*muValues[index]*xValues[j] + xValues[j]*xValues[j])/(sigmaValue*sigmaValue*sigmaValue);
                }
                sigma->grad[j] += dm->tempp[j];
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
        
            for(size_t j = 0; j < dim; j++){
                double sigmaValue = sigmaValues[index];
                if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                    sigmaValue = sqrt(1.0/sigmaValue);
                }
                dm->tempp[j] = (muValues[index] - xValues[j])/(sigmaValue*sigmaValue);
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
    return 0;
}

double DistributionModel_normal_dlogP(DistributionModel* dm, const Parameter* p){
    //TODO: implement
    // // derivative wrt x
    // if (p == dm->x) {
    //     Parameter* mu = Parameters_at(dm->parameters, 0);
    //     Parameter* sigma = Parameters_at(dm->parameters, 1);
    //     size_t dim = Parameter_size(dm->x);

    //     double muValue = -1;
    //     double sigmaValue = -1

    //     if(Parameter_size(mu) > 1){
    //         muValue = Parameter_value_at(mu, i);
    //         sigmaValue = Parameter_value_at(sigma, i);
    //     }
    //     else{
    //         muValue = Parameter_value(mu);
    //         sigmaValue = Parameter_value(sigma);
    //     }
        
    //     if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
    //         sigmaValue = sqrt(1.0/sigmaValue);
    //     }
    //     return (mu - Parameter_value(p))/sigma/sigma;
    // }
    return 0;
}

double DistributionModel_normal_d2logP(DistributionModel* dm, const Parameter* p){
    //TODO: implement
    // // derivative wrt x
    // if (p == dm->x) {
    //     Parameter* mu = Parameters_at(dm->parameters, 0);
    //     Parameter* sigma = Parameters_at(dm->parameters, 1);
    //     size_t dim = Parameter_size(dm->x);

    //     double muValue = -1;
    //     double sigmaValue = -1

    //     if(Parameter_size(mu) > 1){
    //         muValue = Parameter_value_at(mu, i);
    //         sigmaValue = Parameter_value_at(sigma, i);
    //     }
    //     else{
    //         muValue = Parameter_value(mu);
    //         sigmaValue = Parameter_value(sigma);
    //     }
        
    //     if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
    //         sigmaValue = sqrt(1.0/sigmaValue);
    //     }
    //     return -1.0/sigma/sigma;
    // }
    return 0;
}


void DistributionModel_normal_sample(DistributionModel* dm){
    const double* mu = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigma = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double sigmaValue = *sigma;
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigmaValue = sqrt(1.0/sigmaValue);
        }
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_gaussian(dm->rng, sigmaValue) + *mu;
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
                double sigmaValue = sigma[index];
                if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                    sigmaValue = sqrt(1.0/sigmaValue);
                }
                // dm->lp += log(gsl_ran_gaussian_pdf(values[i] - mu[index], sigmaValue));
                dm->tempx[i] = gsl_ran_gaussian(dm->rng, sigmaValue) + mu[index];
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

static void DistributionModel_normal_rsample(DistributionModel* dm){
    const double* mu = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigma = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double sigmaValue = *sigma;
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigmaValue = sqrt(1.0/sigmaValue);
        }
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t sizeX = Parameter_size(x);
            const double* values = Parameter_values(x);
            double* temp = dvector(sizeX);
            for (size_t i = 0; i < sizeX; i++) {
                dm->tempx[i] = rnorm();
                temp[i] = dm->tempx[i] * sigmaValue + *mu;
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
                double sigmaValue = sigma[index];
                if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                    sigmaValue = sqrt(1.0/sigmaValue);
                }
                dm->tempx[index] = rnorm();
                temp[i] = dm->tempx[index] * sigmaValue + mu[index];
                index++;
            }
            Parameter_set_values(x, temp);
            free(temp);
        }
    }
}

// Calculate the gradient of the log PDF at x wrt its parameters where x was sampled using the reparmetrization trick
// eta ~ Normal(0, 1)
// x = μ + σ*eta
// dL(x)/dμ = dL(x)/dx * dx/dμ = dL(x)/dx * 1
// dL/dσ = dL/dx * dx/dσ = dL/dx * eta
static void DistributionModel_normal_rgradient(DistributionModel* dm) {
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* sigma = Parameters_at(dm->parameters, 1);
    const double* muValues = Parameter_values(mu);
    const double* sigmaValues = Parameter_values(sigma);
    size_t dimX = Parameters_count(dm->x);
    size_t index = 0;
    for (size_t j = 0; j < dimX; j++) {
        Parameter* x = Parameters_at(dm->x, j);
        size_t dim = Parameter_size(x);
        for (size_t i = 0; i < dim; i++) {
            double sigmaValue = sigmaValues[index];
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigmaValue = sqrt(1.0 / sigmaValue);
            }
            // dm->tempp[i] = x->grad[i] * (Parameter_value_at(x,
            // i)-muValues[i])/sigmaValue;
            dm->tempp[i] = x->grad[i] * dm->tempx[index];  // eta
            sigma->grad[index] += dm->tempp[i];
            mu->grad[index] += x->grad[i];
            index++;
        }
        if (sigma->transform != NULL) {
            sigma->transform->backward(sigma->transform, dm->tempp);
        }
    }
}

static double DistributionModel_normal_entropy(DistributionModel* dm) {
    // Entropy: 0.5(1 + log(2 \pi s^2))
    // 0.5(1+log(2\pi) + sum(log(sigma))
    const double* sigmaValues = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimP = Parameter_size(Parameters_at(dm->parameters, 1));
    double entropy = 0;
    for (size_t i = 0; i < dimP; i++) {
        double sigmaValue = sigmaValues[i];
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigmaValue = sqrt(1.0 / sigmaValue);
        }
        entropy += log(sigmaValue);
    }
    entropy += 0.5 * dimP * (1.0 + LOG_TWO_PI);
    return entropy;
}

static void DistributionModel_normal_entropy_gradient(DistributionModel* dm,
                                                      const Parameters* parameters) {
    Parameter* sigma = Parameters_at(dm->parameters, 1);
    Parameter* sigmax = Parameters_depends(parameters, sigma);
    const double* sigmaValues = Parameter_values(Parameters_at(dm->parameters, 1));
    if (sigmax != NULL) {
        size_t dimP = Parameter_size(sigma);
        for (size_t i = 0; i < dimP; i++) {
            // mu->grad[i] += 0.0;
            dm->tempp[i] = 1.0 / sigmaValues[i];
            sigma->grad[i] += dm->tempp[i];
        }
        if (sigma != sigmax) {
            sigma->transform->backward(sigma->transform, dm->tempp);
        }
    }
}

DistributionModel* new_NormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_NORMAL;
    dm->parameterization = parameterization;
    dm->logP = DistributionModel_normal_logP;
    // dm->logP_with_values = DistributionModel_normal_logP_with_values;
    dm->gradient2 = DistributionModel_normal_gradient2;
    dm->rgradient = DistributionModel_normal_rgradient;
    dm->dlogP = DistributionModel_normal_dlogP;
    dm->d2logP = DistributionModel_normal_d2logP;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_normal_sample;
    dm->rsample = DistributionModel_normal_rsample;
    dm->entropy = DistributionModel_normal_entropy;
    dm->gradient_entropy = DistributionModel_normal_entropy_gradient;
    // dm->sample_evaluate = DistributionModel_normal_sample_evaluate;
    dm->shift = -INFINITY; // can't be shifted
    return dm;
}

Model* new_NormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    Parameter* mu = NULL;
    Parameter* sigma = NULL;
    distribution_parameterization parameterization = DISTRIBUTION_NORMAL_MEAN_SIGMA;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* muValues = malloc(sizeof(double)*paramCount);
        double* sigmaValues = malloc(sizeof(double)*paramCount);
        
        for (size_t i = 0; i < paramCount; i++) {
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
        if(sigmaNode == NULL){
            sigmaNode = get_json_node(parametersNode, "tau");
            parameterization = DISTRIBUTION_NORMAL_MEAN_TAU;
        }

        if(muNode == NULL || sigmaNode == NULL){
            fprintf(stderr, "Normal distribution should be parametrized with mean and (sigma or tau)\n");
            exit(13);
        }
        mu = distmodel_parse_parameter(muNode, hash, "", -INFINITY, INFINITY);
        sigma = distmodel_parse_parameter(sigmaNode, hash, "", 0.0, INFINITY);
        // Constraint_set_lower(sigma->cnstr, 0.0);
    //    printf("%e\n", Constraint_lower(mu->cnstr));
    //    printf("%e\n", Constraint_lower(sigma->cnstr));
    //    if(Constraint_lower(sigma->cnstr) == 0)  exit(2);
    }

    Parameters_move(parameters, mu);
    Parameters_move(parameters, sigma);

    DistributionModel*dm = new_NormalDistributionModel_with_parameters(parameters, x, parameterization);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}

#pragma mark Halfnormal

double DistributionModel_half_normal_logP(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* mu = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* sigma = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double sigmaValue = sigma[0];
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigmaValue = sqrt(1.0/sigmaValue);
        }
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(sqrt(2)/(sigmaValue*sqrt(M_PI))) - pow(values[i], 2.0)/(2.0*sigmaValue*sigmaValue);
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
                double sigmaValue = sigma[index];
                if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                    sigmaValue = sqrt(1.0/sigmaValue);
                }
                dm->lp += log(sqrt(2)/(sigmaValue*sqrt(M_PI))) - pow(values[i], 2.0)/(2.0*sigmaValue*sigmaValue);
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_half_normal_gradient2(DistributionModel* dm, const Parameters* parameters){
    // size_t offset = 0;
    // for(size_t i = 0; i < Parameters_count(parameters); i++){
    //     Parameter* p = Parameters_at(parameters, i);
    //     // derivative wrt x
    //     if (p == dm->x) {
    //         Parameter* sigma = Parameters_at(dm->parameters, 0);
    //         size_t dim = Parameter_size(dm->x);

    //         double muValue = -1;
    //         double sigmaValue = -1
    //         size_t dim = Parameter_size(sigma);

    //         if(dim > 1){
    //             for(size_t j = 0; j < dim; j++){
    //                 sigmaValue = Parameter_value_at(sigma, j);
    //                 if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
    //                     sigmaValue = sqrt(1.0/sigmaValue);
    //                 }
    //                 grad[offset] = -Parameter_value_at(p, j)/sigmaValue/sigmaValue;
    //                 offset++;
    //             }
    //         }
    //         else{
    //             sigmaValue = Parameter_value(sigma);
    //             if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
    //                 sigmaValue = sqrt(1.0/sigmaValue);
    //             }
    //             for (size_t j = 0; j < dim; j++) {
    //                 grad[offset] = -Parameter_value_at(p, j)/sigmaValue/sigmaValue;
    //                 offset++;
    //             }
    //         }
    //     }
    //     //TODO: wrt sigma
    //     else if(p == Parameters_at(dm->parameters, 0)){
    //         error("DistributionModel_half+normal_gradient2 not implemented for mu and sigma\n");
    //     }
    //     else{
    //         offset += Parameter_size(p);
    //     }
    // }
    return 0;
}

double DistributionModel_half_normal_dlogP(DistributionModel* dm, const Parameter* p){
    // // derivative wrt x
    // for (int i = 0; i < Parameters_count(dm->x); i++) {
    //     if (p == Parameters_at(dm->x,i)) {
    //         double sigma = Parameters_value(dm->parameters[0], 0);
    //         if(Parameters_count(dm->parameters[0]) > 1){
    //             sigma = Parameters_value(dm->parameters[0], i);
    //         }
    //         return -Parameters_value(dm->x, 0)/sigma/sigma;
    //     }
    // }
    return 0;
}

double DistributionModel_half_normal_d2logP(DistributionModel* dm, const Parameter* p){
    // for (int i = 0; i < Parameters_count(dm->x); i++) {
    //     if (p == Parameters_at(dm->x,i)) {
    //         double sigma = Parameters_value(dm->parameters[0], 0);
    //         if(Parameters_count(dm->parameters[0]) > 1){
    //             sigma = Parameters_value(dm->parameters[0], i);
    //         }
    //         return -1.0/sigma/sigma;
    //     }
    // }
    return 0;
}

static void DistributionModel_half_normal_sample(DistributionModel* dm){
    fprintf(stderr, "DistributionModel_half_normal_sample not yet implemented\n");
    exit(2);
}


// static double DistributionModel_half_normal_sample_evaluate(DistributionModel* dm){
//     fprintf(stderr, "DistributionModel_half_normal_sample not DistributionModel_half_normal_sample_evaluate implemented\n");
//     exit(2);
//     return 0;
// }

DistributionModel* new_HalfNormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization param){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_HALFNORMAL;
    dm->parameterization = param;
    dm->logP = DistributionModel_half_normal_logP;
    // dm->logP_with_values = DistributionModel_half_normal_logP_with_values;
    dm->gradient2 = DistributionModel_half_normal_gradient2;
    dm->dlogP = DistributionModel_half_normal_dlogP;
    dm->d2logP = DistributionModel_half_normal_d2logP;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_half_normal_sample;
    // dm->sample_evaluate = DistributionModel_half_normal_sample_evaluate;
    dm->shift = 0;
    return dm;
}


Model* new_HalfNormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(1);
    Parameter* sigma = NULL;
    distribution_parameterization parameterization = DISTRIBUTION_NORMAL_MEAN_SIGMA;
    bool tau = false;

    char* parameterization_string = get_json_node_value_string(node, "parameterization");
    if(parameterization_string!= NULL && strcasecmp(parameterization_string, "tau") == 0){
        tau = true;
    }

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* sigmaValues = malloc(sizeof(double)*paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            const double* vec = Vector_data(samples[i]);
            double sum2 = 0;
            size_t n = Vector_length(samples[i]);
            for(size_t j = 0; j < n; j++){
                sum2 += vec[j]*vec[j];
            }
            sigmaValues[i] = sqrt(sum2/n); // ML estimate
            if(tau){
                sigmaValues[i] = 1.0/sigmaValues[i];
            }
            free_Vector(samples[i]);
        }
        
        sigma = new_Parameter2("sigma", sigmaValues, paramCount, new_Constraint(0, INFINITY));
        free(sigmaValues);
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* sigmaValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            sigmaValues[i] = 1.0;
        }
        sigma = new_Parameter2("sigma", sigmaValues, paramCount, new_Constraint(0, INFINITY));
        free(sigmaValues);
    }
    else{
        json_node* parametersNode = get_json_node(node, "parameters");
        json_node* sigmaNode = get_json_node(parametersNode, "sigma");
        if(sigmaNode == NULL){
            sigmaNode = get_json_node(parametersNode, "tau");
            parameterization = DISTRIBUTION_NORMAL_MEAN_TAU;
        }

        if(sigmaNode == NULL){
            fprintf(stderr, "Half Normal distribution should be parametrized with sigma or tau\n");
            exit(13);
        }

        sigma = distmodel_parse_parameter(sigmaNode, hash, "", 0, INFINITY);
    }

    Parameters_move(parameters, sigma);
        
    DistributionModel*dm = new_HalfNormalDistributionModel_with_parameters(parameters, x, parameterization);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
