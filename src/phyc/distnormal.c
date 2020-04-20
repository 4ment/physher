//
//  distnormal.c
//  physher
//
//  Created by mathieu on 24/1/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distnormal.h"

#include <gsl/gsl_randist.h>
#include <strings.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"


double DistributionModel_normal_logP(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double x = Parameters_value(dm->x, i);
            dm->lp += log(gsl_ran_gaussian_pdf(x - mu, sigma));
        }
    }
    else{
        double mu = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double x = Parameters_value(dm->x, i);
            dm->lp += log(gsl_ran_gaussian_pdf(x - mu, sigma));
        }
    }
    return dm->lp;
}

double DistributionModel_normal_logP_with_values(DistributionModel* dm, const double* values){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            logP += log(gsl_ran_gaussian_pdf(values[i] - mu, sigma));
        }
    }
    else{
        double mu = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            logP += log(gsl_ran_gaussian_pdf(values[i] - mu, sigma));
        }
    }
    return logP;
}

double DistributionModel_normal_dlogP(DistributionModel* dm, const Parameter* p){
    // derivative wrt x
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters[0], 0);
            double sigma = Parameters_value(dm->parameters[1], 0);
            if(Parameters_count(dm->parameters[0]) > 1){
                mu = Parameters_value(dm->parameters[0], i);
                sigma = Parameters_value(dm->parameters[1], i);
            }
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            return (mu - Parameter_value(p))/sigma/sigma;
        }
    }exit(2);
    return 0;
}

double DistributionModel_normal_d2logP(DistributionModel* dm, const Parameter* p){
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters[0], 0);
            double sigma = Parameters_value(dm->parameters[1], 0);
            if(Parameters_count(dm->parameters[0]) > 1){
                mu = Parameters_value(dm->parameters[0], i);
                sigma = Parameters_value(dm->parameters[1], i);
            }
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            return -1.0/sigma/sigma;
        }
    }
    return 0;
}

static void DistributionModel_normal_sample(DistributionModel* dm, double* samples){
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters[1], i);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            samples[i] = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters[0], i);
        }
    }
    else{
        double mean = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            samples[i] = gsl_ran_gaussian(dm->rng, sigma) + mean;
        }
    }
}


static double DistributionModel_normal_sample_evaluate(DistributionModel* dm){
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters[1], i);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double sample = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters[0], i);
            Parameters_set_value(dm->x, i, sample);
        }
    }
    else{
        double mean = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sample = gsl_ran_gaussian(dm->rng, sigma) + mean;
            Parameters_set_value(dm->x, i, sample);
        }
    }
    return DistributionModel_normal_logP(dm);
}

DistributionModel* new_NormalDistributionModel_with_parameters(Parameters** parameters, Parameters* x, distribution_parameterization parameterization){
    DistributionModel* dm = new_DistributionModel(parameters, 2, x);
    dm->type = DISTRIBUTION_NORMAL;
    dm->parameterization = parameterization;
    dm->logP = DistributionModel_normal_logP;
    dm->logP_with_values = DistributionModel_normal_logP_with_values;
    dm->dlogP = DistributionModel_normal_dlogP;
    dm->d2logP = DistributionModel_normal_d2logP;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_normal_sample;
    dm->sample_evaluate = DistributionModel_normal_sample_evaluate;
    dm->shift = -INFINITY; // can't be shifted
    return dm;
}

Model* new_NormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters** parameters = NULL;
    size_t parameters_dim = 2;
    bool tau = false;
    // empirical
    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        size_t paramCount = Parameters_count(x);
        parameters = malloc(sizeof(Parameters*)*2);
        parameters[0] = new_Parameters(paramCount);
        parameters[1] = new_Parameters(paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            Parameters_move(parameters[0], new_Parameter("mu", m, new_Constraint(-INFINITY, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("sigma", sqrt(v), new_Constraint(0, INFINITY)));
        }
    }
    else if(get_json_node(node, "parameters") == NULL){
        parameters = malloc(sizeof(Parameters*)*2);
        parameters[0] = new_Parameters(Parameters_count(x));
        parameters[1] = new_Parameters(Parameters_count(x));
        for (int i = 0; i < Parameters_count(x); i++) {
            Parameters_move(parameters[0], new_Parameter("mu", 0, new_Constraint(-INFINITY, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("sigma", 1, new_Constraint(0, INFINITY)));
        }
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "tau") == 0) {
                tau = true;
            }
            else if (strcasecmp(parameters_node->children[i]->key, "mean") != 0 && strcasecmp(parameters_node->children[i]->key, "sigma") != 0) {
                fprintf(stderr, "Normal distribution should be parametrized with mean and (sigma or tau)\n");
                exit(13);
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
        
        if (strcasecmp(Parameters_name2(parameters[0]), "mean") != 0) {
            Parameters* temp = parameters[0];
            parameters[0] = parameters[1];
            parameters[1] = temp;
        }
    }
    distribution_parameterization parameterization = tau ? DISTRIBUTION_NORMAL_MEAN_TAU: DISTRIBUTION_NORMAL_MEAN_SIGMA;
    DistributionModel*dm = new_NormalDistributionModel_with_parameters(parameters, x, parameterization);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
    
    return model;
}

#pragma mark Halfnormal

double DistributionModel_half_normal_logP(DistributionModel* dm){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters[0], i);
            if (dm->parameterization == DISTRIBUTION_HALFNORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double x = Parameters_value(dm->x, i);
            logP += log(sqrt(2)/(sigma*sqrt(M_PI))) - x*x/(2.0*sigma*sigma);
        }
    }
    else{
        double sigma = Parameters_value(dm->parameters[0], 0);
        if (dm->parameterization == DISTRIBUTION_HALFNORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double x = Parameters_value(dm->x, i);
            logP += log(sqrt(2)/(sigma*sqrt(M_PI))) - x*x/(2.0*sigma*sigma);
        }
    }
    return logP;
}

double DistributionModel_half_normal_logP_with_values(DistributionModel* dm, const double* values){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters[0], i);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            logP += log(sqrt(2)/(sigma*sqrt(M_PI))) - values[i]*values[i]/(2.0*sigma*sigma);
        }
    }
    else{
        double sigma = Parameters_value(dm->parameters[0], 0);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            logP += log(sqrt(2)/(sigma*sqrt(M_PI))) - values[i]*values[i]/(2.0*sigma*sigma);
        }
    }
    return logP;
}

double DistributionModel_half_normal_dlogP(DistributionModel* dm, const Parameter* p){
    // derivative wrt x
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double sigma = Parameters_value(dm->parameters[0], 0);
            if(Parameters_count(dm->parameters[0]) > 1){
                sigma = Parameters_value(dm->parameters[0], i);
            }
            return -Parameters_value(dm->x, 0)/sigma/sigma;
        }
    }
    return 0;
}

double DistributionModel_half_normal_d2logP(DistributionModel* dm, const Parameter* p){
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double sigma = Parameters_value(dm->parameters[0], 0);
            if(Parameters_count(dm->parameters[0]) > 1){
                sigma = Parameters_value(dm->parameters[0], i);
            }
            return -1.0/sigma/sigma;
        }
    }
    return 0;
}

static void DistributionModel_half_normal_sample(DistributionModel* dm, double* samples){
    fprintf(stderr, "DistributionModel_half_normal_sample not yet implemented\n");
    exit(2);
}


static double DistributionModel_half_normal_sample_evaluate(DistributionModel* dm){
    fprintf(stderr, "DistributionModel_half_normal_sample not DistributionModel_half_normal_sample_evaluate implemented\n");
    exit(2);
    return 0;
}

DistributionModel* new_HalfNormalDistributionModel_with_parameters(Parameters** parameters, Parameters* x, distribution_parameterization param){
    DistributionModel* dm = new_DistributionModel(parameters, 1, x);
    dm->type = DISTRIBUTION_HALFNORMAL;
    dm->parameterization = param;
    dm->logP = DistributionModel_half_normal_logP;
    dm->logP_with_values = DistributionModel_half_normal_logP_with_values;
    dm->dlogP = DistributionModel_half_normal_dlogP;
    dm->d2logP = DistributionModel_half_normal_d2logP;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_half_normal_sample;
    dm->sample_evaluate = DistributionModel_half_normal_sample_evaluate;
    dm->shift = 0;
    return dm;
}


Model* new_HalfNormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters** parameters = NULL;
    size_t parameters_dim = 1;
    bool tau = false;
    
    char* parameterization_string = get_json_node_value_string(node, "parameterization");
    
    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        size_t paramCount = Parameters_count(x);
        parameters = malloc(sizeof(Parameters*));
        parameters[0] = new_Parameters(paramCount);
        
        if(parameterization_string!= NULL && strcasecmp(parameterization_string, "tau") == 0){
            tau = true;
        }
        
        for (int i = 0; i < paramCount; i++) {
            double* vec = Vector_data(samples[i]);
            double sum2 = 0;
            size_t n = Vector_length(samples[i]);
            for(size_t j = 0; j < n; j++){
                sum2 += vec[j]*vec[j];
            }
            double sigma = sqrt(sum2/n); // ML estimate
            if(tau){
                Parameters_move(parameters[0], new_Parameter("tau", 1.0/sigma, NULL));
            }
            else{
                Parameters_move(parameters[0], new_Parameter("sigma", sigma, NULL));
            }
            free_Vector(samples[i]);
        }
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        parameters = malloc(sizeof(Parameters*));
        parameters[0] = new_Parameters(Parameters_count(x));
        
        if(parameterization_string!= NULL && strcasecmp(parameterization_string, "tau") == 0){
            tau = true;
            for (int i = 0; i < Parameters_count(x); i++) {
                Parameters_move(parameters[0], new_Parameter("tau", 1.0, new_Constraint(0, INFINITY)));
            }
        }
        else{
            for (int i = 0; i < Parameters_count(x); i++) {
                Parameters_move(parameters[0], new_Parameter("sigma", 1.0, new_Constraint(0, INFINITY)));
            }
        }
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "tau") == 0) {
                tau = true;
            }
            else if (strcasecmp(parameters_node->children[i]->key, "sigma") != 0) {
                fprintf(stderr, "Normal distribution should be parametrized with scale or precision (sigma or tau)\n");
                exit(13);
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
    }
    
    distribution_parameterization parameterization = tau ? DISTRIBUTION_HALFNORMAL_MEAN_TAU: DISTRIBUTION_HALFNORMAL_MEAN_SIGMA;
    DistributionModel*dm = new_HalfNormalDistributionModel_with_parameters(parameters, x, parameterization);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free(parameters);
    
    return model;
}
