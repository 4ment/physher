//
//  distnormal.c
//  physher
//
//  Created by mathieu on 24/1/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distnormal.h"

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"


double DistributionModel_normal_logP(DistributionModel* dm){
    double logP = 0;
    if(Parameters_count(dm->parameters) > 2){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double mu = Parameters_value(dm->parameters, i*2);
            double sigma = Parameters_value(dm->parameters, i*2+1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double x = Parameters_value(dm->x, i);
            logP += log(gsl_ran_gaussian_pdf(x - mu, sigma));
        }
    }
    else{
        double mu = Parameters_value(dm->parameters, 0);
        double sigma = Parameters_value(dm->parameters, 1);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double x = Parameters_value(dm->x, i);
            logP += log(gsl_ran_gaussian_pdf(x - mu, sigma));
        }
    }
    return logP;
}

double DistributionModel_normal_logP_with_values(DistributionModel* dm, const double* values){
    double logP = 0;
    if(Parameters_count(dm->parameters) > 2){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double mu = Parameters_value(dm->parameters, i*2);
            double sigma = Parameters_value(dm->parameters, i*2+1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            logP += log(gsl_ran_gaussian_pdf(values[i] - mu, sigma));
        }
    }
    else{
        double mu = Parameters_value(dm->parameters, 0);
        double sigma = Parameters_value(dm->parameters, 1);
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
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters, 0);
            double sigma = Parameters_value(dm->parameters, 1);
            if(Parameters_count(dm->parameters) > 2){
                mu = Parameters_value(dm->parameters, i*2);
                sigma = Parameters_value(dm->parameters, i*2+1);
            }
            return (mu - Parameters_value(dm->x, 0))/sigma/sigma;
        }
    }
    return 0;
}

double DistributionModel_normal_d2logP(DistributionModel* dm, const Parameter* p){
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters, 0);
            double sigma = Parameters_value(dm->parameters, 1);
            if(Parameters_count(dm->parameters) > 2){
                mu = Parameters_value(dm->parameters, i*2);
                sigma = Parameters_value(dm->parameters, i*2+1);
            }
            return -1.0/sigma/sigma;
        }
    }
    return 0;
}

static void DistributionModel_normal_sample(DistributionModel* dm, double* samples){
    if(Parameters_count(dm->parameters) > 2){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters, i*2+1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            samples[i] = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters, i*2);
        }
    }
    else{
        double mean = Parameters_value(dm->parameters, 0);
        double sigma = Parameters_value(dm->parameters, 1);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            samples[i] = gsl_ran_gaussian(dm->rng, sigma) + mean;
        }
    }
}


static double DistributionModel_normal_sample_evaluate(DistributionModel* dm){
    if(Parameters_count(dm->parameters) > 2){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters, i*2+1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double sample = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters, i*2);
            Parameters_set_value(dm->x, i, sample);
        }
    }
    else{
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double mean = Parameters_value(dm->parameters, 0);
            double sigma = Parameters_value(dm->parameters, 1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double sample = gsl_ran_gaussian(dm->rng, sigma) + mean;
            Parameters_set_value(dm->x, i, sample);
        }
    }
    return DistributionModel_normal_logP(dm);
}

DistributionModel* new_NormalDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_NORMAL;
    dm->logP = DistributionModel_normal_logP;
    dm->logP_with_values = DistributionModel_normal_logP_with_values;
    dm->dlogP = DistributionModel_normal_dlogP;
    dm->d2logP = DistributionModel_normal_d2logP;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_normal_sample;
    dm->sample_evaluate = DistributionModel_normal_sample_evaluate;
    return dm;
}

DistributionModel* new_NormalDistributionModel(const double mean, const double sigma, const Parameters* x){
    Parameters* ps = new_Parameters(2);
    Parameters_move(ps, new_Parameter("normal.mean", mean, NULL));
    Parameters_move(ps, new_Parameter("normal.sigma", sigma, NULL));
    DistributionModel* dm = new_NormalDistributionModel_with_parameters(ps, x);
    free_Parameters(ps);
    return dm;
}

Model* new_NormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    Parameters* x = new_Parameters(1);
    get_parameters_references2(node, hash, x, "x");
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    bool tau = false;
    // empirical
    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        size_t paramCount = Parameters_count(x);
        
        for (int i = 0; i < paramCount; i++) {
            double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            Parameters_move(parameters, new_Parameter("mu", m, NULL));
            Parameters_move(parameters, new_Parameter("sigma", sqrt(v), NULL));
        }
    }
    else if(get_json_node(node, "parameters") == NULL){
        get_parameters_references2(node, hash, x, "x");
        for (int i = 0; i < Parameters_count(x); i++) {
            Parameters_move(parameters, new_Parameter("mu", 0, new_Constraint(-INFINITY, INFINITY)));
            Parameters_move(parameters, new_Parameter("sigma", 1, new_Constraint(0, INFINITY)));
        }
    }
    else{
        get_parameters_references2(node, hash, x, "x");
        json_node* x_node = get_json_node(node, "parameters");
        for (int i = 0; i < x_node->child_count; i++) {
            if (strcasecmp(x_node->children[i]->key, "tau") == 0) {
                tau = true;
            }
            else if (strcasecmp(x_node->children[i]->key, "mean") != 0 && strcasecmp(x_node->children[i]->key, "sigma") != 0) {
                fprintf(stderr, "Normal distribution should be parametrized with mean and (sigma ot tau)\n");
                exit(13);
            }
        }
        get_parameters_references(node, hash, parameters);
        if (strcasecmp(x_node->children[0]->key, "mean") != 0) {
            Parameters_swap_index(parameters, 0, 1);
        }
        for (int i = 0; i < Parameters_count(parameters); i++) {
            Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
        }
    }
    DistributionModel*dm = new_NormalDistributionModel_with_parameters(parameters, x);
    
    dm->parameterization = tau ? DISTRIBUTION_NORMAL_MEAN_TAU: DISTRIBUTION_NORMAL_MEAN_SIGMA;
    dm->shift = -INFINITY; // can't be shifted
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters);
    free_Parameters(x);
    
    return model;
}

#pragma mark Halfnormal

double DistributionModel_half_normal_logP(DistributionModel* dm){
    double logP = 0;
    if(Parameters_count(dm->parameters) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters, i);
            if (dm->parameterization == DISTRIBUTION_HALFNORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double x = Parameters_value(dm->x, i);
            logP += log(sqrt(2)/(sigma*sqrt(M_PI))) - x*x/(2.0*sigma*sigma);
        }
    }
    else{
        double sigma = Parameters_value(dm->parameters, 0);
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
    if(Parameters_count(dm->parameters) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters, i);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            logP += log(sqrt(2)/(sigma*sqrt(M_PI))) - values[i]*values[i]/(2.0*sigma*sigma);
        }
    }
    else{
        double sigma = Parameters_value(dm->parameters, 0);
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
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double sigma = Parameters_value(dm->parameters, 0);
            if(Parameters_count(dm->parameters) > 1){
                sigma = Parameters_value(dm->parameters, i);
            }
            return -Parameters_value(dm->x, 0)/sigma/sigma;
        }
    }
    return 0;
}

double DistributionModel_half_normal_d2logP(DistributionModel* dm, const Parameter* p){
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double sigma = Parameters_value(dm->parameters, 0);
            if(Parameters_count(dm->parameters) > 1){
                sigma = Parameters_value(dm->parameters, i);
            }
            return -1.0/sigma/sigma;
        }
    }
    return 0;
}

static void DistributionModel_half_normal_sample(DistributionModel* dm, double* samples){
    fprintf(stderr, "DistributionModel_half_normal_sample not yet implemented\n");
    exit(2);
    if(Parameters_count(dm->parameters) > 2){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters, i*2+1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            samples[i] = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters, i*2);
        }
    }
    else{
        double mean = Parameters_value(dm->parameters, 0);
        double sigma = Parameters_value(dm->parameters, 1);
        if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
            sigma = sqrt(1.0/sigma);
        }
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            samples[i] = gsl_ran_gaussian(dm->rng, sigma) + mean;
        }
    }
}


static double DistributionModel_half_normal_sample_evaluate(DistributionModel* dm){
    fprintf(stderr, "DistributionModel_half_normal_sample not DistributionModel_half_normal_sample_evaluate implemented\n");
    exit(2);
    if(Parameters_count(dm->parameters) > 2){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sigma = Parameters_value(dm->parameters, i*2+1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double sample = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters, i*2);
            Parameters_set_value(dm->x, i, sample);
        }
    }
    else{
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double mean = Parameters_value(dm->parameters, 0);
            double sigma = Parameters_value(dm->parameters, 1);
            if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
                sigma = sqrt(1.0/sigma);
            }
            double sample = gsl_ran_gaussian(dm->rng, sigma) + mean;
            Parameters_set_value(dm->x, i, sample);
        }
    }
    return DistributionModel_normal_logP(dm);
}

DistributionModel* new_HalfNormalDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_HALFNORMAL;
    dm->logP = DistributionModel_half_normal_logP;
    dm->logP_with_values = DistributionModel_half_normal_logP_with_values;
    dm->dlogP = DistributionModel_half_normal_dlogP;
    dm->d2logP = DistributionModel_half_normal_d2logP;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_half_normal_sample;
    dm->sample_evaluate = DistributionModel_half_normal_sample_evaluate;
    return dm;
}

DistributionModel* new_HalfNormalDistributionModel(const double mean, const double sigma, const Parameters* x){
    Parameters* ps = new_Parameters(1);
    Parameters_move(ps, new_Parameter("halfnormal.sigma", sigma, NULL));
    DistributionModel* dm = new_HalfNormalDistributionModel_with_parameters(ps, x);
    free_Parameters(ps);
    return dm;
}

Model* new_HalfNormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    Parameters* x = new_Parameters(1);
    get_parameters_references2(node, hash, x, "x");
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(1);
    bool tau = false;
    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        size_t paramCount = Parameters_count(x);
        
        for (int i = 0; i < paramCount; i++) {
            double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            Parameters_move(parameters, new_Parameter("sigma", sqrt(v), NULL));
        }
    }
    else if(get_json_node(node, "parameters") == NULL){
        get_parameters_references2(node, hash, x, "x");
        for (int i = 0; i < Parameters_count(x); i++) {
            Parameters_move(parameters, new_Parameter("sigma", 1.0, new_Constraint(0, INFINITY)));
        }
    }
    else{
        get_parameters_references2(node, hash, x, "x");
        json_node* x_node = get_json_node(node, "parameters");
        for (int i = 0; i < x_node->child_count; i++) {
            if (strcasecmp(x_node->children[i]->key, "tau") == 0) {
                tau = true;
            }
            else if (strcasecmp(x_node->children[i]->key, "sigma") != 0) {
                fprintf(stderr, "Normal distribution should be parametrized with scale or precision (sigma or tau)\n");
                exit(13);
            }
        }
        get_parameters_references(node, hash, parameters);
        Hashtable_add(hash, Parameters_name(parameters, 0), Parameters_at(parameters, 0));
    }
    DistributionModel*dm = new_HalfNormalDistributionModel_with_parameters(parameters, x);
    
    dm->parameterization = tau ? DISTRIBUTION_HALFNORMAL_MEAN_TAU: DISTRIBUTION_HALFNORMAL_MEAN_SIGMA;
    dm->shift = 0;
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters);
    free_Parameters(x);
    
    return model;
}
