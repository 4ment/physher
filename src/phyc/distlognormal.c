//
//  distlognormal.c
//  physher
//
//  Created by mathieu on 17/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distlognormal.h"

#include <strings.h>

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"

double DistributionModel_lognormal_logP(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        size_t dim = Parameters_count(dm->x);
        for (int i = 0; i < dim; i++) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            double x = Parameters_value(dm->x, i);
            dm->lp += log(gsl_ran_lognormal_pdf(x, mu, sigma));
        }
    }
    else{
        double mu = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double x = Parameters_value(dm->x, i);
            dm->lp += log(gsl_ran_lognormal_pdf(x, mu, sigma));
        }
    }
    return dm->lp;
}

double DistributionModel_lognormal_logP_with_values(DistributionModel* dm, const double* values){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            logP += log(gsl_ran_lognormal_pdf(values[i], mu, sigma));
        }
    }
    else{
        double mu = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            logP += log(gsl_ran_lognormal_pdf(values[i], mu, sigma));
        }
    }
    return logP;
}

// multiple Xs one distribution
double DistributionModel_lognormal_dlogP(DistributionModel* dm, const Parameter* p){
    // multiple Xs one distribution
    // derivative wrt mu
    if (p == Parameters_at(dm->parameters[0], 0)) {
        double mu = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        double dlogf = 0;
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            dlogf += -(mu - log(Parameters_value(dm->x, i)))/(sigma*sigma);
        }
        return dlogf;
    }
    // derivative wrt sigma
    else if (p == Parameters_at(dm->parameters[1], 0)) {
        double mu = Parameters_value(dm->parameters[0], 0);
        double sigma = Parameters_value(dm->parameters[1], 0);
        double dlogf = 0;
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            dlogf += (pow(mu - log(Parameters_value(dm->x, i)), 2.0) - sigma*sigma)/(sigma*sigma*sigma);
        }
        return dlogf;
    }
    
    // derivative wrt x
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters[0], 0);
            double sigma = Parameters_value(dm->parameters[1], 0);
            double x = Parameter_value(p);
            return -1.0/x - (log(x) - mu)/(sigma*sigma*x);
        }
    }
    return 0;
}

// multiple Xs multiple distribution
double DistributionModel_lognormal_dlogP_multi(DistributionModel* dm, const Parameter* p){
    // derivative wrt mu
    size_t pdim = Parameters_count(dm->parameters[0]);
    for(size_t i = 0; i < pdim; i++){
        if (p == Parameters_at(dm->parameters[0], i)) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            return -(mu - log(Parameters_value(dm->x, i)))/(sigma*sigma);
        }
        // derivative wrt sigma
        else if (p == Parameters_at(dm->parameters[1], i)) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            return (pow(mu - log(Parameters_value(dm->x, i)), 2.0) - sigma*sigma)/(sigma*sigma*sigma);
        }
    }
    
    // derivative wrt x
    size_t xdim = Parameters_count(dm->x);
    for (int i = 0; i < xdim; i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            double x = Parameter_value(p);
            return -1.0/x - (log(x) - mu)/(sigma*sigma*x);
        }
    }
    return 0;
}


double DistributionModel_lognormal_d2logP(DistributionModel* dm, const Parameter* p){
    // TODO: derivative wrt mu and sigma
    if(Parameters_count(dm->parameters[0]) > 1){
        fprintf(stderr, "derivative wrt mu and sigma not implemented DistributionModel_lognormal_d2logP\n");
    }
    
    // derivative wrt x
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters[0], 0);
            double sigma = Parameters_value(dm->parameters[1], 0);
            double x = Parameter_value(p);
            return -1.0/(x*x) + (log(x) - mu)/(sigma*sigma*x*x);
        }
    }
    return 0;
}

// multiple Xs multiple distribution
double DistributionModel_lognormal_d2logP_multi(DistributionModel* dm, const Parameter* p){
    // TODO: derivative wrt mu and sigma
    if(Parameters_count(dm->parameters[0]) > 1){
        fprintf(stderr, "derivative wrt mu and sigma not implemented DistributionModel_lognormal_d2logP_multi\n");
    }
    
    // derivative wrt x
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x,i)) {
            double mu = Parameters_value(dm->parameters[0], i);
            double sigma = Parameters_value(dm->parameters[1], i);
            double x = Parameter_value(p);
            return -1.0/(x*x) + (log(x) - mu)/(sigma*sigma*x*x);
        }
    }
    return 0;
}

static void DistributionModel_lognormal_sample(DistributionModel* dm, double* samples){
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            samples[i] = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i));
        }
    }
    else{
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            samples[i] = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters[0], 0), Parameters_value(dm->parameters[1], 0));
        }
    }
}


static double DistributionModel_lognormal_sample_evaluate(DistributionModel* dm){
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sample = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i));
            Parameters_set_value(dm->x, i, sample);
        }
    }
    else{
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sample = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters[0], 0), Parameters_value(dm->parameters[1], 0));
            Parameters_set_value(dm->x, i, sample);
        }
    }
    return DistributionModel_lognormal_logP(dm);
}

DistributionModel* new_LogNormalDistributionModel_with_parameters(Parameters** parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, 2, x);
    dm->type = DISTRIBUTION_LOGNORMAL;
    dm->logP = DistributionModel_lognormal_logP;
    dm->logP_with_values = DistributionModel_lognormal_logP_with_values;
    dm->dlogP = DistributionModel_lognormal_dlogP;
    dm->d2logP = DistributionModel_lognormal_d2logP;
    dm->ddlogP = DistributionModel_ddlog_0;
    if(Parameters_count(parameters[0]) > 1){
        dm->dlogP = DistributionModel_lognormal_dlogP_multi;
        dm->d2logP = DistributionModel_lognormal_d2logP_multi;
    }
    dm->sample = DistributionModel_lognormal_sample;
    dm->sample_evaluate = DistributionModel_lognormal_sample_evaluate;
    dm->shift = 0;
    return dm;
}

Model* new_LogNormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters** parameters = NULL;
    size_t parameters_dim = 2;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        size_t paramCount = Parameters_count(x);
        parameters = malloc(sizeof(Parameters*)*2);
        parameters[0] = new_Parameters(paramCount);
        parameters[1] = new_Parameters(paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            //TODO: that coming from normal
            const double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            Parameters_move(parameters[0], new_Parameter("mu", m, new_Constraint(-INFINITY, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("sigma", sqrt(v), new_Constraint(0, INFINITY)));
            free_Vector(samples[i]);
        }
        free(samples);
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
            if (strcasecmp(parameters_node->children[i]->key, "mu") != 0 && strcasecmp(parameters_node->children[i]->key, "sigma") != 0) {
                fprintf(stderr, "LogNormal distribution should be parametrized with mean and (sigma ot tau)\n");
                exit(13);
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
        
        if (strcasecmp(Parameters_name2(parameters[0]), "mu") != 0) {
            Parameters* temp = parameters[0];
            parameters[0] = parameters[1];
            parameters[1] = temp;
        }
    }
    
    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
//    free_Parameters(x);
    
    return model;
}
