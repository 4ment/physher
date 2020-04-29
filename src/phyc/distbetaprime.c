//
//  distbetaprime.c
//  physher
//
//  Created by mathieu on 29/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distbetaprime.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"

double DistributionModel_log_betaprime(DistributionModel* dm){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double alpha = Parameters_value(dm->parameters[0], i);
            double beta = Parameters_value(dm->parameters[1], i);
            double x = Parameters_value(dm->x, i);
            logP += (alpha - 1.0)*log(x) - (alpha + beta)*log(1.0 + x) - gsl_sf_lnbeta(alpha, beta);
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double x = Parameters_value(dm->x, i);
            logP += (alpha - 1.0)*log(x) - (alpha + beta)*log(1.0 + x) - gsl_sf_lnbeta(alpha, beta);
        }
    }
    return logP;
}

double DistributionModel_log_betaprime_with_values(DistributionModel* dm, const double* values){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double alpha = Parameters_value(dm->parameters[0], i);
            double beta = Parameters_value(dm->parameters[1], i);
            logP += (alpha - 1.0)*log(values[i]) - (alpha + beta)*log(1.0 + values[i]) - gsl_sf_lnbeta(alpha, beta);
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            logP += (alpha - 1.0)*log(values[i]) - (alpha + beta)*log(1.0 + values[i]) - gsl_sf_lnbeta(alpha, beta);
        }
    }
    return logP;
}

double DistributionModel_dlog_betaprime(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_dlog_betaprime not implemented\n");
    return 0;
}

double DistributionModel_d2log_betaprime(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_d2log_betaprime not implemented\n");
    return 0;
}

static void DistributionModel_betaprime_sample(DistributionModel* dm, double* samples){
    double sample;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            sample = gsl_ran_beta(dm->rng, Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i));
            samples[i] = sample/(1.0 - sample);
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            sample = gsl_ran_beta(dm->rng, alpha, beta);
            samples[i] = sample/(1.0 - sample);
        }
    }
}

static double DistributionModel_betaprime_sample_evaluate(DistributionModel* dm){
    double sample;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            sample = gsl_ran_beta(dm->rng, Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i));
            Parameters_set_value(dm->x, i, sample/(1.0 - sample));
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            sample = gsl_ran_beta(dm->rng, alpha, beta);
            Parameters_set_value(dm->x, i, sample/(1.0 - sample));
        }
    }
    return DistributionModel_log_betaprime(dm);
}

DistributionModel* new_BetaPrimeDistributionModel_with_parameters(Parameters** parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, 2, x);
    dm->type = DISTRIBUTION_BETA_PRIME;
    dm->logP = DistributionModel_log_betaprime;
    dm->logP_with_values = DistributionModel_log_betaprime_with_values;
    dm->dlogP = DistributionModel_dlog_betaprime;
    dm->d2logP = DistributionModel_d2log_betaprime;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_betaprime_sample;
    dm->sample_evaluate = DistributionModel_betaprime_sample_evaluate;
    dm->shift = 0;
    return dm;
}

Model* new_BetaPrimeDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters** parameters = NULL;
    size_t parameters_dim = 2;

    // empirical
    if (file != NULL) {
        fprintf(stderr, "Empirical estimation of beta prime parameters not implemented yet\n");
        exit(1);
    }
    else if(get_json_node(node, "parameters") == NULL){
        parameters = malloc(sizeof(Parameters*)*2);
        parameters[0] = new_Parameters(Parameters_count(x));
        parameters[1] = new_Parameters(Parameters_count(x));
        for (int i = 0; i < Parameters_count(x); i++) {
            Parameters_move(parameters[0], new_Parameter("alpha", 0, new_Constraint(0, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("beta", 1, new_Constraint(0, INFINITY)));
        }
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "alpha") != 0 && strcasecmp(parameters_node->children[i]->key, "beta") != 0) {
                fprintf(stderr, "Beta distribution should be parametrized with alpha and beta\n");
                exit(13);
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
        
        if (strcasecmp(Parameters_name2(parameters[0]), "alpha") != 0) {
            Parameters* temp = parameters[0];
            parameters[0] = parameters[1];
            parameters[1] = temp;
        }
    }
    
    DistributionModel* dm = new_BetaPrimeDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
    
    return model;
}
