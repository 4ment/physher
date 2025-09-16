//
//  distbetaprime.c
//  physher
//
//  Created by mathieu on 29/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distbetaprime.h"

#include <strings.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"


// log of density
double ldbetaprime(double x, double alpha, double beta){
    return (alpha - 1.0)*log(x) - (alpha + beta)*log(1.0 + x) - gsl_sf_lnbeta(alpha, beta);
    // return gsl_ran_beta_pdf(x/(1.0 + x), alpha, beta)/(1.0 + x);
}

double rbetaprime(const gsl_rng *rng, double alpha, double beta){
    double sample = gsl_ran_beta(rng, alpha, beta);
    return sample/(1.0 - sample);
}

double DistributionModel_log_betaprime(DistributionModel* dm){
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
                // dm->lp += log(dbeta(values[i], *alpha, *beta));
                dm->lp += ldbetaprime(values[i], *alpha, *beta);
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
                // dm->lp += log(dbeta(values[i], alpha[index], beta[index]));
                dm->lp += ldbetaprime(values[i], alpha[index], beta[index]);
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}


double DistributionModel_gradient2_betaprime(DistributionModel* dm, const Parameters* p){
    error("DistributionModel_gradient2_betaprime not implemented\n");
    return 0.0;
}

double DistributionModel_dlog_betaprime(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_dlog_betaprime not implemented\n");
    return 0;
}

double DistributionModel_d2log_betaprime(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_d2log_betaprime not implemented\n");
    return 0;
}

static void DistributionModel_betaprime_sample(DistributionModel* dm){
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
                dm->tempx[i] = rbetaprime(dm->rng, *alpha, *beta);
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
                dm->tempx[i] = rbetaprime(dm->rng, alpha[index], beta[index]);
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

DistributionModel* new_BetaPrimeDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_BETA_PRIME;
    dm->logP = DistributionModel_log_betaprime;
    // dm->logP_with_values = DistributionModel_log_betaprime_with_values;
    dm->gradient2 = DistributionModel_gradient2_betaprime;
    dm->dlogP = DistributionModel_dlog_betaprime;
    dm->d2logP = DistributionModel_d2log_betaprime;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_betaprime_sample;
    // dm->sample_evaluate = DistributionModel_betaprime_sample_evaluate;
    dm->shift = 0;
    return dm;
}

Model* new_BetaPrimeDistributionModel_from_json(json_node* node, Hashtable* hash){
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
        fprintf(stderr, "Empirical estimation of beta prime parameters not implemented yet\n");
        exit(1);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* alphaValues = malloc(sizeof(double)*paramCount);
        double* betaValues = malloc(sizeof(double)*paramCount);
        for (int i = 0; i < paramCount; i++) {
            alphaValues[i] = 0;
            betaValues[i] = 1;
        }
        alpha = new_Parameter2("alpha", alphaValues, paramCount, new_Constraint(0, INFINITY));
        beta = new_Parameter2("beta", betaValues, paramCount, new_Constraint(0, INFINITY));
    }
    else{
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
    
    DistributionModel* dm = new_BetaPrimeDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
