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

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"

double DistributionModel_log_beta(DistributionModel* dm){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double alpha = Parameters_value(dm->parameters[0], i);
            double beta = Parameters_value(dm->parameters[1], i);
            double x = Parameters_value(dm->x, i);
            logP += log(gsl_ran_beta_pdf(x, alpha, beta));
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double x = Parameters_value(dm->x, i);
            logP += log(gsl_ran_beta_pdf(x, alpha, beta));
        }
    }
    return logP;
}

double DistributionModel_log_beta_with_values(DistributionModel* dm, const double* values){
    double logP = 0;
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double alpha = Parameters_value(dm->parameters[0], i);
            double beta = Parameters_value(dm->parameters[1], i);
            logP += log(gsl_ran_beta_pdf(values[i], alpha, beta));
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            logP += log(gsl_ran_beta_pdf(values[i], alpha, beta));
        }
    }
    return logP;
}

double DistributionModel_dlog_beta(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_dlog_beta not implemented\n");
    return 0;
}

double DistributionModel_d2log_beta(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_d2log_beta not implemented\n");
    return 0;
}

static void DistributionModel_beta_sample(DistributionModel* dm, double* samples){
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            samples[i] = gsl_ran_beta(dm->rng, Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i));
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            samples[i] = gsl_ran_beta(dm->rng, alpha, beta);
        }
    }
}

static double DistributionModel_beta_sample_evaluate(DistributionModel* dm){
    if(Parameters_count(dm->parameters[0]) > 1){
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sample = gsl_ran_beta(dm->rng, Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i));
            Parameters_set_value(dm->x, i, sample);
        }
    }
    else{
        double alpha = Parameters_value(dm->parameters[0], 0);
        double beta = Parameters_value(dm->parameters[1], 0);
        for (int i = 0; i < Parameters_count(dm->x); i++) {
            double sample = gsl_ran_beta(dm->rng, alpha, beta);
            Parameters_set_value(dm->x, i, sample);
        }
    }
    return DistributionModel_log_beta(dm);
}

DistributionModel* new_BetaDistributionModel_with_parameters(Parameters** parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, 2, x);
    dm->type = DISTRIBUTION_BETA;
    dm->logP = DistributionModel_log_beta;
    dm->logP_with_values = DistributionModel_log_beta_with_values;
    dm->dlogP = DistributionModel_dlog_beta;
    dm->d2logP = DistributionModel_d2log_beta;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_beta_sample;
    dm->sample_evaluate = DistributionModel_beta_sample_evaluate;
    dm->shift = 0;
    return dm;
}

Model* new_BetaDistributionModel_from_json(json_node* node, Hashtable* hash){
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
            double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            double alpha = m*(m*(1.0 - m)/v - 1.0);
            double beta = (1.0 - m)*(m*(1.0 - m)/v - 1.0);
            Parameters_move(parameters[0], new_Parameter("alpha", alpha, new_Constraint(0, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("beta", beta, new_Constraint(0, INFINITY)));
            free_Vector(samples[i]);
        }
        free(samples);
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
    
    DistributionModel* dm = new_BetaDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
    
    return model;
}
