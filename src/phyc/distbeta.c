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

double dbeta(double x, double alpha, double beta){
    return gsl_ran_beta_pdf(x, alpha, beta);
}

double rbeta(const gsl_rng *rng, double alpha, double beta){
    return gsl_ran_beta(rng, alpha, beta);
}

double DistributionModel_log_beta(DistributionModel* dm){
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
                dm->lp += log(dbeta(values[i], *alpha, *beta));
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
                dm->lp += log(dbeta(values[i], alpha[index], beta[index]));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_gradient2_beta(DistributionModel* dm, const Parameters* p){
    error("DistributionModel_gradient2_beta not implemented\n");
    return 0.0;
}

double DistributionModel_dlog_beta(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_dlog_beta not implemented\n");
    return 0;
}

double DistributionModel_d2log_beta(DistributionModel* dm, const Parameter* p){
    error("DistributionModel_d2log_beta not implemented\n");
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
                dm->tempx[i] = rbeta(dm->rng, *alpha, *beta);
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
                dm->tempx[i] = rbeta(dm->rng, alpha[index], beta[index]);
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

// static void DistributionModel_beta_sample(DistributionModel* dm, double* samples){
//     const double* alpha = Parameter_values(Parameters_at(dm->parameters, 0));
//     const double* beta = Parameter_values(Parameters_at(dm->parameters, 1));
//     size_t dim = Parameter_size(dm->x);

//     double* out = samples;
//     if(samples == NULL){
//         out = dm->x->value;
//     }

//     if(Parameter_size(Parameters_at(dm->parameters, 0)) > 1){
//         for (size_t i = 0; i < dim; i++) {
//             out[i] = rbeta(dm->rng, alpha[i], beta[i]);
//         }
//     }
//     else{
//         for (size_t i = 0; i < dim; i++) {
//             out[i] = rbeta(dm->rng, *alpha, *beta);
//         }
//     }

//     if(samples == NULL){
//         Parameter_fire(dm->x);
//     }
// }

// static double DistributionModel_beta_sample_evaluate(DistributionModel* dm){
//     DistributionModel_beta_sample(dm, NULL);
//     return DistributionModel_log_beta(dm);
// }

DistributionModel* new_BetaDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_BETA;
    dm->logP = DistributionModel_log_beta;
    // dm->logP_with_values = DistributionModel_log_beta_with_values;
    dm->dlogP = DistributionModel_dlog_beta;
    dm->gradient2 = DistributionModel_gradient2_beta;
    dm->d2logP = DistributionModel_d2log_beta;
    dm->ddlogP = DistributionModel_ddlog_0;
    dm->sample = DistributionModel_beta_sample;
    // dm->sample_evaluate = DistributionModel_beta_sample_evaluate;
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
