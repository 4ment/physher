//
//  distmultinormal.c
//  physher
//
//  Created by mathieu on 17/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "distmultinormal.h"

#include <strings.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"


void _update_gsl_parameters(DistributionModel* dm){
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* cov = Parameters_at(dm->parameters, 1);
    Parameter* x = Parameters_at(dm->x, 0);
    size_t dim = Parameter_size(x);
    size_t row = 0;
    for (size_t i = 0; i < dim; i++) {
        gsl_vector_set(wrapper->x, i,  Parameter_value_at(x, i));
        gsl_vector_set(wrapper->mu, i, Parameter_value_at(mu, i));
        for (size_t j = 0; j <= i; j++) {
            gsl_matrix_set(wrapper->L, i, j, Parameter_value_at(cov, row));
            row++;
        }
    }
    dm->need_update = false;
}

static double DistributionModel_multivariate_normal_log_prob(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    Parameter* x = Parameters_at(dm->x, 0);
    size_t dim = Parameter_size(x);
    _update_gsl_parameters(dm);
    
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    
    double logJocobian = 0;
    const double* values = Parameter_values(x);
    if(wrapper->transform){
        for (size_t j = 0; j < dim; j++) {
            double y = log(values[j]);
            logJocobian -= y;
            gsl_vector_set(wrapper->x, j,  y);
        }
    }
    double logP;
    dm->lp = logP + logJocobian;
    gsl_ran_multivariate_gaussian_log_pdf (wrapper->x,  wrapper->mu, wrapper->L, &logP, wrapper->work);
    dm->need_update = false;
    return dm->lp;
}

// if dst is NULL we assign directly the sampled values to dm->x
static void DistributionModel_multivariate_normal_sample(DistributionModel* dm){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    gsl_ran_multivariate_gaussian(wrapper->rng, wrapper->mu, wrapper->L, wrapper->x);
    
    Parameter* x = Parameters_at(dm->x, 0);
    size_t dim = Parameter_size(x);

    if(wrapper->transform){
        for (size_t i = 0; i < dim; i++) {
            dm->tempx[i] = exp(gsl_vector_get(wrapper->x, i));
        }
    }
    else{
        for (size_t i = 0; i < dim; i++) {
            dm->tempx[i] = gsl_vector_get(wrapper->x, i);
        }
    }
    Parameter_set_values(x, dm->tempx);
}

static void _free_dist_gsl_multivariate_normal(DistributionModel*dm){
    free_Parameters(dm->x);
    free_Parameters(dm->parameters);
    if(dm->tempx != NULL) free(dm->tempx);
    if(dm->tempp != NULL) free(dm->tempp);
    // if(dm->simplex != NULL) free_Simplex(dm->simplex);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    gsl_vector_free(wrapper->mu);
    gsl_vector_free(wrapper->x);
    gsl_vector_free(wrapper->work);
    gsl_matrix_free(wrapper->L);
    free(wrapper);
    free(dm);
}

DistributionModel* new_MultivariateNormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, x);
    dm->type = DISTRIBUTION_NORMAL_MULTIVARIATE;
    dm->free = _free_dist_gsl_multivariate_normal;
    dm->log_prob = DistributionModel_multivariate_normal_log_prob;
    dm->sample = DistributionModel_multivariate_normal_sample;
    dm->need_update = true;
    
    size_t dim = Parameter_size(Parameters_at(x, 0));
    gsl_multivariate_normal_wrapper_t* wrapper = (gsl_multivariate_normal_wrapper_t*)malloc(sizeof(gsl_multivariate_normal_wrapper_t));
    wrapper->mu = gsl_vector_calloc(dim);
    wrapper->L = gsl_matrix_calloc(dim, dim);
    wrapper->x = gsl_vector_calloc(dim);
    wrapper->work = gsl_vector_calloc(dim);
//
//    for (int i = 0; i < dim; i++) {
//        gsl_vector_set(wrapper->mu, i, Parameters_value(parameters[0], i));
//        for (int j = 0; j < dim; j++) {
//            gsl_matrix_set(wrapper->L, i, j, sigma[i*dim+j]);
//        }
//    }
//    gsl_linalg_cholesky_decomp1(wrapper->L);
    dm->data = wrapper;
    return dm;
}

static Model* _dist_model_clone_mvn( Model *self, Hashtable* hash ){
    Model* clone = self->clone(self, hash);
    DistributionModel* clonedm = clone->obj;
    DistributionModel* dm = self->obj;
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* x = Parameters_at(dm->x, 0);
    size_t dim = Parameter_size(x);
    gsl_multivariate_normal_wrapper_t* clonewrapper = (gsl_multivariate_normal_wrapper_t*)malloc(sizeof(gsl_multivariate_normal_wrapper_t));
    clonewrapper->mu = gsl_vector_calloc(dim);
    clonewrapper->L = gsl_matrix_calloc(dim, dim);
    clonewrapper->x = gsl_vector_calloc(dim);
    clonewrapper->work = gsl_vector_calloc(dim);
    clonewrapper->transform = wrapper->transform;
    clonewrapper->rng = wrapper->rng;
    gsl_vector_memcpy(clonewrapper->mu, wrapper->mu);
    gsl_vector_memcpy(clonewrapper->x, wrapper->x);
    gsl_matrix_memcpy(clonewrapper->L, wrapper->L);
    clonedm->data = clonewrapper;
    return clone;
}

Model* new_MultivariateNormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = new_Parameters(1);
    Parameters_add(x, new_Parameter_from_json(x_node, hash));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));;
    Parameter* mu = NULL;
    Parameter* sigma = NULL;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);

        double* muValues = dvector(paramCount);
        double* sigmaValues = dvector(paramCount*paramCount);
        int n = Vector_length(samples[0]);
        
        for (size_t i = 0; i < paramCount; i++) {
            double* vv = Vector_mutable_data(samples[i]);
            for (int j = 0; j < n; j++) {
                vv[j] = log(vv[j]);
            }
            muValues[i] = mean(vv, n);
        }
        
        // Calculate covariance matrix
        for (int i = 0; i < paramCount; i++) {
            const double* pp = Vector_data(samples[i]);
            sigmaValues[i*paramCount+i] = variance(pp, n, muValues[i]);
            for (int j = i+1; j < paramCount; j++) {
                const double* pp2 = Vector_data(samples[j]);
                sigmaValues[i*paramCount+j] = sigmaValues[j*paramCount+i] = covariance(pp, pp2, muValues[i], muValues[j], n);
            }
            free_Vector(samples[i]);
        }
        mu = new_Parameter2("mu", muValues, paramCount, new_Constraint(-INFINITY, INFINITY));
        sigma = new_Parameter2("sigma", sigmaValues, paramCount*paramCount, new_Constraint(0, INFINITY));

        free(muValues);
        free(sigmaValues);
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* muValues = malloc(sizeof(double)*paramCount);
        double* sigmaValues = malloc(sizeof(double)*paramCount*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            muValues[i] = 0;
            for (size_t j = 0; j < paramCount; j++) {
                sigmaValues[i*paramCount+j] = 1;
            }
        }

        mu = new_Parameter2("mu", muValues, paramCount, new_Constraint(-INFINITY, INFINITY));
        sigma = new_Parameter2("sigma", sigmaValues, paramCount*paramCount, new_Constraint(0, INFINITY));

        free(muValues);
        free(sigmaValues);
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "mu") != 0 &&
                strcasecmp(parameters_node->children[i]->key, "sigma") != 0 &&
                strcasecmp(parameters_node->children[i]->key, "L") != 0) {
                fprintf(stderr, "Multivariate Normal distribution should be parameterized with mean and (sigma or L)\n");
                exit(13);
            }
        }
        
        json_node* mu_node = get_json_node(parameters_node, "mu");
        json_node* sigma_node = get_json_node(parameters_node, "sigma");
        json_node* L_node = get_json_node(parameters_node, "L");
        mu = new_Parameter_from_json(mu_node, hash);
        if(sigma_node != NULL){
            sigma = new_Parameter_from_json(sigma_node, hash);
            size_t p_count = Parameter_size(mu);
            if(Parameter_size(sigma) != p_count*p_count){
                fprintf(stderr, "MultivariateNormal distribution - Dimension of covariance matrix (dim: %zu) does not match mu (dim: %zu). Expected: %zu\n", Parameter_size(sigma), p_count, p_count*p_count);
            }
        }
        else {
            sigma = new_Parameter_from_json(L_node, hash);
            size_t p_count = Parameter_size(mu);
            if(Parameter_size(sigma) != p_count*(p_count-1)/2 + p_count){
                fprintf(stderr, "MultivariateNormal distribution - Dimension of cholesky decomposition (dim: %zu) does not match mu (dim: %zu). Expected: %zu\n", Parameter_size(sigma), p_count, p_count*(p_count-1)/2 + p_count);
            }
        }
    }

    Parameters_add(parameters, mu);
    Parameters_add(parameters, sigma);
    
    DistributionModel* dm = new_MultivariateNormalDistributionModel_with_parameters(parameters, x);
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    model->clone = _dist_model_clone_mvn;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");

    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
