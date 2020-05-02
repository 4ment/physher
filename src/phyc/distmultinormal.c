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
    size_t dim = Parameters_count(dm->x);
    size_t row = 0;
    for (int i = 0; i < dim; i++) {
        gsl_vector_set(wrapper->x, i,  Parameters_value(dm->x, i));
        gsl_vector_set(wrapper->mu, i, Parameters_value(dm->parameters[0], i));
        for (int j = 0; j <= i; j++) {
            gsl_matrix_set(wrapper->L, i, j, Parameters_value(dm->parameters[1], row));
            row++;
        }
    }
    dm->need_update = false;
}

static double _multivariate_normal_logP(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    
    _update_gsl_parameters(dm);
    
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    
    double logJocobian = 0;
    if(wrapper->transform){
        for (int j = 0; j < Parameters_count(dm->x); j++) {
            double y = log(Parameters_value(dm->x, j));
            logJocobian -= log(Parameters_value(dm->x, j));
            gsl_vector_set(wrapper->x, j,  y);
        }
    }
    double logP;
    gsl_ran_multivariate_gaussian_log_pdf (wrapper->x,  wrapper->mu, wrapper->L, &logP, wrapper->work);

    return logP + logJocobian;
}

static double _multivariate_normal_logP_with_values(DistributionModel* dm, const double* values){
    if(!dm->need_update) return dm->lp;
    
    _update_gsl_parameters(dm);
    
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    
    double logJocobian = 0;
    if(wrapper->transform){
        for (int j = 0; j < Parameters_count(dm->x); j++) {
            double y = log(values[j]);
            logJocobian -= y;
            gsl_vector_set(wrapper->x, j,  y);
        }
    }
    double logP;
    gsl_ran_multivariate_gaussian_log_pdf (wrapper->x,  wrapper->mu, wrapper->L, &logP, wrapper->work);

    return logP + logJocobian;
}


static double _DistributionModel_dlog_mvn(DistributionModel* dm, const Parameter* p){
    fprintf(stderr, "%s is not implemented (line %d of file %s)\n", __func__, __LINE__, __FILE__);
    exit(1);
    return 0.0;
}

static double _DistributionModel_d2log_mvn(DistributionModel* dm, const Parameter* p){
    fprintf(stderr, "%s is not implemented (line %d of file %s)\n", __func__, __LINE__, __FILE__);
    exit(1);
    return 0.0;
}

static double _DistributionModel_ddlog_mvn(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
    fprintf(stderr, "%s is not implemented (line %d of file %s)\n", __func__, __LINE__, __FILE__);
    exit(1);
    return 0.0;
}

// if dst is NULL we assign directly the sampled values to dm->x
static void _sample_multivariate_normal(DistributionModel* dm, double* dst){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    gsl_ran_multivariate_gaussian(wrapper->rng, wrapper->mu, wrapper->L, wrapper->x);
    
    if(wrapper->transform){
        if(dst == NULL){
            for (int j = 0; j < Parameters_count(dm->x); j++) {
                Parameters_set_value(dm->x, j, exp(gsl_vector_get(wrapper->x, j)));
            }
        }
        else{
            for (int j = 0; j < Parameters_count(dm->x); j++) {
                dst[j] = exp(gsl_vector_get(wrapper->x, j));
            }
        }
    }
    else{
        if(dst == NULL){
            for (int j = 0; j < Parameters_count(dm->x); j++) {
                Parameters_set_value(dm->x, j, gsl_vector_get(wrapper->x, j));
            }
        }
        else{
            for (int j = 0; j < Parameters_count(dm->x); j++) {
                dst[j] = gsl_vector_get(wrapper->x, j);
            }
        }
    }
}

static double _multivariate_normal_sample_evaluate(DistributionModel* dm){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    gsl_ran_multivariate_gaussian(wrapper->rng, wrapper->mu, wrapper->L, wrapper->x);
    
    double logJocobian = 0;
    if(wrapper->transform){
        for (int j = 0; j < Parameters_count(dm->x); j++) {
            double y = exp(gsl_vector_get(wrapper->x, j));
            logJocobian -= gsl_vector_get(wrapper->x, j);
            Parameters_set_value(dm->x, j, y);
        }
    }
    else{
        for (int j = 0; j < Parameters_count(dm->x); j++) {
            Parameters_set_value(dm->x, j, gsl_vector_get(wrapper->x, j));
        }
    }
    double logP;
    gsl_ran_multivariate_gaussian_log_pdf (wrapper->x,  wrapper->mu, wrapper->L, &logP, wrapper->work);

    return logP + logJocobian;
}

static void _free_dist_gsl_multivariate_normal(DistributionModel*dm){
    if(dm->x != NULL) free_Parameters(dm->x);
    if(dm->parameters != NULL){
        for(size_t i = 0; i < dm->parameter_count; i++){
            free_Parameters(dm->parameters[i]);
        }
        free(dm->parameters);
    }
    if(dm->tempx != NULL) free(dm->tempx);
    if(dm->tempp != NULL) free(dm->tempp);
    if(dm->simplex != NULL) free_Simplex(dm->simplex);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    gsl_vector_free(wrapper->mu);
    gsl_vector_free(wrapper->x);
    gsl_vector_free(wrapper->work);
    gsl_matrix_free(wrapper->L);
    free(wrapper);
    free(dm);
}

DistributionModel* new_MultivariateNormalDistributionModel_with_parameters(Parameters** parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters, 2, x);
    dm->type = DISTRIBUTION_NORMAL_MULTIVARIATE;
    dm->free = _free_dist_gsl_multivariate_normal;
    dm->logP = _multivariate_normal_logP;
    dm->logP_with_values = _multivariate_normal_logP_with_values;
    dm->dlogP = _DistributionModel_dlog_mvn;
    dm->d2logP = _DistributionModel_d2log_mvn;
    dm->ddlogP = _DistributionModel_ddlog_mvn;
    dm->sample = _sample_multivariate_normal;
    dm->sample_evaluate = _multivariate_normal_sample_evaluate;
    dm->need_update = true;
    
    size_t dim = Parameters_count(x);
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
    size_t dim = Parameters_count(dm->x);
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
        parameters[1] = new_Parameters(paramCount*paramCount);
        double* mu = dvector(paramCount);
        double* sigma = dvector(paramCount*paramCount);
        int n = Vector_length(samples[0]);
        
        for (int i = 0; i < paramCount; i++) {
            double* vv = Vector_data(samples[i]);
            for (int j = 0; j < n; j++) {
                vv[j] = log(vv[j]);
            }
            mu[i] = mean(vv, n);
            Parameters_move(parameters[0], new_Parameter("mu", mu[i], new_Constraint(-INFINITY, INFINITY)));
        }
        
        // Calculate covariance matrix
        for (int i = 0; i < paramCount; i++) {
            double* pp = Vector_data(samples[i]);
            sigma[i*paramCount+i] = variance(pp, n, mu[i]);
            for (int j = i+1; j < paramCount; j++) {
                double* pp2 = Vector_data(samples[j]);
                sigma[i*paramCount+j] = sigma[j*paramCount+i] = covariance(pp, pp2, mu[i], mu[j], n);
            }
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
            if (strcasecmp(parameters_node->children[i]->key, "mu") != 0 &&
                strcasecmp(parameters_node->children[i]->key, "sigma") != 0 &&
                strcasecmp(parameters_node->children[i]->key, "L") != 0) {
                fprintf(stderr, "Multivariate Normal distribution should be parameterized with mean and (sigma or L)\n");
                exit(13);
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
        
        if (strcasecmp(parameters_node->children[0]->key, "mu") != 0) {
            Parameters* temp = parameters[0];
            parameters[0] = parameters[1];
            parameters[1] = temp;
        }
        if (strcasecmp(parameters_node->children[0]->key, "L") == 0 || strcasecmp(parameters_node->children[1]->key, "L") == 0) {
            size_t p_count = Parameters_count(parameters[0]);
            if(Parameters_count(parameters[1]) != p_count*(p_count-1)/2 + p_count){
                fprintf(stderr, "MultivariateNormal distribution - Dimension of cholesky decomposition (dim: %zu) does not match mu (dim: %zu). Expected: %zu\n", Parameters_count(parameters[1]), Parameters_count(parameters[0]), p_count*(p_count-1)/2 + p_count);
            }
        }
        if (strcasecmp(parameters_node->children[0]->key, "sigma") == 0 || strcasecmp(parameters_node->children[1]->key, "sigma") == 0) {
            size_t p_count = Parameters_count(parameters[0]);
            if(Parameters_count(parameters[1]) != p_count*p_count){
                fprintf(stderr, "MultivariateNormal distribution - Dimension of covariance matrix (dim: %zu) does not match mu (dim: %zu). Expected: %zu\n", Parameters_count(parameters[1]), Parameters_count(parameters[0]), p_count*p_count);
            }
        }
    }
    
    DistributionModel* dm = new_MultivariateNormalDistributionModel_with_parameters(parameters, x);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    model->clone = _dist_model_clone_mvn;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
//    free_Parameters(x);
    
    return model;
}
