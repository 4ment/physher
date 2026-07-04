// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "distmultinormal.h"

#include <math.h>
#include <strings.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"

#define LOG_TWO_PI (log(2.0) + log(M_PI))


// Load mu and the Cholesky factor L. These parameterise the distribution and are
// shared across every observation in dm->x; the per-observation vector is loaded
// separately by _load_gsl_x. The MVN dimension is Parameter_size(mu).
void _update_gsl_parameters(DistributionModel* dm){
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* cov = Parameters_at(dm->parameters, 1);
    size_t dim = Parameter_size(mu);
    for (size_t i = 0; i < dim; i++) {
        gsl_vector_set(wrapper->mu, i, Parameter_value_at(mu, i));
    }
    if (wrapper->cholesky) {
        // cov is a full (row-major) covariance matrix: copy it and factor it.
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                gsl_matrix_set(wrapper->L, i, j, Parameter_value_at(cov, i * dim + j));
            }
        }
        gsl_linalg_cholesky_decomp1(wrapper->L);
    } else {
        // cov holds the packed lower-triangular Cholesky factor L (row-major).
        size_t row = 0;
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j <= i; j++) {
                gsl_matrix_set(wrapper->L, i, j, Parameter_value_at(cov, row));
                row++;
            }
        }
    }
    dm->need_update = false;
}

// Load one observation (a dim-vector) from dm->x into the gsl work vector.
static void _load_gsl_x(gsl_multivariate_normal_wrapper_t* wrapper, Parameter* x,
                        size_t dim){
    for (size_t i = 0; i < dim; i++) {
        gsl_vector_set(wrapper->x, i, Parameter_value_at(x, i));
    }
}

// Sum of the MVN log density over every observation in dm->x (i.i.d. draws sharing
// the same mu and L).
static double _multivariate_normal_logP(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    _update_gsl_parameters(dm);

    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* mu = Parameters_at(dm->parameters, 0);
    size_t dim = Parameter_size(mu);
    dm->lp = 0;
    for (size_t k = 0; k < Parameters_count(dm->x); k++) {
        _load_gsl_x(wrapper, Parameters_at(dm->x, k), dim);
        double logP;
        gsl_ran_multivariate_gaussian_log_pdf(wrapper->x, wrapper->mu, wrapper->L,
                                              &logP, wrapper->work);
        dm->lp += logP;
    }
    return dm->lp;
}

// Solve for w = L^{-1} r and alpha = Sigma^{-1} r = L^{-T} w, where r = x - mu.
// Both are needed by the gradient wrt mu, x and the Cholesky factor L.
static void _multivariate_normal_solve(gsl_multivariate_normal_wrapper_t* wrapper,
                                       size_t dim, gsl_vector* w, gsl_vector* alpha){
    for (size_t i = 0; i < dim; i++) {
        gsl_vector_set(w, i, gsl_vector_get(wrapper->x, i) - gsl_vector_get(wrapper->mu, i));
    }
    gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, wrapper->L, w);  // w = L^{-1} r
    gsl_vector_memcpy(alpha, w);
    gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, wrapper->L, alpha);  // alpha = L^{-T} w
}

// Gradient of the log density wrt x, mu and (when it holds L directly) the covariance
// parameter. The internal log-transform is not differentiated here — use an external
// Parameter transform for constrained variables.
static void _multivariate_normal_gradient(DistributionModel* dm, Parameters* parameters){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* cov = Parameters_at(dm->parameters, 1);
    size_t dim = Parameter_size(mu);
    size_t covSize = Parameter_size(cov);

    // mu and L are shared across observations, so their gradients sum over dm->x.
    // The gradient wrt a full covariance matrix is not implemented (cholesky case).
    Parameter* mux = Parameters_depends(parameters, mu);
    Parameter* covx = wrapper->cholesky ? NULL : Parameters_depends(parameters, cov);
    double* muGrad = mux != NULL ? dvector(dim) : NULL;          // summed over x
    double* covGrad = covx != NULL ? dvector(covSize) : NULL;    // summed over x

    gsl_vector* w = gsl_vector_alloc(dim);
    gsl_vector* alpha = gsl_vector_alloc(dim);

    for (size_t k = 0; k < Parameters_count(dm->x); k++) {
        Parameter* x = Parameters_at(dm->x, k);
        _load_gsl_x(wrapper, x, dim);
        _multivariate_normal_solve(wrapper, dim, w, alpha);

        // d logP / d x_k = -Sigma^{-1} r = -alpha (each observation is independent)
        Parameter* xx = Parameters_depends(parameters, x);
        if(xx != NULL){
            for (size_t i = 0; i < dim; i++) {
                dm->tempp[i] = -gsl_vector_get(alpha, i);
                x->grad[i] += dm->tempp[i];
            }
            if(x != xx) x->transform->backward(x->transform, dm->tempp);
        }

        // d logP / d mu += Sigma^{-1} r = alpha
        if(mux != NULL){
            for (size_t i = 0; i < dim; i++) {
                double g = gsl_vector_get(alpha, i);
                muGrad[i] += g;
                mu->grad[i] += g;
            }
        }

        // d logP / d L_ab += alpha_a * w_b - [a==b]/L_aa (packed lower-triangular).
        // The -1/L_aa diagonal term contributes once per observation.
        if(covx != NULL){
            size_t idx = 0;
            for (size_t a = 0; a < dim; a++) {
                double alpha_a = gsl_vector_get(alpha, a);
                for (size_t b = 0; b <= a; b++) {
                    double g = alpha_a * gsl_vector_get(w, b);
                    if(a == b) g -= 1.0 / gsl_matrix_get(wrapper->L, a, a);
                    covGrad[idx] += g;
                    cov->grad[idx] += g;
                    idx++;
                }
            }
        }
    }

    if(mux != NULL){
        if(mu != mux) mu->transform->backward(mu->transform, muGrad);
        free(muGrad);
    }
    if(covx != NULL){
        if(cov != covx) cov->transform->backward(cov->transform, covGrad);
        free(covGrad);
    }
    gsl_vector_free(w);
    gsl_vector_free(alpha);
}

static void _sample_multivariate_normal(DistributionModel* dm){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* mu = Parameters_at(dm->parameters, 0);
    size_t dim = Parameter_size(mu);

    for (size_t k = 0; k < Parameters_count(dm->x); k++) {
        Parameter* x = Parameters_at(dm->x, k);
        gsl_ran_multivariate_gaussian(dm->rng, wrapper->mu, wrapper->L, wrapper->x);
        for (size_t i = 0; i < dim; i++) {
            dm->tempx[i] = gsl_vector_get(wrapper->x, i);
        }
        Parameter_set_values(x, dm->tempx);
    }
}

// Reparameterized sample: z ~ N(0, I), x = mu + L z, one draw per observation. Each
// observation's z is stashed at dm->tempx + k*dim so _multivariate_normal_rgradient
// can reuse it.
static void _rsample_multivariate_normal(DistributionModel* dm){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* mu = Parameters_at(dm->parameters, 0);
    size_t dim = Parameter_size(mu);
    double* result = dvector(dim);

    for (size_t k = 0; k < Parameters_count(dm->x); k++) {
        Parameter* x = Parameters_at(dm->x, k);
        double* z = dm->tempx + k * dim;
        for (size_t i = 0; i < dim; i++) {
            z[i] = gsl_ran_ugaussian(dm->rng);
        }
        for (size_t a = 0; a < dim; a++) {
            double v = gsl_vector_get(wrapper->mu, a);
            for (size_t b = 0; b <= a; b++) {
                v += gsl_matrix_get(wrapper->L, a, b) * z[b];
            }
            result[a] = v;
        }
        Parameter_set_values(x, result);
    }
    free(result);
}

// Gradient wrt the variational parameters using the reparameterization trick.
// x_k = mu + L z_k, so given dP/dx_k (in x_k->grad), summed over observations:
//   dP/dmu_a  += dP/dx_k,a
//   dP/dL_ab  += dP/dx_k,a * z_k,b   (packed lower-triangular, a >= b)
static void _multivariate_normal_rgradient(DistributionModel* dm){
    Parameter* mu = Parameters_at(dm->parameters, 0);
    Parameter* cov = Parameters_at(dm->parameters, 1);
    size_t dim = Parameter_size(mu);
    size_t covSize = Parameter_size(cov);
    double* muGrad = dvector(dim);        // summed over x
    double* covGrad = dvector(covSize);   // summed over x

    for (size_t k = 0; k < Parameters_count(dm->x); k++) {
        Parameter* x = Parameters_at(dm->x, k);
        const double* z = dm->tempx + k * dim;
        size_t idx = 0;
        for (size_t a = 0; a < dim; a++) {
            mu->grad[a] += x->grad[a];
            muGrad[a] += x->grad[a];
            for (size_t b = 0; b <= a; b++) {
                double g = x->grad[a] * z[b];
                covGrad[idx] += g;
                cov->grad[idx] += g;
                idx++;
            }
        }
    }
    if(mu->transform != NULL) mu->transform->backward(mu->transform, muGrad);
    if(cov->transform != NULL) cov->transform->backward(cov->transform, covGrad);
    free(muGrad);
    free(covGrad);
}

// Differential entropy of N(mu, Sigma = L L^T):
//   H = k/2 (1 + log 2pi) + 1/2 log|Sigma| = k/2 (1 + log 2pi) + sum_i log L_ii
// summed over the observations in dm->x (independent draws sharing the parameters).
static double _multivariate_normal_entropy(DistributionModel* dm){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    size_t dim = Parameter_size(Parameters_at(dm->parameters, 0));
    size_t nObs = Parameters_count(dm->x);
    double logdet = 0;
    for (size_t i = 0; i < dim; i++) {
        logdet += log(gsl_matrix_get(wrapper->L, i, i));
    }
    return nObs * (0.5 * dim * (1.0 + LOG_TWO_PI) + logdet);
}

// Gradient of the entropy wrt the parameters. It is independent of mu, and depends
// on the covariance only through the diagonal of L:
//   d H / d L_aa = n / L_aa   (off-diagonal and mu contributions are zero)
// Not defined for a full covariance matrix (cholesky case).
static void _multivariate_normal_entropy_gradient(DistributionModel* dm,
                                                  const Parameters* parameters){
    _update_gsl_parameters(dm);
    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    Parameter* cov = Parameters_at(dm->parameters, 1);
    size_t dim = Parameter_size(Parameters_at(dm->parameters, 0));
    size_t nObs = Parameters_count(dm->x);

    Parameter* covx = wrapper->cholesky ? NULL : Parameters_depends(parameters, cov);
    if(covx != NULL){
        size_t covSize = Parameter_size(cov);
        double* covGrad = dvector(covSize);  // zero-initialised; off-diagonals stay 0
        size_t idx = 0;
        for (size_t a = 0; a < dim; a++) {
            for (size_t b = 0; b <= a; b++) {
                if(a == b){
                    double g = nObs / gsl_matrix_get(wrapper->L, a, a);
                    covGrad[idx] = g;
                    cov->grad[idx] += g;
                }
                idx++;
            }
        }
        if(cov != covx) cov->transform->backward(cov->transform, covGrad);
        free(covGrad);
    }
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
    dm->logP = _multivariate_normal_logP;
    dm->gradient = _multivariate_normal_gradient;
    dm->rgradient = _multivariate_normal_rgradient;
    dm->sample = _sample_multivariate_normal;
    dm->rsample = _rsample_multivariate_normal;
    dm->entropy = _multivariate_normal_entropy;
    dm->gradient_entropy = _multivariate_normal_entropy_gradient;
    dm->need_update = true;

    size_t dim = Parameter_size(Parameters_at(x, 0));
    size_t covSize = Parameter_size(Parameters_at(parameters, 1));
    gsl_multivariate_normal_wrapper_t* wrapper = (gsl_multivariate_normal_wrapper_t*)malloc(sizeof(gsl_multivariate_normal_wrapper_t));
    wrapper->mu = gsl_vector_calloc(dim);
    wrapper->L = gsl_matrix_calloc(dim, dim);
    wrapper->x = gsl_vector_calloc(dim);
    wrapper->work = gsl_vector_calloc(dim);
    wrapper->rng = NULL;        // assigned from dm->rng in _from_json
    // A full covariance matrix is dim*dim; a packed Cholesky factor is dim*(dim+1)/2.
    wrapper->cholesky = (covSize == dim * dim);
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
    clonewrapper->cholesky = wrapper->cholesky;
    clonewrapper->rng = wrapper->rng;
    gsl_vector_memcpy(clonewrapper->mu, wrapper->mu);
    gsl_vector_memcpy(clonewrapper->x, wrapper->x);
    gsl_matrix_memcpy(clonewrapper->L, wrapper->L);
    clonedm->data = clonewrapper;
    return clone;
}

Model* new_MultivariateNormalDistributionModel_from_json(json_node* node, Hashtable* hash){
    // The mean and covariance are given under a "parameters" object: "mu" and
    // exactly one of "sigma" (full covariance matrix) or "L" (packed lower-triangular
    // Cholesky factor). Alternatively "file" derives both empirically from a log, or
    // "parameters" is omitted for a standard normal (mu = 0, sigma = I).
    static const json_field schema[] = {
        {"burnin", JSON_OPTIONAL, JSON_NUMBER},
        {"distribution", JSON_REQUIRED, JSON_STRING},
        {"file", JSON_OPTIONAL, JSON_STRING},
        {"parameters", JSON_OPTIONAL, JSON_OBJECT},
        {"x", JSON_REQUIRED, JSON_OBJECT | JSON_STRING | JSON_ARRAY},
    };
    json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

    char* id = get_json_node_value_string(node, "id");

    // x may hold several observations (each a dim-vector) sharing this distribution.
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);

    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));  // MVN dimension
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
        // Exactly one of "sigma" (full covariance) or "L" (Cholesky factor).
        json_validate_xor(parameters_node, "sigma", "L", NULL);

        json_node* mu_node = get_json_node(parameters_node, "mu");
        json_node* sigma_node = get_json_node(parameters_node, "sigma");
        json_node* L_node = get_json_node(parameters_node, "L");
        if(mu_node == NULL){
            json_die(node, "MultivariateNormal - \"parameters\" must contain \"mu\"");
        }
        mu = new_Parameter_from_json(mu_node, hash);
        size_t p_count = Parameter_size(mu);
        if(sigma_node != NULL){
            sigma = new_Parameter_from_json(sigma_node, hash);
            if(Parameter_size(sigma) != p_count*p_count){
                json_die(node, "MultivariateNormal - covariance matrix \"sigma\" has dimension %zu but mu has dimension %zu (expected %zu)",
                         Parameter_size(sigma), p_count, p_count*p_count);
            }
        }
        else{
            sigma = new_Parameter_from_json(L_node, hash);
            if(Parameter_size(sigma) != p_count*(p_count-1)/2 + p_count){
                json_die(node, "MultivariateNormal - Cholesky factor \"L\" has dimension %zu but mu has dimension %zu (expected %zu)",
                         Parameter_size(sigma), p_count, p_count*(p_count-1)/2 + p_count);
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

    gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
    wrapper->rng = dm->rng;

    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
