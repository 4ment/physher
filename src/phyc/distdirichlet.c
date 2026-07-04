// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "distdirichlet.h"

#include <assert.h>
#include <strings.h>

#ifndef GSL_DISABLED
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#endif

#include "matrix.h"
#include "dirichlet.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"

// Flat dirichlet
double DistributionModel_log_flat_dirichlet(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
    dm->lp = 0;
    for(size_t i = 0; i < Parameters_count(dm->x); i++){
        Parameter* x = Parameters_at(dm->x, i);
        dm->lp += log(ddirchlet_flat(Parameter_size(x)));
    }
	dm->need_update = false;
	return dm->lp;
}

void DistributionModel_flat_dirichlet_gradient(DistributionModel* dm, Parameters* parameters){
    // gradient wrt x is 0
    // if it is flat then alpha is constant and set to 1 so gradient should not be calculated
}

double DistributionModel_log_dirichlet(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* x = Parameters_at(dm->x, 0);
#ifndef GSL_DISABLED
    dm->lp = gsl_ran_dirichlet_lnpdf(Parameter_size(x), Parameter_values(alpha), Parameter_values(x));
#else
    dm->lp = ddirchletln(Parameter_values(x), Parameter_size(x), Parameter_values(alpha));
#endif
	dm->need_update = false;
	return dm->lp;
}

static void DistributionModel_dirichlet_sample(DistributionModel* dm){
    Parameter* x = Parameters_at(dm->x, 0);
#ifndef GSL_DISABLED
    gsl_ran_dirichlet(dm->rng, Parameter_size(x), Parameter_values(Parameters_at(dm->parameters, 0)), dm->tempx);
#else
    rdirichlet(dm->tempx, Parameter_size(x), Parameter_values(Parameters_at(dm->parameters, 0)));
#endif
    Parameter_set_values(x, dm->tempx);
}

void DistributionModel_dirichlet_gradient(DistributionModel* dm, Parameters* parameters){
    // The Dirichlet is multivariate: a single simplex draw x parameterized by the
    // whole concentration vector alpha (Parameter_size(x) == Parameter_size(alpha)).
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* x = Parameters_at(dm->x, 0);
    size_t dim = Parameter_size(x);
    const double* xValues = Parameter_values(x);
    const double* alphaValues = Parameter_values(alpha);

    /*
    log pdf(X; alpha) = sum_i (alpha_i - 1) log(x_i) - log B(alpha)
    log B(alpha)      = sum_i log Gamma(alpha_i) - log Gamma(alpha_0),  alpha_0 = sum_j alpha_j

    d log pdf/dalpha_k = log(x_k) - d log B/dalpha_k
                       = log(x_k) - psi(alpha_k) + psi(alpha_0)
    */
    Parameter* alphax = Parameters_depends(parameters, alpha);
    if (alphax != NULL) {
        double alpha0 = 0.0;
        for (size_t i = 0; i < dim; i++) alpha0 += alphaValues[i];
        double psi0 = gsl_sf_psi(alpha0);
        for (size_t i = 0; i < dim; i++) {
            dm->tempp[i] = log(xValues[i]) - gsl_sf_psi(alphaValues[i]) + psi0;
            alpha->grad[i] += dm->tempp[i];
        }
        if (alphax != alpha) {
            alpha->transform->backward(alpha->transform, dm->tempp);
        }
    }

    /*
    d log pdf/dx_k = (alpha_k - 1)/x_k
    d log pdf/dy_k = sum_j d log pdf/dx_j dx_j/dy_k   (simplex transform backward)
    */
    Parameter* xx = Parameters_depends(parameters, x);
    if (xx != NULL) {
        double* grad = dvector(dim);
        for (size_t i = 0; i < dim; i++) {
            grad[i] = (alphaValues[i] - 1.0)/xValues[i];
            x->grad[i] += grad[i];
        }
        if (xx != x) {
            x->transform->backward(x->transform, grad);
        }
        free(grad);
    }
}

DistributionModel* new_FlatDirichletDistributionModel(Parameters* x){
	DistributionModel* dm = new_DistributionModel(NULL, x);
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_flat_dirichlet;
    dm->gradient = DistributionModel_flat_dirichlet_gradient;
	dm->sample = DistributionModel_dirichlet_sample;
    dm->shift = 0;
	return dm;
}

DistributionModel* new_DirichletDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters,  x);
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_dirichlet;
    dm->gradient = DistributionModel_dirichlet_gradient;
	dm->sample = DistributionModel_dirichlet_sample;
    dm->shift = 0;
	return dm;
}


Model* new_DirichletDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = new_Parameters(1);
    char* ref = (char*)x_node->value;
    if (safe_is_reference(ref, id)) {
        Parameter* xx = safe_get_reference_parameter(ref, hash, id);
        Parameters_add(x, xx);
    }
    else{
        char* x_type = get_json_node_value_string(x_node, "type");
        Parameter* xx = NULL;
        if(strcasecmp(x_type, "simplex") != 0){
            // xx = new_SimplexParameter_from_json(x_node, hash);
            xx = new_Parameter_from_json(x_node, hash);
        }
       else{
            xx = new_Parameter_from_json(x_node, hash);
        }
        Hashtable_add(hash, Parameter_name(xx), xx);
        Parameters_move(x, xx);
    }
    // distmodel_get_parameters(x_node, hash, x);
    if(Parameters_count(x) != 1){
        fprintf(stderr, "%s - Dirichlet distribution should have one x (%zu)\n", id, Parameters_count(x));
        exit(13);
    }
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = NULL;
	DistributionModel* dm = NULL;
	
    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        size_t paramCount = Parameter_size(Parameters_at(x, 0));
        parameters = new_Parameters(1);

        double* means = dvector(paramCount);
        double* variances = dvector(paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            means[i] = dmean(Vector_data(samples[i]), Vector_length(samples[i]));
            variances[i] = variance(Vector_data(samples[i]), Vector_length(samples[i]), means[i]);
        }
        double num = 0;
        double denom = 0;
        for (size_t i = 0; i < paramCount; i++) {
            num += means[i]*pow(1.0 - means[i], 2);
            denom += means[i]*variances[i]*(1.0 - means[i]);
        }
        double mhat = num/denom - 1.0;
        for (size_t i = 0; i < paramCount; i++) {
            means[i] *= mhat;
        }
        Parameters_move(parameters, new_Parameter2("alpha", means, paramCount, new_Constraint(0, INFINITY)));
        for (size_t i = 0; i < paramCount; i++) {
            free_Vector(samples[i]);
        }
        free(samples);
        free(means);
        free(variances);

		dm = new_DirichletDistributionModel_with_parameters(parameters, x);
    }
    // Flat dirichlet
    else if(get_json_node(node, "parameters") == NULL){
        dm = new_FlatDirichletDistributionModel(x);
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        
        if (parameters_node->child_count != 1 && strcasecmp(parameters_node->children[0]->key, "concentration") != 0 && strcasecmp(parameters_node->children[0]->key, "alpha") != 0) {
            fprintf(stderr, "Dirichlet distribution should be parametrized with concentration parameter\n");
            exit(13);
        }

        json_node* alpha_node = get_json_node(parameters_node, "alpha");
        Parameter* alpha = new_Parameter_from_json(alpha_node, hash);
        parameters = new_Parameters(1);
        Parameters_add(parameters, alpha);

		dm = new_DirichletDistributionModel_with_parameters(parameters, x);
    }
    
    dm->parameterization = 0;
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
#ifndef GSL_DISABLED
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
#endif
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
