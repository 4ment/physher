//
//  distdirichlet.c
//  physher
//
//  Created by Mathieu Fourment on 2/04/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distdirichlet.h"

#include <assert.h>
#include <strings.h>

#ifndef GSL_DISABLED
#include <gsl/gsl_randist.h>
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

double DistributionModel_flat_dirichlet_gradient2(DistributionModel* dm, const Parameters* parameters){
    // gradient wrt x is 0
    // if it is flat then alpha is constant and set to 1 so gradient should not be calculated
    return 0;
}

double DistributionModel_d2log_flat_dirichlet(DistributionModel* dm, const Parameter* p){
	return 0.0;
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

// static double DistributionModel_dirichlet_sample_evaluate(DistributionModel* dm){
//     DistributionModel_dirichlet_sample(dm, NULL);
//     return DistributionModel_log_dirichlet(dm);
// }

double DistributionModel_dirichlet_gradient2(DistributionModel* dm, const Parameters* parameters){
    // TODO: implement gradient wrt alpha
    for(size_t k = 0; k < Parameters_count(dm->x); k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        
        if (xx != NULL) {
            /*
            log pdf(X; \alpha) = \sum_i (\alpha_i - 1)log(x_i) - log B(\alpha)
            d log pdf(X)/dx_k  = \sum_i (\alpha_i - 1) d log(x_i)/dx_k
                              &= (\alpha_k - 1)/x_k
            
            d log pdf(X)/dy_k  = sum_j d log pdf(X)/dx_j dx_j/dy_k
            
            dx_i/dy_k = 0 for i < k
            */
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            const double* alphaValues = Parameter_values(Parameters_at(dm->parameters, 0));
            double* grad = dvector(dim);
            for(size_t i = 0; i < dim; i++){
                grad[i] = (alphaValues[i] - 1.0)/xValues[i];
                x->grad[i] += grad[i];
            }
            if(xx != x){
                x->transform->backward(x->transform, grad);
            }
            free(grad);
        }
    }
    return 0;
}

// IMPORTANT: The derivative is wrt unconstrained parameter of the simplex
double DistributionModel_dlog_dirichlet(DistributionModel* dm, const Parameter* p){
    /*
    log pdf(X; \alpha) = \sum_i \alpha_i log(x_i) - log B(\alpha)
    d log pdf(X)/dz_k = \sum_i \alpha_i d log(x_i)/dz_k
                     &= \sum_i \alpha_i d log(x_i)/dx_i dx_i/dz_k
                     &= \sum_i \alpha_i/x_i dx_i/dz_k
    
    dx_i/dz_k = 0 for i < k
    */
   //TODO: implement
	// for (int i = 0; i < Parameters_count(dm->x); i++) {
    //     if( p == Parameters_at(dm->x, i)){
    //         const double* values = dm->simplex->get_values(dm->simplex);
    //         double dlogp = 0;
    //         for (size_t j = i; j < Parameters_count(dm->x); j++) {
    //             dm->simplex->gradient(dm->simplex, i, dm->tempp);
    //             dlogp += (Parameters_value(dm->parameters[0], j)-1.0)/values[j] * dm->tempp[j];
    //         }
    //         return dlogp;
    //     }
    // }
	return 0;
}

//TODO: implement
double DistributionModel_d2log_dirichlet(DistributionModel* dm, const Parameter* p){
	// find corresponding alpha
	// return -(alpha-1.0)/(Parameter_value(p)*Parameter_value(p));
	exit(1);
	return 0;
}

static void _DistributionModel_error_sample_dirichlet(DistributionModel* dm, double* samples){
    fprintf(stderr, "_DistributionModel_error_sample_dirichlet not implemented\n");
    exit(1);
}

// static double _DistributionModel_error_sample_evaluate_dirichlet(DistributionModel* dm){
//     fprintf(stderr, "_DistributionModel_error_sample_evaluate_dirichlet not implemented\n");
//     exit(1);
// }

DistributionModel* new_FlatDirichletDistributionModel(Parameters* x){
	DistributionModel* dm = new_DistributionModel(NULL, x);
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_flat_dirichlet;
    dm->gradient2 = DistributionModel_flat_dirichlet_gradient2;
	dm->d2logP = DistributionModel_d2log_flat_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_dirichlet_sample;
    // dm->logP_with_values = DistributionModel_log_flat_dirichlet_with_values;
	// dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
	// dm->tempp = dvector(Parameter_size(x));
	// for (int i = 0; i < Parameter_size(x); i++) {
	// 	dm->tempp[i] = 1;
	// }
    dm->shift = 0;
	return dm;
}

DistributionModel* new_DirichletDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
    DistributionModel* dm = new_DistributionModel(parameters,  x);
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_dirichlet;
    dm->gradient2 = DistributionModel_dirichlet_gradient2;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_dirichlet_sample;
	// dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	// dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
    // dm->tempp = dvector(Parameter_size(x));
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
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
#ifndef GSL_DISABLED
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
#endif
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
