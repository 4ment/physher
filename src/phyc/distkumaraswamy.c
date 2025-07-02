//
//  distkumaraswamy.c
//  physher
//
//  Created by Mathieu Fourment on 19/8/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distkumaraswamy.h"

#include <strings.h>

#include <gsl/gsl_randist.h>

#include "parametersio.h"


double DistributionModel_log_kumaraswamy(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* a = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* b = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(a[0]*b[0]*pow(values[i], a[0] - 1.0) * pow( 1.0 - pow(values[i], a[0]), b[0] - 1.0));
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
                dm->lp += log(a[index]*b[index]*pow(values[i], a[index] - 1.0) * pow( 1.0 - pow(values[i], a[index]), b[index] - 1.0));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_kumaraswamy_log_prob_gradient(DistributionModel* dm, const Parameters* parameters){
    size_t dimX = Parameters_count(dm->x);
    Parameter* a = Parameters_at(dm->parameters, 0);
    Parameter* b = Parameters_at(dm->parameters, 1);
    double* aValues = NULL;
    double* bValues = NULL;

    distmodel_expand_2parameters(dm->x, dm->parameters, &aValues, &bValues);

    Parameter* ax = Parameters_depends(parameters, a);
    Parameter* bx = Parameters_depends(parameters, b);

    if(ax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                double aValue = aValues[index];
                double bValue = bValues[index];
                dm->tempp[index] = 1.0/aValue + log(xValues[j]) + (bValue - 1.0) * (-pow(xValues[j], aValue)*log(xValues[j]))/(1.0 - pow(xValues[j], aValue));
                a->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(a != ax){
            a->transform->backward(a->transform, dm->tempp);
        }
    }

    if(bx != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                dm->tempp[index] =  1.0/bValues[index] + log(1.0-pow(xValues[j], aValues[index]));
                b->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(b != bx){
            b->transform->backward(b->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for(size_t k = 0; k < dimX; k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        if (xx != NULL) {
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
        
            for(size_t j = 0; j < dim; j++){
                double aValue = aValues[index];
                double bValue = bValues[index];
                dm->tempp[j] = (aValue - 1.0)/xValues[j] - (aValue*bValue*pow(xValues[j], 2*aValue - 1.0))/(1.0-pow(xValues[j], aValue));
                x->grad[j] += dm->tempp[j];
                index++;
            }
            if(x != xx){
                x->transform->backward(x->transform, dm->tempp);
            }
        }
        else{
            index += Parameter_size(x);
        }
    }

    if(Parameter_values(a) != aValues){
        free(aValues);
        free(bValues);
    }
    return 0;
}

// F_X(x) = p(X <= x)
double DistributionModel_kumaraswamy_inverse_CDF(double p, double a, double b){
	return pow(1.0 - pow(1.0 - p, 1.0/b), 1.0/a);
}

static void DistributionModel_kumaraswamy_sample(DistributionModel* dm){
    const double* a = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* b = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), *a, *b);
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
                dm->tempx[i] = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a[index], b[index]);
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

DistributionModel* new_KumaraswamyDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_KUMARASWAMY;
	dm->parameterization = 0;
	
	dm->log_prob = DistributionModel_log_kumaraswamy;
	dm->log_prob_grad = DistributionModel_kumaraswamy_log_prob_gradient;
	dm->sample = DistributionModel_kumaraswamy_sample;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = 1;
	return dm;
}

Model* new_KumaraswamyDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = new_Parameters(1);
    distmodel_get_parameters(x_node, hash, x);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
	size_t paramCount = Parameter_size(Parameters_at(x, 0));
	Parameter* a = NULL;
    Parameter* b = NULL;

    if(file != NULL){
           error("Cannot estimate parameters of the Kumaraswamy distribution from a sample\n");
           exit(3);
       }
    

    if(get_json_node(node, "parameters") == NULL){
        double* aValues = malloc(sizeof(double)*paramCount);
        double* bValues = malloc(sizeof(double)*paramCount);

		for (size_t i = 0; i < paramCount; i++) {
            aValues[i] = 0;
            bValues[i] = 1;
        }

        a = new_Parameter2("a", aValues, paramCount, new_Constraint(0, INFINITY));
        b = new_Parameter2("b", bValues, paramCount, new_Constraint(0, INFINITY));

		free(aValues);
        free(bValues);
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "a") != 0 && strcasecmp(parameters_node->children[i]->key, "b") != 0) {
                fprintf(stderr, "Kumaraswamy distribution should be parametrized with a and b\n");
                exit(13);
            }
        }
        json_node* a_node = get_json_node(parameters_node, "a");
        json_node* b_node = get_json_node(parameters_node, "b");
        a = new_Parameter_from_json(a_node, hash);
        b = new_Parameter_from_json(b_node, hash);
    }

	Parameters_add(parameters, a);
    Parameters_add(parameters, b);
    
    DistributionModel* dm = new_KumaraswamyDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
