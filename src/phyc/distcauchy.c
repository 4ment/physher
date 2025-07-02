//
//  distcauchy.c
//  physher
//
//  Created by Mathieu Fourment on 8/6/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distcauchy.h"

#include <strings.h>

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"
#include "mathconstant.h"


double dcauchy(double x, double location, double alpha){
    return gsl_ran_cauchy_pdf(x - location, alpha);
}

double rcauchy(const gsl_rng *rng, double location, double alpha){
    return gsl_ran_cauchy(rng, alpha) + location;
}

double DistributionModel_cauchy_log_prob(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* location = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(dcauchy(values[i], *location, *alpha));
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
                dm->lp += log(dcauchy(values[i], location[index], alpha[index]));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

double DistributionModel_cauchy_log_prob_gradient(DistributionModel* dm, const Parameters* parameters){
    size_t dimX = Parameters_count(dm->x);
    Parameter* location = Parameters_at(dm->parameters, 0);
    Parameter* alpha = Parameters_at(dm->parameters, 1);
    double* locationValues = NULL;
    double* alphaValues = NULL;

    distmodel_expand_2parameters(dm->x, dm->parameters, &locationValues, &alphaValues);

    Parameter* locationx = Parameters_depends(parameters, location);
    Parameter* alphax = Parameters_depends(parameters, alpha);

    if(locationx != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                dm->tempp[index] = 2.0*(xValues[j] -locationValues[index])/(alphaValues[index]*alphaValues[index] + pow(xValues[j] - locationValues[index], 2.0));
                location->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(location != locationx){
            location->transform->backward(location->transform, dm->tempp);
        }
    }

    if(alphax != NULL){
        size_t index = 0;
        for(size_t k = 0; k < dimX; k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            double z = 0.0;
            for(size_t j = 0; j < dim; j++){
                double locationValue = locationValues[index];
                double alphaValue = alphaValues[index];
                z = (xValues[j] - locationValue) / alphaValue;
                dm->tempp[index] =  -1.0/alphaValue + 2.0*z*z/(alphaValue*(1.0+z*z));
                alpha->grad[index] += dm->tempp[index];
                index++;
            }
        }
        if(alpha != alphax){
            alpha->transform->backward(alpha->transform, dm->tempp);
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
                dm->tempp[j] = -2.0*(xValues[j] -locationValues[index])/(alphaValues[index]*alphaValues[index] + pow(xValues[j] - locationValues[index], 2.0));
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

    if(Parameter_values(location) != locationValues){
        free(locationValues);
        free(alphaValues);
    }
    return 0;
}

static void DistributionModel_cauchy_sample(DistributionModel* dm){
    const double* location = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = rcauchy(dm->rng, *location, *alpha);
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
                dm->tempx[i] = rcauchy(dm->rng, location[index], alpha[index]);
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

DistributionModel* new_CauchyDistributionModel_with_parameters(Parameters* parameters, Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_CAUCHY;
	dm->log_prob = DistributionModel_cauchy_log_prob;
    dm->log_prob_grad = DistributionModel_cauchy_log_prob_gradient;
	dm->sample = DistributionModel_cauchy_sample;
    dm->shift = -INFINITY;
	return dm;
}

Model* new_CauchyDistributionModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = new_Parameters(1);
    distmodel_get_parameters(x_node, hash, x);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));
    Parameter* location = NULL;
    Parameter* scale = NULL;

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* locationValues = malloc(sizeof(double)*paramCount);
        double* alphaValues = malloc(sizeof(double)*paramCount);
        
        for (size_t i = 0; i < paramCount; i++) {
            //TODO: that coming from normal
            const double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            locationValues[i] = m;
            alphaValues[i] = sqrt(v);
            free_Vector(samples[i]);
        }

        location = new_Parameter2("location", locationValues, paramCount, new_Constraint(-INFINITY, INFINITY));
        scale = new_Parameter2("scale", alphaValues, paramCount, new_Constraint(0, INFINITY));

        free(locationValues);
        free(alphaValues);
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* locationValues = malloc(sizeof(double)*paramCount);
        double* alphaValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            locationValues[i] = 0;
            alphaValues[i] = 1;
        }

        location = new_Parameter("location", 0, new_Constraint(-INFINITY, INFINITY));
        scale = new_Parameter("scale", 1, new_Constraint(0, INFINITY));
        
        free(locationValues);
        free(alphaValues);
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "location") != 0 && strcasecmp(parameters_node->children[i]->key, "scale") != 0) {
                fprintf(stderr, "Normal distribution should be parametrized with location and scale\n");
                exit(13);
            }
        }

        json_node* alpha_node = get_json_node(parameters_node, "location");
        json_node* beta_node = get_json_node(parameters_node, "scale");
        location = new_Parameter_from_json(alpha_node, hash);
        scale = new_Parameter_from_json(beta_node, hash);
    }

    Parameters_add(parameters, location);
    Parameters_add(parameters, scale);
    
    DistributionModel* dm = new_CauchyDistributionModel_with_parameters(parameters, x);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
