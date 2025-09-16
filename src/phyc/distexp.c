//
//  distexp.c
//  physher
//
//  Created by Mathieu Fourment on 15/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distexp.h"

#include <string.h>
#include <strings.h>

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "statistics.h"


// parameter is lambda or mean
static void _expand_parameters(distribution_parameterization parameterization, Parameter* x, Parameter* parameter, double* out){
    size_t parameterSize = Parameter_size(parameter);
    memcpy(out, Parameter_values(parameter), sizeof(double)*parameterSize);
    if(parameterization == DISTRIBUTION_EXPONENTIAL_MEAN){
        for(size_t k = 0; k < parameterSize; k++){
            out[k] = 1.0/out[k];
        }
    }

    if(parameterSize == 1){
        size_t dim = Parameter_size(x);
        if(dim > 1){
            for(size_t k = 1; k < dim; k++){
                out[k] = out[0];
            }
        }
    }
}

double DistributionModel_log_exp(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    size_t dimX = Parameters_count(dm->x);
    for(size_t j = 0; j < dimX; j++){
        Parameter* x = Parameters_at(dm->x, j);
        Parameter* parameter = Parameters_at(dm->parameters, 0);
        _expand_parameters(dm->parameterization, x, parameter, dm->tempp);
        size_t dim = Parameter_size(x);
        const double* values = Parameter_values(x);
        for (size_t i = 0; i < dim; i++) {
            dm->lp += log(dm->tempp[i]) - dm->tempp[i] * values[i];
        }
    }
    
    dm->need_update = false;
    return dm->lp;
}


// double DistributionModel_log_exp_with_values(DistributionModel* dm, const double* values){
// 	double logP = 0;
// 	if(Parameters_count(dm->parameters[0]) > 1){
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double lambda = Parameters_value(dm->parameters[0], i);
// 			logP += log(lambda) - lambda * values[i];
// 		}
// 	}
// 	else{
// 		double lambda = Parameters_value(dm->parameters[0], 0);
// 		logP = log(lambda) * Parameters_count(dm->x);
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			logP -= lambda * values[i];
// 		}
// 	}
// 	return logP;
// }

// double DistributionModel_log_exp_with_values_mean(DistributionModel* dm, const double* values){
// 	double logP = 0;
// 	if(Parameters_count(dm->parameters[0]) > 1){
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double mean = Parameters_value(dm->parameters[0], i);
// 			logP += -log(mean) - values[i]/mean;
// 		}
// 	}
// 	else{
// 		double mean = Parameters_value(dm->parameters[0], 0);
// 		logP = -log(mean) * Parameters_count(dm->x);
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			logP -= values[i]/mean;
// 		}
// 	}
// 	return logP;
// }

double DistributionModel_gradient_exp(DistributionModel* dm, const Parameters* parameters){
    size_t dimX = Parameters_count(dm->x);
    Parameter* parameter = Parameters_at(dm->parameters, 0); // lambda or mean

	for (size_t i = 0; i < dimX; i++) {
        Parameter* x = Parameters_at(dm->x, i);
        Parameter* xx = Parameters_depends(parameters, x); // check if it is x or x->transform->parameter
        // x is constrained and matches support of distribution (e.g. branch lengths)
        // xx is unconstrained (e.g. transformed branch lengths)
		if (xx != NULL) {
            // tempp contains lambda
            _expand_parameters(dm->parameterization, x, parameter, dm->tempp);
            size_t sizeX = Parameter_size(x);
            for (size_t j = 0; j < sizeX; j++) {
                dm->tempx[j] = -dm->tempp[j];
                x->grad[j] += dm->tempx[j];
            }
            if(xx != x){
                Transform* t = x->transform;
                t->backward(t, dm->tempx);
            }
		}
	}
    Parameter* px = Parameters_depends(parameters, parameter);
    if(px != NULL){
        for (size_t i = 0; i < dimX; i++) {
            Parameter* x = Parameters_at(dm->x, i);
            const double* xValues = Parameter_values(x);
            size_t sizeX = Parameter_size(x);
            if(dm->parameterization == DISTRIBUTION_EXPONENTIAL_RATE){
                _expand_parameters(dm->parameterization, x, parameter, dm->tempp);
                for (size_t j = 0; j < sizeX; j++) {
                    dm->tempx[j] = -xValues[j] - 1.0/dm->tempp[j];
                    x->grad[j] += dm->tempx[j];
                }
            }
            else {
                // tempp contains mean
                _expand_parameters(DISTRIBUTION_EXPONENTIAL_RATE, x, parameter, dm->tempp);
                for (size_t j = 0; j < sizeX; j++) {
                    dm->tempx[j] = xValues[j]/(pow(dm->tempp[j], 2)) + 1.0/dm->tempp[j];
                    x->grad[j] += dm->tempx[j];
                }
            }
            if(px != parameter){
                Transform* t = parameter->transform;
                t->backward(t, dm->tempx);
            }
        }
    }
	return 0;
}

double DistributionModel_dlog_exp(DistributionModel* dm, const Parameter* p){
	//TODO: implement
    // if (p == Parameters_at(dm->parameters[0], 0)) {
    //     double dlogf = Parameters_count(dm->x)/Parameters_value(dm->parameters[0], 0);
    //     for(int i = 0; i < Parameters_count(dm->x); i++){
    //         dlogf -= Parameters_value(dm->x, i);
    //     }
    //     return dlogf;
    // }
	// for (int i = 0; i < Parameters_count(dm->x); i++) {
	// 	if (p == Parameters_at(dm->x, i)) {
    //         if(Parameters_count(dm->parameters[0]) == 1){
    //             return -Parameters_value(dm->parameters[0], 0);
    //         }
    //         return -Parameters_value(dm->parameters[0], i);
	// 	}
	// }
	return 0;
}

double DistributionModel_dlog_exp_mean(DistributionModel* dm, const Parameter* p){
	//TODO: implement
    // if (p == Parameters_at(dm->parameters[0], 0)) {
    //     double mean = Parameters_value(dm->parameters[0], 0);
    //     double dlogf = -Parameters_count(dm->x)/mean;
    //     double mean2 = mean;
    //     for(int i = 0; i < Parameters_count(dm->x); i++){
    //         dlogf += Parameters_value(dm->x, i)/mean2;
    //     }
    //     return dlogf;
    // }
    
	// for (int i = 0; i < Parameters_count(dm->x); i++) {
	// 	if (p == Parameters_at(dm->x, i)) {
    //         if(Parameters_count(dm->parameters[0]) == 1){
    //             return -1.0/Parameters_value(dm->parameters[0], 0);
    //         }
    //         return -1.0/Parameters_value(dm->parameters[0], i);
	// 	}
	// }
	return 0;
}

static void DistributionModel_exp_sample(DistributionModel* dm){
    const double* lambda = Parameter_values(Parameters_at(dm->parameters, 0));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
				dm->tempx[i] = gsl_ran_exponential(dm->rng, 1.0/ lambda[0]);
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
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_exponential(dm->rng, 1.0/ lambda[index]);
                index++;
            }
			Parameter_set_values(x, dm->tempx);
        }
    }
}

static void DistributionModel_exp_sample_mean(DistributionModel* dm){
 	const double* mean = Parameter_values(Parameters_at(dm->parameters, 0));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            for (size_t i = 0; i < dim; i++) {
				dm->tempx[i] = gsl_ran_exponential(dm->rng, mean[0]);
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
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_exponential(dm->rng, mean[index]);
                index++;
            }
			Parameter_set_values(x, dm->tempx);
        }
    }
}

// static double DistributionModel_exp_sample_evaluate(DistributionModel* dm){
// 	if(Parameters_count(dm->parameters[0]) > 1){
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double sample = gsl_ran_exponential(dm->rng, 1.0/Parameters_value(dm->parameters[0], i));
// 			Parameters_set_value(dm->x, i, sample);
// 		}
// 	}
// 	else{
// 		double mean = 1.0/Parameters_value(dm->parameters[0], 0);
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double sample = gsl_ran_exponential(dm->rng, mean);
// 			Parameters_set_value(dm->x, i, sample);
// 		}
// 	}
// 	return DistributionModel_log_exp(dm);
// }

// static double DistributionModel_exp_sample_evaluate_mean(DistributionModel* dm){
// 	if(Parameters_count(dm->parameters[0]) > 1){
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double sample = gsl_ran_exponential(dm->rng, Parameters_value(dm->parameters[0], i));
// 			Parameters_set_value(dm->x, i, sample);
// 		}
// 	}
// 	else{
// 		double mean = Parameters_value(dm->parameters[0], 0);
// 		for (int i = 0; i < Parameters_count(dm->x); i++) {
// 			double sample = gsl_ran_exponential(dm->rng, mean);
// 			Parameters_set_value(dm->x, i, sample);
// 		}
// 	}
// 	return DistributionModel_log_exp_mean(dm);
// }

DistributionModel* new_ExponentialDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_EXPONENTIAL;
	dm->parameterization = parameterization;
    dm->logP = DistributionModel_log_exp;
    dm->gradient2 = DistributionModel_gradient_exp;
    dm->sample = DistributionModel_exp_sample;
	if(parameterization == DISTRIBUTION_EXPONENTIAL_RATE){
		// dm->logP_with_values = DistributionModel_log_exp_with_values;
		dm->dlogP = DistributionModel_dlog_exp;
		// dm->sample_evaluate = DistributionModel_exp_sample_evaluate;
	}
	else{
		// dm->logP_with_values = DistributionModel_log_exp_with_values_mean;
		dm->dlogP = DistributionModel_dlog_exp_mean;
		dm->sample = DistributionModel_exp_sample_mean;
		// dm->sample_evaluate = DistributionModel_exp_sample_evaluate_mean;
	}
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = INFINITY;
	return dm;
}

Model* new_ExponentialDistributionModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
	size_t paramCount = Parameter_size(Parameters_at(x, 0));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(1);
    Parameters_set_name2(parameters, "parameters.exp");
	Parameter* parameter = NULL;
    distribution_parameterization parameterization = DISTRIBUTION_EXPONENTIAL_RATE;

    // empirical
    if (file != NULL) {
        parameterization = DISTRIBUTION_EXPONENTIAL_MEAN;
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* pValues = malloc(sizeof(double)*paramCount);
        
        for (size_t i = 0; i < paramCount; i++) {
            const double* vec = Vector_data(samples[i]);
            pValues[i] = mean(vec, Vector_length(samples[i]));
            // Parameters_move(parameters[0], new_Parameter("mean", m, new_Constraint(0, INFINITY)));
            free_Vector(samples[i]);
        }

		parameter = new_Parameter2("mean", pValues, paramCount, new_Constraint(0, INFINITY));

		free(pValues);
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        parameterization = DISTRIBUTION_EXPONENTIAL_MEAN;
		double* pValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            pValues[i] = 1;
    	}
		parameter = new_Parameter2("mean", pValues, paramCount, new_Constraint(0, INFINITY));
		free(pValues);
	}
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        json_node* parameter_node = get_json_node(parameters_node, "lambda");
        if(parameter_node == NULL){
            parameter_node = get_json_node(parameters_node, "mean");
            if(parameter_node != NULL){
                parameterization = DISTRIBUTION_EXPONENTIAL_MEAN;
            }
        }

        if(parameter_node == NULL){
            fprintf(stderr, "Exponential distribution should be parametrized with mean or lambda\n");
            exit(13);
        }

        parameter = distmodel_parse_parameter(parameter_node, hash, "", 0, INFINITY);
    }
    
	Parameters_move(parameters, parameter);

    DistributionModel* dm = new_ExponentialDistributionModel_with_parameters(parameters, x, parameterization);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(x);
    free_Parameters(parameters);
	
	return model;
}
