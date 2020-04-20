//
//  distgamma.c
//  physher
//
//  Created by Mathieu Fourment on 30/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distgamma.h"

#include <strings.h>

#include <gsl/gsl_randist.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"

double DistributionModel_log_gamma(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	dm->lp = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters[0], i);
			double beta = Parameters_value(dm->parameters[1], i);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(gsl_ran_gamma_pdf(x, alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters[0], 0);
		double beta = Parameters_value(dm->parameters[1], 0);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(gsl_ran_gamma_pdf(x, alpha, beta));
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_gamma_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters[0], i);
			double beta = Parameters_value(dm->parameters[1], i);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			logP += log(gsl_ran_gamma_pdf(values[i] - dm->shift, alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters[0], 0);
		double beta = Parameters_value(dm->parameters[1], 0);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(gsl_ran_gamma_pdf(values[i] - dm->shift, alpha, beta));
		}
	}
	return logP;
}

double DistributionModel_dlog_gamma(DistributionModel* dm, const Parameter* p){
    double dlogP = 0;
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (p == Parameters_at(dm->x, i)) {
			double alpha = Parameters_value(dm->parameters[0], 0);
			double beta = Parameters_value(dm->parameters[1], 0);
			double x = Parameter_value(p) - dm->shift;
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
				beta = 1.0/beta;
			}
			return (alpha-1.0)/x - beta;
			
		}
	}
	return 0;
}

double DistributionModel_dlog_gamma_multi(DistributionModel* dm, const Parameter* p){
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (p == Parameters_at(dm->x, i)) {
            double alpha = Parameters_value(dm->parameters[0], i);
            double beta = Parameters_value(dm->parameters[1], i);
            double x = Parameter_value(p) - dm->shift;
            if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
                beta = 1.0/beta;
            }
            return (alpha-1.0)/x - beta;
        }
    }
    return 0;
}

double DistributionModel_d2log_gamma(DistributionModel* dm, const Parameter* p){
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
            double alpha = Parameters_value(dm->parameters[0], 0);
            double x = Parameter_value(p) - dm->shift;
            return -(alpha-1.0)/x/x;
        }
    }
    return 0;
}

double DistributionModel_d2log_gamma_multi(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters[0], i);
			double x = Parameter_value(p) - dm->shift;
			return -(alpha-1.0)/x/x;
		}
	}
	return 0;
}

static void DistributionModel_gamma_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters[0], i);
			double beta = Parameters_value(dm->parameters[1], i);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			samples[i] = gsl_ran_gamma(dm->rng, alpha, beta);
			samples[i] += dm->shift;
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters[0], 0);
		double beta = Parameters_value(dm->parameters[1], 0);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_gamma(dm->rng, alpha, beta);
			samples[i] += dm->shift;
		}
	}
}

static double DistributionModel_gamma_sample_evaluate(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters[0], i);
			double beta = Parameters_value(dm->parameters[1], i);
			if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
				beta = 1.0/beta;
			}
			double sample = gsl_ran_gamma(dm->rng, alpha, beta);
			Parameters_set_value(dm->x, i, sample + dm->shift);
			logP += log(gsl_ran_gamma_pdf(sample, alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters[0], 0);
		double beta = Parameters_value(dm->parameters[1], 1);
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			beta = 1.0/beta;
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_gamma(dm->rng, alpha, beta);
			Parameters_set_value(dm->x, i, sample + dm->shift);
			logP += log(gsl_ran_gamma_pdf(sample, alpha, beta));
		}
	}
	return logP;
}

DistributionModel* new_GammaDistributionModel_with_parameters(Parameters** parameters, Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, 2, x);
	dm->type = DISTRIBUTION_GAMMA;
	dm->logP = DistributionModel_log_gamma;
	dm->logP_with_values = DistributionModel_log_gamma_with_values;
    if(Parameters_count(parameters[0]) > 1){
        dm->dlogP = DistributionModel_dlog_gamma_multi;
        dm->d2logP = DistributionModel_d2log_gamma_multi;
    }
    else{
        dm->dlogP = DistributionModel_dlog_gamma;
        dm->d2logP = DistributionModel_d2log_gamma;
    }
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_gamma_sample;
	dm->sample_evaluate = DistributionModel_gamma_sample_evaluate;
	dm->parameterization = parameterization;
    dm->shift = 0;
	return dm;
}

Model* new_GammaDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    
    char* file = get_json_node_value_string(node, "file");
    Parameters** parameters = NULL;
    size_t parameters_dim = 2;
    
    bool scale = false;
    
    char* parameterization_string = get_json_node_value_string(node, "parameterization");
    if(parameterization_string!= NULL && strcasecmp(parameterization_string, "scale") == 0){
        scale = true;
    }

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
            Parameters_move(parameters[0], new_Parameter("shape", m*m/v, new_Constraint(0, INFINITY)));
            if (scale) {
                Parameters_move(parameters[1], new_Parameter("scale", v/m, new_Constraint(0, INFINITY)));
            }
            else{
                Parameters_move(parameters[1], new_Parameter("rate", m/v, new_Constraint(0, INFINITY)));
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
            Parameters_move(parameters[0], new_Parameter("shape", 0, new_Constraint(-INFINITY, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("rate", 1, new_Constraint(0, INFINITY)));
        }
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "shape") != 0 && strcasecmp(parameters_node->children[i]->key, "rate") != 0 && strcasecmp(parameters_node->children[i]->key, "scale") != 0) {
                fprintf(stderr, "Gamma distribution should be parametrized with shape and (rate or scale)\n");
                exit(13);
            }
            if(strcasecmp(parameters_node->children[i]->key, "scale") == 0){
                scale = true;
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
        
        if (strcasecmp(Parameters_name2(parameters[0]), "shape") != 0) {
            Parameters* temp = parameters[0];
            parameters[0] = parameters[1];
            parameters[1] = temp;
        }
    }
    distribution_parameterization param = scale ? DISTRIBUTION_GAMMA_SHAPE_SCALE: DISTRIBUTION_GAMMA_SHAPE_RATE;
    DistributionModel* dm = new_GammaDistributionModel_with_parameters(parameters, x, param);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
    
    return model;
}
