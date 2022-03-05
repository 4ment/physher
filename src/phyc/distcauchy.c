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

double DistributionModel_log_cauchy(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	dm->lp = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
            double location = Parameters_value(dm->parameters[0], i);
            double alpha = Parameters_value(dm->parameters[1], i);
			double x = Parameters_value(dm->x, i);
			dm->lp += log(gsl_ran_cauchy_pdf(x - location, alpha));
		}
	}
	else{
        double location = Parameters_value(dm->parameters[0], 0);
		double alpha = Parameters_value(dm->parameters[1], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i);
			dm->lp += log(gsl_ran_cauchy_pdf(x - location, alpha));
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_cauchy_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
            double location = Parameters_value(dm->parameters[0], i);
			double alpha = Parameters_value(dm->parameters[1], i);
			logP += log(gsl_ran_cauchy_pdf(values[i] - location, alpha));
		}
	}
	else{
        double location = Parameters_value(dm->parameters[0], 0);
		double alpha = Parameters_value(dm->parameters[1], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(gsl_ran_cauchy_pdf(values[i] - location, alpha));
		}
	}
	return logP;
}

double DistributionModel_dlog_cauchy(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (p == Parameters_at(dm->x, i)) {
            double alpha;
            double location;
            if(Parameters_count(dm->parameters[0]) > 1){
                location = Parameters_value(dm->parameters[0], i);
                alpha = Parameters_value(dm->parameters[1], i);
            }
            else{
                location = Parameters_value(dm->parameters[0], 0);
                alpha = Parameters_value(dm->parameters[1], 0);
            }
            double x = Parameter_value(p) - location;
            return -2.0*x/(alpha*alpha*(x*x/alpha/alpha + 1.0));
		}
	}
	return 0;
}

double DistributionModel_d2log_cauchy(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
            double location = Parameters_value(dm->parameters[0], 0);
			double alpha = Parameters_value(dm->parameters[1], 0);
			if(Parameters_count(dm->parameters[0]) > 1){
                location = Parameters_value(dm->parameters[0], i);
                alpha = Parameters_value(dm->parameters[1], i);
			}
			double x = Parameter_value(p) - location;
			double alpha2 = alpha*alpha;
			return -2.0*(alpha2 - x*x)/(alpha2*alpha2 + 2.0*alpha2*x*x + x*x*x*x);
		}
	}
	return 0;
}

static void DistributionModel_cauchy_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
            double location = Parameters_value(dm->parameters[0], i);
            double alpha = Parameters_value(dm->parameters[1], i);
			samples[i] = gsl_ran_cauchy(dm->rng, alpha) + location;
		}
	}
	else{
        double location = Parameters_value(dm->parameters[0], 0);
        double alpha = Parameters_value(dm->parameters[1], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_cauchy(dm->rng, alpha) + location;
		}
	}
}

static double DistributionModel_cauchy_sample_evaluate(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
            double location = Parameters_value(dm->parameters[0], i);
            double alpha = Parameters_value(dm->parameters[1], i);
			double sample = gsl_ran_cauchy(dm->rng, alpha);
			Parameters_set_value(dm->x, i, sample + location);
			logP += log(gsl_ran_cauchy_pdf(sample - location, alpha));
		}
	}
	else{
        double location = Parameters_value(dm->parameters[0], 0);
        double alpha = Parameters_value(dm->parameters[1], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_cauchy(dm->rng, alpha);
			Parameters_set_value(dm->x, i, sample + location);
			logP += log(gsl_ran_cauchy_pdf(sample - location, alpha));
		}
	}
	return logP;
}

DistributionModel* new_CauchyDistributionModel_with_parameters(Parameters** parameters, Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, 2, x);
	dm->type = DISTRIBUTION_CAUCHY;
	dm->logP = DistributionModel_log_cauchy;
	dm->logP_with_values = DistributionModel_log_cauchy_with_values;
	dm->dlogP = DistributionModel_dlog_cauchy;
	dm->d2logP = DistributionModel_d2log_cauchy;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_cauchy_sample;
	dm->sample_evaluate = DistributionModel_cauchy_sample_evaluate;
    dm->shift = -INFINITY;
	return dm;
}

Model* new_CauchyDistributionModel_from_json(json_node* node, Hashtable* hash){
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
        parameters[1] = new_Parameters(paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            //TODO: that coming from normal
            const double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            Parameters_move(parameters[0], new_Parameter("location", m, new_Constraint(-INFINITY, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("scale", sqrt(v), new_Constraint(0, INFINITY)));
            free_Vector(samples[i]);
        }
        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        parameters = malloc(sizeof(Parameters*)*2);
        parameters[0] = new_Parameters(Parameters_count(x));
        parameters[1] = new_Parameters(Parameters_count(x));
        for (int i = 0; i < Parameters_count(x); i++) {
            Parameters_move(parameters[0], new_Parameter("location", 0, new_Constraint(-INFINITY, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("scale", 1, new_Constraint(0, INFINITY)));
        }
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "location") != 0 && strcasecmp(parameters_node->children[i]->key, "scale") != 0) {
                fprintf(stderr, "Normal distribution should be parametrized with location and scale\n");
                exit(13);
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
        
        if (strcasecmp(Parameters_name2(parameters[0]), "location") != 0) {
            Parameters* temp = parameters[0];
            parameters[0] = parameters[1];
            parameters[1] = temp;
        }
    }
    
    DistributionModel* dm = new_CauchyDistributionModel_with_parameters(parameters, x);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
    
    return model;
}
