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
	
	dm->lp = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters[0], i);
			double b = Parameters_value(dm->parameters[1], i);
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(a*b*pow(x, a - 1.0) * pow( 1.0 - pow(x, a), b - 1.0));
		}
	}
	else{
		double a = Parameters_value(dm->parameters[0], 0);
		double b = Parameters_value(dm->parameters[1], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i) - dm->shift;
			dm->lp += log(a*b*pow(x, a - 1.0) * pow( 1.0 - pow(x, a), b - 1.0));
		}
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_kumaraswamy_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters[0], i);
			double b = Parameters_value(dm->parameters[1], i);
			logP += log(a*b*pow(values[i], a - 1.0) * pow( 1.0 - pow(values[i], a), b - 1.0));
		}
	}
	else{
		double a = Parameters_value(dm->parameters[0], 0);
		double b = Parameters_value(dm->parameters[1], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(a*b*pow(values[i], a - 1.0) * pow( 1.0 - pow(values[i], a), b - 1.0));
		}
	}
	return logP;
}


double DistributionModel_dlog_kumaraswamy(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (p == Parameters_at(dm->x,i)) {
			double a;
			double b;
			if(Parameters_count(dm->parameters[0]) > 1){
				a = Parameters_value(dm->parameters[0], i);
				b = Parameters_value(dm->parameters[1], i);
			}
			else{
				a = Parameters_value(dm->parameters[0], 0);
				b = Parameters_value(dm->parameters[1], 0);
			}
			double x = Parameters_value(dm->x, i);
			return (a - 1.0)/x + (b - 1.0)*a*pow(x, a - 1.0)/(pow(x, a) - 1.0);
		}
	}
	return 0;
}

// F_X(x) = p(X <= x)
double DistributionModel_kumaraswamy_inverse_CDF(double p, double a, double b){
	return pow(1.0 - pow(1.0 - p, 1.0/b), 1.0/a);
}

static void DistributionModel_kumaraswamy_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters[0], i);
			double b = Parameters_value(dm->parameters[1], i);
			samples[i] = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
		}
	}
	else{
		double a = Parameters_value(dm->parameters[0], 0);
		double b = Parameters_value(dm->parameters[1], 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
		}
	}
}

static double DistributionModel_kumaraswamy_sample_evaluate(DistributionModel* dm){
	if(Parameters_count(dm->parameters[0]) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double a = Parameters_value(dm->parameters[0], i);
			double b = Parameters_value(dm->parameters[1], i);
			double sample = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		double a = Parameters_value(dm->parameters[0], 0);
		double b = Parameters_value(dm->parameters[1], 0);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = DistributionModel_kumaraswamy_inverse_CDF(gsl_ran_flat(dm->rng, 0, 1), a, b);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_log_kumaraswamy(dm);
}

DistributionModel* new_KumaraswamyDistributionModel_with_parameters(Parameters** parameters, Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, 2, x);
	dm->type = DISTRIBUTION_KUMARASWAMY;
	dm->parameterization = 0;
	
	dm->logP = DistributionModel_log_kumaraswamy;
	dm->logP_with_values = DistributionModel_log_kumaraswamy_with_values;
	dm->dlogP = DistributionModel_dlog_kumaraswamy;
	dm->sample = DistributionModel_kumaraswamy_sample;
	dm->sample_evaluate = DistributionModel_kumaraswamy_sample_evaluate;
	
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
    dm->shift = 0;
	return dm;
}

Model* new_KumaraswamyDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    
    char* file = get_json_node_value_string(node, "file");
    if(file != NULL){
           error("Cannot estimate parameters of the Kumaraswamy distribution from a sample\n");
           exit(3);
       }
    
    Parameters** parameters = NULL;
    size_t parameters_dim = 2;

    if(get_json_node(node, "parameters") == NULL){
        parameters = malloc(sizeof(Parameters*)*2);
        parameters[0] = new_Parameters(Parameters_count(x));
        parameters[1] = new_Parameters(Parameters_count(x));
        for (int i = 0; i < Parameters_count(x); i++) {
            Parameters_move(parameters[0], new_Parameter("a", 0, new_Constraint(0, INFINITY)));
            Parameters_move(parameters[1], new_Parameter("b", 1, new_Constraint(0, INFINITY)));
        }
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        for (int i = 0; i < parameters_node->child_count; i++) {
            if (strcasecmp(parameters_node->children[i]->key, "a") != 0 && strcasecmp(parameters_node->children[i]->key, "b") != 0) {
                fprintf(stderr, "Kumaraswamy distribution should be parametrized with a and b\n");
                exit(13);
            }
        }
        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
        
        if (strcasecmp(Parameters_name2(parameters[0]), "a") != 0) {
            Parameters* temp = parameters[0];
            parameters[0] = parameters[1];
            parameters[1] = temp;
        }
    }
    
    DistributionModel* dm = new_KumaraswamyDistributionModel_with_parameters(parameters, x);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    free_Parameters(parameters[0]);
    free_Parameters(parameters[1]);
    free(parameters);
//    free_Parameters(x);
    
    return model;
}
