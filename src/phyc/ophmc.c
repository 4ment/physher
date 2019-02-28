//
//  ophmc.c
//  physher
//
//  Created by Mathieu Fourment on 11/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "ophmc.h"

#include "simplex.h"
#include "matrix.h"
#include "utilsgsl.h"
#include "gaussian.h"

void operator_hmc_store(Operator* op){
	Parameters* ps = op->x;
	if(op->model_count > 1){
		ps = ((Simplex*)op->models[1]->obj)->parameters;
	}
	for(int i = 0; i < Parameters_count(ps); i++){
		Parameter_store(Parameters_at(ps, i));
	}
}


void operator_hmc_restore(Operator* op){
	Parameters* ps = op->x;
	if(op->model_count > 1){
		ps = ((Simplex*)op->models[1]->obj)->parameters;
	}
	for(int i = 0; i < Parameters_count(ps); i++){
		Parameter_restore(Parameters_at(ps, i));
	}
}

bool operator_hmc(Operator* op, double* logHR){
	op->indexes[0] = gsl_rng_uniform_int(op->rng, Parameters_count(op->x));
	Parameter* p = Parameters_at(op->x, op->indexes[0]);
	Model* posterior = op->models[0];
	double p0 = rnorm();
	double x0 = Parameter_value(p);
	double delta = op->parameters[0];
	int steps = op->parameters[1];
	
	double dU = -posterior->dlogP(posterior, p);
	double pStep = p0 - delta/2.0* dU;
 
	// Full step
	double xStep = x0 + delta*pStep;
	
	for (int i = 0; i < steps; i++) {
		// Update momentum
		Parameter_set_value(p, xStep);
		double dU  = -posterior->dlogP(posterior, p);
		pStep = pStep - delta*dU;
		
		// Update position
		xStep = xStep + delta*pStep;
		if (xStep < 0) {
			Parameter_set_value(p, x0);
			return false;
		}
//			printf("%f %f\n", x0, xStep);
	}
//		printf("\n");
	Parameter_set_value(p, xStep);
	dU  = -posterior->dlogP(posterior, p);
	pStep = pStep - delta/2*dU;
	Parameter_set_value(p, xStep);
	*logHR = p0*p0/2 - pStep*pStep/2;
	return true;
}

Operator* new_HMCOperator_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"algorithm",
		"coalescent",
		"delay",
		"model",
		"parameters",
		"tree",
		"weight",
		"x"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	Operator* op = malloc(sizeof(Operator));
	const char* id_string = get_json_node_value_string(node, "id");
	op->weight = get_json_node_value_double(node, "weight", 1);
	op->name = String_clone(id_string);
	
	op->x = new_Parameters(1);
	get_parameters_references2(node, hash, op->x, "x");
	char* ref = get_json_node_value_string(node, "model");
	op->models = malloc(sizeof(Model*));
	// posterior model
	op->models[0] = Hashtable_get(hash, ref+1);
	op->models[0]->ref_count++;
	op->model_count = 1;
	
	op->propose = operator_hmc;
	op->store = operator_hmc_store;
	op->restore = operator_hmc_restore;
	op->optimize = NULL;
	op->parameters = dvector(2);
	op->parameters[0] = get_json_node_value_double(node, "delta", 0.0001);
	op->parameters[1] = get_json_node_value_double(node, "steps", 35);
	op->indexes = ivector(3);
	
	op->rejected_count = 0;
	op->accepted_count = 0;
	op->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return op;
}
