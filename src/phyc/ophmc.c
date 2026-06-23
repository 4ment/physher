//
//  ophmc.c
//  physher
//
//  Created by Mathieu Fourment on 11/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "ophmc.h"

#include <string.h>

#include "matrix.h"
#include "utilsgsl.h"
#include "gaussian.h"

void operator_hmc_optimize(Operator* op, double logAlpha){
	long count, accepted;
	bool useAcceptanceRate = true;
	if(operator_tuning_stats(op, &count, &accepted) && (!useAcceptanceRate || count >= 10)){
		double prob = useAcceptanceRate ? (double)accepted / count : exp(logAlpha);
		double oldStepSize = op->parameters[0];
		double newStepSize = log(oldStepSize) + (prob - op->target) / (2 + count);
		op->parameters[0] = exp(newStepSize);
	}
}

bool operator_hmc(Operator* op, double* logHR){
	size_t paramCount = Parameters_count(op->x);
	size_t dim = Parameters_size(op->x);
	Model* posterior = op->models[0];
	double stepSize = op->parameters[0];
	int steps = op->parameters[1];
	
	double* momentum0 = dvector(dim);
	double* momentum = dvector(dim);
	double* position = dvector(dim);

	for(size_t i = 0; i < dim; i++){
		momentum0[i] = rnorm();
	}
	memcpy(momentum, momentum0, dim*sizeof(double));
	size_t offset = 0;
	for(size_t i = 0; i < paramCount; i++){
		Parameter* p = Parameters_at(op->x, i);
		const double* values = Parameter_values(p);
		memcpy(position + offset, values, Parameter_size(p)*sizeof(double));
		offset += Parameter_size(p);
	}

	// const double U0 = -posterior->logP(posterior);
	Parameters_zero_grad(op->x);
	posterior->gradient(posterior, op->x);
	offset = 0;
	for(size_t i = 0; i < paramCount; i++){
		Parameter* p = Parameters_at(op->x, i);
		for(size_t j = 0; j < Parameter_size(p); j++){
			double dU = -p->grad[j];
			momentum[offset] -= stepSize/2.0 * dU;
			offset++;
		}
	}

	for (size_t s = 0; s < steps; s++) {
		for(size_t i = 0; i < dim; i++){
			position[i] += stepSize * momentum[i];
		}
		Parameters_set_values(op->x, position);

		Parameters_zero_grad(op->x);
		posterior->gradient(posterior, op->x);
		offset = 0;
		for(size_t i = 0; i < paramCount; i++){
			Parameter* p = Parameters_at(op->x, i);
			for(size_t j = 0; j < Parameter_size(p); j++){
				double dU = -p->grad[j];
				momentum[offset] -= stepSize * dU;
				offset++;
			}
		}
	}

	offset = 0;
	for(size_t i = 0; i < paramCount; i++){
		Parameter* p = Parameters_at(op->x, i);
		for(size_t j = 0; j < Parameter_size(p); j++){
			double dU = -p->grad[j];
			momentum[offset] += stepSize/2.0 * dU;
			offset++;
		}
	}

	double K0 = 0.0;
	double K1 = 0.0;
    for(size_t i = 0; i < dim; i++){
        K0 += momentum0[i] * momentum0[i];
        K1 += momentum[i] * momentum[i];
    }
    K0 /= 2.0;
    K1 /= 2.0;
	// const double U1 = -posterior->logP(posterior);
	// *logHR = (U0 + K0) - (U1 + K1);
	*logHR = K0 - K1;

	free(momentum0);
	free(momentum);
	free(position);

	return true;
}

Operator* new_HMCOperator_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"algorithm",
		"coalescent",
		"delay",
		"model",
		"parameters",
		"stepsize",
		"steps",
		"target",
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
	op->optimize = operator_hmc_optimize;
	op->parameters = dvector(2);
	op->parameters[0] = get_json_node_value_double(node, "stepsize", 0.01);
	op->parameters[1] = get_json_node_value_double(node, "steps", 5);
	
	op->rejected_count = 0;
	op->accepted_count = 0;
	op->failure_count = 0;
	op->tuning_delay = get_json_node_value_size_t(node, "delay", 0);
	op->accepted_at_delay = 0;
	op->count_at_delay = 0;
	op->tuning_started = false;
	op->all = false;
	op->target = get_json_node_value_double(node, "target", 0.8);
	op->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return op;
}
