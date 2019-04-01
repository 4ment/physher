//
//  discreteoperator.c
//  physher
//
//  Created by Mathieu Fourment on 1/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "discreteoperator.h"

#include "discreteparameter.h"

bool operator_discrete_bitflip(Operator* op, double* logHR){
	Model* mdp = op->models[0];
	DiscreteParameter* dp = mdp->obj;
	
	size_t index = gsl_rng_uniform_int(op->rng, dp->length);
	unsigned value = dp->values[index];
	dp->set_value(dp, index, 1-value);
	*logHR = 0;
	return true;
}

bool operator_discrete_exchange(Operator* op, double* logHR){
	Model* mdp = op->models[0];
	DiscreteParameter* dp = mdp->obj;
	size_t index1 = gsl_rng_uniform_int(op->rng, dp->length);
	size_t index2 = index1;
	while (index1 == index2) {
		index2 = gsl_rng_uniform_int(op->rng, dp->length);
	}
	
	int v1 = dp->values[index1];
	int v2 = dp->values[index2];
	
	int d = gsl_rng_uniform_int(op->rng, (int)round(op->parameters[0])) + 1;
	
	v1 = round(v1 - d);
	v2 = round(v2 + d);
	if (v1  <= 0 || v2 <=0) {
		return false;
	}

	if (v1+v2 != dp->values[index1]+dp->values[index2]) {
		fprintf(stderr, "discrete exchange failed from %d %d to %d %d\n", dp->values[index1], dp->values[index2], v1, v2);
		exit(2);
	}
	dp->set_value(dp, index1, v1);
	dp->set_value(dp, index2, v2);
	
	*logHR = 0;
	return true;
}