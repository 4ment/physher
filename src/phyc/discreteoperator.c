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