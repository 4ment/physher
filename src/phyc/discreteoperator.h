//
//  discreteoperator.h
//  physher
//
//  Created by Mathieu Fourment on 1/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef discreteoperator_h
#define discreteoperator_h

#include <stdio.h>

#include "operator.h"

#include <gsl/gsl_rng.h>

bool operator_discrete_bitflip(Operator* op, double* logHR);

bool operator_discrete_exchange(Operator* op, double* logHR);

#endif /* discreteoperator_h */
