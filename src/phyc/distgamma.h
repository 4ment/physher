//
//  distgamma.h
//  physher
//
//  Created by Mathieu Fourment on 30/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef distgamma_h
#define distgamma_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

// Derivative of the regularized lower incomplete gamma function P(a, x) with
// respect to the shape parameter a (exposed for testing). Used for the implicit
// reparameterization gradient of the Gamma distribution.
double gamma_p_grad_a(double a, double x);

DistributionModel* new_GammaDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization);

Model* new_GammaDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distgamma_h */
