// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
