// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distexp_h
#define distexp_h

#include <stdio.h>

#include "parameters.h"
#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_ExponentialDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization);

Model* new_ExponentialDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distexp_h */
