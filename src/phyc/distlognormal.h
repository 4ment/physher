// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distlognormal_h
#define distlognormal_h

#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"
#include "distmodel.h"

DistributionModel* new_LogNormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization);

Model* new_LogNormalDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distlognormal_h */
