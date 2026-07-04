// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distnormal_h
#define distnormal_h

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_NormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization);

DistributionModel* new_HalfNormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization);

Model* new_NormalDistributionModel_from_json(json_node* node, Hashtable* hash);

Model* new_HalfNormalDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distnormal_h */
