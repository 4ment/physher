// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distcauchy_h
#define distcauchy_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_CauchyDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_CauchyDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distcauchy_h */
