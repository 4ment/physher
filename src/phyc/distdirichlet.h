// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distdirichlet_h
#define distdirichlet_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_DirichletDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_DirichletDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distdirichlet_h */
