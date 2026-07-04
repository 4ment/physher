// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distbeta_h
#define distbeta_h

#include "mjson.h"
#include "hashtable.h"
#include "parameters.h"
#include "distmodel.h"

DistributionModel* new_BetaDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_BetaDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distbeta_h */
