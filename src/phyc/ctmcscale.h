// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef ctmcscale_h
#define ctmcscale_h

#include "distmodel.h"

Model* new_CTMCScaleModel_from_json(json_node* node, Hashtable* hash);

DistributionModel* new_CTMCScale_with_parameters(Parameters* x, Tree* tree);

Model* new_CTMCScaleModel(const char* name, DistributionModel* dm, Model* tree);

static void _calculate_height_gradient(Tree* tree, double rate, double shape, double totalTreeTime, double* gradient);

void CTMCModel_gradient(Model *self, int flags, double* gradient);

#endif /* ctmcscale_h */
