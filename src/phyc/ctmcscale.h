//
//  ctmcscale.h
//  physher
//
//  Created by mathieu on 31/3/21.
//  Copyright © 2021 Mathieu Fourment. All rights reserved.
//

#ifndef ctmcscale_h
#define ctmcscale_h

#include "distmodel.h"

Model* new_CTMCScaleModel_from_json(json_node* node, Hashtable* hash);

DistributionModel* new_CTMCScale_with_parameters(Parameters* x, Tree* tree);

Model* new_CTMCScaleModel(const char* name, DistributionModel* dm, Model* tree);

size_t DistributionModel_initialize_gradient(Model *self, int flags);
static void _calculate_height_gradient(Tree* tree, double rate, double shape, double totalTreeTime, double* gradient);

double* DistributionModel_gradient(Model *self);

#endif /* ctmcscale_h */
