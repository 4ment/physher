//
//  distlognormal.h
//  physher
//
//  Created by mathieu on 17/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef distlognormal_h
#define distlognormal_h

#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"
#include "distmodel.h"

DistributionModel* new_LogNormalDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_LogNormalDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distlognormal_h */
