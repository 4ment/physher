//
//  distoneonx.h
//  physher
//
//  Created by mathieu on 9/7/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef distoneonx_h
#define distoneonx_h

#include "distmodel.h"
#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"

DistributionModel* new_OneOnXDistributionModel(Parameters* x);

Model* new_OneOnXDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distoneonx_h */
