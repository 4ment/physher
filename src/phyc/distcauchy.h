//
//  distcauchy.h
//  physher
//
//  Created by Mathieu Fourment on 8/6/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef distcauchy_h
#define distcauchy_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_CauchyDistributionModel_with_parameters(Parameters** parameters, Parameters* x);

Model* new_CauchyDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distcauchy_h */
