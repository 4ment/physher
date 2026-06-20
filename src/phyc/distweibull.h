//
//  distweibull.h
//  physher
//
//  Created by mathieu on 20/6/26.
//  Copyright © 2026 Mathieu Fourment. All rights reserved.
//

#ifndef distweibull_h
#define distweibull_h

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_WeibullDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_WeibullDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distweibull_h */
