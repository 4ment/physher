//
//  distbeta.h
//  physher
//
//  Created by mathieu on 17/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef distbeta_h
#define distbeta_h

#include "mjson.h"
#include "hashtable.h"
#include "parameters.h"
#include "distmodel.h"

DistributionModel* new_BetaDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_BetaDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distbeta_h */
