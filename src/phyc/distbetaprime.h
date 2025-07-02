//
//  distbetaprime.h
//  physher
//
//  Created by mathieu on 29/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef distbetaprime_h
#define distbetaprime_h

#include "mjson.h"
#include "hashtable.h"
#include "parameters.h"
#include "distmodel.h"

DistributionModel* new_BetaPrimeDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_BetaPrimeDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distbetaprime_h */
