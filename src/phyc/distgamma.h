//
//  distgamma.h
//  physher
//
//  Created by Mathieu Fourment on 30/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef distgamma_h
#define distgamma_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

Model* new_GammaDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distgamma_h */
