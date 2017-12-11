//
//  ophmc.h
//  physher
//
//  Created by Mathieu Fourment on 11/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef ophmc_h
#define ophmc_h

#include <stdio.h>

#include "operator.h"

Operator* new_HMCOperator_from_json(json_node* node, Hashtable* hash);

#endif /* ophmc_h */
