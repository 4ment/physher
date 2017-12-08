//
//  opvb.h
//  physher
//
//  Created by Mathieu Fourment on 7/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef opvb_h
#define opvb_h

#include <stdio.h>

#include "operator.h"

Operator* new_VariationalOperator_from_json(json_node* node, Hashtable* hash);

#endif /* opvb_h */
