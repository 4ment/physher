//
//  klpq.h
//  physher
//
//  Created by Mathieu Fourment on 8/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef klpq_h
#define klpq_h

#include <stdio.h>

#include "mjson.h"
#include "hashtable.h"
#include "parameters.h"

Model* new_KLpqBound_from_json(json_node* node, Hashtable* hash);

#endif /* klpq_h */
