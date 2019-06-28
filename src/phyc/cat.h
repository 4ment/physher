//
//  cat.h
//  physher
//
//  Created by Mathieu Fourment on 4/6/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef cat_h
#define cat_h

#include <stdio.h>

#include "mjson.h"
#include "hashtable.h"

void cat_estimator_from_json(json_node* node, Hashtable* hash);

#endif /* cat_h */
