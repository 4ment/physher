//
//  gmrf.h
//  physher
//
//  Created by Mathieu Fourment on 18/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef gmrf_h
#define gmrf_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"

Model* new_GMRFModel_from_json(json_node* node, Hashtable* hash);

#endif /* gmrf_h */
