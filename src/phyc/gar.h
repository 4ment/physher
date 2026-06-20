//
//  gar.h
//  physher
//
//  Created by Mathieu Fourment on 3/11/2025.
//  Copyright © 2025 Mathieu Fourment. All rights reserved.
//

#ifndef gar_h
#define gar_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"

Model* new_GARModel_from_json(json_node* node, Hashtable* hash);

#endif /* gar_h */
