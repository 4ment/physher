//
//  distnormal.h
//  physher
//
//  Created by mathieu on 24/1/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef distnormal_h
#define distnormal_h

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

Model* new_NormalDistributionModel_from_json(json_node* node, Hashtable* hash);

Model* new_HalfNormalDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distnormal_h */
