#ifndef _DIST_MODEL_FACTORY_H_
#define _DIST_MODEL_FACTORY_H_

#include "parameters.h"

Model* new_DistributionModel_from_json(json_node* node, Hashtable* hash);

#endif