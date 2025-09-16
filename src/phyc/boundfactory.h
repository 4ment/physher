#include "parameters.h"
#include "mjson.h"
#include "hashtable.h"

Model* new_BoundModel_from_json(json_node* node, Hashtable* hash);