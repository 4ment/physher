#include "boundfactory.h"

#include <strings.h>

#include "elbo.h"
#include "klpq.h"

Model* new_BoundModel_from_json(json_node* node, Hashtable* hash){
    char* boundString = get_json_node_value_string(node, "bound");
    if(strcasecmp(boundString, "elbo") == 0){
        return new_ELBO_from_json(node, hash);
    }
    else if(strcasecmp(boundString, "klpq") == 0){
        return new_KLpqBound_from_json(node, hash);
    }
    fprintf(stderr, "Unknown bound type %s\n", boundString);
    exit(1);
    return NULL;
}