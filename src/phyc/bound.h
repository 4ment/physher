//
//  bound.h
//  physher
//
//  Created by Mathieu Fourment on 18/6/24.
//  Copyright © 2024 Mathieu Fourment. All rights reserved.
//

#ifndef bound_h
#define bound_h

#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"

typedef struct Bound {
    Model* joint;
    Model* variational;  // DistributionModel or CompoundModel of DistributionModels
    size_t samples;
    size_t kSamples;
    bool entropy;
    double (*logP)(struct Bound*);
    double (*gradient)(struct Bound*, Parameters* parameters);
    Parameters* parameters;  // parameters of the joint model approximated by
                             // variational distribution
} Bound;

Model* new_BoundModel(const char* name, Bound* bound);

Model* new_AbstractBoundModel_from_json(json_node* node, Hashtable* hash);

#endif /* bound_h */