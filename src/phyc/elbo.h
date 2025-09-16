//
//  klqp.h
//  physher
//
//  Created by Mathieu Fourment on 18/6/24.
//  Copyright © 2024 Mathieu Fourment. All rights reserved.
//

#ifndef elbo_h
#define elbo_h

#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"

typedef struct ELBO {
    Model* joint;
    Model* variational;  // DistributionModel or CompoundModel of DistributionModels
    size_t samples;
    size_t kSamples;
    bool entropy;
    double (*logP)(struct ELBO*);
    double (*gradient)(struct ELBO*, Parameters* parameters);
    Parameters* parameters;  // parameters of the joint model approximated by
                             // variational distribution
} ELBO;

Model* new_ELBO_from_json(json_node* node, Hashtable* hash);

#endif /* elbo_h */