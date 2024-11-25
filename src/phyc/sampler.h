//
//  sampler.h
//  physher
//
//  Created by Mathieu Fourment on 19/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef sampler_h
#define sampler_h

#include <stdio.h>

#include "logmcmc.h"
#include "parameters.h"

typedef struct Sampler {
    Model* model;
    Log** loggers;
    size_t logger_count;
    size_t samples;
    void (*sample)(struct Sampler*);
    void (*initialize)(struct Sampler*);
    void (*finalize)(struct Sampler*);
    void (*free)(struct Sampler*);
} Sampler;

Sampler* new_Sampler_from_json(json_node* node, Hashtable* hash);

#endif /* sampler_h */