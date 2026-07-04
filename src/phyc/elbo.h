// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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