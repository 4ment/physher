// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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