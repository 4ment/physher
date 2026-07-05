// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef sampler_h
#define sampler_h

#include <stdio.h>

#include "tracelogger.h"
#include "parameters.h"

typedef struct Sampler {
    Model* model;
    Trace** loggers;
    size_t logger_count;
    size_t samples;
    void (*sample)(struct Sampler*);
    void (*initialize)(struct Sampler*);
    void (*finalize)(struct Sampler*);
    void (*free)(struct Sampler*);
} Sampler;

Sampler* new_Sampler_from_json(json_node* node, Hashtable* hash);

#endif /* sampler_h */