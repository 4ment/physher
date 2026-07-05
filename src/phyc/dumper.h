// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef dumper_h
#define dumper_h

#include <stdio.h>

#include "parameters.h"

// A Dumper serializes one or more Models back to JSON (via each Model's
// `jsonize`), unlike a logger which streams parameter/model values.
struct Dumper{
    Parameters** parameters;
    size_t parameter_count;
    Model** models;
    size_t model_count;
    void (*dump)(struct Dumper*);
    void (*free)(struct Dumper*);
    FILE* file;
    char* filename;
};

struct Dumper* new_Dumper_from_json(json_node* node, Hashtable* hash);

#endif /* dumper_h */
