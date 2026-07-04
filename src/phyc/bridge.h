// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef bridge_h
#define bridge_h

#include <stdio.h>

#include "parameters.h"

#include <gsl/gsl_rng.h>

typedef struct BridgeSampling{
	Model* model;
	Parameters* x;
	size_t burnin;
	char* file;
	char* likelihood_tag;
	gsl_rng* rng;
	void (*run)(struct BridgeSampling*);
	void (*free)(struct BridgeSampling*);
}BridgeSampling;

BridgeSampling* new_BridgeSampling_from_json(json_node* node, Hashtable* hash);

#endif /* bridge_h */
