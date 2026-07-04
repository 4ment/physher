// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef is_h
#define is_h

#include <stdio.h>
#include "parameters.h"

#include <gsl/gsl_rng.h>

typedef struct ImportanceSampler{
	Model* model; // distribution of interest
	Model** distribution; // importance distribution or mixture of distribution
	double* weights;
	size_t distribution_count;
	size_t samples;
	Parameters* parameters;
	bool normalize;
	double(*calculate)(struct ImportanceSampler*);
	void(*free)(struct ImportanceSampler*);
	gsl_rng* rng;
}ImportanceSampler;

ImportanceSampler* new_ImportanceSampler_from_json(json_node* node, Hashtable* hash);

#endif /* is_h */
