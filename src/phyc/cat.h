// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef cat_h
#define cat_h

#include <stdio.h>

#include "mjson.h"
#include "hashtable.h"
#include "treelikelihood.h"

// How a pattern picks its rate category out of the per-rate likelihood profile.
typedef enum {
	// FastTree's rule: the category maximizing the profile (times the prior).
	// Fast, and what every published CAT result was produced with.
	CAT_ASSIGNMENT_ARGMAX,
	// Posterior-mean rate per pattern, then the categories are the K clusters
	// of those rates. The arg-max of a noisy profile is biased away from the
	// mean of the rate distribution -- the ordinary selection bias of taking a
	// maximum -- while the posterior mean shrinks by exactly the amount the
	// column is uninformative. See docs/methods/cat-improvements.md.
	CAT_ASSIGNMENT_POSTERIOR_MEAN
} cat_assignment_t;

typedef struct {
	cat_assignment_t assignment;
	// Shape of a mean-one Gamma prior on the category rate; <= 0 disables it.
	double prior_shape;
	// Number of rates the profile is built on. Only used by the posterior mean,
	// which needs a finer grid than the K categories it ends up reporting.
	int probe_count;
	int verbosity;
} CatOptions;

CatOptions cat_options_default(cat_assignment_t assignment);

void fasttree_cat(SingleTreeLikelihood* tlk, const CatOptions* options);

void cat_estimator_from_json(json_node* node, Hashtable* hash);

#endif /* cat_h */
