// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef cat_h
#define cat_h

#include <stdbool.h>
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

// Where the prior over the rate grid comes from. Both assignment rules weight the
// per-pattern likelihood profile by it, so it is the one thing they share.
typedef enum {
	// Picked in advance and the same for every alignment: a mean-one Gamma of
	// shape `prior_shape` when that is positive, flat over the grid when it is
	// not. FastTree's Gamma(3, 1/3) is this with prior_shape = 3.
	CAT_PRIOR_FIXED,
	// Estimated from the alignment instead of assumed: the distribution over the
	// grid rates that maximizes the marginal likelihood of the whole profile, with
	// no parametric family imposed on it (Kiefer & Wolfowitz 1956). Costs no tree
	// traversals -- the profile the assignment already needs is the only input --
	// and `prior_shape` is unused.
	CAT_PRIOR_NPMLE
} cat_prior_t;

typedef struct {
	cat_assignment_t assignment;
	cat_prior_t prior;
	// Shape of a mean-one Gamma prior on the category rate; <= 0 disables it.
	// Ignored unless prior is CAT_PRIOR_FIXED.
	double prior_shape;
	// Number of rates the profile is built on. Only used by the posterior mean,
	// which needs a finer grid than the K categories it ends up reporting.
	int probe_count;
	// Cap on the EM passes the NPMLE is allowed, and the Kiefer--Wolfowitz gap it
	// stops at. Unused unless prior is CAT_PRIOR_NPMLE.
	int npmle_iterations;
	double npmle_tolerance;
	int verbosity;
} CatOptions;

CatOptions cat_options_default(cat_assignment_t assignment);

// Outcome of one fasttree_cat call.
typedef struct CatResult {
	double logP;      // log-likelihood the call left behind
	// Full likelihood traversals performed, so a caller keeping an evaluation
	// budget can charge for the work.
	int evaluations;
	// The new assignment scored worse than the one that came in and the old one
	// was put back: the call changed nothing.
	bool reverted;
	// How the estimated prior turned out, when there was one (zeroed otherwise).
	// `npmle_atoms` is the number of grid rates it left carrying real mass: the
	// nonparametric MLE of a mixing distribution is discrete, but it is approached
	// from a dense start, so this counts what a bounded iteration got round to
	// discarding as much as it counts the true support. The reading that does not
	// depend on that is a collapse onto one rate, which says no mixture beats
	// rescaling the tree. `npmle_gap` is the Kiefer--Wolfowitz optimality gap the
	// iteration stopped at -- zero at the maximum -- and `npmle_passes` how many
	// passes that took.
	int npmle_atoms;
	int npmle_passes;
	double npmle_gap;
} CatResult;

// Read the CAT keys ("assignment", "prior", "probe", "npmle_iterations",
// "npmle_tolerance", "verbosity") off `node`, on top of the defaults the assignment
// rule implies. Shared by the standalone "cat" estimator and by "algorithm": "cat"
// in the optimizer, so the two spell their options identically.
CatOptions cat_options_from_json(json_node* node);

// The estimator rewrites sm->site_category, which only the empirical CAT site
// model has. Dies at `node` if `tlk` has not got one.
void cat_check_sitemodel(json_node* node, const SingleTreeLikelihood* tlk);

// Reassign every pattern to a rate category and reset the category rates from the
// per-pattern likelihood profile. Guarded: the likelihood is taken before and
// after, and the incoming assignment is restored if the new one scores worse, so
// the call can only improve the model or leave it alone.
CatResult fasttree_cat(SingleTreeLikelihood* tlk, const CatOptions* options);

void cat_estimator_from_json(json_node* node, Hashtable* hash);

#endif /* cat_h */
