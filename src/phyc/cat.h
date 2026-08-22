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
	// column is uninformative. See docs/methods/cat.md.
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
	// Log-likelihood of the mixture the estimated prior defines: the grid rates
	// weighted by the estimate, with the category label summed out rather than
	// selected. Unlike `logP` this is a likelihood and comparable with one,
	// though of a richer model -- `npmle_atoms` components rather than K -- and it
	// is taken at the tree the assignment saw. Zero unless prior is CAT_PRIOR_NPMLE.
	double npmle_logP;
	// Read off the per-pattern posterior over the grid the profile was built on --
	// the K categories under CAT_ASSIGNMENT_ARGMAX, the probe grid under
	// CAT_ASSIGNMENT_POSTERIOR_MEAN -- and weighted by pattern site counts.
	// `grid` says which, so `entropy` can be read against its ceiling log(grid).
	// Together they say whether the grid is finer than the columns can resolve
	// before anything else is fitted: a `confidence` near 1/grid and an `entropy`
	// near its ceiling means no column knows its own rate.
	int grid;
	// Site-weighted mean of max_j pi_ij. Under the arg-max rule that maximum is
	// the selected category, so this is pi_bar(khat) exactly.
	double confidence;
	// Site-weighted mean of H(pi_i), in nats, between 0 and log(grid).
	double entropy;
	// Plug-in mutual information between a column and its rate: H(pi_bar) minus
	// the mean entropy above, where pi_bar is the site-weighted aggregate
	// posterior. Equals the log(grid) - entropy of docs/methods/cat.md when that
	// aggregate is flat.
	double information;
} CatResult;

// What a CAT fit is worth as a mixture, all taken at the tree, rates and
// assignment the model currently holds. See docs/methods/cat.md, "What you may
// compare, and what you may not".
//
// The label that is summed out here is the per-pattern rate category, not a
// parameter: nothing is integrated over a prior and this is not the model evidence
// the `bridgesampling`, `nest`, `is` and `laplace` actions estimate. What comes out
// is an ordinary log P(D | theta) -- the same quantity every non-CAT site model
// reports as "treelikelihood", and marginal over its categories in exactly the way
// a discrete Gamma's is.
typedef struct CatMixture {
	// The number that is comparable with a Gamma, free-rate or single-rate score:
	// the log-likelihood of the K-component free-rate mixture whose rates are the
	// category rates and whose weights are the site-weighted assignment
	// frequencies. That mixture has mean rate one by construction, because the
	// category rates were normalised against the same assignment, so it scores the
	// same branch lengths on the same scale -- it is the model CAT fitted, with
	// the label summed out instead of selected.
	double logP_mixture;
	// The CAT score, sum_i n_i log f_{i khat_i}: what "treelikelihood" reports and
	// what must not be compared with anything.
	double logP_cat;
	// sum_i n_i (log f_{i khat_i} + log p_{khat_i}), the mixture's ELBO at a point
	// mass on the selected category. A lower bound on logP_mixture for *any*
	// assignment, so it costs nothing and cannot be wrong.
	double logP_bound;
	// logP_cat - logP_mixture: the pointwise mutual information the CAT score is
	// credited with and the mixture is not. Positive whenever the assignment is
	// the mixture's own arg-max, which is not guaranteed under a prior or under
	// the posterior-mean rule.
	double gap;
	// Categories carrying at least one site: the mixture has this many components,
	// the rest having weight zero and dropping out of it.
	//
	// 2*used - 2 (rates plus weights, less the simplex and the unit-mean
	// constraints) is an *upper bound* on the free parameters, not the k an
	// information criterion wants, and on real data it is a loose one. Two reasons,
	// both measured in docs/methods/cat.md, "Are they free parameters?":
	//
	// - The values are not maximum-likelihood estimates of this mixture. The rates
	//   are cluster centres of posterior means and the weights are site counts of a
	//   hard assignment; neither was chosen to maximise the number reported here.
	//   Fitting the same mixture by ML gains 6.8 log units on fluA, so k is
	//   charging for optimisation that did not happen.
	// - The components are not distinct. The quantizer stacks categories on the
	//   same rate, and the likelihood is flat along the directions that move weight
	//   between two of them, so the Fisher information is singular and those
	//   directions are not parameters at all. On fluA at K = 20, 15 occupied
	//   categories sit on 6 distinct rates, and the free-rate family saturates at
	//   3 -- an identified dimension of 4, against a nominal 28.
	//
	// The identified dimension is found by fitting "distribution": "discrete" at
	// increasing K until the likelihood stops moving. There is no way to read it
	// off this struct, which is why the number here is reported as a bound.
	int used;
	// Full likelihood traversals performed: one per category, plus the two that
	// bracket them.
	int evaluations;
} CatMixture;

// Sum the category label out of the current CAT fit, so its likelihood can be
// compared with any other site model's. Costs K + 2 traversals and leaves the model
// exactly as it found it -- assignment, category rates and partials all restored --
// so it can be run between optimizer actions. `tlk` must carry a CAT site model.
//
// At `verbosity > 0` the mixture is also printed component by component: rate and
// weight per category, which is what another program -- or physher's own
// "distribution": "discrete" site model -- needs to be handed the same model.
CatMixture cat_mixture(SingleTreeLikelihood* tlk, int verbosity);

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
