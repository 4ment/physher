// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "cat.h"

#include "parameters.h"
#include "treelikelihood.h"
#include "utils.h"
#include "matrix.h"


#define CAT_PROBE_DEFAULT 20
#define CAT_LLOYD_ITERATIONS 100

CatOptions cat_options_default(cat_assignment_t assignment){
	CatOptions options;
	options.assignment = assignment;
	// FastTree hard-codes a Gamma(3, 1/3) prior and the arg-max rule was tuned
	// against it, so that pairing is kept. The posterior mean already shrinks by
	// the right amount on its own, and stacking the prior on top of it shrinks
	// twice: it costs 10 % of the tree length on large trees, so the default
	// there is a flat prior over the probe grid.
	options.prior_shape = (assignment == CAT_ASSIGNMENT_ARGMAX ? 3.0 : 0.0);
	options.probe_count = (assignment == CAT_ASSIGNMENT_ARGMAX ? 0 : CAT_PROBE_DEFAULT);
	options.verbosity = 0;
	return options;
}

// log of a mean-one Gamma(shape, 1/shape) density, up to a constant. The constant
// is the same for every category and only differences are ever used.
static void _cat_log_prior(double* log_prior, const double* rates, int count,
                           double shape){
	for (int i = 0; i < count; i++) {
		log_prior[i] = (shape - 1.0)*log(rates[i]) - shape*rates[i];
	}
}

// Fill likelihoods[count x P] with the per-pattern log-likelihood at each rate.
// Called with sm->cat_count forced to 1, so every pattern is scored against
// matrix block 0 and one traversal per rate is enough.
static void _cat_probe(SingleTreeLikelihood* tlk, const double* rates, int count,
                       double* likelihoods){
	SiteModel* sm = tlk->sm;
	for (int i = 0; i < count; i++) {
		sm->cat_rates[0] = rates[i];
		SingleTreeLikelihood_update_all_nodes(tlk);
		tlk->calculate(tlk);
		memcpy(likelihoods + (i*tlk->sp->count), tlk->pattern_lk,
		       sizeof(double)*tlk->sp->count);
	}
}

typedef struct {
	double value;
	double weight;
} _cat_point;

static int _cat_compare_points(const void* a, const void* b){
	double x = ((const _cat_point*)a)->value;
	double y = ((const _cat_point*)b)->value;
	return (x > y) - (x < y);
}

// Weighted Lloyd (k-means) in one dimension, seeded at the weighted quantiles.
// The rate estimates are already sorted along a single axis, so this converges in
// a handful of passes and its optimum is the least-squares K-point quantizer of
// the estimated rate distribution -- the empirical counterpart of the fixed
// 1/K..K grid it replaces.
static void _cat_quantize(const double* x, const double* weights, int n,
                          double* centres, int* labels, int count){
	_cat_point* sorted = malloc(sizeof(_cat_point)*n);
	for (int i = 0; i < n; i++) {
		sorted[i].value = x[i];
		sorted[i].weight = weights[i];
	}
	qsort(sorted, n, sizeof(_cat_point), _cat_compare_points);

	// Seed at the weighted quantiles of the estimates, so the categories start out
	// where the sites are. Several seeds can land on the same value -- identical
	// columns share a profile and therefore an estimate, and a real alignment puts
	// a large atom on one of them -- and the duplicates are nudged apart below
	// rather than spread over the support. Seeding one category per *distinct*
	// value instead is the obvious alternative and measures worse: it costs 6 % of
	// the tree length at 8-16 taxa, because it spends categories on a sparse tail
	// that carries almost no sites (docs/methods/cat-improvements.md).
	double total = 0;
	for (int i = 0; i < n; i++) total += sorted[i].weight;
	double accumulated = 0;
	int placed = 0;
	for (int i = 0; i < n && placed < count; i++) {
		accumulated += sorted[i].weight;
		if (accumulated >= total*(placed + 0.5)/count) centres[placed++] = sorted[i].value;
	}
	// Fewer usable quantiles than categories: the spares sit just above the last one
	// so the centres stay ordered, stay distinct enough for Lloyd to move them
	// apart, and no category is left uninitialised.
	while (placed < count) {
		centres[placed] = centres[placed - 1] + 1e-3;
		placed++;
	}
	free(sorted);

	for (int i = 0; i < n; i++) labels[i] = -1;
	for (int iteration = 0; iteration < CAT_LLOYD_ITERATIONS; iteration++) {
		bool moved = false;
		for (int i = 0; i < n; i++) {
			int best = 0;
			double best_distance = fabs(x[i] - centres[0]);
			for (int j = 1; j < count; j++) {
				double distance = fabs(x[i] - centres[j]);
				if (distance < best_distance) {
					best_distance = distance;
					best = j;
				}
			}
			if (labels[i] != best) {
				labels[i] = best;
				moved = true;
			}
		}
		if (!moved) break;
		for (int j = 0; j < count; j++) {
			double sum_weight = 0;
			double sum_value = 0;
			for (int i = 0; i < n; i++) {
				if (labels[i] == j) {
					sum_weight += weights[i];
					sum_value += weights[i]*x[i];
				}
			}
			// An empty cluster keeps its centre: dropping it would renumber the
			// categories and the rate vector is indexed by category.
			if (sum_weight > 0) centres[j] = sum_value/sum_weight;
		}
	}
}

// Posterior-mean rate per pattern, then K categories at the cluster centres.
// Both halves are standard: the per-site posterior mean is the empirical-Bayes
// rate estimator (Mayrose et al. 2004), and clustering it is the least-squares
// quantizer (Lloyd 1982). What they buy together is that the categories are
// placed where the data are instead of on a fixed grid, and that each pattern is
// shrunk toward the mean by exactly as much as its column is uninformative.
// Returns the number of traversals _cat_probe performed.
static int _cat_assign_posterior_mean(SingleTreeLikelihood* tlk, double* rates,
                                      int cat_count, const CatOptions* options){
	int pattern_count = tlk->sp->count;
	int probe_count = options->probe_count;
	if (probe_count < cat_count) probe_count = cat_count;

	// A wider grid than the categories will end up spanning: the estimate has to
	// be free to sit outside 1/K..K, which is what caps the dispersion the fixed
	// grid can represent.
	double* probe_rates = dvector(probe_count);
	log_spaced_spaced_vector2(probe_rates, 1.0/probe_count, probe_count, probe_count);
	double* log_prior = dvector(probe_count);
	if (options->prior_shape > 0) {
		_cat_log_prior(log_prior, probe_rates, probe_count, options->prior_shape);
	}

	double* likelihoods = malloc(sizeof(double)*probe_count*pattern_count);
	_cat_probe(tlk, probe_rates, probe_count, likelihoods);

	// The quantizer works on log rates: the profile is close to symmetric there,
	// and the rate distributions this approximates (gamma, lognormal) put their
	// structure on the log scale.
	double* estimates = dvector(pattern_count);
	double* weights = dvector(pattern_count);
	for (int i = 0; i < pattern_count; i++) {
		weights[i] = tlk->sp->weights[i];
		double max = -DBL_MAX;
		for (int j = 0; j < probe_count; j++) {
			double logP = likelihoods[j*pattern_count + i] + log_prior[j];
			if (logP > max) max = logP;
		}
		double numerator = 0;
		double denominator = 0;
		for (int j = 0; j < probe_count; j++) {
			double posterior = exp(likelihoods[j*pattern_count + i] + log_prior[j] - max);
			numerator += posterior*probe_rates[j];
			denominator += posterior;
		}
		estimates[i] = log(numerator/denominator);
	}

	double* centres = dvector(cat_count);
	int* labels = ivector(pattern_count);
	_cat_quantize(estimates, weights, pattern_count, centres, labels, cat_count);
	for (int i = 0; i < cat_count; i++) rates[i] = exp(centres[i]);
	for (int i = 0; i < pattern_count; i++) tlk->sm->site_category[i] = labels[i];

	free(probe_rates);
	free(log_prior);
	free(likelihoods);
	free(estimates);
	free(weights);
	free(centres);
	free(labels);
	return probe_count;
}

// FastTree's rule: score every pattern at each of the K grid rates and keep the
// best, under the Gamma prior when there is one.
// Returns the number of traversals _cat_probe performed.
static int _cat_assign_argmax(SingleTreeLikelihood* tlk, const double* rates,
                              int cat_count, const CatOptions* options,
                              int* counts){
	int pattern_count = tlk->sp->count;
	double* log_prior = dvector(cat_count);
	if (options->prior_shape > 0) {
		_cat_log_prior(log_prior, rates, cat_count, options->prior_shape);
	}
	double* likelihoods = malloc(sizeof(double)*pattern_count*cat_count);
	_cat_probe(tlk, rates, cat_count, likelihoods);

	for (int i = 0; i < pattern_count; i++) {
		int best = 0;
		double bestLnl = -DBL_MAX;
		for (int j = 0; j < cat_count; j++) {
			double logP = likelihoods[j*pattern_count + i] + log_prior[j];
			if(logP > bestLnl){
				bestLnl = logP;
				best = j;
			}
		}
		tlk->sm->site_category[i] = best;
		if(options->verbosity > 0){
			printf("CAT pattern %d rate index %d rate %f [%f\n", i, best, rates[best],
			       bestLnl);
			counts[best]++;
		}
	}
	free(log_prior);
	free(likelihoods);
	return cat_count;
}

// Assign every pattern one rate category, by the rule in options->assignment.
//
// The probe scores every pattern against the block the pattern is currently
// assigned to, so it is only meaningful while every pattern sits in category 0:
// sm->cat_count is forced to 1 for the duration, and the transition-matrix refresh
// loop (treelikelihood.c, SingleTreeLikelihood_update_Q) therefore rewrites block 0
// alone. A pattern left in a category >= 1 by an earlier call would be scored against
// a block that never changes between probes, its profile would be flat, and the
// arg-max would hand it back to category 0 -- collapsing the whole assignment on the
// second call. Clearing the assignment first restores the precondition the first call
// gets for free, and is what makes repeated calls (an assign/optimize loop) possible.
// FastTree does the same, structurally, in AllocRateCategories.
int fasttree_cat(SingleTreeLikelihood* tlk, const CatOptions* options){
	SiteModel* sm = tlk->sm;
	int cat_count = sm->cat_count;
	memset(sm->site_category, 0, sizeof(int)*tlk->sp->count);
	sm->cat_count = 1;
	int* counts = NULL;
	if(options->verbosity > 0) counts = ivector(cat_count);
	double* rates = dvector(cat_count);
	log_spaced_spaced_vector2(rates, 1.0/cat_count, cat_count, cat_count);
//	log_spaced_spaced_vector2(rates+1, 1.0/(cat_count-1), cat_count-1, cat_count-1); // invariant
	if(options->verbosity > 0) print_dvector(rates, cat_count);
	tlk->calculate(tlk);// make sure everything is up-to-date
	sm->need_update = false;
	int evaluations = 1;

	if (options->assignment == CAT_ASSIGNMENT_POSTERIOR_MEAN) {
		evaluations += _cat_assign_posterior_mean(tlk, rates, cat_count, options);
		if(options->verbosity > 0){
			for (int i = 0; i < tlk->sp->count; i++) counts[sm->site_category[i]]++;
			print_dvector(rates, cat_count);
		}
	}
	else {
		evaluations += _cat_assign_argmax(tlk, rates, cat_count, options, counts);
	}

	// The category rates are one vector parameter. Every element but the last goes
	// in quietly and the last one loudly, so the listeners fire once, after the
	// whole vector is in place rather than on the first element of it.
	Parameter* cat_rates = Parameters_at(sm->rates, 0);
	for (int i = 0; i < cat_count - 1; i++) {
		Parameter_set_value_at_quietly(cat_rates, rates[i], i);
	}
	Parameter_set_value_at(cat_rates, rates[cat_count - 1], cat_count - 1);
	if(options->verbosity > 0){
		for (int i = 0; i < cat_count; i++) {
			if(counts[i] == 0) fprintf(stdout, "CAT rate %f not used\n", rates[i]);
		}
		free(counts);
	}
	SingleTreeLikelihood_use_rescaling(tlk, false);
	SingleTreeLikelihood_update_all_nodes(tlk);
	
	sm->cat_count = cat_count;
	free(rates);
	return evaluations;
}

void cat_check_sitemodel(json_node* node, const SingleTreeLikelihood* tlk){
	if (tlk->sm == NULL || tlk->sm->site_category == NULL) {
		json_die(node, "the tree likelihood does not carry an empirical CAT site "
		               "model: give its \"sitemodel\" a \"categories\" count and a "
		               "rate vector of that size, with no rate distribution");
	}
}

CatOptions cat_options_from_json(json_node* node){
	const char* assignment = get_json_node_value_string(node, "assignment");
	cat_assignment_t rule = CAT_ASSIGNMENT_ARGMAX;
	if (assignment != NULL) {
		if (strcasecmp(assignment, "argmax") == 0) {
			rule = CAT_ASSIGNMENT_ARGMAX;
		}
		else if (strcasecmp(assignment, "posterior_mean") == 0) {
			rule = CAT_ASSIGNMENT_POSTERIOR_MEAN;
		}
		else {
			json_die(node, "\"assignment\" is \"argmax\" or \"posterior_mean\", "
			               "not \"%s\"", assignment);
		}
	}

	CatOptions options = cat_options_default(rule);
	// Shape of the mean-one Gamma prior on the category rate. FastTree's value, and
	// on by default for the arg-max: without it the selected rate is unshrunk, the
	// assigned rates are over-dispersed and the branches stretch -- badly so on
	// small trees. Set to 0 to select on the likelihood alone. The posterior mean
	// defaults to 0 for the opposite reason: it shrinks on its own.
	options.prior_shape = get_json_node_value_double(node, "prior", options.prior_shape);
	options.probe_count = get_json_node_value_int(node, "probe", options.probe_count);
	options.verbosity = get_json_node_value_int(node, "verbosity", 0);
	if (rule == CAT_ASSIGNMENT_ARGMAX && get_json_node(node, "probe") != NULL) {
		json_die(node, "\"probe\" sets the size of the grid the posterior mean is "
		               "built on; \"assignment\": \"argmax\" scores the categories "
		               "themselves and has no separate grid");
	}
	if (options.probe_count < 0) {
		json_die(node, "\"probe\" must be positive, not %d", options.probe_count);
	}
	return options;
}

void cat_estimator_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"assignment", JSON_OPTIONAL, JSON_STRING},
	    {"id", JSON_OPTIONAL, JSON_STRING},
	    {"model", JSON_REQUIRED, JSON_STRING},
	    {"prior", JSON_OPTIONAL, JSON_NUMBER},
	    {"probe", JSON_OPTIONAL, JSON_NUMBER},
	    {"type", JSON_OPTIONAL, JSON_STRING},
	    {"verbosity", JSON_OPTIONAL, JSON_NUMBER},
	};
	json_validate(node, schema, sizeof(schema)/sizeof(schema[0]));

	char* ref = get_json_node_value_string(node, "model");
	Model* mtlk = Hashtable_get(hash, ref+1);
	if (mtlk == NULL || mtlk->type != MODEL_TREELIKELIHOOD) {
		json_die(node, "\"model\" must reference a tree likelihood: '%s'", ref);
	}
	SingleTreeLikelihood *tlk = mtlk->obj;
	cat_check_sitemodel(node, tlk);

	CatOptions options = cat_options_from_json(node);
	fasttree_cat(tlk, &options);
}
