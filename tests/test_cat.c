// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

// The empirical CAT assignment rules (src/phyc/cat.c). These run the whole tree
// likelihood: fasttree_cat only makes sense against real per-pattern profiles, and
// the profile it builds internally is not otherwise reachable.

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "minunit.h"
#include "phyc/cat.h"
#include "phyc/filereader.h"
#include "phyc/hashtable.h"
#include "phyc/matrix.h"
#include "phyc/optimizer.h"
#include "phyc/node.h"
#include "phyc/tree.h"
#include "phyc/parameters.h"
#include "phyc/sitemodel.h"
#include "phyc/treelikelihood.h"

#define TOL 1.e-10

static Hashtable* _new_hash() {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);
    return hash;
}

static Model* _treelikelihood_from_file(const char* file, Hashtable* hash) {
    char* content = load_file(file);
    json_node* json = create_json_tree(content);
    free(content);
    Model* model = new_TreeLikelihoodModel_from_json(json->children[0], hash);
    json_free_tree(json);
    return model;
}

// The grid fasttree_cat builds for K categories, recomputed here so the tests do
// not have to trust the copy under test.
static double* _grid(int count) {
    double* rates = dvector(count);
    log_spaced_spaced_vector2(rates, 1.0 / count, count, count);
    return rates;
}

// The K x P profile, computed the slow honest way: a single-rate model evaluated
// once per rate. Independent of everything in cat.c.
static double* _profile(const char* file, const double* rates, int count,
                        int* pattern_count) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file(file, hash);
    SingleTreeLikelihood* tlk = model->obj;
    Parameter* mu = Hashtable_get(hash, "mu");
    *pattern_count = tlk->sp->count;

    double* profile = dvector(count * tlk->sp->count);
    for (int c = 0; c < count; c++) {
        Parameter_set_value(mu, rates[c]);
        model->logP(model);
        memcpy(profile + c * tlk->sp->count, tlk->pattern_lk,
               sizeof(double) * tlk->sp->count);
    }
    model->free(model);
    free_Hashtable(hash);
    return profile;
}

// The nonparametric MLE of the prior over the grid, recomputed here so the tests do
// not have to trust the copy in cat.c. Plain EM from the flat start,
// g_j <- g_j D_j / N with D_j = sum_i n_i L_ij / sum_k g_k L_ik, run for a fixed
// number of passes so the comparison is deterministic. The floor mirrors
// CAT_NPMLE_FLOOR in cat.c: without it the two iterations part company as soon as a
// rate is driven out of the estimate. Writes the Kiefer--Wolfowitz gap of the
// distribution it returns to `gap`.
static double* _npmle(const double* profile, const double* weights,
                      int pattern_count, int count, int passes, double* gap) {
    double* scaled = dvector(pattern_count * count);
    for (int i = 0; i < pattern_count; i++) {
        double max = -INFINITY;
        for (int j = 0; j < count; j++) {
            double value = profile[j * pattern_count + i];
            if (value > max) max = value;
        }
        for (int j = 0; j < count; j++) {
            scaled[i * count + j] = exp(profile[j * pattern_count + i] - max);
        }
    }
    double sites = 0;
    for (int i = 0; i < pattern_count; i++) sites += weights[i];

    double* g = dvector(count);
    double* gradient = dvector(count);
    for (int j = 0; j < count; j++) g[j] = 1.0 / count;
    for (int pass = 0;; pass++) {
        for (int j = 0; j < count; j++) gradient[j] = 0;
        for (int i = 0; i < pattern_count; i++) {
            const double* row = scaled + i * count;
            double mixture = 0;
            for (int j = 0; j < count; j++) mixture += g[j] * row[j];
            double weight = weights[i] / mixture;
            for (int j = 0; j < count; j++) gradient[j] += weight * row[j];
        }
        *gap = 0;
        for (int j = 0; j < count; j++) {
            double excess = gradient[j] / sites - 1.0;
            if (excess > *gap) *gap = excess;
        }
        if (pass >= passes) break;
        for (int j = 0; j < count; j++) {
            g[j] = fmax(g[j] * gradient[j] / sites, 1.e-100);
        }
    }
    free(scaled);
    free(gradient);
    return g;
}

static char* test_options_default(void) {
    CatOptions argmax = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    mu_assert(argmax.assignment == CAT_ASSIGNMENT_ARGMAX, "CAT: rule not kept");
    mu_assert(argmax.prior_shape == 3.0, "CAT: the arg-max defaults to FastTree's prior");

    CatOptions mean = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    mu_assert(mean.prior_shape == 0.0,
              "CAT: the posterior mean shrinks on its own and defaults to no prior");
    mu_assert(mean.probe_count > 0, "CAT: the posterior mean needs a probe grid");
    return NULL;
}

// The shipped rule is an arg-max over the profile plus the log prior, so
// recomputing the profile independently has to reproduce the assignment pattern for
// pattern. Run at both prior settings: with the prior on this is also what pins its
// sign, since the wrong one would move every pattern the other way.
static char* _check_argmax_matches_reference(double prior_shape) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;
    int count = (int)sm->cat_count;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    options.prior_shape = prior_shape;
    fasttree_cat(tlk, &options);

    double* rates = _grid(count);
    int pattern_count = 0;
    double* profile = _profile("jc69-cat-ref.json", rates, count, &pattern_count);
    mu_assert(pattern_count == tlk->sp->count, "CAT: reference has other patterns");

    for (int i = 0; i < pattern_count; i++) {
        int best = 0;
        double best_value = -INFINITY;
        for (int c = 0; c < count; c++) {
            double value = profile[c * pattern_count + i];
            if (prior_shape > 0) {
                value += (prior_shape - 1.0) * log(rates[c]) - prior_shape * rates[c];
            }
            if (value > best_value) {
                best_value = value;
                best = c;
            }
        }
        mu_assert(sm->site_category[i] == best, "CAT: arg-max disagrees with the profile");
    }

    free(rates);
    free(profile);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

static char* test_argmax_matches_reference(void) {
    return _check_argmax_matches_reference(0.0);
}

static char* test_argmax_prior_matches_reference(void) {
    return _check_argmax_matches_reference(3.0);
}

// The estimated prior, end to end. The assignment under "prior": "npmle" has to be
// the one an independently computed nonparametric MLE implies, pattern for pattern,
// and the optimality gap the call reports has to be that estimate's gap -- which
// together pin both halves: that the profile reaching the estimator is the right
// one, and that the estimator is solving the problem it claims to.
//
// The pass cap is small and the tolerance zero so the cap is what stops both
// iterations and the comparison is against a defined iterate rather than against
// whatever two convergence tests happened to accept.
static char* test_npmle_matches_reference(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat-short.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;
    int count = (int)sm->cat_count;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    options.prior = CAT_PRIOR_NPMLE;
    options.npmle_iterations = 200;
    options.npmle_tolerance = 0.0;
    CatResult result = fasttree_cat(tlk, &options);

    double* rates = _grid(count);
    int pattern_count = 0;
    double* profile = _profile("jc69-cat-short-ref.json", rates, count,
                               &pattern_count);
    mu_assert(pattern_count == tlk->sp->count, "CAT: reference has other patterns");

    double gap = 0;
    double* g = _npmle(profile, tlk->sp->weights, pattern_count, count, 200, &gap);

    char* failure = NULL;
    if (result.npmle_passes != 200) {
        failure = (char*)"CAT: the estimated prior did not run the passes it was given";
    }
    if (failure == NULL && fabs(result.npmle_gap - gap) > 1.e-8 * fmax(1.0, gap)) {
        failure = (char*)"CAT: the reported optimality gap is not the gap of the "
                         "nonparametric MLE of the prior";
    }
    int atoms = 0;
    for (int j = 0; j < count; j++) {
        if (g[j] > 1.e-6) atoms++;
    }
    if (failure == NULL && result.npmle_atoms != atoms) {
        failure = (char*)"CAT: the reported support of the estimated prior is not "
                         "the support of the nonparametric MLE";
    }
    for (int i = 0; failure == NULL && i < pattern_count; i++) {
        int best = 0;
        double best_value = -INFINITY;
        for (int c = 0; c < count; c++) {
            double value = profile[c * pattern_count + i] + log(g[c]);
            if (value > best_value) {
                best_value = value;
                best = c;
            }
        }
        if (sm->site_category[i] != best) {
            failure = (char*)"CAT: the arg-max under the estimated prior disagrees "
                             "with the reference";
        }
    }

    free(g);
    free(profile);
    free(rates);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// The gap is a certificate, not a progress meter: when the loop stops on it the
// answer really is that close to the maximum. Given enough passes and a tolerance
// it can reach, the run has to end on the tolerance rather than on the cap.
static char* test_npmle_stops_on_its_certificate(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    options.prior = CAT_PRIOR_NPMLE;
    options.npmle_iterations = 100000;
    options.npmle_tolerance = 1.e-3;
    CatResult result = fasttree_cat(tlk, &options);

    char* failure = NULL;
    if (result.npmle_passes >= options.npmle_iterations) {
        failure = (char*)"CAT: the estimated prior never met its own tolerance";
    }
    if (failure == NULL && !(result.npmle_gap <= options.npmle_tolerance)) {
        failure = (char*)"CAT: the estimated prior stopped above the gap it was "
                         "asked for";
    }
    // A distribution has to put its mass somewhere, and it cannot put it on more
    // rates than the grid offers.
    if (failure == NULL
        && (result.npmle_atoms < 1 || result.npmle_atoms > options.probe_count)) {
        failure = (char*)"CAT: the estimated prior is supported on an impossible "
                         "number of rates";
    }

    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// One category is one grid rate, and the only distribution on one point is the one
// that puts everything there. That is already optimal, so the iteration has nothing
// to do and must say so rather than spending its budget.
static char* test_npmle_single_category(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat1.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    options.prior = CAT_PRIOR_NPMLE;
    CatResult result = fasttree_cat(tlk, &options);

    char* failure = NULL;
    if (result.npmle_atoms != 1 || result.npmle_passes != 0
        || fabs(result.npmle_gap) > TOL) {
        failure = (char*)"CAT: the estimated prior on a single rate is not the point "
                         "mass, already optimal";
    }

    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// The diagnostics belong to the estimated prior and mean nothing without one, so a
// fixed prior has to leave them alone rather than report a stale or invented
// number.
static char* test_npmle_diagnostics_only_when_estimated(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    CatResult result = fasttree_cat(tlk, &options);

    char* failure = NULL;
    if (result.npmle_atoms != 0 || result.npmle_passes != 0 || result.npmle_gap != 0) {
        failure = (char*)"CAT: a fixed prior reported an estimated one's diagnostics";
    }

    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// Estimating the prior has to change something, or the option is decoration. The
// flat prior is the one it starts from, so the comparison is against that rather
// than against FastTree's Gamma, and it is made on the category rates because those
// are what the estimate is meant to move.
static char* test_npmle_moves_the_rates(void) {
    Hashtable* hash = _new_hash();
    Model* flat_model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* flat_tlk = flat_model->obj;
    CatOptions flat = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    fasttree_cat(flat_tlk, &flat);
    double* flat_rates = clone_dvector(flat_tlk->sm->cat_rates,
                                       flat_tlk->sm->cat_count);
    flat_model->free(flat_model);
    free_Hashtable(hash);

    hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    options.prior = CAT_PRIOR_NPMLE;
    fasttree_cat(tlk, &options);

    bool moved = false;
    for (size_t c = 0; c < tlk->sm->cat_count; c++) {
        if (fabs(tlk->sm->cat_rates[c] - flat_rates[c]) > TOL) moved = true;
    }
    char* failure = NULL;
    if (!moved) {
        failure = (char*)"CAT: estimating the prior left the categories where the "
                         "flat prior put them";
    }

    free(flat_rates);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// A prior this strong swamps every likelihood difference, so the choice stops
// depending on the data: every pattern lands in the same category, and a mean-one
// prior cannot pick an end of the grid. Catches a prior applied with the wrong sign
// or on the wrong scale, which the reference check above would also catch but only
// while the data happen to disagree with it.
static char* test_strong_prior_ignores_data(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    options.prior_shape = 1.e4;
    fasttree_cat(tlk, &options);

    int chosen = sm->site_category[0];
    for (int i = 1; i < tlk->sp->count; i++) {
        mu_assert(sm->site_category[i] == chosen,
                  "CAT: an overwhelming prior should ignore the data");
    }
    mu_assert(chosen > 0 && chosen < (int)sm->cat_count - 1,
              "CAT: a mean-one prior should not pick an extreme of the grid");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// Running the action twice in one process used to collapse the assignment onto
// category 0, because a pattern already in a category >= 1 was scored against a
// matrix block that never changed between probes. That is what this pins: a second
// call still has to produce a usable assignment over more than one category, and
// -- the guard -- it can never leave the model worse than it found it.
//
// What it deliberately does not pin is that the second call changes nothing. The
// posterior mean rebuilds its probe grid around the previous call's centres
// (docs/methods/cat-vs-raxml.md, section 2), so a second call is a refinement of
// the first rather than a repeat of it.
static char* _check_second_call(cat_assignment_t rule) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(rule);
    options.prior_shape = 0.0;
    CatResult first = fasttree_cat(tlk, &options);
    int* categories = clone_ivector(sm->site_category, tlk->sp->count);

    CatResult second = fasttree_cat(tlk, &options);

    char* failure = NULL;
    bool spread_before = false;
    bool spread_after = false;
    for (int i = 0; i < tlk->sp->count; i++) {
        if (categories[i] != 0) spread_before = true;
        if (sm->site_category[i] != 0) spread_after = true;
    }
    if (!spread_before) {
        failure = (char*)"CAT: the first call put every pattern in one category";
    }
    if (failure == NULL && !spread_after) {
        failure = (char*)"CAT: a second call collapsed the assignment onto category 0";
    }
    if (failure == NULL && second.logP < first.logP - TOL*fabs(first.logP)) {
        failure = (char*)"CAT: a second call lost likelihood";
    }
    if (failure == NULL && fabs(second.logP - model->logP(model)) > TOL) {
        failure = (char*)"CAT: a second call reported a likelihood the model does "
                         "not have";
    }

    free(categories);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

static char* test_argmax_second_call(void) {
    return _check_second_call(CAT_ASSIGNMENT_ARGMAX);
}

static char* test_posterior_mean_second_call(void) {
    return _check_second_call(CAT_ASSIGNMENT_POSTERIOR_MEAN);
}

// The arg-max scores the categories themselves, on a grid rebuilt identically every
// call, so its fixed point is exact: a second call reproduces the first to the last
// bit. It is the one rule the warm start leaves alone, and this is the assertion the
// posterior mean can no longer make.
static char* test_argmax_idempotent(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    options.prior_shape = 0.0;
    fasttree_cat(tlk, &options);
    int* first = clone_ivector(sm->site_category, tlk->sp->count);
    double* first_rates = clone_dvector(sm->cat_rates, sm->cat_count);

    fasttree_cat(tlk, &options);
    char* failure = NULL;
    for (int i = 0; i < tlk->sp->count && failure == NULL; i++) {
        if (sm->site_category[i] != first[i]) {
            failure = (char*)"CAT: a second arg-max call changed the assignment";
        }
    }
    for (size_t c = 0; c < sm->cat_count && failure == NULL; c++) {
        if (fabs(sm->cat_rates[c] - first_rates[c]) > TOL) {
            failure = (char*)"CAT: a second arg-max call changed the rates";
        }
    }

    free(first);
    free(first_rates);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// The warm start itself. A second posterior-mean call rebuilds its probe grid
// around the centres the first one produced, so it lands somewhere the fixed grid
// could not: on this fixture it moves the rates and the guard keeps the move, which
// is the whole point of section 2 -- successive rounds refine instead of repeating.
// A second call that came back identical would mean the grid had stopped following
// the centres.
static char* test_posterior_mean_refines(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    options.prior_shape = 0.0;
    CatResult first = fasttree_cat(tlk, &options);
    double* first_rates = clone_dvector(sm->cat_rates, sm->cat_count);

    CatResult second = fasttree_cat(tlk, &options);

    char* failure = NULL;
    bool moved = false;
    for (size_t c = 0; c < sm->cat_count; c++) {
        if (fabs(sm->cat_rates[c] - first_rates[c]) > TOL) moved = true;
    }
    if (!moved) {
        failure = (char*)"CAT: a second posterior-mean call reproduced the first, so "
                         "the probe grid is not following the centres";
    }
    if (failure == NULL && second.reverted) {
        failure = (char*)"CAT: the refined grid scored worse than the fixed one";
    }
    if (failure == NULL && second.logP <= first.logP) {
        failure = (char*)"CAT: the refined grid did not improve on the fixed one";
    }

    free(first_rates);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// Whatever the rule, the assignment has to be usable: every pattern in range, and
// every rate a positive finite number the site model can normalise.
static char* _check_wellformed(cat_assignment_t rule, cat_prior_t prior) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(rule);
    options.prior = prior;
    fasttree_cat(tlk, &options);

    for (int i = 0; i < tlk->sp->count; i++) {
        mu_assert(sm->site_category[i] >= 0 && sm->site_category[i] < (int)sm->cat_count,
                  "CAT: a pattern was assigned a category that does not exist");
    }
    Parameter* rates = Parameters_at(sm->rates, 0);
    for (size_t c = 0; c < sm->cat_count; c++) {
        double rate = Parameter_value_at(rates, c);
        mu_assert(isfinite(rate) && rate > 0.0, "CAT: a category rate is not positive");
    }
    // The rate vector must have reached the parameter, or the likelihood would go
    // on using whatever the config supplied.
    double logP = model->logP(model);
    mu_assert(isfinite(logP), "CAT: the likelihood is not finite after assignment");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

static char* test_argmax_wellformed(void) {
    return _check_wellformed(CAT_ASSIGNMENT_ARGMAX, CAT_PRIOR_FIXED);
}

static char* test_posterior_mean_wellformed(void) {
    return _check_wellformed(CAT_ASSIGNMENT_POSTERIOR_MEAN, CAT_PRIOR_FIXED);
}

// The estimated prior can be sparse enough to leave a category with no rate near
// it, so an assignment made under it is worth the same check as one made under a
// fixed prior: every pattern in range, every rate positive and finite.
static char* test_argmax_npmle_wellformed(void) {
    return _check_wellformed(CAT_ASSIGNMENT_ARGMAX, CAT_PRIOR_NPMLE);
}

static char* test_posterior_mean_npmle_wellformed(void) {
    return _check_wellformed(CAT_ASSIGNMENT_POSTERIOR_MEAN, CAT_PRIOR_NPMLE);
}

// The point of the posterior mean: it shrinks, and it shrinks by construction. A
// posterior mean over the probe grid is a convex combination of its rates, so it
// lies strictly inside the grid's range however extreme the profile is -- where the
// arg-max sits on the endpoints as soon as a column looks constant or saturated.
// That is the whole difference between the two rules, and unlike the size of the
// shrinkage it does not depend on which tree the assignment was run on.
static char* _check_posterior_mean_stays_inside_the_grid(cat_prior_t prior) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    options.prior = prior;
    fasttree_cat(tlk, &options);

    double lower = 1.0 / options.probe_count;
    double upper = options.probe_count;
    Parameter* rates = Parameters_at(sm->rates, 0);
    for (size_t c = 0; c < sm->cat_count; c++) {
        double rate = Parameter_value_at(rates, c);
        mu_assert(rate > lower && rate < upper,
                  "CAT: a posterior mean escaped the grid it was averaged over");
    }
    model->free(model);
    free_Hashtable(hash);

    // The arg-max on the same alignment does reach its endpoints, so the bound
    // above is a real constraint rather than one no rule could violate.
    hash = _new_hash();
    model = _treelikelihood_from_file("jc69-cat.json", hash);
    tlk = model->obj;
    sm = tlk->sm;
    CatOptions argmax = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    argmax.prior_shape = 0.0;
    fasttree_cat(tlk, &argmax);
    bool extreme = false;
    for (int i = 0; i < tlk->sp->count; i++) {
        if (sm->site_category[i] == 0 ||
            sm->site_category[i] == (int)sm->cat_count - 1) {
            extreme = true;
        }
    }
    mu_assert(extreme, "CAT: expected the arg-max to use the ends of its grid");
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

static char* test_posterior_mean_stays_inside_the_grid(void) {
    return _check_posterior_mean_stays_inside_the_grid(CAT_PRIOR_FIXED);
}

// Convexity is a property of the average, not of the prior it averages against, so
// estimating the prior cannot let an estimate out of the grid however much mass the
// estimate moves to one end of it.
static char* test_posterior_mean_npmle_stays_inside_the_grid(void) {
    return _check_posterior_mean_stays_inside_the_grid(CAT_PRIOR_NPMLE);
}

// A single category cannot express any heterogeneity: every pattern goes to 0 and
// the model has to collapse onto the single-rate one. Exercises the quantizer's
// degenerate case, where the seeding has no quantiles to spread over.
static char* test_posterior_mean_single_category(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat1.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;
    mu_assert(sm->cat_count == 1, "CAT: fixture should have one category");

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    fasttree_cat(tlk, &options);
    for (int i = 0; i < tlk->sp->count; i++) {
        mu_assert(sm->site_category[i] == 0, "CAT: one category leaves no choice");
    }
    double logP = model->logP(model);
    // Normalised to mean one over the alignment, a single category is rate one.
    mu_assert(fabs(sm->get_rate(sm, 0) - 1.0) < 1.e-8,
              "CAT: the only category should normalise to rate one");
    mu_assert(isfinite(logP), "CAT: single-category likelihood is not finite");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// More categories than the profile can distinguish: the spare ones must not be
// handed a rate of zero or a duplicate that breaks the normalisation.
static char* test_posterior_mean_more_categories_than_patterns(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    options.probe_count = 2;  // a two-point profile cannot fill four categories
    fasttree_cat(tlk, &options);

    Parameter* rates = Parameters_at(sm->rates, 0);
    for (size_t c = 0; c < sm->cat_count; c++) {
        double rate = Parameter_value_at(rates, c);
        mu_assert(isfinite(rate) && rate > 0.0, "CAT: a spare category has no valid rate");
    }
    mu_assert(isfinite(model->logP(model)), "CAT: likelihood is not finite");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The guard: the call reports the likelihood it left behind, and that likelihood is
// never below the one it started from. On a fresh model the assignment is the first
// one there has ever been, so it has only the single-rate likelihood to beat and the
// step must be accepted.
static char* _check_guard_does_not_lose_ground(cat_assignment_t rule) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    double before = model->logP(model);
    CatOptions options = cat_options_default(rule);
    CatResult result = fasttree_cat(tlk, &options);

    mu_assert(!result.reverted, "CAT: the first assignment should beat a single rate");
    mu_assert(result.logP >= before, "CAT: the guard let the likelihood drop");
    // What it reports has to be the model's own value, or a schedule budgeting or
    // converging on it would be reading a number from somewhere else.
    mu_assert(fabs(result.logP - model->logP(model)) < TOL,
              "CAT: reported a likelihood the model does not have");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

static char* test_argmax_guard_does_not_lose_ground(void) {
    return _check_guard_does_not_lose_ground(CAT_ASSIGNMENT_ARGMAX);
}

static char* test_posterior_mean_guard_does_not_lose_ground(void) {
    return _check_guard_does_not_lose_ground(CAT_ASSIGNMENT_POSTERIOR_MEAN);
}

// The revert path. A prior strong enough to swamp the data collapses every pattern
// into one category, which is the single-rate model and so strictly worse than a
// fitted assignment -- run after one, it has to be refused, and refused exactly: the
// categories, the rates and the likelihood all back where they were.
static char* test_guard_reverts_a_worse_assignment(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;
    Parameter* rates = Parameters_at(sm->rates, 0);

    CatOptions fit = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    double fitted = fasttree_cat(tlk, &fit).logP;
    int* categories = clone_ivector(sm->site_category, tlk->sp->count);
    double* fitted_rates = dvector(sm->cat_count);
    for (size_t c = 0; c < sm->cat_count; c++) {
        fitted_rates[c] = Parameter_value_at(rates, c);
    }

    CatOptions collapse = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    collapse.prior_shape = 1.e4;
    CatResult result = fasttree_cat(tlk, &collapse);

    char* failure = NULL;
    if (!result.reverted) {
        failure = (char*)"CAT: collapsing onto one rate should have been refused";
    }
    if (failure == NULL && fabs(result.logP - fitted) > TOL) {
        failure = (char*)"CAT: the revert did not recover the likelihood";
    }
    // One traversal to bring the likelihood up to date, one per category probed,
    // one to score what the probe chose and one more to put the old assignment
    // back: the restore is a full evaluation and has to be charged for.
    if (failure == NULL && result.evaluations != 3 + (int)sm->cat_count) {
        failure = (char*)"CAT: the revert was not charged for";
    }
    for (int i = 0; i < tlk->sp->count && failure == NULL; i++) {
        if (sm->site_category[i] != categories[i]) {
            failure = (char*)"CAT: the revert did not restore the assignment";
        }
    }
    for (size_t c = 0; c < sm->cat_count && failure == NULL; c++) {
        if (fabs(Parameter_value_at(rates, c) - fitted_rates[c]) > TOL) {
            failure = (char*)"CAT: the revert did not restore the rates";
        }
    }
    // The parameter is only half of it: the model has to see the restored value too,
    // which it does not if the revert wrote the rates without firing the listeners.
    if (failure == NULL && fabs(model->logP(model) - fitted) > TOL) {
        failure = (char*)"CAT: the model did not pick the restored rates up";
    }

    free(categories);
    free(fitted_rates);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// What fasttree_cat charges for: the traversal that brings the likelihood up to
// date, one per rate it probes at, and the one the guard spends scoring the
// assignment it just made. A schedule budgets on this number, so a probe loop that
// grew or shrank without the count following it would spend an evaluation allowance
// it never reports.
//
// Each count is taken on its own model: a call that reverts pays for one more
// traversal putting the old assignment back, which is charged for in the revert test
// and would otherwise make these numbers depend on what ran before them.
static char* _check_evaluation_count(const CatOptions* options, int probes,
                                     const char* message) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    CatResult result = fasttree_cat(tlk, options);
    char* failure = NULL;
    if (result.reverted || result.evaluations != 2 + probes) {
        failure = (char*)message;
    }

    model->free(model);
    free_Hashtable(hash);
    return failure;
}

static char* test_evaluation_count(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    int cat_count = (int)((SingleTreeLikelihood*)model->obj)->sm->cat_count;
    model->free(model);
    free_Hashtable(hash);

    CatOptions argmax = cat_options_default(CAT_ASSIGNMENT_ARGMAX);
    char* failure = _check_evaluation_count(&argmax, cat_count,
                                            "CAT: the arg-max probes once per category");
    if (failure != NULL) return failure;

    CatOptions mean = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    mean.probe_count = 7;
    failure = _check_evaluation_count(&mean, 7,
                                      "CAT: the posterior mean probes once per grid point");
    if (failure != NULL) return failure;

    // Fewer probes than categories is widened to the categories, and the count
    // has to follow it there too.
    mean.probe_count = 1;
    return _check_evaluation_count(&mean, cat_count,
                                   "CAT: a grid narrower than the categories is widened to them");
}

// The same reassignment reached through "algorithm": "cat" in an optimizer. What
// is under test is the wiring -- the model lookup, the option keys and the value
// reported back -- not the assignment, which the tests above cover; so the
// optimizer is checked against a direct fasttree_cat call on a second copy.
// `model_key` is the JSON key naming the tree likelihood: "model" and the older
// "treelikelihood" have to be interchangeable.
static char* _check_optimizer(const char* model_key, const char* prior_json,
                              cat_prior_t prior) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    Hashtable_add(hash, "treelikelihood", model);
    SingleTreeLikelihood* tlk = model->obj;

    char config[512];
    snprintf(config, sizeof(config),
             "{\"id\": \"catopt\", \"type\": \"optimizer\","
             " \"algorithm\": \"cat\", \"target\": \"@treelikelihood\","
             " \"%s\": \"@treelikelihood\","
             " \"assignment\": \"posterior_mean\", \"prior\": %s}",
             model_key, prior_json);
    json_node* node = create_json_tree(config);
    Optimizer* opt = new_Optimizer_from_json(node, hash);

    double fmin = 0;
    mu_assert(opt_optimize(opt, &fmin) == OPT_SUCCESS,
              "CAT optimizer: did not report success");
    // The objective is the target's negative log-likelihood, so a meta schedule
    // can compare what this entry reports with what every other one does.
    mu_assert(fabs(fmin + model->logP(model)) < TOL,
              "CAT optimizer: reported something other than the target's value");

    int pattern_count = tlk->sp->count;
    int* assigned = ivector(pattern_count);
    memcpy(assigned, tlk->sm->site_category, sizeof(int) * pattern_count);
    double* assigned_rates = dvector(tlk->sm->cat_count);
    for (size_t c = 0; c < tlk->sm->cat_count; c++) {
        assigned_rates[c] = Parameter_value_at(Parameters_at(tlk->sm->rates, 0), c);
    }

    free_Optimizer(opt);
    json_free_tree(node);
    model->free(model);
    free_Hashtable(hash);

    Hashtable* reference_hash = _new_hash();
    Model* reference = _treelikelihood_from_file("jc69-cat.json", reference_hash);
    SingleTreeLikelihood* reference_tlk = reference->obj;
    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    options.prior = prior;
    options.prior_shape = 3.0;
    fasttree_cat(reference_tlk, &options);

    char* failure = NULL;
    for (int i = 0; i < pattern_count && failure == NULL; i++) {
        if (assigned[i] != reference_tlk->sm->site_category[i]) {
            failure = (char*)"CAT optimizer: assigned a different category than "
                             "fasttree_cat";
        }
    }
    for (size_t c = 0; c < reference_tlk->sm->cat_count && failure == NULL; c++) {
        double rate = Parameter_value_at(Parameters_at(reference_tlk->sm->rates, 0), c);
        if (fabs(assigned_rates[c] - rate) > TOL) {
            failure = (char*)"CAT optimizer: set a different rate than fasttree_cat";
        }
    }

    free(assigned);
    free(assigned_rates);
    reference->free(reference);
    free_Hashtable(reference_hash);
    return failure;
}

static char* test_optimizer_model_key(void) {
    return _check_optimizer("model", "3.0", CAT_PRIOR_FIXED);
}

// "prior" carries two different kinds of value, so the string form needs its own
// pass through the JSON: a number is a Gamma shape and the word "npmle" asks for
// the prior to be estimated instead. Checked the same way as the shape -- against a
// direct call that asks for it in C -- so this is about the key reaching the
// estimator, not about what the estimator then does.
static char* test_optimizer_npmle_key(void) {
    return _check_optimizer("model", "\"npmle\"", CAT_PRIOR_NPMLE);
}

// The key serial Brent and EM have always used for the same thing, kept working.
static char* test_optimizer_treelikelihood_key(void) {
    return _check_optimizer("treelikelihood", "3.0", CAT_PRIOR_FIXED);
}

// The branch sweep (serial_brent_optimize_tree) never evaluates the whole tree: it
// combines the upper partials of one node with the lower partials below it and reads
// the likelihood off that. At unchanged branch lengths the result has to equal the
// ordinary lower likelihood, at every node -- that identity is the whole justification
// for the sweep, and it is what the CAT kernels in treelikelihood4CAT.c have to
// preserve.
//
// Walking the nodes in postorder is what the sweep does, so this also exercises
// _calculate_uppper's incremental path (the sibling's uppers are rebuilt from the
// previous node) and not just its from-scratch one.
static char* _check_upper_matches_lower(const char* file, bool cat, bool sse) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file(file, hash);
    SingleTreeLikelihood* tlk = model->obj;
    SingleTreeLikelihood_enable_SSE(tlk, sse);
    if (cat) {
        CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
        fasttree_cat(tlk, &options);
    }

    double lower = model->logP(model);
    mu_assert(isfinite(lower), "upper: the lower likelihood is not finite");

    tlk->node_upper = NULL;
    tlk->use_upper = true;
    tlk->update_upper = true;
    SingleTreeLikelihood_update_uppers(tlk);

    Node** nodes = Tree_get_nodes(tlk->tree, POSTORDER);
    char* failure = NULL;
    for (int i = 0; i < Tree_node_count(tlk->tree) && failure == NULL; i++) {
        Node* node = nodes[i];
        // The root, and the child of the root pinned to a zero branch, are the two the
        // sweep skips: they have no branch of their own to optimize.
        if (!Node_has_distance(node)) continue;
        double upper = tlk->calculate_upper(tlk, node);
        if (!isfinite(upper) || fabs(upper - lower) > 1.e-6) {
            failure = (char*)"upper: a node's upper likelihood disagrees with the "
                             "lower likelihood at the same branch lengths";
        }
        // _calculate does this for the optimizer; _calculate_uppper reads it to decide
        // whether it can update incrementally.
        tlk->node_upper = node;
    }
    tlk->use_upper = false;

    model->free(model);
    free_Hashtable(hash);
    return failure;
}

static char* test_upper_matches_lower_SSE(void) {
    return _check_upper_matches_lower("jc69-cat.json", true, true);
}

static char* test_upper_matches_lower(void) {
    return _check_upper_matches_lower("jc69-cat.json", true, false);
}

// A single category leaves no room for the per-pattern matrix offset to be wrong, so
// it separates a broken offset from a broken recursion.
static char* test_upper_matches_lower_single_category(void) {
    return _check_upper_matches_lower("jc69-cat1.json", true, true);
}

// The same identity on a model that is not CAT. If this failed too, the check itself
// would be wrong rather than the CAT kernels.
static char* test_upper_matches_lower_without_cat(void) {
    return _check_upper_matches_lower("jc69-freerate.json", false, true);
}

static char* all_tests() {
    mu_suite_start();
    mu_run_test(test_options_default);
    mu_run_test(test_argmax_matches_reference);
    mu_run_test(test_argmax_prior_matches_reference);
    mu_run_test(test_strong_prior_ignores_data);
    mu_run_test(test_argmax_second_call);
    mu_run_test(test_posterior_mean_second_call);
    mu_run_test(test_argmax_idempotent);
    mu_run_test(test_posterior_mean_refines);
    mu_run_test(test_argmax_wellformed);
    mu_run_test(test_posterior_mean_wellformed);
    mu_run_test(test_posterior_mean_stays_inside_the_grid);
    mu_run_test(test_npmle_matches_reference);
    mu_run_test(test_npmle_stops_on_its_certificate);
    mu_run_test(test_npmle_single_category);
    mu_run_test(test_npmle_diagnostics_only_when_estimated);
    mu_run_test(test_npmle_moves_the_rates);
    mu_run_test(test_argmax_npmle_wellformed);
    mu_run_test(test_posterior_mean_npmle_wellformed);
    mu_run_test(test_posterior_mean_npmle_stays_inside_the_grid);
    mu_run_test(test_posterior_mean_single_category);
    mu_run_test(test_posterior_mean_more_categories_than_patterns);
    mu_run_test(test_argmax_guard_does_not_lose_ground);
    mu_run_test(test_posterior_mean_guard_does_not_lose_ground);
    mu_run_test(test_guard_reverts_a_worse_assignment);
    mu_run_test(test_evaluation_count);
    mu_run_test(test_optimizer_model_key);
    mu_run_test(test_optimizer_npmle_key);
    mu_run_test(test_optimizer_treelikelihood_key);
    mu_run_test(test_upper_matches_lower_SSE);
    mu_run_test(test_upper_matches_lower);
    mu_run_test(test_upper_matches_lower_single_category);
    mu_run_test(test_upper_matches_lower_without_cat);
    return NULL;
}

RUN_TESTS(all_tests);
