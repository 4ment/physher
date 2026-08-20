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
    Model* model = _treelikelihood_from_file("jc69-cat-ref.json", hash);
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
// matrix block that never changed between probes.
static char* _check_idempotent(cat_assignment_t rule) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(rule);
    options.prior_shape = 0.0;
    fasttree_cat(tlk, &options);
    int* first = clone_ivector(sm->site_category, tlk->sp->count);
    double* first_rates = clone_dvector(sm->cat_rates, sm->cat_count);

    fasttree_cat(tlk, &options);
    for (int i = 0; i < tlk->sp->count; i++) {
        mu_assert(sm->site_category[i] == first[i],
                  "CAT: a second call changed the assignment");
    }
    for (size_t c = 0; c < sm->cat_count; c++) {
        mu_assert(fabs(sm->cat_rates[c] - first_rates[c]) < TOL,
                  "CAT: a second call changed the rates");
    }

    free(first);
    free(first_rates);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

static char* test_argmax_idempotent(void) {
    return _check_idempotent(CAT_ASSIGNMENT_ARGMAX);
}

static char* test_posterior_mean_idempotent(void) {
    return _check_idempotent(CAT_ASSIGNMENT_POSTERIOR_MEAN);
}

// Whatever the rule, the assignment has to be usable: every pattern in range, and
// every rate a positive finite number the site model can normalise.
static char* _check_wellformed(cat_assignment_t rule) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(rule);
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
    return _check_wellformed(CAT_ASSIGNMENT_ARGMAX);
}

static char* test_posterior_mean_wellformed(void) {
    return _check_wellformed(CAT_ASSIGNMENT_POSTERIOR_MEAN);
}

// The point of the posterior mean: it shrinks, and it shrinks by construction. A
// posterior mean over the probe grid is a convex combination of its rates, so it
// lies strictly inside the grid's range however extreme the profile is -- where the
// arg-max sits on the endpoints as soon as a column looks constant or saturated.
// That is the whole difference between the two rules, and unlike the size of the
// shrinkage it does not depend on which tree the assignment was run on.
static char* test_posterior_mean_stays_inside_the_grid(void) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
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

static char* all_tests() {
    mu_suite_start();
    mu_run_test(test_options_default);
    mu_run_test(test_argmax_matches_reference);
    mu_run_test(test_argmax_prior_matches_reference);
    mu_run_test(test_strong_prior_ignores_data);
    mu_run_test(test_argmax_idempotent);
    mu_run_test(test_posterior_mean_idempotent);
    mu_run_test(test_argmax_wellformed);
    mu_run_test(test_posterior_mean_wellformed);
    mu_run_test(test_posterior_mean_stays_inside_the_grid);
    mu_run_test(test_posterior_mean_single_category);
    mu_run_test(test_posterior_mean_more_categories_than_patterns);
    return NULL;
}

RUN_TESTS(all_tests);
