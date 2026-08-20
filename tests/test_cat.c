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
static char* _check_optimizer(const char* model_key) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-cat.json", hash);
    Hashtable_add(hash, "treelikelihood", model);
    SingleTreeLikelihood* tlk = model->obj;

    char config[512];
    snprintf(config, sizeof(config),
             "{\"id\": \"catopt\", \"type\": \"optimizer\","
             " \"algorithm\": \"cat\", \"target\": \"@treelikelihood\","
             " \"%s\": \"@treelikelihood\","
             " \"assignment\": \"posterior_mean\", \"prior\": 3.0}",
             model_key);
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

static char* test_optimizer_model_key(void) { return _check_optimizer("model"); }

// The key serial Brent and EM have always used for the same thing, kept working.
static char* test_optimizer_treelikelihood_key(void) {
    return _check_optimizer("treelikelihood");
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
    mu_run_test(test_argmax_idempotent);
    mu_run_test(test_posterior_mean_idempotent);
    mu_run_test(test_argmax_wellformed);
    mu_run_test(test_posterior_mean_wellformed);
    mu_run_test(test_posterior_mean_stays_inside_the_grid);
    mu_run_test(test_posterior_mean_single_category);
    mu_run_test(test_posterior_mean_more_categories_than_patterns);
    mu_run_test(test_argmax_guard_does_not_lose_ground);
    mu_run_test(test_posterior_mean_guard_does_not_lose_ground);
    mu_run_test(test_guard_reverts_a_worse_assignment);
    mu_run_test(test_evaluation_count);
    mu_run_test(test_optimizer_model_key);
    mu_run_test(test_optimizer_treelikelihood_key);
    mu_run_test(test_upper_matches_lower_SSE);
    mu_run_test(test_upper_matches_lower);
    mu_run_test(test_upper_matches_lower_single_category);
    mu_run_test(test_upper_matches_lower_without_cat);
    return NULL;
}

RUN_TESTS(all_tests);
