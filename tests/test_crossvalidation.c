// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

// The identities the "crossvalidation" action reduces held-out sites with
// (src/phyc/crossvalidation.c). The action itself needs an optimizer and a JSON
// run list, and is exercised end to end by examples/fluA/regression; what is
// checked here is the arithmetic underneath it, which is where a wrong answer
// would be silent rather than loud.

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "minunit.h"
#include "phyc/cat.h"
#include "phyc/filereader.h"
#include "phyc/hashtable.h"
#include "phyc/matrix.h"
#include "phyc/model.h"
#include "phyc/sitemodel.h"
#include "phyc/sitepattern.h"
#include "phyc/treelikelihood.h"

#define TOL 1.e-10

static Hashtable* _new_hash(void) {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);
    return hash;
}

static Model* _load(const char* file, Hashtable* hash) {
    char* content = load_file(file);
    json_node* json = create_json_tree(content);
    free(content);
    Model* model = new_TreeLikelihoodModel_from_json(json->children[0], hash);
    json_free_tree(json);
    return model;
}

// A held-out fold, deterministically: every third site of every pattern. Nothing
// here depends on the real splitter, only on the weights being a partition.
static double* _held_out(const SitePattern* sp) {
    double* test = dvector(sp->count);
    for (int i = 0; i < sp->count; i++) {
        test[i] = floor(sp->weights[i] / 3.0);
    }
    return test;
}

static void _install(SingleTreeLikelihood* tlk, const double* weights) {
    double total = 0;
    for (int i = 0; i < tlk->sp->count; i++) total += weights[i];
    memcpy(tlk->sp->weights, weights, sizeof(double) * tlk->sp->count);
    tlk->sp->nsites = (int)(total + 0.5);
    SingleTreeLikelihood_update_weights(tlk);
    if (tlk->sm != NULL && tlk->sm->site_category != NULL) tlk->sm->need_update = true;
}

// The reduction the action performs on a non-CAT model: the per-pattern array is
// dotted with the held-out multiplicities. That is only a log-likelihood over the
// held-out sites if the same array dotted with the *observed* multiplicities is
// the observed log-likelihood, which is the identity tested here.
static char* test_pattern_lk_is_the_likelihood(void) {
    Hashtable* hash = _new_hash();
    Model* model = _load("jc69-freerate.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    double logP = model->logP(model);
    SingleTreeLikelihood_update_all_nodes(tlk);
    tlk->calculate(tlk);

    double dot = 0;
    for (int i = 0; i < tlk->sp->count; i++) {
        dot += tlk->pattern_lk[i] * tlk->sp->weights[i];
    }

    char* failure = NULL;
    if (fabs(dot - logP) > 1.e-8) {
        failure = (char*)"crossvalidation: the per-pattern log likelihoods do not "
                         "sum to the log likelihood";
    }
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// The same identity for the CAT path, which cannot use pattern_lk: under a hard
// assignment that array is the *selected*-category density, so the action takes
// the marginalized mixture instead. The array cat_mixture_patterns writes has to
// be the one whose weighted sum cat_mixture already reports.
static char* test_cat_mixture_patterns_agree(void) {
    Hashtable* hash = _new_hash();
    Model* model = _load("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    cat_assign(tlk, &options);

    double* patterns = dvector(tlk->sp->count);
    CatMixture with = cat_mixture_patterns(tlk, 0, patterns);
    CatMixture without = cat_mixture(tlk, 0);

    char* failure = NULL;
    // Adding the out-parameter must not have moved the number itself.
    if (fabs(with.logP_mixture - without.logP_mixture) > TOL) {
        failure = (char*)"crossvalidation: cat_mixture_patterns disagrees with "
                         "cat_mixture";
    }
    double dot = 0;
    for (int i = 0; i < tlk->sp->count; i++) {
        dot += patterns[i] * tlk->sp->weights[i];
    }
    if (failure == NULL && fabs(dot - with.logP_mixture) > 1.e-8) {
        failure = (char*)"crossvalidation: the per-pattern mixture densities do "
                         "not sum to the mixture log likelihood";
    }
    // The score is the conditional density and the mixture is not; a fixture with
    // any rate variation at all must tell them apart, or the test is vacuous.
    if (failure == NULL && fabs(with.gap) < 1.e-6) {
        failure = (char*)"crossvalidation: the fixture cannot distinguish the CAT "
                         "score from its mixture";
    }

    free(patterns);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

// An empirical CAT site model normalises its category rates to mean one over
// sp->weights and sp->nsites. Installing a training fold moves both, so the
// normaliser has to be recomputed against the fold actually installed -- if it is
// not, every rate is off by the ratio of the two totals and the branch lengths the
// fold chose are read at the wrong scale. This is the invariant that makes it safe
// to fit on a subset of the sites.
static char* test_cat_rates_stay_mean_one_on_a_fold(void) {
    Hashtable* hash = _new_hash();
    Model* model = _load("jc69-cat.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SitePattern* sp = tlk->sp;

    CatOptions options = cat_options_default(CAT_ASSIGNMENT_POSTERIOR_MEAN);
    cat_assign(tlk, &options);

    double observed_logP = model->logP(model);
    double* weights = clone_dvector(sp->weights, sp->count);
    int observed_nsites = sp->nsites;

    double* test = _held_out(sp);
    double* train = dvector(sp->count);
    for (int i = 0; i < sp->count; i++) train[i] = weights[i] - test[i];

    _install(tlk, train);
    model->logP(model);
    SiteModel* sm = tlk->sm;
    sm->update(sm);

    char* failure = NULL;
    double mean = 0;
    for (int i = 0; i < sp->count; i++) {
        mean += sm->get_rate(sm, sm->get_site_category(sm, i)) * sp->weights[i];
    }
    mean /= sp->nsites;
    if (fabs(mean - 1.0) > 1.e-8) {
        failure = (char*)"crossvalidation: the CAT rates are not mean one over the "
                         "training fold";
    }
    // A fold really is a different dataset: if the likelihood did not move, the
    // installation did nothing and the invariant above was never tested.
    if (failure == NULL && fabs(model->logP(model) - observed_logP) < 1.e-6) {
        failure = (char*)"crossvalidation: installing a training fold left the "
                         "likelihood unchanged";
    }

    // And it has to be reversible, or whatever runs after cross-validation sees a
    // mutilated alignment.
    _install(tlk, weights);
    sp->nsites = observed_nsites;
    if (failure == NULL && fabs(model->logP(model) - observed_logP) > TOL) {
        failure = (char*)"crossvalidation: the observed likelihood was not restored";
    }

    free(test);
    free(train);
    free(weights);
    model->free(model);
    free_Hashtable(hash);
    return failure;
}

char* all_tests(void) {
    mu_suite_start();
    mu_run_test(test_pattern_lk_is_the_likelihood);
    mu_run_test(test_cat_mixture_patterns_agree);
    mu_run_test(test_cat_rates_stay_mean_one_on_a_fold);
    return NULL;
}

RUN_TESTS(all_tests);
