// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <stdio.h>

#ifndef GSL_DISABLED
#include <gsl/gsl_rng.h>
#endif

#include "minunit.h"
#include "phyc/bootstrap.h"
#include "phyc/filereader.h"
#include "phyc/hashtable.h"
#include "phyc/matrix.h"
#include "phyc/sitepattern.h"
#include "phyc/treelikelihood.h"

static Model* _load(const char* file, Hashtable* hash) {
    char* content = load_file(file);
    json_node* json = create_json_tree(content);
    free(content);
    Model* model = new_TreeLikelihoodModel_from_json(json->children[0], hash);
    json_free_tree(json);
    return model;
}

// A replicate is the observed pattern set with new multiplicities. Deterministic
// stand-in for a multinomial draw: drop every third pattern and move its sites
// onto its neighbour, so the total number of sites is unchanged, some patterns
// have weight 0 and others are inflated -- the two features a real replicate has
// and the observed data does not.
static double* _fake_replicate(const SitePattern* sp) {
    double* counts = clone_dvector(sp->weights, sp->count);
    for (int i = 0; i < sp->count; i += 3) {
        int neighbour = (i + 1) % sp->count;
        counts[neighbour] += counts[i];
        counts[i] = 0;
    }
    return counts;
}

// The two paths of the bootstrap -- overwriting the weights in place (leaving
// zeros in the pattern set) and swapping in a compacted pattern set -- must
// agree exactly, and both must leave the likelihood restorable to the observed
// value. Run over both tip representations: with tip partials the compacted
// path has to rebuild the leaf partials, which the tip-states path does not.
static char* _test_bootstrap_paths(const char* file) {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);

    Model* model = _load(file, hash);
    SingleTreeLikelihood* tlk = model->obj;
    SitePattern* sp = tlk->sp;

    double observed_logP = model->logP(model);
    double* weights = clone_dvector(sp->weights, sp->count);
    int count = sp->count;

    // Compacting the observed weights keeps every pattern and must reproduce
    // the observed likelihood: the identity case of the swap machinery.
    SitePattern* view = new_SitePattern_view(sp);
    SitePattern_compact_into(sp, weights, view);
    mu_assert(view->count == count, "compacting the observed weights dropped patterns");
    SingleTreeLikelihood_set_sitepattern(tlk, view);
    mu_assert(fabs(model->logP(model) - observed_logP) < 1.e-10,
              "compacted observed data does not reproduce the observed logP");

    double* counts = _fake_replicate(sp);

    // Compacted replicate.
    SitePattern_compact_into(sp, counts, view);
    mu_assert(view->count < count, "the replicate should have dropped patterns");
    double total = 0;
    for (int i = 0; i < view->count; i++) total += view->weights[i];
    mu_assert(fabs(total - sp->nsites) < 1.e-10, "the replicate lost sites");
    SingleTreeLikelihood_set_sitepattern(tlk, view);
    double compact_logP = model->logP(model);

    // Same replicate as weights over the full pattern set, zeros included.
    SingleTreeLikelihood_set_sitepattern(tlk, sp);
    memcpy(sp->weights, counts, sizeof(double) * count);
    SingleTreeLikelihood_update_weights(tlk);
    double weighted_logP = model->logP(model);

    mu_assert(fabs(compact_logP - weighted_logP) < 1.e-10,
              "compacted and weighted replicates disagree");
    mu_assert(fabs(compact_logP - observed_logP) > 1.e-6,
              "the replicate is indistinguishable from the observed data");

    // Back to the observed data.
    memcpy(sp->weights, weights, sizeof(double) * count);
    SingleTreeLikelihood_update_weights(tlk);
    mu_assert(fabs(model->logP(model) - observed_logP) < 1.e-10,
              "the observed likelihood was not restored");

    free(counts);
    free(weights);
    free_SitePattern_view(view);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

char* test_bootstrap_tipstates() { return _test_bootstrap_paths("jc69-distance.json"); }

char* test_bootstrap_tippartials() {
    return _test_bootstrap_paths("jc69-distance-partials.json");
}

// Run a bootstrap of the jc69.json likelihood with `threads` workers, logging
// every replicate to `file`.
static void _run_bootstrap(size_t threads, const char* file) {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);
#ifndef GSL_DISABLED
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, 3);
    Hashtable_add(hash, "RANDOM_GENERATOR!@", rng);
#endif
    // Kept until the end, as physher does: the hashtable keys point into it.
    char* content = load_file("jc69.json");
    json_node* config = create_json_tree(content);
    free(content);
    Model* model = new_TreeLikelihoodModel_from_json(config->children[0], hash);
    Hashtable_add(hash, model->name, model);

    char node[1024];
    snprintf(node, sizeof(node),
             "{\"id\": \"boot\", \"type\": \"bootstrap\", \"replicates\": 6,"
             " \"include_observed\": true, \"threads\": %zu,"
             " \"treelikelihood\": \"@treelikelihood\","
             " \"optimizer\": {\"id\": \"opt\", \"type\": \"optimizer\","
             "   \"algorithm\": \"bfgs\", \"target\": \"@treelikelihood\","
             "   \"parameters\": \"&tree.distances\", \"verbosity\": 0},"
             " \"loggers\": [{\"id\": \"log\", \"type\": \"logger\","
             "   \"file\": \"%s\", \"every\": 1,"
             "   \"columns\": [\"@treelikelihood\", \"&tree.distances\"]}]}",
             threads, file);
    json_node* json = create_json_tree(node);
    Bootstrap* bootstrap = new_Bootstrap_from_json(json, hash);
    bootstrap->run(bootstrap);
    bootstrap->free(bootstrap);
    json_free_tree(json);

    model->free(model);
    free_Hashtable(hash);
    json_free_tree(config);
#ifndef GSL_DISABLED
    gsl_rng_free(rng);
#endif
}

// Replicates are drawn in order and logged in order from the original model,
// so the log must not depend on how many workers fitted them.
char* test_bootstrap_threads() {
    _run_bootstrap(1, "bootstrap-threads1.csv");
    _run_bootstrap(3, "bootstrap-threads3.csv");
    char* serial = load_file("bootstrap-threads1.csv");
    char* threaded = load_file("bootstrap-threads3.csv");
    mu_assert(strlen(serial) > 0, "nothing logged");
    mu_assert(strcmp(serial, threaded) == 0, "threaded log differs from the serial one");
    free(serial);
    free(threaded);
    remove("bootstrap-threads1.csv");
    remove("bootstrap-threads3.csv");
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_bootstrap_tipstates);
    mu_run_test(test_bootstrap_tippartials);
    mu_run_test(test_bootstrap_threads);
    return NULL;
}

RUN_TESTS(all_tests);
