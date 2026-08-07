// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>

#include "minunit.h"
#include "phyc/hashtable.h"
#include "phyc/parameters.h"
#include "phyc/sitemodel.h"

#define TOL 1.e-10

static Hashtable* _new_hash() {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);
    return hash;
}

static Model* _sitemodel_from_string(const char* json, Hashtable* hash) {
    json_node* root = create_json_tree(json);
    Model* model = new_SiteModel_from_json(root, hash);
    json_free_tree(root);
    return model;
}

// sum_k p_k r_k, the quantity the unit-mean constraint pins at 1.
static double _weighted_mean(SiteModel* sm) {
    double mean = 0;
    for (size_t i = 0; i < sm->cat_count; i++) {
        mean += sm->get_proportion(sm, i) * sm->get_rate(sm, i);
    }
    return mean;
}

// Mean-contribution simplex ("mean_contribution"): the free
// parameter is s_k = p_k r_k, so the rates are recovered by r_k = s_k / p_k and
// the unit mean holds by construction -- no weighted-mean denominator anywhere.
char* test_mean_contribution() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"mean_contribution\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    const double s[4] = {0.4, 0.3, 0.2, 0.1};
    const double p[4] = {0.1, 0.2, 0.3, 0.4};
    const double expected[4] = {4.0, 1.5, 2.0 / 3.0, 0.25};

    mu_assert(sm->cat_count == 4, "mean contribution: wrong number of categories");
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_proportion(sm, i) - p[i]) < TOL,
                  "mean contribution: proportions not matching");
        mu_assert(fabs(sm->get_rate(sm, i) - expected[i]) < TOL,
                  "mean contribution: rates not matching");
        mu_assert(fabs(sm->get_rate(sm, i) - s[i] / p[i]) < TOL,
                  "mean contribution: rate is not s_k/p_k");
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "mean contribution: rates do not have unit weighted mean");

    // A change to the simplex must invalidate the cached categories and the new
    // rates must again satisfy the constraint. An equal-contribution simplex makes
    // every category contribute 1/K to the mean, i.e. r_k = 1/(K p_k).
    const double s2[4] = {0.25, 0.25, 0.25, 0.25};
    Parameter_set_values(Hashtable_get(hash, "s"), s2);
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - 1.0 / (4.0 * p[i])) < TOL,
                  "mean contribution: rates not updated after a change");
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "mean contribution: unit weighted mean lost after a change");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// With an invariant class the constraint only runs over the variable categories,
// sum_{k>0} p_k r_k = 1, so s has one element fewer than the proportions simplex
// and category 0 is pinned at rate 0.
char* test_mean_contribution_invariant() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"mean_contribution\","
        "\"categories\":3,"
        "\"rates\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.2,0.3,0.5]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.25,0.15,0.25,0.35]}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    const double s[3] = {0.2, 0.3, 0.5};
    const double p[4] = {0.25, 0.15, 0.25, 0.35};

    // The proportions simplex being one longer than "categories" turns on +I.
    mu_assert(sm->invariant, "mean contribution +I: invariant class not detected");
    mu_assert(sm->cat_count == 4, "mean contribution +I: wrong number of categories");
    mu_assert(sm->get_rate(sm, 0) == 0.0,
              "mean contribution +I: invariant category is not rate 0");
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_proportion(sm, i) - p[i]) < TOL,
                  "mean contribution +I: proportions not matching");
    }
    for (size_t i = 1; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - s[i - 1] / p[i]) < TOL,
                  "mean contribution +I: rates not matching");
    }
    // The invariant category contributes nothing, so the mean over all categories
    // is still 1.
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "mean contribution +I: rates do not have unit weighted mean");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The mean-contribution and rate-shape parameterizations are two coordinate
// systems on the same constraint surface: with x proportional to s_k/p_k they
// must produce identical category rates.
char* test_mean_contribution_dual() {
    const char* mean_contribution =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"mean_contribution\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    // x_k = (s_k/p_k) / sum_j (s_j/p_j), i.e. the same rates up to scale.
    const char* rate_shape =
        "{\"id\":\"sitemodel2\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_shape\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"x\",\"type\":\"simplex\","
        "\"x\":[0.6233766233766234,0.23376623376623376,0.1038961038961039,"
        "0.03896103896103896]},"
        "\"proportions\":{\"id\":\"p2\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";

    Hashtable* hash = _new_hash();
    Model* model_b = _sitemodel_from_string(mean_contribution, hash);
    Model* model_a = _sitemodel_from_string(rate_shape, hash);
    SiteModel* smb = model_b->obj;
    SiteModel* sma = model_a->obj;

    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sma->get_rate(sma, i) - smb->get_rate(smb, i)) < 1.e-9,
                  "mean contribution: the two simplex parameterizations disagree");
    }
    mu_assert(fabs(_weighted_mean(sma) - 1.0) < TOL,
              "rate shape: rates do not have unit weighted mean");

    model_a->free(model_a);
    model_b->free(model_b);
    free_Hashtable(hash);
    return NULL;
}

// Ordered increments: the free vector holds the gaps between consecutive
// categories, so the raw rates are its running sum, normalised afterwards by the
// weighted mean. Equal gaps give evenly spaced rates.
char* test_rate_increments() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_increments\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0,1.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    // Raw rates 1,2,3,4 with weighted mean 3, so r = (1,2,3,4)/3.
    const double expected[4] = {1.0 / 3.0, 2.0 / 3.0, 1.0, 4.0 / 3.0};
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - expected[i]) < TOL,
                  "increments: rates not matching");
    }
    for (size_t i = 1; i < 4; i++) {
        mu_assert(sm->get_rate(sm, i) > sm->get_rate(sm, i - 1),
                  "increments: positive gaps must give increasing rates");
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "increments: rates do not have unit weighted mean");

    // The normalisation is homogeneous of degree zero, so scaling every gap by a
    // constant leaves the rates untouched: the parameterization carries one
    // redundant degree
    // of freedom and the likelihood is exactly flat along it.
    const double doubled[4] = {2.0, 2.0, 2.0, 2.0};
    Parameter_set_values(Hashtable_get(hash, "gaps"), doubled);
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - expected[i]) < TOL,
                  "increments: rates are not invariant to the scale of the gaps");
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// Ordered ratios: each raw rate is a fraction of the next one up, the chain
// hanging off a free rate for the last category.
char* test_rate_ratios() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_ratios\","
        "\"categories\":4,"
        "\"rates\":[{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5],\"lower\":1.e-8,\"upper\":0.99},"
        "{\"id\":\"top\",\"type\":\"parameter\",\"x\":8.0,\"lower\":0}],"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    // Raw rates 1,2,4,8 with weighted mean 4.9.
    const double raw[4] = {1.0, 2.0, 4.0, 8.0};
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - raw[i] / 4.9) < TOL,
                  "ratios: rates not matching");
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "ratios: rates do not have unit weighted mean");

    // The rate of the last category is redundant for the same reason the
    // increments carry a spare degree of freedom: normalisation cancels it.
    const double top2 = 16.0;
    Parameter_set_values(Hashtable_get(hash, "top"), &top2);
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - raw[i] / 4.9) < TOL,
                  "ratios: rates are not invariant to the rate of the last category");
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// Naming a parameterization that the parser would have inferred anyway must not
// change what is computed: a plain vector infers the increments, a pair of
// parameters the ratios.
char* test_rate_parameterization_inferred() {
    const char* increments =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rates\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[0.5,1.0,1.5,2.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    const char* named =
        "{\"id\":\"sitemodel2\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_increments\",\"categories\":4,"
        "\"rates\":{\"id\":\"gaps2\",\"type\":\"parameter\","
        "\"x\":[0.5,1.0,1.5,2.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p2\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";

    Hashtable* hash = _new_hash();
    Model* model_inferred = _sitemodel_from_string(increments, hash);
    Model* model_named = _sitemodel_from_string(named, hash);
    SiteModel* sm_inferred = model_inferred->obj;
    SiteModel* sm_named = model_named->obj;

    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm_inferred->get_rate(sm_inferred, i) -
                       sm_named->get_rate(sm_named, i)) < TOL,
                  "increments: inferred and named parameterizations disagree");
    }
    mu_assert(fabs(_weighted_mean(sm_inferred) - 1.0) < TOL,
              "increments: inferred rates do not have unit weighted mean");

    model_named->free(model_named);
    model_inferred->free(model_inferred);
    free_Hashtable(hash);
    return NULL;
}

// Parse a site model in a forked child and report its exit status, so the tests
// can check that a malformed "parameterization" is rejected instead of silently
// reading the wrong number of simplex elements.
static int _parse_status(const char* json) {
    fflush(stdout);
    fflush(stderr);
    pid_t pid = fork();
    if (pid == 0) {
        freopen("/dev/null", "w", stderr);
        freopen("/dev/null", "w", stdout);
        Hashtable* hash = _new_hash();
        Model* model = _sitemodel_from_string(json, hash);
        model->free(model);
        free_Hashtable(hash);
        _exit(0);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

char* test_mean_contribution_rejects() {
    const char* good =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"mean_contribution\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(good) == 0, "mean contribution: valid model should parse");

    // Unknown parameterization name.
    const char* unknown =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"parameterization\":\"cumprod\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(unknown) == 13,
              "mean contribution: unknown parameterization should die");

    // s must be a simplex: a plain positive vector does not sum to one, so the
    // unit mean would not hold.
    const char* not_simplex =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"mean_contribution\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"s\",\"type\":\"parameter\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(not_simplex) == 2,
              "mean contribution: non-simplex rates should die");

    // One contribution per variable category, no more, no less.
    const char* wrong_dim =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"mean_contribution\","
        "\"categories\":4,"
        "\"rates\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.5,0.3,0.2]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(wrong_dim) == 2,
              "mean contribution: mismatched rates dimension should die");

    // It parameterises free rates; it means nothing for a rate distribution
    // whose categories come from a quantile function.
    const char* not_discrete =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\","
        "\"parameterization\":\"mean_contribution\","
        "\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0}}";
    mu_assert(_parse_status(not_discrete) == 2,
              "mean contribution: non-discrete distribution should die");

    return NULL;
}

// Neither ordered parameterization checked the size of "rates" before: the
// category count
// comes from the proportions simplex, so a short vector was read off the end.
char* test_rate_parameterization_rejects() {
    // One increment per category, whether the parameterization was named...
    const char* short_increments =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_increments\",\"categories\":4,"
        "\"rates\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(short_increments) == 2,
              "increments: mismatched rates dimension should die");

    // ...or inferred from the shape of "rates".
    const char* short_inferred =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rates\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(short_inferred) == 2,
              "increments: mismatched rates dimension should die when inferred");

    // A simplex is a rate shape, not a sequence of increments.
    const char* simplex_increments =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_increments\",\"categories\":4,"
        "\"rates\":{\"id\":\"gaps\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(simplex_increments) == 2,
              "increments: a simplex should die");

    // This parameterization needs the ratios *and* the rate they hang off.
    const char* lone_ratios =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_ratios\",\"categories\":4,"
        "\"rates\":{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(lone_ratios) == 2,
              "ratios: a single rates parameter should die");

    // K-1 ratios for K categories.
    const char* wrong_ratios =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"parameterization\":\"rate_ratios\",\"categories\":4,"
        "\"rates\":[{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5,0.5],\"lower\":0},"
        "{\"id\":\"top\",\"type\":\"parameter\",\"x\":8.0,\"lower\":0}],"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(wrong_ratios) == 2,
              "ratios: mismatched ratio count should die");

    // Both ordered parameterizations write their first rate into category 0,
    // which an
    // invariant class would need pinned at 0.
    const char* increments_invariant =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"invariant\":true,"
        "\"parameterization\":\"rate_increments\",\"categories\":3,"
        "\"rates\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(increments_invariant) == 2,
              "increments: an invariant category should die");

    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_mean_contribution);
    mu_run_test(test_mean_contribution_invariant);
    mu_run_test(test_mean_contribution_dual);
    mu_run_test(test_mean_contribution_rejects);
    mu_run_test(test_rate_increments);
    mu_run_test(test_rate_ratios);
    mu_run_test(test_rate_parameterization_inferred);
    mu_run_test(test_rate_parameterization_rejects);
    return NULL;
}

RUN_TESTS(all_tests);
