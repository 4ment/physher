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
        "\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
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
        "\"distribution\":\"discrete\",\"invariant\":true,"
        "\"categories\":3,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.2,0.3,0.5]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.25,0.15,0.25,0.35]}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    const double s[3] = {0.2, 0.3, 0.5};
    const double p[4] = {0.25, 0.15, 0.25, 0.35};

    // "invariant" turns on +I, and "categories" counts the variable categories
    // only, so the proportions simplex is one longer than it.
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
        "\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    // x_k = (s_k/p_k) / sum_j (s_j/p_j), i.e. the same rates up to scale.
    const char* rate_shape =
        "{\"id\":\"sitemodel2\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"categories\":4,"
        "\"rate_shape\":{\"id\":\"x\",\"type\":\"simplex\","
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
        "\"categories\":4,"
        "\"rate_increments\":{\"id\":\"gaps\",\"type\":\"parameter\","
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
        "\"categories\":4,"
        "\"rate_ratios\":{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5],\"lower\":1.e-8,\"upper\":0.99},"
        "\"top_rate\":{\"id\":\"top\",\"type\":\"parameter\",\"x\":8.0,\"lower\":0},"
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

// The key selects the parameterization, so the same simplex means different things
// under "rate_shape" and "mean_contribution": x_k/sum_j p_j x_j against x_k/p_k.
// These two are indistinguishable by shape, which is exactly what the key resolves.
char* test_rate_parameterization_key_selects() {
    const char* shape =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rate_shape\":{\"id\":\"x\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    const char* contribution =
        "{\"id\":\"sitemodel2\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p2\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";

    Hashtable* hash = _new_hash();
    Model* model_shape = _sitemodel_from_string(shape, hash);
    Model* model_contribution = _sitemodel_from_string(contribution, hash);
    SiteModel* sm_shape = model_shape->obj;
    SiteModel* sm_contribution = model_contribution->obj;

    const double x[4] = {0.4, 0.3, 0.2, 0.1};
    const double p[4] = {0.1, 0.2, 0.3, 0.4};
    // sum_j p_j x_j = 0.04+0.06+0.06+0.04 = 0.2
    bool differ = false;
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm_shape->get_rate(sm_shape, i) - x[i] / 0.2) < TOL,
                  "rate shape: rates are not x_k normalised by the weighted mean");
        mu_assert(
            fabs(sm_contribution->get_rate(sm_contribution, i) - x[i] / p[i]) < TOL,
                  "mean contribution: rates are not s_k/p_k");
        differ |= fabs(sm_shape->get_rate(sm_shape, i) -
                       sm_contribution->get_rate(sm_contribution, i)) > TOL;
    }
    mu_assert(differ,
              "the two simplex parameterizations should not agree on this simplex");
    mu_assert(fabs(_weighted_mean(sm_shape) - 1.0) < TOL,
              "rate shape: rates do not have unit weighted mean");
    mu_assert(fabs(_weighted_mean(sm_contribution) - 1.0) < TOL,
              "mean contribution: rates do not have unit weighted mean");

    model_contribution->free(model_contribution);
    model_shape->free(model_shape);
    free_Hashtable(hash);
    return NULL;
}

// Parse a site model in a forked child and report its exit status, so the tests can
// check that a malformed rate parameterization is rejected instead of silently
// reading the wrong number of simplex elements. json_die exits 12; the semantic
// checks in _check_rate_parameterization exit 2.
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
        "\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(good) == 0, "mean contribution: valid model should parse");

    // An unknown parameterization is now an unknown key rather than an unknown
    // string, so the schema rejects it before anything is built.
    const char* unknown =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"cumprod\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(unknown) == 12,
              "mean contribution: unknown parameterization key should die");

    // The keys are alternative coordinate systems on the same rates, not layers.
    const char* two_keys =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
        "\"rate_shape\":{\"id\":\"x\",\"type\":\"simplex\",\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(two_keys) == 12,
              "two parameterization keys at once should die");

    // s must be a simplex: a plain positive vector does not sum to one, so the
    // unit mean would not hold.
    const char* not_simplex =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"parameter\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(not_simplex) == 2,
              "mean contribution: non-simplex rates should die");

    // One contribution per variable category, no more, no less.
    const char* wrong_dim =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\","
        "\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\",\"x\":[0.5,0.3,0.2]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(wrong_dim) == 2,
              "mean contribution: mismatched rates dimension should die");

    // It parameterises free rates; it means nothing for a rate distribution
    // whose categories come from a quantile function.
    const char* not_discrete =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0}}";
    mu_assert(_parse_status(not_discrete) == 2,
              "mean contribution: non-discrete distribution should die");

    // Nor is it silently ignored on a model that declares no distribution at all.
    const char* no_distribution =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\",\"categories\":4,"
        "\"mean_contribution\":{\"id\":\"s\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]}}";
    mu_assert(_parse_status(no_distribution) == 12,
              "mean contribution: no \"distribution\" should die");

    return NULL;
}

// "proportions" is checked against "categories" wherever it is accepted: it is
// copied straight into cat_proportions, which is sized from "categories" (plus the
// invariant class), so a mis-sized simplex runs off the end of that buffer.
char* test_proportions_rejects() {
    // One weight per category...
    const char* long_proportions =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"quadrature\":\"discrete\","
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.2,0.2,0.2,0.2,0.2]}}";
    mu_assert(_parse_status(long_proportions) == 2,
              "proportions: a simplex longer than \"categories\" should die");

    // ...and one more for the invariant class when there is one, no fewer.
    const char* short_proportions_invariant =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"invariant\":true,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"quadrature\":\"discrete\","
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.25,0.25,0.25,0.25]}}";
    mu_assert(_parse_status(short_proportions_invariant) == 2,
              "proportions +I: a simplex of \"categories\" elements should die");

    const char* good_invariant =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"invariant\":true,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"quadrature\":\"discrete\","
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.2,0.2,0.2,0.2,0.2]}}";
    mu_assert(_parse_status(good_invariant) == 0,
              "proportions +I: \"categories\"+1 elements should parse");

    // The weights weight the categories of a rate distribution. With no
    // distribution there is one category of weight 1 and nothing to weight.
    const char* no_distribution =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\",\"categories\":4,"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(no_distribution) == 12,
              "proportions: no \"distribution\" should die");

    // Two spellings of the same weights.
    const char* both =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"proportion_invariant\":0.2,"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\",\"x\":[0.2,0.8]}}";
    mu_assert(_parse_status(both) == 2,
              "\"proportions\" and \"proportion_invariant\" together should die");

    return NULL;
}

// The category count comes from the proportions simplex, so a parameter of the
// wrong size would be read off the end. Every shape is checked against the key
// that named it.
char* test_rate_parameterization_rejects() {
    // One increment per category, no more, no less.
    const char* short_increments =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rate_increments\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(short_increments) == 2,
              "increments: mismatched dimension should die");

    // A simplex is a rate shape, not a sequence of increments. Under the old
    // shape-inferred scheme this silently built a rate-shape model instead.
    const char* simplex_increments =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rate_increments\":{\"id\":\"gaps\",\"type\":\"simplex\","
        "\"x\":[0.4,0.3,0.2,0.1]},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(simplex_increments) == 2,
              "increments: a simplex should die");

    // This parameterization needs the ratios *and* the rate they hang off.
    const char* lone_ratios =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rate_ratios\":{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(lone_ratios) == 12,
              "ratios: \"rate_ratios\" without \"top_rate\" should die");

    // ...and the rate they hang off means nothing on its own.
    const char* lone_top =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"top_rate\":{\"id\":\"top\",\"type\":\"parameter\",\"x\":8.0,\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(lone_top) == 12,
              "ratios: \"top_rate\" without \"rate_ratios\" should die");

    // K-1 ratios for K categories.
    const char* wrong_ratios =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rate_ratios\":{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5,0.5],\"lower\":0},"
        "\"top_rate\":{\"id\":\"top\",\"type\":\"parameter\",\"x\":8.0,\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(wrong_ratios) == 2,
              "ratios: mismatched ratio count should die");

    // The rate the chain hangs off is a single rate, not a vector.
    const char* vector_top =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rate_ratios\":{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5],\"lower\":0},"
        "\"top_rate\":{\"id\":\"top\",\"type\":\"parameter\","
        "\"x\":[8.0,8.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(vector_top) == 2,
              "ratios: a vector \"top_rate\" should die");

    // With an invariant class there is one increment per *variable* category, one
    // fewer than the proportions simplex.
    const char* increments_invariant_wrong_dim =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"invariant\":true,\"categories\":3,"
        "\"rate_increments\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0,1.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(increments_invariant_wrong_dim) == 2,
              "increments +I: an increment for the invariant category should die");

    // The ratio chain still runs down to category 0, which +I needs pinned at 0.
    const char* ratios_invariant =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"invariant\":true,\"categories\":3,"
        "\"rate_ratios\":{\"id\":\"ratios\",\"type\":\"parameter\","
        "\"x\":[0.5,0.5,0.5],\"lower\":0},"
        "\"top_rate\":{\"id\":\"top\",\"type\":\"parameter\",\"x\":8.0,\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    mu_assert(_parse_status(ratios_invariant) == 2,
              "ratios: an invariant category should die");

    return NULL;
}

// An invariant class can be prepended to the increments, as to the mean-contribution
// simplex: category 0 is pinned at rate 0 and contributes nothing to the mean, so the
// running sum starts at category 1 and theta carries one element per variable
// category.
char* test_rate_increments_invariant() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"invariant\":true,\"categories\":3,"
        "\"rate_increments\":{\"id\":\"gaps\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0],\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.1,0.2,0.3,0.4]}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(sm->invariant, "increments +I: invariant class not detected");
    mu_assert(sm->cat_count == 4, "increments +I: wrong number of categories");
    mu_assert(sm->get_rate(sm, 0) == 0.0,
              "increments +I: invariant category is not rate 0");

    // Raw rates 0,1,2,3 with weighted mean 0.2*1 + 0.3*2 + 0.4*3 = 2.
    const double expected[4] = {0.0, 0.5, 1.0, 1.5};
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - expected[i]) < TOL,
                  "increments +I: rates not matching");
    }
    for (size_t i = 2; i < 4; i++) {
        mu_assert(sm->get_rate(sm, i) > sm->get_rate(sm, i - 1),
                  "increments +I: positive gaps must give increasing rates");
    }
    // The invariant category contributes nothing, so the mean over all categories
    // is still 1.
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "increments +I: rates do not have unit weighted mean");

    // Still scale-free, and the invariant category stays pinned across a change.
    const double doubled[3] = {2.0, 2.0, 2.0};
    Parameter_set_values(Hashtable_get(hash, "gaps"), doubled);
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - expected[i]) < TOL,
                  "increments +I: rates are not invariant to the scale of the gaps");
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_mean_contribution);
    mu_run_test(test_mean_contribution_invariant);
    mu_run_test(test_mean_contribution_dual);
    mu_run_test(test_mean_contribution_rejects);
    mu_run_test(test_rate_increments);
    mu_run_test(test_rate_increments_invariant);
    mu_run_test(test_rate_ratios);
    mu_run_test(test_rate_parameterization_key_selects);
    mu_run_test(test_rate_parameterization_rejects);
    mu_run_test(test_proportions_rejects);
    return NULL;
}

RUN_TESTS(all_tests);
