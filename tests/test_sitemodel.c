// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>

#include "minunit.h"
#include "phyc/em.h"
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

    // Weighting the categories individually only means something to the
    // "discrete" quadrature; every other rule derives the weights from the
    // discretisation. Left accepted, the median rule matches no branch of the rate
    // construction and every rate comes out infinite, and the beta rule reads the
    // first two elements of the simplex as if it were the invariant pair.
    const char* median_proportions =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"invariant\":true,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.2,0.2,0.2,0.2,0.2]}}";
    mu_assert(_parse_status(median_proportions) == 12,
              "proportions: the default median quadrature should die");

    const char* mean_proportions =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.25,0.25,0.25,0.25]}}";
    mu_assert(_parse_status(mean_proportions) == 12,
              "proportions: the mean quadrature should die");

    const char* beta_proportions =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"quadrature\":\"beta\","
        "\"shape\":{\"id\":\"shape\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"alpha\":{\"id\":\"a\",\"type\":\"parameter\",\"x\":1.0,\"lower\":0},"
        "\"beta\":{\"id\":\"b\",\"type\":\"parameter\",\"x\":1.0,\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.25,0.25,0.25,0.25]}}";
    mu_assert(_parse_status(beta_proportions) == 12,
              "proportions: the beta quadrature should die");

    // The free-rates model is not affected: "discrete" places no quantiles, so it
    // never reads "quadrature" and its weights are the simplex by definition.
    const char* free_rates =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"rate_increments\":{\"id\":\"r\",\"type\":\"parameter\","
        "\"x\":[1.0,1.0,1.0,1.0],\"lower\":0.0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.25,0.25,0.25,0.25]}}";
    mu_assert(_parse_status(free_rates) == 0,
              "proportions: +R should not need a \"quadrature\"");

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

// Conditional mean of a unit-mean exponential (a Weibull of shape 1) over the k-th
// equal-probability bin, times K. With e^{-u_j} = 1 - j/K at the bin boundaries
// u_j = -log(1 - j/K), the partial expectation collapses to a closed form, which
// gives an independent check on the incomplete-gamma construction.
static double _exponential_bin_mean(size_t k, size_t K) {
    // int_a^b r e^{-r} dr = (1 + a) e^{-a} - (1 + b) e^{-b}, and 1 + u_j is
    // 1 - log(e^{-u_j}) with e^{-u_j} = 1 - j/K. The last bin runs to infinity.
    double lower_mass = 1.0 - (double)(k - 1) / K;
    double upper_mass = 1.0 - (double)k / K;
    double lower = (1.0 - log(lower_mass)) * lower_mass;
    double upper = (k == K ? 0.0 : (1.0 - log(upper_mass)) * upper_mass);
    return K * (lower - upper);
}

// Mean quadrature places category k at the conditional mean of its equal-probability
// bin rather than at the bin's median. For the Weibull the substitution u = r^alpha
// turns that partial expectation into an incomplete gamma of order 1 + 1/alpha
// evaluated at fixed boundaries, so the rates come out with unit mean by
// construction and are never renormalized.
char* test_weibull_mean_quadrature() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"weibull\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":1.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(sm->cat_count == 4, "weibull mean: wrong number of categories");
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_proportion(sm, i) - 0.25) < TOL,
                  "weibull mean: bins are not equiprobable");
        mu_assert(fabs(sm->get_rate(sm, i) - _exponential_bin_mean(i + 1, 4)) < TOL,
                  "weibull mean: shape 1 rates are not the exponential bin means");
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "weibull mean: rates do not have unit weighted mean");

    // Away from the exponential there is no closed form; these are the same
    // integrals K * int r f(r) dr evaluated numerically outside physher.
    const double expected_half[4] = {0.012812320429430, 0.120440304532478,
                                     0.519546986081804, 3.347200388955839};
    const double quad_tol = 1.e-9;  // accuracy of the numerical reference
    sm->set_rate(sm, 0, 0.5);
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - expected_half[i]) < quad_tol,
                  "weibull mean: rates not matching after a change of shape");
        if (i > 0) {
            mu_assert(sm->get_rate(sm, i) > sm->get_rate(sm, i - 1),
                      "weibull mean: rates are not increasing");
        }
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "weibull mean: unit weighted mean lost after a change of shape");

    // The boundaries are constants in u and the shape only moves the order of the
    // incomplete gamma, so nothing can overflow the way a Weibull *quantile* does at
    // a small shape: the discretization degenerates to a single fast category
    // instead of returning NaN.
    sm->set_rate(sm, 0, 0.01);
    for (size_t i = 0; i < 4; i++) {
        mu_assert(isfinite(sm->get_rate(sm, i)),
                  "weibull mean: rate is not finite at the smallest allowed shape");
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "weibull mean: unit weighted mean lost at the smallest allowed shape");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The invariant class contributes nothing to the weighted mean, so the constraint
// still reads sum_{k>0} p_k r_k = 1 with p_k = (1 - p_inv)/K: the conditional means
// are scaled up by the variable proportion.
char* test_weibull_mean_quadrature_invariant() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"weibull\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":1.0,\"lower\":0},"
        "\"proportion_invariant\":{\"id\":\"pinv\",\"type\":\"parameter\","
        "\"x\":0.2,\"lower\":0.0,\"upper\":1.0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(sm->invariant, "weibull mean +I: invariant class not detected");
    mu_assert(sm->cat_count == 5, "weibull mean +I: wrong number of categories");
    mu_assert(sm->get_rate(sm, 0) == 0.0,
              "weibull mean +I: invariant category is not rate 0");
    mu_assert(fabs(sm->get_proportion(sm, 0) - 0.2) < TOL,
              "weibull mean +I: invariant proportion not matching");
    for (size_t i = 1; i < 5; i++) {
        mu_assert(fabs(sm->get_proportion(sm, i) - 0.8 / 4.0) < TOL,
                  "weibull mean +I: variable categories do not share 1 - p_inv");
        mu_assert(fabs(sm->get_rate(sm, i) - _exponential_bin_mean(i, 4) / 0.8) < TOL,
                  "weibull mean +I: rates are not scaled by the variable proportion");
    }
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "weibull mean +I: rates do not have unit weighted mean");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// A synthetic objective L = sum_k p_k c_k r_k^2 over the category rates, in the
// convention treelikelihood.c passes gradients back in: ingrad[k] is dL/dr_k with the
// category weight p_k factored out, and ingrad[0] is the part of dL/dp_inv that runs
// through the proportions rather than through the rates.
static const double _grad_weights[5] = {0.3, 1.1, -0.7, 2.3, 0.5};

static double _quadratic_objective(SiteModel* sm) {
    double L = 0;
    for (size_t i = 0; i < sm->cat_count; i++) {
        double r = sm->get_rate(sm, i);
        L += sm->get_proportion(sm, i) * _grad_weights[i] * r * r;
    }
    return L;
}

static void _fill_rate_gradient(SiteModel* sm, double* ingrad) {
    for (size_t i = 0; i < sm->cat_count; i++) {
        ingrad[i] = _grad_weights[i] * 2.0 * sm->get_rate(sm, i);
    }
}

// The mean-quadrature rates are partial expectations, not quantiles, so they need a
// derivative of their own -- the median-quadrature formulas would return a plausible
// but wrong number. Check both parameters against central differences of the same
// objective.
static char* _check_mean_quadrature_gradient(const char* json, bool invariant) {
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;
    double ingrad[5];

    Parameter* shape = Parameters_at(sm->rates, 0);
    double alpha = Parameter_value(shape);
    _fill_rate_gradient(sm, ingrad);
    double analytic = sm->derivative(sm, ingrad, shape);

    double h = 1.e-6 * alpha;
    sm->set_rate(sm, 0, alpha + h);
    double plus = _quadratic_objective(sm);
    sm->set_rate(sm, 0, alpha - h);
    double minus = _quadratic_objective(sm);
    sm->set_rate(sm, 0, alpha);
    mu_assert(fabs(analytic - (plus - minus) / (2.0 * h)) < 1.e-6,
              "mean quadrature: shape derivative does not match finite differences");

    if (invariant) {
        Parameter* pinv = sm->proportions->transform->parameter;
        double p0 = Parameter_value(pinv);
        // dp_0/dp_inv = 1 and dp_k/dp_inv = -1/K for the variable categories.
        double through_proportions = 0;
        for (size_t i = 0; i < sm->cat_count; i++) {
            double r = sm->get_rate(sm, i);
            double dp = (i == 0 ? 1.0 : -1.0 / (sm->cat_count - 1));
            through_proportions += dp * _grad_weights[i] * r * r;
        }
        _fill_rate_gradient(sm, ingrad);
        ingrad[0] = through_proportions;
        analytic = sm->derivative(sm, ingrad, pinv);

        h = 1.e-6;
        Parameter_set_value(pinv, p0 + h);
        plus = _quadratic_objective(sm);
        Parameter_set_value(pinv, p0 - h);
        minus = _quadratic_objective(sm);
        Parameter_set_value(pinv, p0);
        mu_assert(fabs(analytic - (plus - minus) / (2.0 * h)) < 1.e-6,
                  "mean quadrature: p_inv derivative does not match finite differences");
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

char* test_mean_quadrature_gradient() {
    char* result = _check_mean_quadrature_gradient(
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"weibull\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.7,\"lower\":0}}",
        false);
    if (result != NULL) return result;

    result = _check_mean_quadrature_gradient(
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"weibull\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.7,\"lower\":0},"
        "\"proportion_invariant\":{\"id\":\"pinv\",\"type\":\"parameter\","
        "\"x\":0.2,\"lower\":0.0,\"upper\":1.0}}",
        true);
    if (result != NULL) return result;

    // The gamma shares the construction, differing only in that its bin boundaries
    // move with the shape.
    result = _check_mean_quadrature_gradient(
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.9,\"lower\":0},"
        "\"proportion_invariant\":{\"id\":\"pinv\",\"type\":\"parameter\","
        "\"x\":0.3,\"lower\":0.0,\"upper\":1.0}}",
        true);
    return result;
}

// The conditional mean of a bin is derived per family; only the gamma and the Weibull
// have one here, and the others must be rejected rather than silently discretized as
// a gamma.
char* test_mean_quadrature_rejects() {
    const char* weibull =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"weibull\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0}}";
    mu_assert(_parse_status(weibull) == 0, "mean quadrature: weibull should parse");

    const char* lognormal =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"lognormal\",\"categories\":4,\"quadrature\":\"mean\","
        "\"scale\":{\"id\":\"sigma\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0}}";
    mu_assert(_parse_status(lognormal) == 12,
              "mean quadrature: lognormal should die");

    return NULL;
}

// ---------------------------------------------------------------------------
// Gauss-Laguerre quadrature

// The gamma has mean 1, so its density is alpha^alpha r^(alpha-1) e^(-alpha r) /
// Gamma(alpha). Substituting z = alpha*r turns E[g(r)] into the Gauss-Laguerre
// weight function z^(alpha-1) e^(-z) integrated against g(z/alpha), divided by
// Gamma(alpha). So the n-point rule with alf = alpha - 1 gives the categories
// directly: r_k = z_k / alpha and p_k = w_k / Gamma(alpha), which is exactly what
// _gamma_approx_laguerre does. Unlike the quantile discretizations the rates are
// not equiprobable and are never renormalized -- the unit mean falls out of the
// rule being exact for polynomials.
char* test_gausslaguerre_quadrature() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"quadrature\":\"gausslaguerre\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":1.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(sm->cat_count == 4, "gausslaguerre: wrong number of categories");

    // Shape 1 is alf = 0 and Gamma(1) = 1, so the categories must be the textbook
    // 4-point Gauss-Laguerre abscissas and weights unchanged.
    const double nodes[4] = {0.322547689619392, 1.745761101158347,
                             4.536620296921128, 9.395070912301133};
    const double weights[4] = {0.603154104341634, 0.357418692437800,
                               0.038887908515005, 0.000539294705561};
    const double quad_tol = 1.e-12;  // accuracy of the tabulated reference
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - nodes[i]) < quad_tol,
                  "gausslaguerre: shape 1 rates are not the Laguerre abscissas");
        mu_assert(fabs(sm->get_proportion(sm, i) - weights[i]) < quad_tol,
                  "gausslaguerre: shape 1 proportions are not the Laguerre weights");
    }

    // An n-point rule integrates polynomials up to degree 2n-1 exactly, so with 4
    // categories the first 7 moments of the gamma must come out to machine
    // precision: E[r^k] = Gamma(alpha+k) / (Gamma(alpha) alpha^k). k = 0 is the
    // proportions summing to 1 and k = 1 is the unit mean.
    const double shapes[3] = {1.0, 0.1, 7.5};
    for (size_t s = 0; s < 3; s++) {
        const double alpha = shapes[s];
        sm->set_rate(sm, 0, alpha);
        for (size_t k = 0; k <= 7; k++) {
            double moment = 0;
            for (size_t i = 0; i < 4; i++) {
                moment += sm->get_proportion(sm, i) * pow(sm->get_rate(sm, i), k);
            }
            // Gamma(alpha+k)/(Gamma(alpha) alpha^k) built up as a product to keep
            // it accurate at a large shape.
            double expected = 1.0;
            for (size_t j = 0; j < k; j++) {
                expected *= (alpha + j) / alpha;
            }
            mu_assert(fabs(moment - expected) < 1.e-9 * expected,
                      "gausslaguerre: moment does not match the gamma");
        }
        for (size_t i = 0; i < 4; i++) {
            mu_assert(sm->get_rate(sm, i) > 0,
                      "gausslaguerre: rate is not positive");
            mu_assert(sm->get_proportion(sm, i) > 0,
                      "gausslaguerre: proportion is not positive");
            if (i > 0) {
                mu_assert(sm->get_rate(sm, i) > sm->get_rate(sm, i - 1),
                          "gausslaguerre: rates are not increasing");
            }
        }
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// ---------------------------------------------------------------------------
// EM for free rates (+R)

static Model* _treelikelihood_from_file(const char* file, Hashtable* hash) {
    char* content = load_file(file);
    json_node* json = create_json_tree(content);
    free(content);
    Model* model = new_TreeLikelihoodModel_from_json(json->children[0], hash);
    json_free_tree(json);
    return model;
}

// Everything the parameterization guarantees and EM must not break: rates
// increasing (so the increments stay positive), weights a simplex, unit weighted
// mean.
static char* _check_freerate_invariants(SiteModel* sm) {
    double sum = 0;
    for (size_t i = 0; i < sm->cat_count; i++) {
        double p = sm->get_proportion(sm, i);
        mu_assert(p > 0.0 && p < 1.0, "EM: a weight left the simplex");
        sum += p;
        if (i > 0) {
            mu_assert(sm->get_rate(sm, i) >= sm->get_rate(sm, i - 1),
                      "EM: the rates are no longer increasing");
        }
    }
    mu_assert(fabs(sum - 1.0) < 1.e-8, "EM: the weights do not sum to 1");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < 1.e-8,
              "EM: the rates lost their unit weighted mean");
    return NULL;
}

// Each step is a full E/M sweep and must not decrease the likelihood -- the one
// property EM is bought for. Run them one at a time so every step is checked,
// rather than only the endpoints.
static char* _check_freerate_em(const char* file) {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file(file, hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    double logP = model->logP(model);
    char* result = _check_freerate_invariants(sm);
    if (result != NULL) return result;

    double start = logP;
    for (int step = 0; step < 6; step++) {
        SiteModelEM em = SiteModel_optimize_freerate_EM(tlk, 1, 1.e-8);
        mu_assert(!isnan(em.logP), "EM refused a model it is meant to handle");
        mu_assert(em.steps == 1, "EM did not run the step it was asked for");
        // What it reports is what it leaves behind.
        mu_assert(fabs(model->logP(model) - em.logP) < 1.e-8,
                  "EM: the reported logP is not the model's");
        mu_assert(em.logP > logP - 1.e-6, "EM: a step decreased the log-likelihood");
        logP = em.logP;
        result = _check_freerate_invariants(sm);
        if (result != NULL) return result;
    }
    mu_assert(logP > start + 1.0, "EM: six steps bought nothing");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

char* test_freerate_em() { return _check_freerate_em("jc69-freerate.json"); }

// Weights-only EM (§5.3): with the increments fixed, the closed-form weight
// update still moves the normalized rates, so the branch-length rescaling has to
// happen here too or the step is not monotone.
char* test_freerate_em_fixed_rates() {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-freerate.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;
    Parameter* increments = Parameters_at(sm->rates, 0);
    Parameter_set_estimate(increments, false);

    double* before = clone_dvector(Parameter_values(increments),
                                  Parameter_size(increments));
    double logP = model->logP(model);
    for (int step = 0; step < 4; step++) {
        SiteModelEM em = SiteModel_optimize_freerate_EM(tlk, 1, 1.e-8);
        mu_assert(!isnan(em.logP), "EM refused a weights-only model");
        mu_assert(em.logP > logP - 1.e-6,
                  "EM: a weights-only step decreased the log-likelihood");
        logP = em.logP;
    }
    const double* after = Parameter_values(increments);
    for (size_t i = 0; i < Parameter_size(increments); i++) {
        mu_assert(after[i] == before[i], "EM moved increments it was told to fix");
    }
    char* result = _check_freerate_invariants(sm);
    if (result != NULL) return result;

    free(before);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// With an invariant class, category 0 is pinned at rate 0 and outside the
// running sum of the increments; EM estimates its weight and nothing else.
char* test_freerate_em_invariant() {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-freerate-invariant.json", hash);
    SingleTreeLikelihood* tlk = model->obj;
    SiteModel* sm = tlk->sm;

    double before = model->logP(model);
    SiteModelEM em = SiteModel_optimize_freerate_EM(tlk, 20, 1.e-5);
    mu_assert(!isnan(em.logP), "EM refused a +R+I model");
    mu_assert(em.logP > before, "EM did not improve a +R+I model");
    mu_assert(sm->get_rate(sm, 0) == 0.0, "EM moved the invariant category off 0");
    char* result = _check_freerate_invariants(sm);
    if (result != NULL) return result;

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// Anything that is not a free-rate mixture on ordered increments has to come
// back untouched rather than have some other parameter optimized in its place.
char* test_freerate_em_rejects() {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-distance.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    double before = model->logP(model);
    SiteModelEM em = SiteModel_optimize_freerate_EM(tlk, 10, 1.e-5);
    mu_assert(isnan(em.logP), "EM accepted a single-category site model");
    mu_assert(em.steps == 0, "EM did work on a model it rejected");
    mu_assert(model->logP(model) == before, "EM changed a model it rejected");

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The M-step evaluates through the mixture's own likelihood with a one-category
// site model swapped in for tlk->sm and the posterior counts w_ic swapped in for
// the pattern weights. Both have to be put back: what EM borrows it does not own,
// and the caller's next calculate() must be the mixture's, not the last class's.
char* test_freerate_em_restores_likelihood() {
    Hashtable* hash = _new_hash();
    Model* model = _treelikelihood_from_file("jc69-freerate.json", hash);
    SingleTreeLikelihood* tlk = model->obj;

    model->logP(model);
    SiteModel* sm = tlk->sm;
    int cat_count = tlk->cat_count;
    int cat_capacity = tlk->cat_capacity;
    bool use_upper = tlk->use_upper;
    int node_id = tlk->node_id;
    double* weights = clone_dvector(tlk->sp->weights, tlk->sp->count);

    SiteModelEM em = SiteModel_optimize_freerate_EM(tlk, 4, 1.e-8);
    mu_assert(!isnan(em.logP), "EM refused a model it is meant to handle");
    mu_assert(em.evaluations > 0, "EM ran no M-step, so nothing was swapped");

    mu_assert(tlk->sm == sm, "EM left its one-category site model installed");
    mu_assert(tlk->cat_count == cat_count, "EM did not restore cat_count");
    mu_assert(tlk->cat_capacity == cat_capacity, "EM moved the buffer capacity");
    mu_assert(tlk->use_upper == use_upper, "EM did not restore use_upper");
    mu_assert(tlk->node_id == node_id, "EM did not restore node_id");
    for (int i = 0; i < tlk->sp->count; i++) {
        mu_assert(tlk->sp->weights[i] == weights[i],
                  "EM left a class's posterior counts in the pattern weights");
    }
    // The partials the M-step overwrote are the mixture's, so a plain recompute
    // has to agree with what EM reported rather than with the last class.
    mu_assert(fabs(tlk->calculate(tlk) - em.logP) < 1.e-8,
              "EM left the likelihood holding a one-category value");

    free(weights);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// A discretised gamma plus a free rate class (+G+F): category 0 is a point mass
// of weight p_f whose rate is estimated instead of being pinned at 0. The rate is
// given here as the class's contribution to the unit mean, c = p_f r_f, so that
// r_f = c/p_f and the gamma categories -- which the quadrature normalised to
// carry the whole mean -- are scaled by 1-c to leave room for it.
#define FREE_CLASS_P 0.25
#define FREE_CLASS_C 0.4
static const char* FREE_CLASS_G4 =
    "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
    "\"distribution\":\"gamma\",\"categories\":4,"
    "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
    "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
    "\"x\":0.25,\"lower\":0,\"upper\":1},"
    "\"free_class_contribution\":{\"id\":\"cf\",\"type\":\"parameter\","
    "\"x\":0.4,\"lower\":0,\"upper\":1}}";

// The same shape and weight, with category 0 pinned at rate 0: the +G+I model the
// free class generalises. Different ids so the two can share a hashtable.
static const char* INVARIANT_G4 =
    "{\"id\":\"sitemodel2\",\"type\":\"sitemodel\","
    "\"distribution\":\"gamma\",\"categories\":4,"
    "\"shape\":{\"id\":\"alpha2\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
    "\"proportion_invariant\":{\"id\":\"pinv\",\"type\":\"parameter\","
    "\"x\":0.25,\"lower\":0,\"upper\":1}}";

char* test_free_class() {
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(FREE_CLASS_G4, hash);
    SiteModel* sm = model->obj;

    // "categories" counts the variable categories, so the point mass is extra.
    mu_assert(sm->cat_count == 5, "+F: wrong number of categories");
    mu_assert(sm->invariant, "+F: category 0 is not a separate class");
    mu_assert(sm->class_rate != NULL, "+F: the class rate was not attached");
    mu_assert(sm->class_rate_parameterization ==
                  FREE_CLASS_PARAMETERIZATION_CONTRIBUTION,
              "+F: \"free_class_contribution\" did not select the contribution chart");

    mu_assert(fabs(sm->get_proportion(sm, 0) - FREE_CLASS_P) < TOL,
              "+F: the free class does not have weight p_f");
    mu_assert(fabs(sm->get_rate(sm, 0) - FREE_CLASS_C / FREE_CLASS_P) < TOL,
              "+F: the free class is not at rate c/p_f");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F: rates do not have unit weighted mean");
    // c > p_f is exactly the condition for the class to be faster than average,
    // which is the point of the model.
    mu_assert(sm->get_rate(sm, 0) > 1.0,
              "+F: c > p_f should put the class above the mean rate");

    // The gamma categories are the +G+I ones scaled by 1-c: the free class takes
    // over that share of the mean, and nothing else about the discretisation moves.
    Model* invariant_model = _sitemodel_from_string(INVARIANT_G4, hash);
    SiteModel* ism = invariant_model->obj;
    for (size_t i = 1; i < 5; i++) {
        mu_assert(fabs(sm->get_proportion(sm, i) - ism->get_proportion(ism, i)) < TOL,
                  "+F: the variable categories do not keep the +I weights");
        mu_assert(fabs(sm->get_rate(sm, i) -
                       (1.0 - FREE_CLASS_C) * ism->get_rate(ism, i)) < TOL,
                  "+F: the variable rates are not the +I rates scaled by 1-c");
    }

    // A change to either new parameter has to invalidate the cached categories,
    // and the constraint has to survive it.
    Parameter_set_value(Hashtable_get(hash, "cf"), 0.1);
    mu_assert(fabs(sm->get_rate(sm, 0) - 0.1 / FREE_CLASS_P) < TOL,
              "+F: rates not updated after a change to the contribution");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F: unit mean lost after a change to the contribution");

    Parameter_set_value(Hashtable_get(hash, "pf"), 0.5);
    mu_assert(fabs(sm->get_proportion(sm, 0) - 0.5) < TOL,
              "+F: weights not updated after a change to the proportion");
    mu_assert(fabs(sm->get_rate(sm, 0) - 0.1 / 0.5) < TOL,
              "+F: rates not updated after a change to the proportion");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F: unit mean lost after a change to the proportion");

    invariant_model->free(invariant_model);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The rate and the contribution are two charts on the same model, related by
// c = p_f r_f. Given matching values they must produce identical categories.
char* test_free_class_charts_agree() {
    const char* rate_chart =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.25,\"lower\":0,\"upper\":1},"
        // r_f = c/p_f = 0.4/0.25
        "\"free_class_rate\":{\"id\":\"rf\",\"type\":\"parameter\","
        "\"x\":1.6,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* rate_model = _sitemodel_from_string(rate_chart, hash);
    SiteModel* rsm = rate_model->obj;
    mu_assert(rsm->class_rate_parameterization == FREE_CLASS_PARAMETERIZATION_RATE,
              "+F: \"free_class_rate\" selected the contribution chart");

    Hashtable* hash2 = _new_hash();
    Model* contribution_model = _sitemodel_from_string(FREE_CLASS_G4, hash2);
    SiteModel* csm = contribution_model->obj;

    mu_assert(rsm->cat_count == csm->cat_count, "+F: charts disagree on the categories");
    for (size_t i = 0; i < rsm->cat_count; i++) {
        mu_assert(fabs(rsm->get_rate(rsm, i) - csm->get_rate(csm, i)) < TOL,
                  "+F: the two charts give different rates");
        mu_assert(fabs(rsm->get_proportion(rsm, i) - csm->get_proportion(csm, i)) < TOL,
                  "+F: the two charts give different weights");
    }
    mu_assert(fabs(_weighted_mean(rsm) - 1.0) < TOL,
              "+F: the rate chart loses the unit mean");

    contribution_model->free(contribution_model);
    free_Hashtable(hash2);
    rate_model->free(rate_model);
    free_Hashtable(hash);
    return NULL;
}

// r_f = 0 is not a limit but a point of the parameter space: it is the invariant
// class, and it must reproduce +G+I exactly rather than approximately.
char* test_free_class_reduces_to_invariant() {
    const char* pinned =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.25,\"lower\":0,\"upper\":1},"
        "\"free_class_rate\":{\"id\":\"rf\",\"type\":\"parameter\","
        "\"x\":0.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(pinned, hash);
    SiteModel* sm = model->obj;
    Model* invariant_model = _sitemodel_from_string(INVARIANT_G4, hash);
    SiteModel* ism = invariant_model->obj;

    mu_assert(sm->get_rate(sm, 0) == 0.0, "+F: r_f = 0 is not rate 0");
    for (size_t i = 0; i < 5; i++) {
        mu_assert(sm->get_rate(sm, i) == ism->get_rate(ism, i),
                  "+F: r_f = 0 does not reproduce +G+I exactly");
        mu_assert(sm->get_proportion(sm, i) == ism->get_proportion(ism, i),
                  "+F: r_f = 0 does not reproduce the +G+I weights exactly");
    }

    invariant_model->free(invariant_model);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The scaling is applied to whatever the quadrature produced, so the mean
// quadrature -- whose rates are conditional means rather than quantiles, and are
// never renormalised -- takes the free class on the same terms.
char* test_free_class_mean_quadrature() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.25,\"lower\":0,\"upper\":1},"
        "\"free_class_contribution\":{\"id\":\"cf\",\"type\":\"parameter\","
        "\"x\":0.4,\"lower\":0,\"upper\":1}}";
    const char* invariant_json =
        "{\"id\":\"sitemodel2\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"quadrature\":\"mean\","
        "\"shape\":{\"id\":\"alpha2\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"proportion_invariant\":{\"id\":\"pinv\",\"type\":\"parameter\","
        "\"x\":0.25,\"lower\":0,\"upper\":1}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;
    Model* invariant_model = _sitemodel_from_string(invariant_json, hash);
    SiteModel* ism = invariant_model->obj;

    mu_assert(fabs(sm->get_rate(sm, 0) - FREE_CLASS_C / FREE_CLASS_P) < TOL,
              "+F mean quadrature: the free class is not at rate c/p_f");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F mean quadrature: rates do not have unit weighted mean");
    for (size_t i = 1; i < 5; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) -
                       (1.0 - FREE_CLASS_C) * ism->get_rate(ism, i)) < TOL,
                  "+F mean quadrature: rates are not the +I rates scaled by 1-c");
    }

    invariant_model->free(invariant_model);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The rate chart's real constraint is the joint p_f r_f < 1, which no box on r_f
// alone can express: beyond it the gamma categories would have to carry a
// negative share of the mean. That is reported as a failed update -- the same
// answer a non-finite quantile gives -- and not as silently negative rates.
char* test_free_class_domain() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.5,\"lower\":0,\"upper\":1},"
        // p_f r_f = 1.5 > 1
        "\"free_class_rate\":{\"id\":\"rf\",\"type\":\"parameter\","
        "\"x\":3.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(sm->update(sm) == false,
              "+F: p_f r_f > 1 should fail the update");

    // p_f r_f = 1 exactly is the same failure: the gamma categories would all
    // collapse onto rate 0.
    Parameter_set_value(Hashtable_get(hash, "rf"), 2.0);
    mu_assert(sm->update(sm) == false, "+F: p_f r_f = 1 should fail the update");

    // ...and back inside the domain it recovers.
    Parameter_set_value(Hashtable_get(hash, "rf"), 1.0);
    mu_assert(sm->update(sm) == true, "+F: p_f r_f < 1 should update");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F: unit mean not restored inside the domain");
    for (size_t i = 1; i < sm->cat_count; i++) {
        mu_assert(sm->get_rate(sm, i) > 0.0, "+F: a variable rate is not positive");
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The fastest of the variable rates, which is what an increment is measured from.
static double _fastest_variable_rate(SiteModel* sm) {
    double g = 0.0;
    for (size_t i = 1; i < sm->cat_count; i++) {
        if (sm->get_rate(sm, i) > g) g = sm->get_rate(sm, i);
    }
    return g;
}

// The increment chart ("fast class"): the free class is placed at d above the
// fastest variable category instead of anywhere on the line. The placement is
// implicit -- the rescaling that makes room for the class moves the very rate it
// is measured from -- so the test is that the defining identity r_f = r_max + d
// holds of the rates the model actually reports.
#define FREE_CLASS_D 2.0
char* test_free_class_increment() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.25,\"lower\":0,\"upper\":1},"
        "\"free_class_increment\":{\"id\":\"df\",\"type\":\"parameter\","
        "\"x\":2.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(sm->class_rate_parameterization == FREE_CLASS_PARAMETERIZATION_INCREMENT,
              "+F: \"free_class_increment\" did not select the increment chart");
    mu_assert(fabs(sm->get_rate(sm, 0) -
                   (_fastest_variable_rate(sm) + FREE_CLASS_D)) < TOL,
              "+F increment: the class is not d above the fastest variable rate");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F increment: rates do not have unit weighted mean");

    // The variable categories are still the +G+I ones scaled by 1-c; the only
    // difference from the other charts is which equation fixes c.
    Model* invariant_model = _sitemodel_from_string(INVARIANT_G4, hash);
    SiteModel* ism = invariant_model->obj;
    const double c = FREE_CLASS_P * sm->get_rate(sm, 0);
    for (size_t i = 1; i < 5; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - (1.0 - c) * ism->get_rate(ism, i)) < TOL,
                  "+F increment: the variable rates are not the +I rates scaled by 1-c");
    }

    // Moving either parameter has to keep the identity, not just the value it
    // happened to have at construction: p_f rescales the bulk the class is
    // measured from, so r_f moves even when d does not.
    Parameter_set_value(Hashtable_get(hash, "pf"), 0.4);
    mu_assert(fabs(sm->get_rate(sm, 0) -
                   (_fastest_variable_rate(sm) + FREE_CLASS_D)) < TOL,
              "+F increment: identity lost after a change to the proportion");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F increment: unit mean lost after a change to the proportion");

    Parameter_set_value(Hashtable_get(hash, "df"), 1.0);
    mu_assert(fabs(sm->get_rate(sm, 0) - (_fastest_variable_rate(sm) + 1.0)) < TOL,
              "+F increment: identity lost after a change to the increment");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F increment: unit mean lost after a change to the increment");

    invariant_model->free(invariant_model);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// d = 0 is the boundary of the fast class, where it sits exactly on the fastest
// variable category: the model is still well defined there (it is a +G with one
// category duplicated), and it is the closest the chart gets to +G+I -- which,
// unlike the other two charts, it cannot reach at all.
char* test_free_class_increment_boundary() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.25,\"lower\":0,\"upper\":1},"
        "\"free_class_increment\":{\"id\":\"df\",\"type\":\"parameter\","
        "\"x\":0.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(fabs(sm->get_rate(sm, 0) - _fastest_variable_rate(sm)) < TOL,
              "+F increment: d = 0 does not put the class on the fastest category");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F increment: d = 0 loses the unit mean");
    // Which is the whole point of the chart: no value of d puts the class below
    // the bulk, so a fit cannot wander into the slow mode.
    for (double d = 0.0; d < 3.5; d += 0.37) {
        Parameter_set_value(Hashtable_get(hash, "df"), d);
        const double rf = sm->get_rate(sm, 0);
        for (size_t i = 1; i < sm->cat_count; i++) {
            mu_assert(rf >= sm->get_rate(sm, i) - TOL,
                      "+F increment: the class fell below a variable rate");
        }
        mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
                  "+F increment: unit mean lost along the increment");
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The increment is a restriction of the free model, not a different one: the
// categories it produces are those of the rate chart at the r_f it lands on. This
// is what checks the closed-form solve of r_f = (1-p_f r_f) g + d.
char* test_free_class_increment_matches_rate() {
    const char* increment =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.25,\"lower\":0,\"upper\":1},"
        "\"free_class_increment\":{\"id\":\"df\",\"type\":\"parameter\","
        "\"x\":2.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(increment, hash);
    SiteModel* sm = model->obj;

    char rate_json[512];
    snprintf(rate_json, sizeof(rate_json),
             "{\"id\":\"sitemodel2\",\"type\":\"sitemodel\","
             "\"distribution\":\"gamma\",\"categories\":4,"
             "\"shape\":{\"id\":\"alpha2\",\"type\":\"parameter\",\"x\":0.5,"
             "\"lower\":0},"
             "\"free_class_proportion\":{\"id\":\"pf2\",\"type\":\"parameter\","
             "\"x\":0.25,\"lower\":0,\"upper\":1},"
             "\"free_class_rate\":{\"id\":\"rf2\",\"type\":\"parameter\","
             "\"x\":%.17g,\"lower\":0}}",
             sm->get_rate(sm, 0));
    Model* rate_model = _sitemodel_from_string(rate_json, hash);
    SiteModel* rsm = rate_model->obj;

    for (size_t i = 0; i < sm->cat_count; i++) {
        mu_assert(fabs(sm->get_rate(sm, i) - rsm->get_rate(rsm, i)) < TOL,
                  "+F increment: rates differ from the rate chart at the same r_f");
        mu_assert(fabs(sm->get_proportion(sm, i) - rsm->get_proportion(rsm, i)) < TOL,
                  "+F increment: weights differ from the rate chart at the same r_f");
    }

    rate_model->free(rate_model);
    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

// The increment's domain reads p_f d < 1 rather than p_f r_f < 1: the two say the
// same thing -- c -> 1 exactly as p_f d -> 1 -- but in the increment chart the
// bulk term cancels, so the bound sits at a fixed 1/p_f instead of moving with
// the shape. Beyond it the variable rates would go negative.
char* test_free_class_increment_domain() {
    const char* json =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":{\"id\":\"pf\",\"type\":\"parameter\","
        "\"x\":0.5,\"lower\":0,\"upper\":1},"
        "\"free_class_increment\":{\"id\":\"df\",\"type\":\"parameter\","
        "\"x\":1.0,\"lower\":0}}";
    Hashtable* hash = _new_hash();
    Model* model = _sitemodel_from_string(json, hash);
    SiteModel* sm = model->obj;

    mu_assert(sm->update(sm) == true, "+F increment: p_f d < 1 should update");
    mu_assert(0.5 * sm->get_rate(sm, 0) < 1.0,
              "+F increment: p_f d < 1 should imply p_f r_f < 1");
    mu_assert(fabs(_weighted_mean(sm) - 1.0) < TOL,
              "+F increment: unit mean lost");

    Parameter_set_value(Hashtable_get(hash, "df"), 2.0);
    mu_assert(sm->update(sm) == false, "+F increment: p_f d = 1 should fail the update");

    Parameter_set_value(Hashtable_get(hash, "df"), 3.0);
    mu_assert(sm->update(sm) == false, "+F increment: p_f d > 1 should fail the update");

    Parameter_set_value(Hashtable_get(hash, "df"), 1.5);
    mu_assert(sm->update(sm) == true, "+F increment: p_f d < 1 should update again");
    for (size_t i = 1; i < sm->cat_count; i++) {
        mu_assert(sm->get_rate(sm, i) > 0.0,
                  "+F increment: a variable rate is not positive");
    }

    model->free(model);
    free_Hashtable(hash);
    return NULL;
}

char* test_free_class_rejects() {
    const char* good =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6}";
    mu_assert(_parse_status(good) == 0, "+F: a valid model should parse");

    // The three keys stand or fall together: a rate with no weight is not a
    // mixture component, and a weight with no rate does not say where it sits.
    const char* no_proportion =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_rate\":1.6}";
    mu_assert(_parse_status(no_proportion) == 12,
              "+F: a rate with no weight should die");

    const char* no_rate =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25}";
    mu_assert(_parse_status(no_rate) == 12, "+F: a weight with no rate should die");

    const char* increment =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_increment\":2.0}";
    mu_assert(_parse_status(increment) == 0, "+F: the increment chart should parse");

    // Charts on the same rate, not layers on top of each other -- any two of the
    // three over-determine it.
    const char* both_charts =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6,"
        "\"free_class_contribution\":0.4}";
    mu_assert(_parse_status(both_charts) == 12,
              "+F: both rate charts together should die");

    const char* rate_and_increment =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6,"
        "\"free_class_increment\":2.0}";
    mu_assert(_parse_status(rate_and_increment) == 12,
              "+F: a rate and an increment together should die");

    // Each of these names the weight of a different category 0.
    const char* with_proportion_invariant =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"proportion_invariant\":0.1,"
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6}";
    mu_assert(_parse_status(with_proportion_invariant) == 12,
              "+F: \"proportion_invariant\" alongside a free class should die");

    const char* with_proportions =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,\"invariant\":true,"
        "\"quadrature\":\"discrete\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"proportions\":{\"id\":\"p\",\"type\":\"simplex\","
        "\"x\":[0.2,0.2,0.2,0.2,0.2]},"
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6}";
    mu_assert(_parse_status(with_proportions) == 12,
              "+F: \"proportions\" alongside a free class should die");

    // The free class sits beside a rate distribution; on its own it would be a
    // two-point mixture, which is what "discrete" already is.
    const char* no_distribution =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6}";
    mu_assert(_parse_status(no_distribution) == 12,
              "+F: no \"distribution\" should die");

    const char* discrete =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"discrete\",\"categories\":4,"
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6}";
    mu_assert(_parse_status(discrete) == 12,
              "+F: \"distribution\": \"discrete\" should die");

    // Gauss-Laguerre takes its weights from the quadrature rule, so there is no
    // separately weighted class to add.
    const char* laguerre =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"quadrature\":\"gausslaguerre\","
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_rate\":1.6}";
    mu_assert(_parse_status(laguerre) == 13,
              "+F: Gauss-Laguerre quadrature should die");

    // A contribution is a share of the unit mean, so it is bounded whatever the
    // weight is; a rate only has to be non-negative.
    const char* contribution_too_big =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_contribution\":1.4}";
    mu_assert(_parse_status(contribution_too_big) == 2,
              "+F: a contribution outside (0,1) should die");

    const char* negative_rate =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_rate\":-1.0}";
    mu_assert(_parse_status(negative_rate) == 2,
              "+F: a negative rate should die");

    // A negative increment would put the class inside the bulk, which is the one
    // thing this chart exists to rule out; ask for the rate chart instead.
    const char* negative_increment =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":0.25,\"free_class_increment\":-1.0}";
    mu_assert(_parse_status(negative_increment) == 2,
              "+F: a negative increment should die");

    const char* weight_out_of_range =
        "{\"id\":\"sitemodel\",\"type\":\"sitemodel\","
        "\"distribution\":\"gamma\",\"categories\":4,"
        "\"shape\":{\"id\":\"alpha\",\"type\":\"parameter\",\"x\":0.5,\"lower\":0},"
        "\"free_class_proportion\":1.0,\"free_class_rate\":1.6}";
    mu_assert(_parse_status(weight_out_of_range) == 2,
              "+F: a weight outside (0,1) should die");

    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_free_class);
    mu_run_test(test_free_class_charts_agree);
    mu_run_test(test_free_class_reduces_to_invariant);
    mu_run_test(test_free_class_mean_quadrature);
    mu_run_test(test_free_class_domain);
    mu_run_test(test_free_class_increment);
    mu_run_test(test_free_class_increment_boundary);
    mu_run_test(test_free_class_increment_matches_rate);
    mu_run_test(test_free_class_increment_domain);
    mu_run_test(test_free_class_rejects);
    mu_run_test(test_freerate_em);
    mu_run_test(test_freerate_em_fixed_rates);
    mu_run_test(test_freerate_em_invariant);
    mu_run_test(test_freerate_em_rejects);
    mu_run_test(test_freerate_em_restores_likelihood);
    mu_run_test(test_weibull_mean_quadrature);
    mu_run_test(test_weibull_mean_quadrature_invariant);
    mu_run_test(test_mean_quadrature_gradient);
    mu_run_test(test_mean_quadrature_rejects);
    mu_run_test(test_gausslaguerre_quadrature);
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
