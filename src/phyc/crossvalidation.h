// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef crossvalidation_h
#define crossvalidation_h

#include <stdint.h>

#include "optimizer.h"
#include "parameters.h"
#include "sitepattern.h"

// K-fold cross-validation over alignment sites, run from the "physher" section
// after a point estimate has been obtained: the parameter values in place when
// run() is called are the estimate every fold is restarted from.
//
// What it estimates is the expected out-of-sample log predictive density,
//
//   elpd = sum_f sum_i n_i^(f) log p(D_i | theta^(-f)),
//
// where fold f is fitted on the complement of its own sites and scores only its
// own. That is a model comparison criterion that needs no parameter count, which
// is the whole reason it is here: a CAT fit has no honest k -- its rates are
// cluster centres and its weights are site counts, so neither was maximized
// over, and its components are not distinct -- so AIC cannot rank it against a discrete
// Gamma, and the in-sample likelihood cannot either because the assignment was chosen
// using the data it is then scored on. Cross-validation charges for exactly that.
//
// Three decisions carry the correctness of the whole thing:
//
// 1. **Folds partition columns, not patterns.** A pattern of multiplicity n_i has
//    its n_i columns dealt across the folds, so the pattern set never changes and
//    only the multiplicities do. Holding out whole patterns would ask a different
//    and much harsher question -- prediction of column types the fit never saw --
//    and would make the folds as unbalanced as the multiplicities are.
//
// 2. **The held-out weights are never installed.** They are only ever a reduction
//    vector applied to per-pattern log densities computed while the *training*
//    weights are in place. Installing them would silently change the model being
//    scored: an empirical CAT site model normalises its category rates to mean one
//    over sp->weights and sp->nsites (_cat_update), so a swap to the held-out
//    weights rescales every rate, and the branch lengths the training fold chose
//    would be read at another rate scale.
//
// 3. **A CAT model is scored by its marginalized mixture, never by its score.**
//    A held-out site has no category label, and taking one off its own column is
//    the leakage this exists to prevent -- note that cat_assign hands even a
//    pattern of zero training weight an arg-max label read from its own profile,
//    so the conditional score leaks whether or not the pattern was in the
//    training fold. The mixture density depends only on the category rates and
//    their weights, which the training sites fixed. Everything else is scored by
//    its ordinary per-pattern likelihood, which is already marginal over its rate
//    categories, so the two are on one scale and directly comparable.
typedef struct CrossValidation {
    // Distinct SitePatterns to split. Indexed by pattern set.
    SitePattern** patterns;  // borrowed from the tree likelihoods
    double** weights;        // owned, the observed multiplicities
    double** train;          // owned, the training multiplicities of one fold
    double** test;           // owned, folds x count held-out multiplicities
    int* observed_nsites;    // owned, sp->nsites as found
    size_t pattern_count;

    // Tree likelihoods reading those pattern sets, and the pattern set each one
    // reads (two partitions may share one alignment).
    Model** treelikelihoods;
    size_t* pattern_index;
    double** logpred;  // owned, per-pattern log predictive density, per likelihood
    double* elpd;      // owned, accumulated per likelihood
    size_t treelikelihood_count;

    // The training schedule, run once per fold on the training weights.
    Optimizer** optimizers;
    size_t optimizer_count;

    Parameters* parameters;  // restored to `estimate` before each fold
    double* estimate;        // owned
    size_t estimate_size;

    size_t folds;
    char* filename;  // owned, per-fold per-pattern output, NULL to skip
    int verbosity;
    // The split runs off its own generator, not the global one, so that two
    // models fitted in two processes hold out the same sites given the same
    // "seed". A paired comparison of their per-site densities is only as good as
    // that: with different splits the pairing survives -- every site is still
    // held out exactly once by each -- but the difference picks up fold-to-fold
    // noise neither model is responsible for.
    uint64_t rng_state;

    void (*run)(struct CrossValidation*);
    void (*free)(struct CrossValidation*);
} CrossValidation;

CrossValidation* new_CrossValidation_from_json(json_node* node, Hashtable* hash);

#endif /* crossvalidation_h */
