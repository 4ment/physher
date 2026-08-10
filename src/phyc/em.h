// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _EM_H_
#define _EM_H_

#include <stdbool.h>
#include <stddef.h>

struct _SingleTreeLikelihood;

// Outcome of one SiteModel_optimize_freerate_EM run.
typedef struct SiteModelEM {
    double logP;   // log-likelihood at the parameters the run left behind
    size_t steps;  // EM steps actually performed
    // Objective evaluations spent in the per-class M-steps -- the cost of the
    // run, near enough: each is a full traversal of the one-category likelihood,
    // against one K-category traversal per step for the E-step.
    size_t evaluations;
    bool converged;  // weights and rates both moved less than `tolerance`
} SiteModelEM;

// EM for a free-rate ("+R") site model, following Wang, Li, Susko & Roger (2008)
//
// The mixture is over the categories of `tlk->sm`, which must be a "discrete"
// site model on the "rate_increments" parameterization (an invariant class is
// allowed; it is category 0 at rate 0 and only its weight is estimated). One
// step is:
//
//   E    gamma_ic = pi_c f_c(D_i) / sum_k pi_k f_k(D_i),  w_ic = n_i gamma_ic
//   M    pi_c  <- (1/N) sum_i w_ic
//        r_c   <- argmax_r sum_i w_ic log f(D_i | r t)      (one Brent per class)
//
// The per-class problem in the second M-step is an ordinary *single*-category
// tree likelihood whose pattern weights are the posterior counts w_ic, so it is
// run on a scratch one-category likelihood over the same tree, substitution
// model and alignment -- the whole point of the method (em_algorithm.md, §5).
//
// The rates come back sorted increasing, which both breaks the label-switching
// symmetry and is what makes the increments positive. Because the M-step fixes
// the branch lengths while the site model re-imposes sum_k p_k r_k = 1 by
// dividing the raw rates through, the tree absorbs the same factor: branch
// lengths are scaled by sum_k p_k r~_k at the end of every step, exactly as
// RateFree::normalizeRates does. That is likelihood-preserving but it does move
// `tree`, and it is why a time tree (whose lengths are heights x clock rate) is
// rejected rather than silently left inconsistent.
//
// `max_steps` of 0 means one step per category; `tolerance`
// of 0 means 1e-4. Returns .logP = NAN, .steps = 0 without touching anything if
// the site model is not one this handles.
SiteModelEM SiteModel_optimize_freerate_EM(struct _SingleTreeLikelihood *tlk,
                                           size_t max_steps, double tolerance);

#endif
