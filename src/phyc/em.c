// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

// EM for finite mixtures over sites.
//
// This lives apart from sitemodel.c because it is the one piece of free-rate
// machinery that needs the *likelihood*: the M-step is a tree traversal per
// class per Brent evaluation. sitemodel.c stays a pure parameters-to-rates map.

#include "em.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "optimizer.h"
#include "parameters.h"
#include "sitemodel.h"
#include "sitepattern.h"
#include "tree.h"
#include "treelikelihood.h"

#define SITEMODEL_EM_MIN_PROP 1.0e-4
#define SITEMODEL_EM_MIN_RATE 1.0e-4
// Slack on the monotonicity check. EM cannot decrease the likelihood, but the
// inner M-steps are Brent runs stopped at a tolerance, so this is generalized
// EM and small backward moves are expected; anything larger is a bug worth
// reporting.
#define SITEMODEL_EM_MONOTONICITY_SLACK 0.1

// Tolerance of the per-class Brent search, deliberately not the caller's EM
// tolerance. The M-step does not need the argmax: monotonicity survives any r_c
// that improves the weighted objective (this is generalized EM), and the next
// step re-solves the same problem from wherever this one stopped. Every extra
// digit costs a handful of *full tree traversals* per class per step, to refine
// a number the outer iteration was going to move anyway. Brent's test is
// relative (tol*|x| + ZEPS), so this is a fraction of the rate.
#define SITEMODEL_EM_RATE_TOL 1.0e-3
// How far either side of its current value a class rate is searched. Left to
// itself find_bracket starts from delta = 0.1*(upper - lower) with upper =
// 1/pi_c, which for a small class is a window tens of times wider than the
// answer -- and every golden step widening or shrinking it is another traversal.
// After the first step or two the rates barely move, so a window this tight is
// almost always already a bracket. It does not cap the step: find_bracket
// expands past it when the objective still falls, and a step that ends pinned to
// the window edge simply moves again on the next one.
#define SITEMODEL_EM_RATE_WINDOW 2.0

// The single-class problem the M-step hands to Brent (em_algorithm.md §5): the
// mixture's own likelihood, with a one-category site model whose only rate is
// `rate` swapped in for `mixture`, and the mixture's pattern weights overwritten
// by the posterior counts w_ic of the class being optimized.
//
// The auxiliary likelihood is the mixture's rather than a second
// SingleTreeLikelihood of its own so that the two cannot disagree. Everything
// the mixture was configured with -- rescaling, root frequencies, whether SSE
// kernels are used, thread count -- is shared by construction instead of being
// copied field by field into a scratch object that a later field would silently
// be left out of.
typedef struct {
    SingleTreeLikelihood* tlk;
    Parameter* rate;
    size_t evaluations;
    // The one-category site model, owned here and installed only for the
    // M-step; and everything the swap has to put back.
    SiteModel* single;
    SiteModel* mixture;
    double* weights;  // the alignment's pattern counts n_i
    int cat_count;
    bool use_upper;
    int node_id;
    size_t recompute_count;
    bool swapped;
} FreeRateMStep;

static double _em_rate_objective(Parameters* x, double* grad, void* data) {
    FreeRateMStep* m = (FreeRateMStep*)data;
    m->evaluations++;
    // The M-step rate is a bare Parameter with no listeners, so the value Brent
    // just wrote marked nothing dirty: every transition matrix has to be redone
    // by hand.
    SingleTreeLikelihood_update_all_nodes(m->tlk);
    return -m->tlk->calculate(m->tlk);
}

// Build the one-category uniform site model the M-step evaluates through, and
// the buffer its swap needs. Done once, outside the EM loop; the swap itself is
// _em_enter/_em_leave.
static void _em_mstep_init(FreeRateMStep* m, SingleTreeLikelihood* tlk,
                           Parameter* rate) {
    Parameters* empty = new_Parameters(1);
    m->single = new_SiteModel_with_parameters(empty, NULL, 1, DISTRIBUTION_UNIFORM,
                                              false, QUADRATURE_QUANTILE_MEDIAN,
                                              RATE_PARAMETERIZATION_AUTO);
    free_Parameters(empty);
    // _get_rate of a uniform model is exactly mu, and its single proportion is 1,
    // so integrating the root partials against it is the plain one-class density.
    SiteModel_set_mu(m->single, rate);

    m->tlk = tlk;
    m->rate = rate;
    m->evaluations = 0;
    m->mixture = tlk->sm;
    m->weights = clone_dvector(tlk->sp->weights, tlk->sp->count);
    m->swapped = false;
}

static void _em_mstep_free(FreeRateMStep* m) {
    free_SiteModel(m->single);
    free(m->weights);
}

// Install the one-category model. Bracketed around the M-step alone rather than
// the whole EM run, so the E-step and _em_refresh always see the mixture and
// which site model tlk carries is answerable line by line.
static void _em_enter(FreeRateMStep* m) {
    SingleTreeLikelihood* tlk = m->tlk;
    m->cat_count = tlk->cat_count;
    m->use_upper = tlk->use_upper;
    m->node_id = tlk->node_id;
    m->recompute_count = tlk->recompute_count;

    // The upper fast path returns a likelihood computed from partials the
    // one-category pass is about to invalidate, and node_id restricts the
    // traversal to a subtree. Both have to be off, exactly as they were on the
    // freshly constructed likelihood this used to allocate. Clearing use_upper
    // also disables the tripod path in _calculate.
    tlk->use_upper = false;
    tlk->node_id = -1;
    // Suppress the partials/matrices ping-pong flip for the whole M-step. A
    // proposal that has already stored leaves partials[1] allocated for good,
    // and _calculate_partials flips into it whenever recompute_count == 1 --
    // which would hand the M-step's throwaway one-category partials to a later
    // restore. Holding the counter above the threshold pins every evaluation to
    // the current slot, and the slot is made good again by the _em_refresh that
    // follows the write-back. Same device as the rescaling re-pass in
    // _calculate_simple.
    tlk->recompute_count = 2;

    SingleTreeLikelihood_set_sitemodel(tlk, m->single);
    m->swapped = true;
}

static void _em_leave(FreeRateMStep* m) {
    if (!m->swapped) return;
    SingleTreeLikelihood* tlk = m->tlk;
    // The pattern weights are the alignment's again before anything can read a
    // weighted sum off them.
    memcpy(tlk->sp->weights, m->weights, sizeof(double) * tlk->sp->count);
    SingleTreeLikelihood_set_sitemodel(tlk, m->mixture);
    // set_sitemodel derives this and should land on the same number; restoring
    // what was actually saved is what makes "left exactly as found" a promise
    // rather than a coincidence.
    tlk->cat_count = m->cat_count;
    tlk->use_upper = m->use_upper;
    tlk->node_id = m->node_id;
    tlk->recompute_count = m->recompute_count;
    m->swapped = false;
}

// Make the weighted objective of class c an ordinary one by injecting the
// posterior counts w_ic as the pattern weights. No node is dirtied: per-pattern
// likelihoods do not depend on the weights, only their sum does, and the rate
// change that follows dirties everything anyway.
static void _em_set_class_weights(FreeRateMStep* m, const double* w) {
    memcpy(m->tlk->sp->weights, w, sizeof(double) * m->tlk->sp->count);
    SingleTreeLikelihood_update_weights(m->tlk);
}

// Fill `out[c*P + i]` with log f_c(D_i), the per-pattern likelihood of category
// c *without* its weight pi_c, and leave tlk->lk holding the mixture
// log-likelihood at the current parameters.
//
// Integrating the root partials against a one-hot weight vector picks out one
// category, so this reuses the very kernels the mixture uses -- rescaling
// included -- rather than reaching into the [category][pattern][state] layout.
// The upper-partials fast path can return without refreshing the root's lower
// partials, hence the forced full pass.
//
// tlk->pattern_lk is the buffer the kernels are allocated against, so it is what
// they write into; the mixture values are put back at the end, leaving tlk
// exactly as a plain calculate() would have.
static void _em_refresh(SingleTreeLikelihood* tlk, double* out) {
    const size_t cat_count = tlk->sm->cat_count;
    const int pattern_count = tlk->sp->count;

    tlk->use_upper = false;
    tlk->sm->need_update = true;
    SingleTreeLikelihood_update_all_nodes(tlk);
    tlk->calculate(tlk);

    int rootID = Node_id(Tree_root(tlk->tree));
    const double* root_partials =
        tlk->partials[tlk->current_partials_indexes[rootID]][rootID];
    const double* freqs = tlk->get_root_frequencies(tlk);

    double* onehot = dvector(cat_count);
    for (size_t c = 0; c < cat_count; c++) {
        memset(onehot, 0, sizeof(double) * cat_count);
        onehot[c] = 1.0;
        tlk->integrate_partials(tlk, root_partials, onehot, tlk->root_partials);
        tlk->node_log_likelihoods(tlk, tlk->root_partials, freqs, tlk->pattern_lk);
        memcpy(out + c * pattern_count, tlk->pattern_lk,
               sizeof(double) * pattern_count);
    }
    free(onehot);

    tlk->integrate_partials(tlk, root_partials, tlk->sm->get_proportions(tlk->sm),
                            tlk->root_partials);
    tlk->node_log_likelihoods(tlk, tlk->root_partials, freqs, tlk->pattern_lk);
}

// Sort `rates[0..count)` increasing, permuting `props` alongside. The mixture is
// invariant under relabelling, so this costs nothing and buys two things: the
// increments r_k - r_{k-1} come out positive, and consecutive EM steps keep
// talking about the same categories. Insertion sort -- `count` is the number of
// rate categories, never more than a couple of dozen.
static void _em_sort_rates(double* rates, double* props, size_t count) {
    for (size_t i = 1; i < count; i++) {
        double rate = rates[i];
        double prop = props[i];
        size_t j = i;
        while (j > 0 && rates[j - 1] > rate) {
            rates[j] = rates[j - 1];
            props[j] = props[j - 1];
            j--;
        }
        rates[j] = rate;
        props[j] = prop;
    }
}

// Reject, with a message naming what is wrong, every site model this cannot
// handle. Getting it wrong is not a numerical near-miss but a silent
// optimization of the wrong parameters, so the checks are exhaustive.
static bool _em_freerate_supported(const SingleTreeLikelihood* tlk) {
    const SiteModel* sm = (tlk == NULL ? NULL : tlk->sm);
    if (sm == NULL) {
        fprintf(stderr, "sitemodel EM: no site model\n");
        return false;
    }
    if (sm->distribution != DISTRIBUTION_DISCRETE) {
        fprintf(stderr,
                "sitemodel EM: only \"distribution\": \"discrete\" (+R) is a "
                "mixture EM can split\n");
        return false;
    }
    if (sm->rate_parameterization != RATE_PARAMETERIZATION_RATE_INCREMENTS) {
        fprintf(stderr,
                "sitemodel EM: implemented for \"rate_increments\" only, not "
                "\"%s\"\n",
                SiteModel_rate_parameterization_name(sm->rate_parameterization));
        return false;
    }
    if (Parameters_count(sm->rates) != 1 || sm->proportions == NULL) {
        fprintf(stderr,
                "sitemodel EM: expected one increment vector and a "
                "proportions simplex\n");
        return false;
    }
    if (sm->cat_count < 2) {
        fprintf(stderr, "sitemodel EM: a %u-category model is not a mixture\n",
                sm->cat_count);
        return false;
    }
    // The weight floor is paid for out of the largest weight, which is at least
    // 1/K; that only works while K categories' worth of deficit still leaves it
    // above the floor itself.
    if ((sm->cat_count + 1.0) * sm->cat_count * SITEMODEL_EM_MIN_PROP >= 1.0) {
        fprintf(stderr,
                "sitemodel EM: %u categories cannot all clear the weight floor "
                "of %g\n",
                sm->cat_count, SITEMODEL_EM_MIN_PROP);
        return false;
    }
    // The M-step maximizes over r_c with the branch lengths held fixed, and the
    // site model then divides the raw rates by sum_k p_k r~_k; the tree has to
    // absorb that factor for the two to describe the same model. On a time tree
    // the lengths are heights x clock rate and there is no single distance
    // vector to scale.
    if (Tree_is_time_mode(tlk->tree)) {
        fprintf(stderr,
                "sitemodel EM: the mean-rate constraint is re-imposed by "
                "scaling branch lengths, which a time tree does not own\n");
        return false;
    }
    return true;
}

SiteModelEM SiteModel_optimize_freerate_EM(SingleTreeLikelihood* tlk, size_t max_steps,
                                           double tolerance) {
    SiteModelEM result = {NAN, 0, 0, false};
    if (!_em_freerate_supported(tlk)) return result;

    SiteModel* sm = tlk->sm;
    Parameter* increments = Parameters_at(sm->rates, 0);
    const bool free_props = Parameter_estimate(sm->proportions);
    const bool free_rates = Parameter_estimate(increments);
    if (!free_props && !free_rates) {
        fprintf(stderr, "sitemodel EM: both the weights and the rates are fixed\n");
        return result;
    }

    const size_t cat_count = sm->cat_count;
    // With an invariant class, category 0 is pinned at rate 0 and sits outside
    // the running sum of the increments, exactly as in
    // _calculate_rates_discrete_increments. Only its *weight* is estimated.
    const size_t first = (sm->invariant ? 1 : 0);
    const int pattern_count = tlk->sp->count;
    // N = sum_i n_i. Summed rather than read off sp->nsites, which is an int and
    // is not refreshed when the weights are overwritten in place (bootstrap).
    double nsites = 0;
    for (int i = 0; i < pattern_count; i++) nsites += tlk->sp->weights[i];
    if (max_steps == 0) max_steps = cat_count;
    if (tolerance <= 0.0) tolerance = 1.0e-4;
    // Increments of exactly 0 (two categories that converged onto the same rate)
    // would sit on, or below, the parameter's own lower bound.
    const double increment_floor = fmax(Parameter_lower(increments),
                                        SITEMODEL_EM_MIN_RATE * SITEMODEL_EM_MIN_RATE);

    // The auxiliary single-class machinery, built once outside the loop.
    Parameter* rate = new_Parameter("sitemodel.em.rate", 1.0,
                                    new_Constraint(SITEMODEL_EM_MIN_RATE, INFINITY));
    FreeRateMStep mstep;
    _em_mstep_init(&mstep, tlk, rate);
    Optimizer* opt = new_Optimizer(OPT_BRENT);
    opt_set_objective_function(opt, _em_rate_objective);
    opt_set_data(opt, &mstep);
    opt_set_tolx(opt, SITEMODEL_EM_RATE_TOL);
    opt_set_verbosity(opt, 0);

    double* logf = dvector(cat_count * pattern_count);     // log f_c(D_i)
    double* weights = dvector(cat_count * pattern_count);  // w_ic = n_i gamma_ic
    double* props = dvector(cat_count);
    double* rates = dvector(cat_count);
    double* old_props = dvector(cat_count);
    double* old_rates = dvector(cat_count);
    double* theta = dvector(cat_count - first);
    double* work = dvector(cat_count);

    _em_refresh(tlk, logf);
    double logP = tlk->lk;

    size_t step = 0;
    for (; step < max_steps; step++) {
        memcpy(old_props, sm->get_proportions(sm), sizeof(double) * cat_count);
        for (size_t c = 0; c < cat_count; c++) old_rates[c] = sm->get_rate(sm, c);
        // Where each Brent search starts.
        memcpy(rates, old_rates, sizeof(double) * cat_count);

        // --- E-step: gamma_ic = pi_c f_c(D_i) / L_i, then w_ic = n_i gamma_ic.
        // Done in log space and normalized by the largest term: log f_c is a
        // sum over thousands of sites and underflows long before the ratio does.
        memset(props, 0, sizeof(double) * cat_count);
        bool degenerate = false;
        for (int i = 0; i < pattern_count; i++) {
            double max = -INFINITY;
            for (size_t c = 0; c < cat_count; c++) {
                work[c] = log(old_props[c]) + logf[c * pattern_count + i];
                if (work[c] > max) max = work[c];
            }
            if (!isfinite(max)) {
                degenerate = true;
                break;
            }
            double sum = 0;
            for (size_t c = 0; c < cat_count; c++) {
                work[c] = exp(work[c] - max);
                sum += work[c];
            }
            const double n = tlk->sp->weights[i];
            for (size_t c = 0; c < cat_count; c++) {
                double w = n * work[c] / sum;
                weights[c * pattern_count + i] = w;
                props[c] += w;
            }
        }
        if (degenerate) {
            fprintf(stderr,
                    "sitemodel EM: every category has zero likelihood at a "
                    "pattern; stopping after %zu step(s)\n",
                    step);
            break;
        }

        // --- M-step part 1: the weights have the closed form pi_c = (1/N) sum_i w_ic.
        // sum_c gamma_ic is 1 by construction, so the simplex constraint holds
        // without renormalizing.
        bool floored = false;
        if (free_props) {
            for (size_t c = 0; c < cat_count; c++) props[c] /= nsites;
            size_t largest = 0;
            for (size_t c = 1; c < cat_count; c++) {
                if (props[c] > props[largest]) largest = c;
            }
            for (size_t c = 0; c < cat_count; c++) {
                if (props[c] < SITEMODEL_EM_MIN_PROP) {
                    // Take the deficit out of the largest weight so the simplex
                    // still sums to 1.
                    props[largest] -= SITEMODEL_EM_MIN_PROP - props[c];
                    props[c] = SITEMODEL_EM_MIN_PROP;
                    floored = true;
                }
            }
        } else {
            memcpy(props, old_props, sizeof(double) * cat_count);
        }

        // --- M-step part 2: one independent 1-D problem per class,
        // r_c = argmax_r sum_i w_ic log f(D_i | r t). Injecting w_ic as the
        // pattern weights of the single-class likelihood turns the weighted
        // objective into an ordinary one, so plain Brent solves it.
        if (free_rates) {
            _em_enter(&mstep);
            for (size_t c = first; c < cat_count; c++) {
                _em_set_class_weights(&mstep, weights + c * pattern_count);
                // sum_c pi_c r_c = 1 with every rate positive, so no single rate
                // can exceed 1/pi_c; within that, search a window around where the
                // class currently sits.
                double ceiling = 1.0 / props[c];
                double lower =
                    fmax(SITEMODEL_EM_MIN_RATE, rates[c] / SITEMODEL_EM_RATE_WINDOW);
                double upper = fmin(ceiling, rates[c] * SITEMODEL_EM_RATE_WINDOW);
                // A rate already sitting on (or past) the ceiling leaves no window
                // to centre; fall back to the full interval.
                if (!(lower < upper)) {
                    lower = SITEMODEL_EM_MIN_RATE;
                    upper = ceiling;
                }
                Parameter_set_bounds(rate, lower, upper);
                Parameter_set_value(rate, fmin(fmax(rates[c], lower), upper));
                double fmax_c;
                opt_maximize_univariate(opt, rate, &fmax_c);
                rates[c] = Parameter_value(rate);
            }
            _em_leave(&mstep);
            _em_sort_rates(rates + first, props + first, cat_count - first);
        }

        // --- Write the step back. `rates` now holds the absolute rates the
        // M-step settled on, against the branch lengths as they stand; the
        // increments are their gaps, and the site model turns those back into
        // rates by running sum and division by norm = sum_k p_k r~_k. That
        // division is no part of the M-step objective, so the tree absorbs the
        // same factor and every r_c t stays exactly where it was.
        //
        // The weights need it just as much as the rates do: with the increments
        // held fixed, moving p still moves the normalized rates -- by exactly
        // norm = sum_k p_k^new r_k^old -- and the weight EM is only valid while
        // the class densities f_c stay put.
        double norm = 0;
        if (free_rates) {
            double running = 0;
            for (size_t c = first; c < cat_count; c++) {
                double increment = (c == first ? rates[c] : rates[c] - rates[c - 1]);
                theta[c - first] = fmax(increment, increment_floor);
                running += theta[c - first];
                rates[c] = running;  // the flooring may have moved it
                norm += props[c] * running;
            }
        } else {
            for (size_t c = first; c < cat_count; c++) norm += props[c] * rates[c];
        }
        if (free_props) Parameter_set_values(sm->proportions, props);
        if (free_rates) Parameter_set_values(increments, theta);
        Tree_scale_distance(tlk->tree, norm);

        _em_refresh(tlk, logf);
        const double new_logP = tlk->lk;
        if (new_logP < logP - SITEMODEL_EM_MONOTONICITY_SLACK) {
            fprintf(stderr,
                    "sitemodel EM: step %zu decreased the log-likelihood, "
                    "%f -> %f\n",
                    step + 1, logP, new_logP);
        }
        logP = new_logP;
        result.steps = step + 1;

        // Read the step back off the model rather than off the arrays above: a
        // weight update alone moves the normalized rates too, and a reordering
        // in the M-step permuted `props`/`rates` but not their `old_`
        // counterparts (that direction is safe -- it can only overstate the
        // step, never declare a false convergence).
        double delta = 0;
        for (size_t c = 0; c < cat_count; c++) {
            delta = fmax(delta, fabs(sm->get_proportion(sm, c) - old_props[c]));
            delta = fmax(delta, fabs(sm->get_rate(sm, c) - old_rates[c]));
        }
        if (delta < tolerance) {
            result.converged = true;
            break;
        }
        // A weight on the floor means the model has more categories than the
        // data supports; further steps only push the others around.
        if (floored) break;
    }
    result.logP = logP;
    result.evaluations = mstep.evaluations;

    free(logf);
    free(weights);
    free(props);
    free(rates);
    free(old_props);
    free(old_rates);
    free(theta);
    free(work);
    free_Optimizer(opt);
    _em_mstep_free(&mstep);
    free_Parameter(rate);
    return result;
}
