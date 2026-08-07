// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

// Tests for the optimizer stopping criteria, in particular the one the meta
// optimizer applies to a whole sweep of its schedule.
//
// The objective is a synthetic separable quadratic rather than a tree
// likelihood: the meta criterion is pure bookkeeping over the values its
// children report, so a closed-form objective pins it down exactly and without
// the numerical slop of a real likelihood.

#include <math.h>
#include <stdlib.h>

#include "minunit.h"
#include "phyc/optimizer.h"
#include "phyc/parameters.h"

// f(x) = offset + sum_i (x_i - target_i)^2, minimized at x = target with value
// `offset`. Separable, so a schedule of one Brent per coordinate reaches the
// minimum in a single sweep. `evaluations` counts calls so a test can assert on
// how much work a criterion caused.
typedef struct {
    Parameters* x;
    double* target;
    double offset;
    // Cross term 2*coupling*(x0-t0)*(x1-t1), for n == 2 only. Positive definite
    // while |coupling| < 1; the closer to 1, the slower coordinate descent
    // crawls along the resulting ridge, which is how a test forces a run to
    // exhaust its iteration budget instead of converging.
    double coupling;
    size_t evaluations;
} Quadratic;

static double quadratic_f(Parameters* params, double* grad, void* data) {
    Quadratic* q = (Quadratic*)data;
    q->evaluations++;
    double sum = q->offset;
    for (size_t i = 0; i < Parameters_count(q->x); i++) {
        double d = Parameters_value(q->x, i) - q->target[i];
        sum += d * d;
    }
    if (q->coupling != 0.0 && Parameters_count(q->x) == 2) {
        sum += 2.0 * q->coupling * (Parameters_value(q->x, 0) - q->target[0]) *
               (Parameters_value(q->x, 1) - q->target[1]);
    }
    return sum;
}

// The same objective with the sign flipped. A child driving this instead of
// quadratic_f is a minimizer working against the meta optimizer's target: it
// makes the sweep end worse than it started, and reports its own (falling)
// value while doing so.
static double negated_quadratic_f(Parameters* params, double* grad, void* data) {
    return -quadratic_f(params, grad, data);
}

static Quadratic* new_Quadratic(size_t n, const double* start, const double* target,
                                double offset) {
    Quadratic* q = (Quadratic*)malloc(sizeof(Quadratic));
    q->x = new_Parameters(n);
    q->target = (double*)malloc(sizeof(double) * n);
    q->offset = offset;
    q->coupling = 0.0;
    q->evaluations = 0;
    for (size_t i = 0; i < n; i++) {
        Parameters_move(q->x, new_Parameter("x", start[i],
                                            new_Constraint(-INFINITY, INFINITY)));
        q->target[i] = target[i];
    }
    return q;
}

static void free_Quadratic(Quadratic* q) {
    free_Parameters(q->x);
    free(q->target);
    free(q);
}

// One Brent child per coordinate, wired to the shared objective.
static Optimizer* new_meta_over_quadratic(Quadratic* q, double precision,
                                          size_t iterations) {
    Optimizer* meta = new_Optimizer(OPT_META);
    opt_set_objective_function(meta, quadratic_f);
    opt_set_data(meta, q);
    opt_set_tolfx(meta, precision);
    opt_set_max_iteration(meta, iterations);
    opt_set_verbosity(meta, 0);

    for (size_t i = 0; i < Parameters_count(q->x); i++) {
        Optimizer* child = new_Optimizer(OPT_BRENT);
        opt_set_objective_function(child, quadratic_f);
        opt_set_data(child, q);
        opt_set_tolx(child, 1.e-8);
        opt_set_max_iteration(child, 200);
        opt_set_verbosity(child, 0);
        Parameters* ps = new_Parameters(1);
        Parameters_add(ps, Parameters_at(q->x, i));
        opt_set_parameters(child, ps);
        free_Parameters(ps);
        opt_add_optimizer(meta, child);
    }
    return meta;
}

// A schedule of coordinate-wise Brent searches must land on the analytic
// minimum, and the value handed back must be the objective at that point.
char* test_meta_reaches_minimum() {
    double start[2] = {-3.0, 5.0};
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, start, target, 7.0);
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-8, 100);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result == OPT_SUCCESS, "meta: did not report convergence");
    mu_assert(fabs(Parameters_value(q->x, 0) - target[0]) < 1.e-4,
              "meta: x0 not at the minimum");
    mu_assert(fabs(Parameters_value(q->x, 1) - target[1]) < 1.e-4,
              "meta: x1 not at the minimum");
    mu_assert(fabs(fmin - q->offset) < 1.e-6, "meta: fmin is not the minimum value");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// Starting at the minimum, the first sweep gains nothing and the run stops
// immediately -- but it must still report the objective it stopped at.
char* test_meta_reports_fmin_when_already_converged() {
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, target, target, 7.0);
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-8, 100);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result == OPT_SUCCESS, "meta: did not report convergence at the optimum");
    mu_assert(!isnan(fmin), "meta: fmin left untouched");
    mu_assert(fabs(fmin - q->offset) < 1.e-6, "meta: fmin is not the minimum value");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// Running out of iterations must still hand back the objective at the point
// the run stopped at. The coupled quadratic keeps every sweep improving by more
// than the tolerance, so the loop can only end on the iteration cap.
char* test_meta_reports_fmin_on_max_iterations() {
    double start[2] = {-3.0, 5.0};
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, start, target, 7.0);
    q->coupling = 0.95;
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-12, 3);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result == OPT_MAXITER, "meta: expected the iteration cap to be hit");
    mu_assert(!isnan(fmin), "meta: fmin left untouched on the OPT_MAXITER path");
    double expected = quadratic_f(NULL, NULL, q);
    mu_assert(fabs(fmin - expected) < 1.e-9,
              "meta: fmin is not the objective at the point it stopped at");
    mu_assert(opt_iterations(meta) == 3, "meta: wrong number of sweeps reported");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// `iter_min` (the JSON `min` key) is a floor on the number of sweeps. Starting
// at the optimum, every sweep gains nothing, so only the floor can keep the run
// going.
char* test_meta_honours_min_iterations() {
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, target, target, 7.0);
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-8, 100);
    opt_set_min_iteration(meta, 4);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result == OPT_SUCCESS, "meta: did not report convergence");
    mu_assert(opt_iterations(meta) >= 4, "meta: stopped before iter_min sweeps");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// `patience` requires several consecutive flat sweeps, not just one.
char* test_meta_honours_patience() {
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, target, target, 7.0);
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-8, 100);
    opt_set_patience(meta, 3);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result == OPT_SUCCESS, "meta: did not report convergence");
    mu_assert(opt_iterations(meta) == 3, "meta: did not wait for patience sweeps");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// A sweep that leaves the objective higher than it started is not convergence.
// The child here also reports a *falling* value of its own objective, so this
// only holds because the criterion re-evaluates the target rather than trusting
// what the last child in the schedule handed back.
char* test_meta_does_not_succeed_on_a_worse_point() {
    double target[1] = {1.0};
    double start[1] = {1.0};
    Quadratic* q = new_Quadratic(1, start, target, 0.0);
    // Bounded, or the child would run the coordinate off to infinity.
    Parameter_set_bounds(Parameters_at(q->x, 0), -5.0, 5.0);

    Optimizer* meta = new_Optimizer(OPT_META);
    opt_set_objective_function(meta, quadratic_f);
    opt_set_data(meta, q);
    opt_set_tolfx(meta, 1.e-8);
    opt_set_max_iteration(meta, 5);
    opt_set_verbosity(meta, 0);

    Optimizer* child = new_Optimizer(OPT_BRENT);
    opt_set_objective_function(child, negated_quadratic_f);
    opt_set_data(child, q);
    opt_set_tolx(child, 1.e-8);
    opt_set_verbosity(child, 0);
    Parameters* ps = new_Parameters(1);
    Parameters_add(ps, Parameters_at(q->x, 0));
    opt_set_parameters(child, ps);
    free_Parameters(ps);
    opt_add_optimizer(meta, child);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result != OPT_SUCCESS, "meta: reported success at a worse point");
    mu_assert(fmin > 1.0, "meta: fmin does not reflect the objective it ended at");
    mu_assert(fabs(fmin - quadratic_f(NULL, NULL, q)) < 1.e-9,
              "meta: fmin was taken from the child rather than the target");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// The schedule is configuration, not scratch space: an Optimizer is reused
// across bootstrap replicates, so a run that rewrites `rounds` gives replicate 0
// a different schedule from every replicate after it.
char* test_meta_does_not_mutate_the_schedule() {
    double start[2] = {-3.0, 5.0};
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, start, target, 7.0);
    // Coupled, so the run takes more than the handful of sweeps a separable
    // objective needs -- the schedule was rewritten on the third sweep.
    q->coupling = 0.95;
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-10, 50);
    OptimizerSchedule* schedule = opt_get_schedule(meta);
    schedule->rounds[0] = 3;
    schedule->rounds[1] = 2;

    double fmin = NAN;
    opt_optimize(meta, &fmin);

    mu_assert(schedule->rounds[0] == 3, "meta: rounds[0] was overwritten by the run");
    mu_assert(schedule->rounds[1] == 2, "meta: rounds[1] was overwritten by the run");

    // A second run from a fresh starting point must behave identically.
    Parameters_set_value(q->x, 0, start[0]);
    Parameters_set_value(q->x, 1, start[1]);
    double fmin2 = NAN;
    opt_optimize(meta, &fmin2);
    mu_assert(fabs(fmin - fmin2) < 1.e-9, "meta: consecutive runs disagree");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// A child that fails without writing fmin -- scaler_optimize returns OPT_FAIL
// without touching it when it has no treelikelihood -- must neither abandon the
// sweep nor leave the criterion comparing a value that entry never produced.
char* test_meta_survives_a_failed_child() {
    double start[1] = {-3.0};
    double target[1] = {1.0};
    Quadratic* q = new_Quadratic(1, start, target, 7.0);

    Optimizer* meta = new_Optimizer(OPT_META);
    opt_set_objective_function(meta, quadratic_f);
    opt_set_data(meta, q);
    opt_set_tolfx(meta, 1.e-8);
    opt_set_max_iteration(meta, 10);
    opt_set_verbosity(meta, 0);

    // The entry that does the actual work must still get there.
    Optimizer* child = new_Optimizer(OPT_BRENT);
    opt_set_objective_function(child, quadratic_f);
    opt_set_data(child, q);
    opt_set_tolx(child, 1.e-8);
    opt_set_verbosity(child, 0);
    Parameters* ps = new_Parameters(1);
    Parameters_add(ps, Parameters_at(q->x, 0));
    opt_set_parameters(child, ps);
    free_Parameters(ps);
    opt_add_optimizer(meta, child);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result == OPT_SUCCESS, "meta: a failed entry aborted the schedule");
    mu_assert(fabs(Parameters_value(q->x, 0) - target[0]) < 1.e-4,
              "meta: the working entry did not reach the minimum");
    mu_assert(fabs(fmin - q->offset) < 1.e-6,
              "meta: fmin does not match the target after a failed child");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// The criterion in isolation, where the sweep values can be dictated exactly.
char* test_check_progress() {
    OptStopCriterion stop = {0};
    stop.tolfx = 1.e-5;
    stop.patience = 1;
    stop.iter_min = 0;
    stop.iter = 1;

    stop.stall = 0;
    mu_assert(opt_check_progress(&stop, 100.0, 90.0) == OPT_PROGRESS_ONGOING,
              "criterion: a real gain must not stop the run");

    stop.stall = 0;
    mu_assert(opt_check_progress(&stop, 100.0, 100.0) == OPT_PROGRESS_CONVERGED,
              "criterion: a flat sweep must stop the run");

    stop.stall = 0;
    mu_assert(opt_check_progress(&stop, 100.0, 110.0) == OPT_PROGRESS_WORSE,
              "criterion: a regression must not be reported as convergence");

    // The tolerance is absolute: a gain of 0.5 log units is a real gain whether
    // the objective is 1 or 1e6. Scaling it by the magnitude would make a run on
    // a large alignment accept a gain of tens of units as convergence.
    stop.stall = 0;
    mu_assert(opt_check_progress(&stop, 1.5, 1.0) == OPT_PROGRESS_ONGOING,
              "criterion: a gain of 0.5 at magnitude 1 is not convergence");
    stop.stall = 0;
    mu_assert(opt_check_progress(&stop, 1000000.5, 1000000.0) == OPT_PROGRESS_ONGOING,
              "criterion: a gain of 0.5 at magnitude 1e6 is not convergence either");

    // ...but a tolerance below the representable resolution of the objective
    // still terminates, rather than asking for a difference doubles cannot hold.
    stop.stall = 0;
    stop.tolfx = 0.0;
    mu_assert(opt_check_progress(&stop, 1.0e12, 1.0e12 - 1.0e-6) == OPT_PROGRESS_CONVERGED,
              "criterion: a sub-ulp gain must not keep the run alive");
    stop.tolfx = 1.e-5;

    // iter_min is a floor even when the sweeps are flat.
    stop.iter_min = 5;
    stop.stall = 0;
    stop.iter = 2;
    mu_assert(opt_check_progress(&stop, 100.0, 100.0) == OPT_PROGRESS_ONGOING,
              "criterion: stopped before iter_min");
    stop.iter = 5;
    stop.stall = 0;
    mu_assert(opt_check_progress(&stop, 100.0, 100.0) == OPT_PROGRESS_CONVERGED,
              "criterion: did not stop once iter_min was reached");

    // patience counts *consecutive* flat sweeps; a real gain resets it.
    stop.iter_min = 0;
    stop.patience = 3;
    stop.stall = 0;
    mu_assert(opt_check_progress(&stop, 100.0, 100.0) == OPT_PROGRESS_ONGOING,
              "criterion: stopped on the first flat sweep despite patience 3");
    mu_assert(opt_check_progress(&stop, 100.0, 100.0) == OPT_PROGRESS_ONGOING,
              "criterion: stopped on the second flat sweep despite patience 3");
    mu_assert(opt_check_progress(&stop, 100.0, 90.0) == OPT_PROGRESS_ONGOING,
              "criterion: a real gain must not stop the run");
    mu_assert(opt_check_progress(&stop, 100.0, 100.0) == OPT_PROGRESS_ONGOING,
              "criterion: a real gain must reset the patience counter");

    return NULL;
}

// The budget check in isolation, where the clock and the counters can be set to
// whatever the case needs.
char* test_check_limits() {
    OptStopCriterion stop = {0};

    // Everything zero means no limit.
    mu_assert(opt_check_limits(&stop) == OPT_KEEP_GOING,
              "limits: an unbudgeted run must not be stopped");

    // A budget that has not been reached yet.
    stop.time_max = 3600.0;
    time(&stop.time_start);
    stop.iter_max = 100;
    stop.iter = 5;
    stop.f_eval_max = 1000;
    stop.f_eval_current = 5;
    mu_assert(opt_check_limits(&stop) == OPT_KEEP_GOING,
              "limits: stopped inside every budget");

    stop.iter = 101;
    mu_assert(opt_check_limits(&stop) == OPT_MAXITER, "limits: iteration budget");
    stop.iter = 5;

    stop.f_eval_current = 1000;
    mu_assert(opt_check_limits(&stop) == OPT_MAXEVAL, "limits: evaluation budget");

    // Time is the most urgent, so it must win over a budget that is also spent
    // -- the version this replaced evaluated all three in sequence and let the
    // last one checked erase an already-detected timeout.
    stop.time_start -= 7200;
    mu_assert(opt_check_limits(&stop) == OPT_MAXTIME,
              "limits: time budget must not be masked by another");
    stop.f_eval_current = 5;
    mu_assert(opt_check_limits(&stop) == OPT_MAXTIME, "limits: time budget");

    return NULL;
}

// An evaluation budget stops the schedule and is reported as OPT_MAXEVAL. The
// coupled quadratic keeps improving, so nothing else can end the run.
char* test_meta_honours_max_evaluations() {
    double start[2] = {-3.0, 5.0};
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, start, target, 7.0);
    q->coupling = 0.95;
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-14, 10000);
    opt_set_max_evaluation(meta, 200);

    double fmin = NAN;
    opt_result result = opt_optimize(meta, &fmin);

    mu_assert(result == OPT_MAXEVAL, "meta: evaluation budget not enforced");
    mu_assert(!isnan(fmin), "meta: fmin left untouched on the OPT_MAXEVAL path");
    // Enforced between sweeps, so the budget is honoured to within one sweep;
    // what must not happen is running the full 10000 sweeps regardless.
    mu_assert(opt_iterations(meta) < 10000, "meta: ran past the evaluation budget");
    // Meta's own tally must be in the same ballpark as the objective's, i.e. the
    // entries of the schedule are actually reporting their work.
    mu_assert(opt_f_evaluations(meta) >= 200, "meta: undercounted its evaluations");
    mu_assert((size_t)opt_f_evaluations(meta) <= q->evaluations,
              "meta: counted more evaluations than the objective saw");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// A wall-clock budget stops the schedule and is reported as OPT_MAXTIME. The
// patience is set beyond reach so only the clock can end the run.
char* test_meta_honours_time_limit() {
    double start[2] = {-3.0, 5.0};
    double target[2] = {1.0, 2.0};
    Quadratic* q = new_Quadratic(2, start, target, 7.0);
    q->coupling = 0.95;
    Optimizer* meta = new_meta_over_quadratic(q, 1.e-14, 100000000);
    opt_set_patience(meta, 100000000);
    opt_set_time_max(meta, 1);

    double fmin = NAN;
    time_t before = time(NULL);
    opt_result result = opt_optimize(meta, &fmin);
    double elapsed = difftime(time(NULL), before);

    mu_assert(result == OPT_MAXTIME, "meta: time budget not enforced");
    mu_assert(!isnan(fmin), "meta: fmin left untouched on the OPT_MAXTIME path");
    mu_assert(elapsed < 30.0, "meta: overshot the time budget badly");

    free_Optimizer(meta);
    free_Quadratic(q);
    return NULL;
}

// opt_check_stop measures each check against the last one that made progress,
// not against the previous check. A run creeping downhill at just under the
// tolerance per step is still descending, and must not be called converged
// merely because no single step cleared the bar.
char* test_check_stop_holds_its_reference() {
    Parameters* x = new_Parameters(1);
    Parameters_move(x, new_Parameter("x", 0.0, new_Constraint(-INFINITY, INFINITY)));

    OptStopCriterion stop = {0};
    stop.tolfx = 1.0;
    stop.patience = 2;
    stop.iter_max = 1000;

    // Seeding call: records the starting value, never stops.
    mu_assert(opt_check_stop(&stop, x, 100.0) == OPT_KEEP_GOING,
              "check_stop: stopped on the seeding call");

    // A steady descent of 0.6 per step against a tolerance of 1.0. No single
    // step clears the bar, so measured against the previous check every one of
    // them is flat and patience 2 would be exhausted by the second. Measured
    // against the last point of real progress, every other step accumulates
    // past the bar and the stall counter keeps resetting.
    double fx = 100.0;
    for (int i = 0; i < 8; i++) {
        fx -= 0.6;
        mu_assert(opt_check_stop(&stop, x, fx) == OPT_KEEP_GOING,
                  "check_stop: called a steady descent converged");
        mu_assert(stop.stall < stop.patience,
                  "check_stop: a descending run exhausted its patience");
    }

    // Genuinely flat now: two checks at the same value exhaust the patience.
    mu_assert(opt_check_stop(&stop, x, fx) == OPT_KEEP_GOING,
              "check_stop: stopped on the first flat step despite patience 2");
    mu_assert(opt_check_stop(&stop, x, fx) == OPT_SUCCESS,
              "check_stop: did not converge on a genuinely flat objective");

    free_Parameters(x);
    return NULL;
}

// A budget still stops the run even while the objective keeps improving.
char* test_check_stop_reports_budgets() {
    Parameters* x = new_Parameters(1);
    Parameters_move(x, new_Parameter("x", 0.0, new_Constraint(-INFINITY, INFINITY)));

    OptStopCriterion stop = {0};
    stop.tolfx = 1.e-8;
    stop.patience = 1;
    stop.iter_max = 10;

    mu_assert(opt_check_stop(&stop, x, 100.0) == OPT_KEEP_GOING,
              "check_stop: stopped on the seeding call");
    stop.iter = 5;
    mu_assert(opt_check_stop(&stop, x, 90.0) == OPT_KEEP_GOING,
              "check_stop: stopped inside the iteration budget");
    stop.iter = 11;
    mu_assert(opt_check_stop(&stop, x, 80.0) == OPT_MAXITER,
              "check_stop: iteration budget not reported");

    free_Parameters(x);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_check_limits);
    mu_run_test(test_check_stop_holds_its_reference);
    mu_run_test(test_check_stop_reports_budgets);
    mu_run_test(test_check_progress);
    mu_run_test(test_meta_reaches_minimum);
    mu_run_test(test_meta_reports_fmin_when_already_converged);
    mu_run_test(test_meta_reports_fmin_on_max_iterations);
    mu_run_test(test_meta_honours_min_iterations);
    mu_run_test(test_meta_honours_patience);
    mu_run_test(test_meta_does_not_succeed_on_a_worse_point);
    mu_run_test(test_meta_does_not_mutate_the_schedule);
    mu_run_test(test_meta_survives_a_failed_child);
    mu_run_test(test_meta_honours_max_evaluations);
    mu_run_test(test_meta_honours_time_limit);
    return NULL;
}

RUN_TESTS(all_tests);
