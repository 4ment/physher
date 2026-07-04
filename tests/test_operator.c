// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>
#include <stdlib.h>

#include "minunit.h"
#include "phyc/operator.h"

// Non-static helpers defined in operator.c but not exposed in operator.h.
extern double tune(long accepted, long count, double target, double tuner,
                   double min, double max, bool inverse);
extern double getDeltaP(long count, double logAlpha, double target);
extern double optimizeScaleFactor(double scaleFactor, long count, double logAlpha,
                                  double target);

// Build a bare Operator with only the fields operator_tuning_stats reads.
static Operator make_op(size_t accepted, size_t rejected, size_t delay) {
    Operator op = {0};  // zeroes tuning_started/accepted_at_delay/count_at_delay
    op.accepted_count = accepted;
    op.rejected_count = rejected;
    op.tuning_delay = delay;
    return op;
}

// With no delay, tuning starts immediately and the snapshot baseline is the
// state at the first call.
char* test_tuning_stats_no_delay() {
    Operator op = make_op(0, 0, 0);
    long count = -1, accepted = -1;

    mu_assert(operator_tuning_stats(&op, &count, &accepted),
              "no delay: tuning should start immediately");
    mu_assert(op.tuning_started, "no delay: tuning_started must be set");
    mu_assert(count == 0 && accepted == 0,
              "no delay: first call must report zero post-delay proposals");

    // Advance the chain: 3 accepts, 2 rejects since the (zero) baseline.
    op.accepted_count = 3;
    op.rejected_count = 2;
    mu_assert(operator_tuning_stats(&op, &count, &accepted),
              "no delay: tuning should remain active");
    mu_assert(count == 5, "no delay: post-delay count must be total proposals");
    mu_assert(accepted == 3, "no delay: post-delay accepted must be total accepts");
    return NULL;
}

// During the delay window the helper reports "not tuning" and never snapshots.
char* test_tuning_stats_within_delay() {
    Operator op = make_op(2, 3, 10);  // total 5 < delay 10
    long count = -1, accepted = -1;

    mu_assert(!operator_tuning_stats(&op, &count, &accepted),
              "within delay: should report not tuning");
    mu_assert(!op.tuning_started, "within delay: must not snapshot yet");

    // Still inside the window at the boundary minus one.
    op.accepted_count = 4;
    op.rejected_count = 5;  // total 9 < 10
    mu_assert(!operator_tuning_stats(&op, &count, &accepted),
              "within delay: total 9 < delay 10 still not tuning");
    mu_assert(!op.tuning_started, "within delay: still no snapshot");
    return NULL;
}

// The bias fix: once tuning starts, accepted/count must measure ONLY the
// post-delay window, not the cumulative totals that include the warm-up.
char* test_tuning_stats_unbiased_after_delay() {
    Operator op = make_op(5, 5, 10);  // total 10 == delay -> tuning begins now
    long count = -1, accepted = -1;

    mu_assert(operator_tuning_stats(&op, &count, &accepted),
              "after delay: tuning should start at the boundary");
    mu_assert(op.tuning_started, "after delay: snapshot must be taken");
    mu_assert(op.accepted_at_delay == 5 && op.count_at_delay == 10,
              "after delay: baseline snapshot must capture current counts");
    mu_assert(count == 0 && accepted == 0,
              "after delay: first post-delay window is empty");

    // Run 10 more proposals: 3 accepts, 7 rejects.
    op.accepted_count = 8;   // 5 warm-up + 3 post-delay
    op.rejected_count = 12;  // 5 warm-up + 7 post-delay
    mu_assert(operator_tuning_stats(&op, &count, &accepted),
              "after delay: tuning remains active");
    mu_assert(count == 10, "after delay: post-delay count excludes warm-up");
    mu_assert(accepted == 3,
              "after delay: post-delay accepted excludes warm-up accepts");

    // The biased (old) ratio would have been 8/20 = 0.4; the correct
    // post-delay ratio is 3/10 = 0.3. They must differ.
    double biased = 8.0 / 20.0;
    double unbiased = (double)accepted / (count + 1);
    mu_assert(fabs(biased - 0.4) < 1e-12, "sanity: biased ratio is 0.4");
    mu_assert(fabs(unbiased - 3.0 / 11.0) < 1e-12,
              "after delay: unbiased ratio uses post-delay window");
    return NULL;
}

// The baseline snapshot is taken exactly once and is not moved by later calls.
char* test_tuning_stats_snapshot_is_sticky() {
    Operator op = make_op(7, 3, 5);  // total 10 >= delay 5
    long count = -1, accepted = -1;

    mu_assert(operator_tuning_stats(&op, &count, &accepted), "should be tuning");
    size_t snap_acc = op.accepted_at_delay;
    size_t snap_cnt = op.count_at_delay;
    mu_assert(snap_acc == 7 && snap_cnt == 10, "snapshot captured at first call");

    op.accepted_count = 100;
    op.rejected_count = 100;
    mu_assert(operator_tuning_stats(&op, &count, &accepted), "still tuning");
    mu_assert(op.accepted_at_delay == snap_acc && op.count_at_delay == snap_cnt,
              "snapshot must not move after it is first taken");
    mu_assert(count == 190 && accepted == 93,
              "post-delay stats are measured from the sticky baseline");
    return NULL;
}

// tune() steps the tuner up when acceptance exceeds target (and down below it),
// with the direction reversed for inverse-scaled parameters.
char* test_tune_direction() {
    // Scaler-like (inverse=false): high acceptance -> larger scale factor.
    double up = tune(80, 100, 0.24, 1.0, 0.0001, 20.0, false);
    mu_assert(up > 1.0, "tune: high acceptance should increase tuner");

    double down = tune(10, 100, 0.24, 1.0, 0.0001, 20.0, false);
    mu_assert(down < 1.0, "tune: low acceptance should decrease tuner");

    // Out-of-range result falls back to the old tuner unchanged.
    double clamped = tune(100, 100, 0.24, 19.999, 0.0001, 20.0, false);
    mu_assert(clamped == 19.999, "tune: out-of-range step keeps old tuner");

    // Inverse (concentration-like): high acceptance -> smaller concentration.
    double inv = tune(80, 100, 0.24, 1.0, 0.01, 10000.0, true);
    mu_assert(inv < 1.0, "tune: inverse high acceptance should decrease tuner");
    return NULL;
}

// getDeltaP and optimizeScaleFactor are the Robbins-Monro pieces shared by the
// scaler/up-down/exchange operators.
char* test_robbins_monro() {
    // At exactly the target acceptance the delta is zero.
    double at_target = getDeltaP(100, log(0.24), 0.24);
    mu_assert(fabs(at_target) < 1e-12, "getDeltaP: zero step at target");

    // Above target -> positive step, below -> negative.
    double above = getDeltaP(100, 0.0, 0.24);  // exp(min(0,0))=1
    mu_assert(fabs(above - (1.0 - 0.24) / 100.0) < 1e-12,
              "getDeltaP: step above target");
    double below = getDeltaP(100, log(0.1), 0.24);
    mu_assert(fabs(below - (0.1 - 0.24) / 100.0) < 1e-12,
              "getDeltaP: step below target");

    // A scale factor of 0.5 at target acceptance is a fixed point.
    double sf = optimizeScaleFactor(0.5, 100, log(0.24), 0.24);
    mu_assert(fabs(sf - 0.5) < 1e-12,
              "optimizeScaleFactor: 0.5 is fixed at target acceptance");
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_tuning_stats_no_delay);
    mu_run_test(test_tuning_stats_within_delay);
    mu_run_test(test_tuning_stats_unbiased_after_delay);
    mu_run_test(test_tuning_stats_snapshot_is_sticky);
    mu_run_test(test_tune_direction);
    mu_run_test(test_robbins_monro);
    return NULL;
}

RUN_TESTS(all_tests);
