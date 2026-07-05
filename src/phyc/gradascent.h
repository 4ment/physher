// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef gradascent_h
#define gradascent_h

#include <stdio.h>
#include "optimizer.h"
#include "tracelogger.h"

opt_result optimize_stochastic_gradient(Parameters* parameters, opt_func f, opt_grad_func grad_f, double eta, void *data, OptStopCriterion *stop, int verbose, double *fmin, OptimizerCheckpoint* checkpointer);
opt_result optimize_stochastic_gradient_adam(bool maximize, Parameters* parameters, opt_func f, opt_grad_func grad_f, double eta, void *data, OptStopCriterion *stop, int verbose, double *fmin, OptimizerCheckpoint* checkpointer, Trace* logger);


opt_result optimize_stochastic_gradient_adapt(Parameters* parameters, opt_func f, opt_grad_func grad_f, void(*reset)(void*),
											  double* etas, size_t eta_count, void *data,
											  OptStopCriterion *stop, int verbose, double *best_eta, size_t nthreads);

#endif /* gradascent_h */
