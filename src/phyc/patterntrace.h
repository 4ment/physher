// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef patterntrace_h
#define patterntrace_h

#include <stddef.h>

#include "model.h"

// Read back the per-pattern log-likelihood block a column logger wrote for
// `model` (a tree likelihood), discarding the first `burnin` rows.
//
// Returns a pattern-major matrix [pattern_count][*sample_count] -- the series
// of every pattern is contiguous, which is what the estimators reduce over --
// to be released with free_dmatrix(trace, pattern_count). `action` only names
// the caller in error messages. Anything malformed is fatal.
double** read_pattern_log_likelihoods(Model* model, const char* filename,
                                      size_t burnin, const char* action,
                                      size_t* sample_count);

#endif /* patterntrace_h */
