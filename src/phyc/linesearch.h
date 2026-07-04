// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _LINESEARCH_H_
#define _LINESEARCH_H_

#include "parameters.h"
#include "optimizer.h"

typedef struct {
    double alpha;
    double fx;
    int nfev;
    int status;
} LineSearchResult;

LineSearchResult strong_wolfe_line_search(Parameters* parameters, opt_func fun,
                                          opt_grad_func grad_f, void* data, const double *x,
                                          const double *p, const double *g, double f0,
                                          double c1, double c2, double alpha0, double amax);

opt_result lnsrch(Parameters *parameters, double* x,  opt_func func, void *data, double fold, double *g, double *p, double *fmin, double stpmax, double alam);

#endif
