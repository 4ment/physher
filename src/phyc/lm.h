// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _LINEAR_MODEL_H_
#define _LINEAR_MODEL_H_


void regression( double *x, double *y, int n, double *slope, double *intercept );

void residuals( const double *x, const double *y, int n, double s, double intcp, double *res );
double residual( double x, double y, double s, double intcp );
double sse( double *x, double *y, int n, double s, double intcp );
double mse( double *x, double *y, int n, double s, double intcp);
double var( double *x, double *y, int n, double s, double intcp );

void CI_xIntercept(double *x, double *y, int n, double slope, double intercept, double *lower, double *upper, double alpha );
double CI_slope( double *x, double *y, int n, double slope, double intercept, double prob );
double CI_intercept( double *x, double *y, int n, double slope, double intercept, double prob );

#endif
