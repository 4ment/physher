// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _EXPONENTIAL_H_
#define _EXPONENTIAL_H_

double rexp( const double lambda );

double dexp( const double x,  const double lambda);

// Cumulative distribution function
double pexp( const double x, const double lambda);


// Inverse cumulative distribution function
double qexp( const double p, const double lambda);

double exp_mean( const double lambda );

void exponential_discretize( const double lambda, double *bins, const int count );


#endif
