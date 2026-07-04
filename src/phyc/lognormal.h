// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _LOGNORMAL_H_
#define _LOGNORMAL_H_

double dlnorm( const double x,  const double mu, const double sigma );

// Cumulative distribution function
double plnorm( const double x, const double mu, const double sigma );


// Inverse cumulative distribution function
double qlnorm( const double p, const double logmu, const double logsigma );

double lognorm_mean( const double logmu, const double logsigma );

void lognormal_discretize( const double mu, const double sigma, double *bins, const int count );

#endif
