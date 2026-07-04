// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "chisq.h"

#include "gamma.h"

// Cumulative distribution function
// return P(X <= x)
double pchisq( const double x, const int df ){
	if (x < 0.0 || df < 0.0) return 1.0;
	return gammp( ((double)df)/2, x/2);
}

// Inverse cumulative distribution function
double qchisq( const double p, const int df ){
	return 2.*invgammp(p,0.5*df);
}
