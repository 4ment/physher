// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "lognormal.h"

#include <math.h>

#include "mathconstant.h"
#include "gaussian.h"

#include <stdio.h>


// Probability density function
double dlnorm( const double x,  const double logmu, const double logsigma ){
	return dnorm( log(x), logmu, logsigma) / x;
}

// Cumulative distribution function
double plnorm( const double x, const double logmu, const double logsigma ){
	return pnorm( log(x), logmu, logsigma);
}


// Inverse cumulative distribution function
double qlnorm( const double p, const double logmu, const double logsigma ){
	return exp(qnorm(p, logmu, logsigma));
}

double lognorm_mean( const double logmu, const double logsigma ){
	return exp( logmu + (logsigma*logsigma*0.5) );
}

void lognormal_discretize( const double logmu, const double logsigma, double *bins, const int count ){
	double s = 1./(double)count;
	double z = s*0.5;
	for ( int i = 0; i < count; i++ ) {
		bins[i] = qlnorm(z, logmu, logsigma );
		z += s;
	}
}
