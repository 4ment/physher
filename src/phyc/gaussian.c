// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "gaussian.h"

#include <math.h>

#include "mathconstant.h"
#include "errorfuntion.h"
#include "random.h"

double dnorm( const double x, const double mu, const double sigma ) {
	double a = 1.0 / ( SQRT_2PI * sigma );
	double b = pow(x - mu, 2) / ( 2.0 *  sigma * sigma );
	return a * exp(-b);
}

double dnorml( const double x, const double mu, const double sigma ) {
    return -log( SQRT_2PI * sigma ) - pow(x - mu, 2) / ( 2.0 *  sigma * sigma );
}

double pnorm( const double x, const double mu, const double sigma ) {
	double a = (x - mu) / ( sigma * sqrt(2.0) );
	return 0.5 * (1.0 + erf(a));
}


double qnorm( const double p, const double mu, const double sigma ) {
	return mu + sigma * sqrt(2.0) * inverf(2.0 * p - 1.0);
}

// Marsaglia polar method
double rnorm(){
    float x1, x2, w, y1;//, y2;
    
    do {
        x1 = 2.0 * random_double() - 1.0;
        x2 = 2.0 * random_double() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    //y2 = x2 * w;
    return y1;
}
