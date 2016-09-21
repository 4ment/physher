/*
 *  lognormal.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 13/8/11.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

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
