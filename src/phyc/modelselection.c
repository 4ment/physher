/*
 *  ic.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/22/11.
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

#include "modelselection.h"

#include <math.h>

#include "utils.h"
#include "chisq.h"

// descending order
static int cmpFunc( const void *a, const void *b ){
	return ( *(int*)b - *(int*)a );
}


// AIC = -2lnL + 2K
double AIC( const double lk, const int k ){
	return -2 * lk + 2 * k;
}

// AICc = -2lnL + 2Kn/(n-K-1)
// AICc = AIC + 2K(K+1)/(n-K-1)
// lk: max likelihood, k: number of free parameters, n: sample size
double AICc( const double lk, const int k, const int n ){
	return -2 * lk + (2*k*n)/(n-k-1);
}

// BIC = -2lnL + Kln(n)
double BIC( const double lk, const int k, const int n ){
	return -2 * lk + k*log(n);
}

// Information criterion of Hannan and Quinn.
// HQ = −2 lnL + 2Kln( ln(n) )
double HQ( const double lk, const int k, const int n ){
	return -2 * lk + 2*k*log(log(n));
}


void model_weight( double *IC, const int n, double *ICw, double *ICcw ) {
	double minv = dmin_vector( IC, n );
	
	double sum = .0;
	int i = 0;
	for ( i = 0; i < n; i++) {
		ICw[i] = exp( -0.5 * (IC[i] - minv) );
		sum += ICw[i];
	}
	
	// normalize over all the models
	for ( i = 0; i < n; i++) {
		ICw[i] /= sum;
	}
	
	qsort(IC, n, sizeof(double), cmpFunc );
	
	
	double cumWeight = 0.0;
	for ( i = 0; i < n; i++) {
		cumWeight += ICw[i];
		ICcw[i] = cumWeight;
	}
		
}

double LRT( const double lk0, const double lk1, const unsigned n0, const unsigned n1  ){
	double delta = 2 * ( lk1 - lk0 );
	if (delta <= 0 ) return 1;
	return 1 - pchisq(delta, (n1 - n0) );
	
}

