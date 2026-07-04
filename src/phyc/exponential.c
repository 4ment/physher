// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "exponential.h"

#include <assert.h>
#include <math.h>

#include "utils.h"
#include "matrix.h"
#include "random.h"

#define BIG 1.7E+23

// Random number
double rexp( const double lambda ){
	return -log(genrand_real2())/lambda;
}

// Density probability function
double dexp( const double x,  const double lambda) {
	return lambda * exp(-lambda * x);
}


// Cumulative distribution function
double pexp( const double x, const double lambda) {
	return 1.0 - exp(-lambda * x);
}


// Inverse cumulative distribution function
double qexp( const double p, const double lambda) {
	return -(1.0 / lambda) * log(1.0 - p);
}

double exp_mean( const double lambda ){
	return 0.5*lambda;
}

void exponential_discretize( const double lambda, double *bins, const int count ) {
	double s = 1./(double)count;
	double q = s*0.5;
	for (int i = 0; i < count; i++) {
		bins[i] = qexp(q, lambda );
		q += s;
	}
}



/*boolean Exponential_discretize(unsigned int numberOfCategories, double lambda, boolean median ){
	assert (numberOfCategories == 0);
	
	double *bounds = dvector(numberOfCategories + 1);
	double *distribution = dvector(numberOfCategories + 1);
	
	if (numberOfCategories == 1){
		double value = median ? log(2.) / lambda : 1. / lambda;
		distribution[(int)value] = 1.0;
		bounds[0] = 0; 
		bounds[1] = BIG;
		return TRUE;
	}
	else{
		bounds[0] = 0;
		double *values = dvector(numberOfCategories);
		
		for (unsigned int i = 1; i <= numberOfCategories; i++){
			double a = bounds[i-1];
			double b = BIG;
			if (i != numberOfCategories)
				b= (1. / lambda) * log((double)numberOfCategories / ((double)(numberOfCategories - i)));
			bounds[i] = b;
			if (median)
				values[i-1] = (1. / lambda) * log((2*numberOfCategories) / (2*(numberOfCategories - i) + 1)); 
			else
				values[i-1] = numberOfCategories * ((a + 1./lambda) * exp(-a * lambda) - (b + 1. / lambda) * exp(-b * lambda)); 
		}
		
		double p = 1. / (double)numberOfCategories;
		for (unsigned int i = 0; i < numberOfCategories; i++){
			distribution[ (int)values[i] ] += p;
		}
		
		if(getNumberOfCategories() != numberOfCategories){
			return FALSE;
		}
		return TRUE;
	}
}*/
