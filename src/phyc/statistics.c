// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <stdio.h>
#include <math.h>

#include "statistics.h"


double correlation( double *x, double *y, int dim ){
	double sumX, sumY, sumXX, sumYY, sumXY;
	sumX = sumY = sumXX = sumYY = sumXY = 0;
	for( int i = 0; i < dim; i++ ){
		sumX  += x[i];
		sumY  += y[i];
		sumXX += pow(x[i], 2);
		sumYY += pow(y[i], 2);
		sumXY += x[i] * y[i];
	}
	double sumSqDevX   = sumXX - pow(sumX, 2) / dim;
	double sumSqDevY   = sumYY - pow(sumY, 2) / dim;
	double sumSqDevXY  = sumXY - sumX * sumY  / dim;
	double denom = sumSqDevX * sumSqDevY;
	return sumSqDevXY / sqrt(denom);
}

float fcorrelation( float *x, float *y, int dim ){
	float sumX, sumY, sumXX, sumYY, sumXY;
	sumX = sumY = sumXX = sumYY = sumXY = 0;
	for( int i = 0; i < dim; i++ ){
		sumX  += x[i];
		sumY  += y[i];
		sumXX += pow(x[i], 2);
		sumYY += pow(y[i], 2);
		sumXY += x[i] * y[i];
	}
	float sumSqDevX   = sumXX - pow(sumX, 2) / dim;
	float sumSqDevY   = sumYY - pow(sumY, 2) / dim;
	float sumSqDevXY  = sumXY - sumX * sumY  / dim;
	float denom = sumSqDevX * sumSqDevY;
	return sumSqDevXY / sqrt(denom);
}

double covariance( const double *x, const double *y, double meanX, double meanY, int dim ){
	double cor = 0;
	
	for ( int i = 0; i < dim; i++ ) {
		cor += (x[i] - meanX)*(y[i] - meanY);
	}
	return cor / (dim-1);
}

double mean( const double *x, int dim ){
	double mean = 0;
	for ( int i = 0; i < dim; i++ ) {
		mean += x[i];
	}
	return mean/dim;
}

double variance( const double *x, int dim, double mean ){
	double var = 0;
	for ( int i = 0; i < dim; i++ ) {
		var += (x[i]-mean) * (x[i]-mean);
	}
	return var/(dim-1); 
}

double standard_deviation( const double *x, int dim, double mean ){
	double stdv = 0;
	for ( int i = 0; i < dim; i++ ) {
		stdv += (x[i]-mean) * (x[i]-mean);
	}
	return sqrt(stdv/(dim-1));
}
