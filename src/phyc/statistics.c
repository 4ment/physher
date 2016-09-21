/*
 *  statistics.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 15/8/12.
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

double covariance( const double *x, const double *y, int dim ){
	double cor = 0;
	double meanX = 0;
	double meanY = 0;
	for( int i = 0; i < dim; i++ ){
		meanX  += x[i];
		meanY  += y[i];
	}
	meanX /= dim;
	meanY /= dim;
	
	for ( int i = 0; i < dim; i++ ) {
		cor += (x[i] - meanX)*(y[i] - meanY);
	}
	return cor / dim;
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
