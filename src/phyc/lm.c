/*
 *  lm.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/29/10.
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


#include "lm.h"
#include <math.h>
//#include "nmath.h"
#include "utils.h"
#include "student.h"
#include "statistics.h"

void residuals( const double *x, const double *y, int n, double slope, double intcp, double *res ){
	for( int i = 0; i < n; i++ ) {
		res[i] = y[i] - ( intcp + (slope * x[i]) );
	}
}

double residual( double x, double y, double slope, double intcp ){
	return  y - ( intcp + (slope * x) );
}

double sse( double *x, double *y, int n, double slope, double intcp ){
	double SSE = 0;
	for( int i = 0; i < n; i++ ) SSE += pow(residual(x[i],y[i],slope,intcp),2);
	return SSE;
}

double mse( double *x, double *y, int n, double slope, double intcp){
	return sse(x,y,n,slope,intcp)/(n-2);
}

double var( double *x, double *y, int n, double slope, double intcp ){
	return sqrt(mse(x,y,n,slope,intcp));
}


double SE_intercept( double *x, double *y, int n, double slope, double intercept ){
	double sumXmeanX = 0;
	double meanX = mean(x,n);
	for ( int j = 0; j < n; j++ ) sumXmeanX += pow( x[j] - meanX, 2 );
	return var(x, y, n, slope, intercept) * sqrt( ( pow(meanX,2)/sumXmeanX )+(1/n) );
}

double SE_slope( double *x, double *y, int n, double slope, double intercept ){
	double sumXmeanX = 0;
	double meanX = mean(x, n);
	for ( int j = 0; j < n; j++) sumXmeanX += pow( x[j] - meanX, 2 );
	return var(x, y, n, slope, intercept) * sqrt( ( 1/sumXmeanX ) );
}

double CI_slope( double *x, double *y, int n, double slope, double intercept, double prob ){
	//return SE_slope(x,y,n,slope,intercept) * qt(prob, n-2, TRUE, FALSE);
	return SE_slope(x,y,n,slope,intercept) * qt(prob, n-2);
}

double CI_intercept( double *x, double *y, int n, double slope, double intercept, double prob ){
	//return SE_intercept(x,y,n,slope,intercept) * qt(prob, n-2, TRUE, FALSE);
	return SE_intercept(x,y,n,slope,intercept) * qt(prob, n-2);
}

void CI_xIntercept(double *x, double *y, int n, double slope, double intercept, double *lower, double *upper, double alpha ) {
	double meanX = mean(x,n);
	double m11 = 0;
	double ss = 0;
	for( int i = 0; i < n; i++ ){
		m11 += pow(x[i]-meanX,2);
		ss += pow(residual(x[i], y[i], slope, intercept), 2);
	}
	ss /= (n-2);
	double tstar = qt(alpha, n-2);//qt(alpha, n-2, TRUE, FALSE);
	//error("Need to code qt again!");
	double g = (pow(tstar, 2) * ss) / (pow(slope,2) * m11 );
	double xIntercept = - intercept / slope;
	double left = (xIntercept - mean(x,n)) * g;
	double right = (tstar * sqrt(ss)/slope) * sqrt( pow((xIntercept - meanX),2)/m11 + (1 - g)/n );
	*lower = xIntercept + (left + right) / (1 - g);
	*upper = xIntercept + (left - right) / (1 - g);
}

/*
 b = sum((x[i] - meanX) * (y[i] - meanY)) / sum((x[i] - meanX)^2 ) 
 a = meanY - b * meanX
 */
void regression( double *x, double *y, int n, double *slope, double *intercept ){
	double sumX, sumY, sumXX, sumYY, sumXY;
	sumX = sumY = sumXX = sumYY = sumXY = 0;
	int i = 0;
	for( ; i < n; i++ ){
		sumX  += x[i];
		sumY  += y[i];
		sumXX += pow(x[i], 2);
		sumYY += pow(y[i], 2);
		sumXY += x[i] * y[i];
	}
	
	double sumSqDevX   = sumXX - pow(sumX, 2) / n;
	double sumSqDevXY  = sumXY - sumX * sumY / n;
	
	*slope     = sumSqDevXY / sumSqDevX;
	*intercept = (sumY - (*slope) * sumX) / n;
}
