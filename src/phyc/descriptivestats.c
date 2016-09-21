/*
 *  descriptivestats.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 26/9/12.
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


#include "descriptivestats.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "matrix.h"


double dmean( const double *v, int len ){
	double tot = 0;
	for ( int i = 0; i < len; i++ ) {
		tot += v[i];
	}
	return tot/len;
}

double dweighted_average( const double *v, const double *weights, int len ){
	double tot = 0;
	for ( int i = 0; i < len; i++ ) {
		tot += v[i] * weights[i];
	}
	return tot;
}


double dmedian( const double *v, int len ){

	double *v2 = clone_dvector( v, len );
	qsort(v2, len, sizeof(double), qsort_asc_dvector);
	double median = dmedian_ordered(v2, len);
	free(v2);
	return median;
}

double dmedian_ordered( const double *v, int len ){
	int rank = len/2;
	// odd
	if ( len & 1 ) {
		return v[rank];
	}
	return (v[rank - 1] + v[rank]) / 2.0;
}

double dpercentile( const double *v, int len, double percentile ){
	double *v2 = clone_dvector( v, len );
    qsort(v2, len, sizeof(double), qsort_asc_dvector);
	double val =  dpercentile_ordered(v2, len, percentile);
	free(v2);
	return val;
}

double dpercentile_ordered( const double *v, int len, double percentile ){
//	if( percentile <= 0 || percentile > 100 ){
//        fprintf(stderr, "Percentile should be > 0 and <= 100 (%f)\n", percentile);
//        exit(1);
//    }
	if ( len == 1 ) {
		return v[0];
	}
	int rank = ceil(percentile * len)-1;
    return v[rank];
}

double dpercentile_ordered2( const double *v, int len, double percentile ){
//    if( percentile <= 0 || percentile > 100 ){
//        fprintf(stderr, "Percentile should be > 0 and <= 100 (%f)\n", percentile);
//        exit(1);
//    }
	if ( len == 1 ) {
		return v[0];
	}
    int pos = percentile*(len+1.0)/100.0;
    double val = 0;
    if( pos < 1 ){
        val = v[0];
    }
    else if( pos >= len ){
        val = v[len-1];
    }
    else {
        double lower = v[(int)floor(pos)];
        double upper = v[(int)floor(pos)+1];
        val = lower + (pos-floor(pos)) * (upper-lower);
    }
	return val;
}


#pragma mark -
#pragma mark Integer

int imedian( const int *v, int len ){
	int *v2 = clone_ivector( v, len );
	int median = imedian_ordered(v2,len);
	free(v2);
	return median;
}

int imedian_ordered( const int *v, int len ){
	int rank = len/2;
	return v[rank];
}

int ipercentile( const int *v, int len, double percentile ){
	int *v2 = clone_ivector( v, len );
	qsort(v2, len, sizeof(int), qsort_asc_ivector);
	int val = ipercentile_ordered(v2, len, percentile);
	free(v2);
	return val;
}

int ipercentile_ordered( const int *v, int len, double percentile ){
	if ( percentile == 0 ) {
		return v[0];
	}
	int rank = ceil(percentile * len)-1;
	return v[rank];
}
