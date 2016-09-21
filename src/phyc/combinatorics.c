/*
 *  combinatorics.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 23/8/12.
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
#include <stdlib.h>

#include "matrix.h"
#include "utils.h"

double bell_number( unsigned int n ){
	if ( n == 0 || n == 1 ) {
		return 1;
	}
	double *prev = NULL;
	int prev_len = 0;
	int cur_len = 0;
	for ( int r = 0; r < n; r++) {
		cur_len = r+1;
		double *cur = dvector(cur_len);
		cur[0] = ( prev_len > 0 ? prev[r - 1] : 1);
		for ( int i = 1; i <= r; i++) {
			cur[i] = prev[i - 1] + cur[i - 1];
		}
		if ( prev != NULL ) {
			free(prev);
			prev = NULL;
		}
		prev = cur;
		prev_len = cur_len;
	}
	double b = prev[prev_len - 1];
	free(prev);
	return b;
}

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int I4_CHOOSE, the number of combinations of N
//    things taken K at a time.

int choose ( unsigned int n, unsigned int k ) {
	int i;
	int mn;
	int mx;
	int value;
	
	mn = imin ( k, n - k );
	
	if ( mn < 0 ) {
		value = 0;
	}
	else if ( mn == 0 ) {
		value = 1;
	}
	else {
		mx = imax ( k, n - k );
		value = mx + 1;
		
		for ( i = 2; i <= mn; i++ ) {
			value = ( value * ( mx + i ) ) / i;
		}
	}
	return value;
}

// return largest value v where v < a and  Choose(v,b) <= x
static long _LargestV(long a, long b, long x) {
    long v = a - 1;
    
    while (choose(v,b) > x){
        --v;
    }
    
    return v;
}

// http://msdn.microsoft.com/en-us/library/aa289166%28v=vs.71%29.aspx
// return the mth lexicographic element of combination C(n,k)
// ans is of length k
void combination_at_index( unsigned int n, unsigned int k, long m, unsigned *ans ) {
    
    long a = n;
    long b = k;
    long x = (choose(n, k) - 1) - m; // x is the "dual" of m
    
    for ( int i = 0; i < k; i++ ) {
        ans[i] = _LargestV(a,b,x); // largest value v, where v < a and vCb < x
        x = x - choose(ans[i], b);
        a = ans[i];
        b--;
    }
    
    for ( long i = 0; i < k; ++i ) {
        ans[i] = (n-1) - ans[i];
    }
    
    //return new Combination(this.n, this.k, ans);
}

