/*
 *  student.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/20/11.
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

#include "student.h"

#include <stdio.h>
#include <math.h>

#include "gaussian.h"

#ifndef M_PI_2
	#define M_PI_2  1.57079632679489661923 /* pi/2 */
#endif

/**
 * Returns the quantile for the two-tailed student's t distribution.
 * This is an implementation of the algorithm in
 * G. W. Hill. "Algorithm 396: Student's t-Quantiles." Communications
 * of the ACM 13(10):619--620.  ACM Press, October, 1970.
 */
double qt(double p, long n) {
	double a, b, c, d, x, y;
	
	if(n < 1) {
		/* you can put your own error handling here */
		fprintf(stderr, "tquantile(%f, %ld): error: second argument must be >= 1 !", p, n);
		return 0.0;
	} else if(p > 1.0 || p <= 0.0) {
		/* you can put your own error handling here */
		fprintf(stderr, "tquantile(%f, %ld): error: first argument must be in (0.0, 1.0] !", p, n);
		return 0.0;
	}
	
	if(n == 1) {
		/* special case */
		p *= M_PI_2;
		return cos(p) / sin(p);
	}
	
	a = 1.0 / (n-0.5);
	b = pow(48.0 / a, 2.0);
	c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
	d = ((94.5 / (b + c) - 3.0) / b + 1.0) * sqrt(a * M_PI_2) * (double)n;
	x = d * p;
	y = pow(x, 2.0/(double)n);
	if(y > 0.05 + a) {
		/* asymptotic inverse expansion about the normal */
		x = qnorm(p*0.5, 0, 1);
		y = x * x;
		if(n < 5) {
			c += 0.3 * ((double)n - 4.5) * (x + 0.6);
			c = (((0.5 * d * x - 0.5) * x - 7.0) * x - 2.0) * x + b + c;
			y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
			y *= a * y;
			if(y > 0.002)
				y = exp(y) - 1.0;
			else
				y += 0.5 * y * y;
		}
	} else
		y = ((1.0/((((double)n + 6.0)/((double)n * y) - 0.089 * d - 0.822) * ((double)n+2.0) * 3.0) + 0.5 / ((double)n+4.0))*y - 1.0) * ((double)n + 1.0) / ((double)n + 2.0) + 1.0 / y;
	
	return sqrt((double)n * y);
}
