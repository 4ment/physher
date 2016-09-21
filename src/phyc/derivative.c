/*
 *  derivative.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 6/20/11.
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

#include "derivative.h"


#include <math.h>

#include "utils.h"
#include "matrix.h"
#include "mathconstant.h"

#define CON 1.4         /* Stepsize is decreased by CON at each iteration. */
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10         /* Sets maximum size of tableau. */
#define SAFE 2.0        /* Return when error is SAFE worse than the best so far. */

/* Returns the derivative of a function func at a point x by Ridders’ method of polynomial extrapolation.
 The value h is input as an estimated initial stepsize; it need not be small, but rather should be
 an increment in x over which func changes substantially. An estimate of the error in the derivative is returned as err. */
double dfridr( Parameters *x, opt_func func, void *data, double h, int index, double *err, bool quick ){
	int i,j; 
	double errt,fac,hh,**a,ans;
	if (h == 0.0){
		error("h must be nonzero in dfridr.");
	}
	
	hh = h;
	//a[0][0] = ((*func)(x+hh)-(*func)(x-hh))/(2.0*hh); 
	
	double xx = Parameters_value( x, index);
	
	*err = BIG;
	
	Parameters_set_value( x, index, xx+hh );
	ans =  func(x, NULL, data);
	
	Parameters_set_value( x, index, xx-hh );
	ans -= (func(x, NULL, data)/(2.0*hh));
	
	if( quick ) return ans;
	
	a = dmatrix(NTAB,NTAB);
	
	a[0][0] = ans;
	
	/* Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation. */
	for ( i = 1;i < NTAB; i++ ){
		hh /= CON; /* Try new, smaller stepsize. */
		//a[0][i] = ((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);	
		xx = Parameters_value( x, index);
		Parameters_set_value( x, index, xx+hh );
		a[0][i] =  func(x, NULL, data);
		
		Parameters_set_value( x, index, xx-hh );
		a[0][i] -= (func(x, NULL, data)/(2.0*hh));
		
		fac = CON2;
		for ( j = 1; j <= i; j++ ) {	/* Compute extrapolations of various orders, requiring no new function evaluations */
			a[j][i] = (a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0); 
			fac = CON2*fac;
			errt = dmax(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			/* The error strategy is to compare each new extrapolation to one order lower, both at the present stepsize and the previous one. */
			if (errt <= *err) {	/* If error is decreased, save the improved answer. */
				*err = errt;
				ans = a[j][i];
			}
		}
		/* If higher order is worse by a significant factor SAFE, then quit early. */
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	
	free_dmatrix(a,NTAB);
	
	return ans;
}

double first_derivative( Parameters *x, opt_func func, void *data ){
	double h = SQRT_EPS*(fabs(Parameters_value( x, 0)) + 1.0);
	
	Parameters_set_value( x, 0, Parameters_value( x, 0) + h);
	double val1 = func(x, NULL, data);
	
	Parameters_set_value( x, 0, Parameters_value( x, 0) - 2*h);
	double val2 = func(x, NULL, data);
	// Centered first derivative
	return (val1 - val2)/(2.0*h);
}
