/*
 *  linesearch.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/11/11.
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

#include "linesearch.h"

#include <math.h>

#include "parameters.h"
#include "utils.h"
#include "optimizer.h"
#include "matrix.h"

#define ALF 1.0e-4	//Ensures sufficient decrease in function value. 
#define TOLX 1.0e-7	//Convergence criterion on ∆x.
/*
 * Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold and g[1..n], and a direction p[1..n],
 * finds a new point x[1..n] along the direction p from xold where the function func has decreased “sufficiently.”
 * The new function value is returned in f.
 * stpmax is an input quantity that limits the length of the steps so that you do not try to evaluate the function in regions
 * where it is undefined or subject to overflow. p is usually the Newton direction. The output quantity check is false (0) on a normal exit.
 * It is true (1) when x is too close to xold. In a minimization algorithm, this usually signals convergence and can be ignored.
 * However, in a zero-finding algorithm the calling program should check whether the convergence is spurious.
 * Some “difficult” problems may require double precision in this routine. 
 */
// opt_result dfpmin( Parameters *p, opt_func f, void *data, const int maxeval, const double gtol, double *fmin ){

//void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[], float *f, float stpmax, int *check, float (*func)(float [])){

opt_result lnsrch( Parameters *xold, Parameters *x,  opt_func func, void *data, double fold, double *g, double *p, double *fmin, double stpmax ){
	int i;
	double a,alam,alamin,b,disc,rhs1,rhs2,slope,temp,test,tmplam;
	double alam2 = 0;
	double f2 = 0;
	
	opt_result status = OPT_ERROR;	
	
	int n = Parameters_count(xold);
	
	double sum = 0.0; 
	for ( i = 0; i < n; i++ ) sum += p[i]*p[i];
	sum = sqrt(sum);
	
	if (sum > stpmax){
		for ( i = 0; i < n; i++ ) p[i] *= stpmax/sum;  //Scale if attempted step is too big.
	}
	
	slope = 0.0;
	for ( i = 0; i < n; i++ ){
		slope += g[i]*p[i];
	}
	
	// Roundoff problem
	if ( slope >= 0.0 ){
		return status;
	}
	
	test = 0.0;	//Compute λmin.
	for ( i = 0; i < n; i++ ) {					
		temp = fabs(p[i])/dmax(fabs( Parameters_value(xold, i) ),1.0); 
		if (temp > test) test = temp;
	} 
	alamin = TOLX/test;
	alam   = 1.0; //Always try full Newton step first	
	
	for (;;) {
		for ( i = 0; i < n; i++ ){
			Parameters_set_value(x, i, Parameters_value(xold, i) + alam*p[i]); //x_i = xold_i + alam*p_i;
		}
		*fmin = func(x, NULL, data);
		
		// Convergence on ∆x. For zero finding, the calling program should verify the convergence
		if (alam < alamin) { 
			for ( i = 0; i < n; i++ ) Parameters_set_value(x, i, Parameters_value(xold,i) ); // x_i = xold_i;
			status = OPT_NEED_CHECK;
			break;
		}
		// Sufficient function decrease
		else if ( *fmin <= fold+ALF*alam*slope ){
			status = OPT_SUCCESS;
			break;
		}
		// Backtrack
		else {
			if (alam == 1.0){
				tmplam = -slope/(2.0*(*fmin-fold-slope));
			}
			else {
				rhs1 = *fmin-fold-alam*slope; 
				rhs2 = f2-fold-alam2*slope;
				a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2); 
				b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				
				if (a == 0.0){
					tmplam = -slope/(2.0*b);
				}
				else {
					disc = b*b-3.0*a*slope; 
					if (disc < 0.0) tmplam = 0.5*alam;
					else if (b <= 0.0) tmplam = (-b+sqrt(disc))/(3.0*a); 
					else tmplam = -slope/(b+sqrt(disc));
					
				}
				//λ ≤ 0.5λ1
				if (tmplam > 0.5*alam) tmplam = 0.5*alam;
			}
			
		}
		alam2 = alam;
		f2    = *fmin;
		alam  = dmax(tmplam,0.1*alam); //λ ≥ 0.1λ1
	}
	
	return status;
}

