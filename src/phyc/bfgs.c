/*
 *  bfgs.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/11/11.
 *  Copyright (C) 2010 Mathieu Fourment. All rights reserved.
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


#include "bfgs.h"

#include <math.h>
#include <string.h>

#include "utils.h"
#include "matrix.h"
#include "optimizer.h"
#include "linesearch.h"
#include "derivative.h"

#define STPMX 100.0  // Scaled maximum step length allowed in line searches
#define EPS 3.0e-8   // Machine accuracy
#define TOLX (4*EPS) // Convergence criterion on x values

static int check_variables( const Parameters *p, const double *grad, bool *active, double **hessian );


/*
 * Given a starting point p[0..n-1] that is a vector of length n, the Broyden-Fletcher-Goldfarb-Shanno variant of Davidon-Fletcher-Powell minimization
 * is performed on a function func, using its gradient as calculated by a routine dfunc. 
 * The convergence requirement on zeroing the gradient is input as gtol. Returned quantities are p[1..n] (the location of the minimum), 
 * iter (the number of iterations that were performed), and fret (the minimum value of the function).
 * The routine lnsrch is called to perform approximate line minimizations.
 */

	
opt_result dfpmin_optimize( Parameters *p, opt_func f, void *data, OptStopCriterion stop, double *fmin ){
	
	int i,its,j; 
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test; 
	
	opt_result status = OPT_SUCCESS;
	
	int maxeval = stop.iter_max;
	double gtol = stop.tolg;
		
	int n = Parameters_count(p);
	
	double *dg      = dvector(n); 
	double *g       = dvector(n);
	double *hdg     = dvector(n);
	double **hessin = dmatrix(n,n);
	double *xi      = dvector(n);
	bool *active    = bvector(n);
	
	Parameters *pnew = clone_Parameters(p);
	
	// Evaluate at point p and get gradients g
	fp = f(p, g, data);
	
	for (i = 0; i < n; i++ ) {
		memset(hessin[i], 0, n*sizeof(double));
		hessin[i][i] = 1.0;
		
		xi[i] = -g[i]; 
		sum  += Parameters_value(p,i)*Parameters_value(p,i);
		
		// force the parameters to be active
		if( Parameters_value(pnew, i) < Parameters_flower(pnew,i) ){
			Parameters_set_value(pnew, i, Parameters_flower(pnew,i) );
		}
		else if( Parameters_value(pnew, i) > Parameters_fupper(pnew,i) ){
			Parameters_set_value(pnew, i, Parameters_fupper(pnew,i) );
		}
		active[i] = true;
	}
	
	stpmax = STPMX * dmax(sqrt(sum),(double)n);

	//Main loop over the iterations. 
	for ( its = 0; its < maxeval; its++ ) {
		
		status = lnsrch( p, pnew, f, data, fp, g, xi, fmin, stpmax );
		
		// The new function evaluation occurs in lnsrch; save the function value in fp for the next line search. It is usually safe to ignore the value of check.
		fp = *fmin;
		
		// Update the line direction, and the current point
		for ( i = 0; i < n; i++ ) {
			xi[i] = Parameters_value(pnew,i) - Parameters_value(p,i);
			Parameters_set_value(p, i, Parameters_value(pnew, i) );
		}
		
		// Test for convergence on ∆x
		test = 0.0;
		for ( i = 0; i < n; i++ ) {
			temp = fabs(xi[i])/dmax(fabs(Parameters_value(p,i)),1.0);
			if (temp > test) test = temp;
		}

		if (test < TOLX) {
			break;
		}
		
		// Save the old gradient
		memcpy(dg, g, n*sizeof(double));

		// Get the new gradient (we only need the gradient here, not the evaluation of the function at this point)
		f(p, g, data);
		
		//Test for convergence on zero gradient
		test = 0.0; 
		den = dmax(*fmin,1.);
		for ( i = 0; i < n; i++ ) {
			temp = fabs(g[i])*dmax(fabs(Parameters_value(p,i)),1.)/den; 
			if (temp > test) test = temp;
		} 
		if (test < gtol) {
			break;
		}
		
		check_variables(p, g, active, hessin);		
		
		// Compute difference of gradients
		for ( i = 0; i < n; i++ ) dg[i] = g[i] - dg[i];
		
		// Compute difference times current Hessian matrix
		for ( i = 0; i < n; i++ ) {
			hdg[i] = 0.0;
			if ( !active[i] ) continue;
			
			for ( j = 0; j < n; j++ ){
				if ( active[j] ) {
					hdg[i] += hessin[i][j]*dg[j];
				}
			}
		}
		
		// Calculate dot products for the denominators.
		fac = fae = sumdg = sumxi = 0.0; 
		for ( i = 0; i < n; i++ ) {
			if ( active[i] ) {
				fac += dg[i]*xi[i];
				fae += dg[i]*hdg[i];
				sumdg += SQR(dg[i]);
				sumxi += SQR(xi[i]);
			}
		}
		
		if (fac > sqrt(EPS*sumdg*sumxi)) {
			fac = 1.0/fac;
			fad = 1.0/fae;
			
			// The vector that makes BFGS different from DFP:
			for ( i = 0; i < n; i++ ){
				if ( active[i] ) {
					dg[i] = fac*xi[i] - fad*hdg[i];
				}
			}
			
			// The BFGS updating formula:
			for ( i = 0; i < n; i++ ) {
				if ( !active[i] ) continue;
				
				for ( j = i; j < n; j++ ) {
					if ( active[j] ) {
						hessin[i][j] += fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j]; 
						hessin[j][i]  = hessin[i][j];
					}
				}
			}
		}
		
		// Calculate the next direction to go
		for (i = 0; i < n; i++ ) {
			xi[i] = 0.0;
			if ( !active[i] ) continue;
			
			for ( j = 0; j < n; j++ ){
				if ( active[j] ) {
					xi[i] -= hessin[i][j]*g[j];
				}
			}
		}
		
	}
	
	free(xi);
	free_Parameters(pnew);
	free_dmatrix(hessin,n);
	free(hdg);
	free(g);
	free(dg);
	free(active);
	
	if ( its == maxeval ) {
		status = OPT_MAXEVAL;
	}
	
	return status;
}

int check_variables( const Parameters *p, const double *grad, bool *active, double **hessian ){
	int numActive = 0;
	//fprintf(stderr, "----------------------------\n");
	//fprintf(stderr, "LineFunction_check_variables\n");
	int n = Parameters_count(p);
	int j = 0;
	for ( int i = 0; i < n; i++ ){
		active[i] = true;
		if ( active[i] && Parameters_value(p, i) <= Parameters_flower(p,i)+EPS ){
			// no search towards lower boundary
			if ( grad[i] > 0 ){
				active[i] = false;
			}
		}
		else if ( active[i] && Parameters_value(p, i) >= Parameters_fupper(p,i)-EPS ){
			// no search towards upper boundary
			if ( grad[i] < 0 ){
				active[i] = false;
			}
		}
		else if( !active[i] ){
			for ( j = 0; j < n; j++ ) {
				hessian[j][i] = 0.0;
			}
			memset(hessian[i], 0.0, n*sizeof(double));
			hessian[i][i] = 1.0;
			numActive++;
		}
		else {
			numActive++;
		}
		//fprintf(stderr, "%s %d\n", Parameters_name(p, i), active[i] );
	}
	//fprintf(stderr, "\n");
	return numActive;
}
