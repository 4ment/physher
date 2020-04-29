/*
 *  powell.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/10/11.
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

#include "powell.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "optimizer.h"
#include "matrix.h"
#include "linefunction.h"
#include "matrix.h"
#include "mathconstant.h"


#define POWELL_DEBUG 0



/*! Powell.
 * Minimization of a function func of n variables. Input consists of an initial starting point p[0..n-1]; an initial matrix xi[0..n-1][0..n-1], 
 * whose columns contain the initial set of directions (usually the n unit vectors); and ftol, the fractional tolerance in the function value
 * such that failure to decrease by more than this amount on one iteration signals doneness. On output, p is set to the best point found, 
 * xi is the then-current direction set, fret is the returned function value at p, and iter is the number of iterations taken. The routine linmin is used.
 */

opt_result powell_optimize( Parameters *p, opt_func f, void *data, OptStopCriterion stop, double *fmin, opt_update_data uf ){
	
	opt_result status = OPT_KEEP_GOING;
	int i, ibig, j;
	double del, fp, fptt, t;
	double fret = t = 0;
	
	int n = Parameters_count(p);
	
	Parameters *ptt = clone_Parameters(p);
	double *pt  = dvector(n);
	double *xit = dvector(n);
	
	fret = f(p, NULL, data);
	
	// init stop condition
	opt_check_stop( &stop, p, *fmin );
	
	// currently active variables
	bool *active = bvector(n);
	
	
	LineFunction *lf = new_LineFunction( p, f, data );

	LineFunction_force_within_bounds( lf, p );
	
	int numActive = LineFunction_set_active_parameters(lf, p, NULL, active); // should bot be NULL
	
	// if no variables are active return
	if ( numActive == 0 ){
		free(active);
		free_Parameters(ptt);
		free(pt);
		free(xit);
		free_LineFunction(lf);
		*fmin = fret;
		return OPT_NEED_CHECK;
	}
	
	double **xi = deye(n);
	
	for ( j = 0; j < n; j++ ) pt[j] = Parameters_value(p, j); // Save the initial point
	
	while( status == OPT_KEEP_GOING ){
		if(POWELL_DEBUG) fprintf(stdout, "\nPowell iteration #%d (lk=%f)\n", (stop.iter+1), fret);
		fp   = fret;
		ibig = 0;   // Index of the biggest function decrease
		del  = 0.0; // Will be the biggest function decrease.
		for ( i = 0; i < n; i++ ) { // In each iteration, loop over all directions in the set.
			for ( j = 0; j < n; j++ ) xit[j] = xi[j][i];	//Copy the direction,
			fptt = fret; 
			
			LineFunction_update(lf, p, xit);
			fret = LineFunction_minimize(lf);
			
			//linemin( p, f, data, xit, &fret,uf, &stop);
			
			if (fptt-fret > del) {
				del  = fptt - fret;
				ibig = i;
			}
			
		} 
		
		//Termination criterion.
		if ( (status = opt_check_stop( &stop, p, fret )) != OPT_KEEP_GOING ){
			break;
		}
		
		//fprintf(stderr, "Extrapolation %f %f\n",(2.0*(fp-fret)), (ftol*(fabs(fp)+fabs(fret))+TINY) );
		
		//Construct the extrapolated point (ptt) and the average direction moved (xit). Save the old starting point (pt).
		for ( j = 0; j < n; j++ ) {
			//Parameter_warn(ptt->list[j], 2.0 * p->list[j]->value - pt[j]);
			Parameters_set_value( ptt, j, 2.0 * Parameters_value(p, j) - pt[j] );
			xit[j] = Parameters_value(p, j) - pt[j];
			pt[j]  = Parameters_value(p, j);
		}
		
		
		fptt = f(ptt, NULL, data); //Function value at extrapolated point.
		
		//fprintf( stderr, "fptt %f fp %f \n",fptt, fp );
		
		t = 0;
		if ( fptt < fp ) {
			t = 2.0 * (fp-2.0*fret+fptt) * SQR(fp-fret-del) - del * SQR(fp-fptt); 
			// Move to the minimum of the new direction and save the new direction
			if (t < 0.0) {
				//linemin( p, f, data, xit, &fret, uf, &stop );
				fret = LineFunction_minimize(lf);
				
				// from this point on xi is not the initial matrix (n unit vectors)
				for ( j = 0; j < n; j++ ) {
					xi[j][ibig] = xi[j][n-1];
					xi[j][n-1]  = xit[j];
				}
			}
		}
		// reset parameters and bounds without recalculating the function (fret has not changed)
		if ( fptt >= fp || t >= 0 ) uf(data, p);
	}
	
	
	free(active);
	free_Parameters(ptt);
	free(pt);
	free(xit);
	free_LineFunction(lf);
	*fmin = fret;
	return status;
}


void powell_reset_xi( Powell *powell ){
	deye2(powell->xi, powell->n);
}
