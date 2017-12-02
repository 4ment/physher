/*
 *  brent.c
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

#include "brent.h"

#include <assert.h>
#include <math.h>

#include "utils.h"
#include "mathconstant.h"

#define ZEPS 1.0e-10


/*
 * 
 */

opt_result serial_brent_optimize( Parameters *ps, opt_func f, void *data, OptStopCriterion *stop, double *fmin ){
	Parameters* temp = new_Parameters(1);
	//for(int j = 0; j < stop->iter_min; j++)
	for(int i = 0; i < Parameters_count(ps); i++){
		Parameters_add(temp, Parameters_at(ps, i));
		stop->iter = 0;
		stop->f_eval_current = 0;
		opt_result status = brent_optimize(temp, f, data, stop, fmin);
		Parameters_pop(temp);
	}
	free_Parameters(temp);
	return OPT_SUCCESS;
}

//opt_result brent_optimize( Parameters *ps, opt_func f, void *data, const int maxeval, const double tol, double *fmin ){
opt_result brent_optimize( Parameters *ps, opt_func f, void *data, OptStopCriterion *stop, double *fmin ){

	double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e = 0.0;
	double d = 0;
	
	double tol = stop->tolx;
    //int *numFun = &stop.f_eval_current;
    size_t *iter = &stop->iter;

	a = ( Parameters_flower(ps,0) < Parameters_fupper(ps,0) ? Parameters_flower(ps,0) : Parameters_fupper(ps,0));
	b = ( Parameters_flower(ps,0) > Parameters_fupper(ps,0) ? Parameters_flower(ps,0) : Parameters_fupper(ps,0));
	x = Parameters_value(ps,0);
    
//    // x is equal to its lower bound
//    if ( x == a) {
//        //ax = a;
//        Parameters_set_value(newps,0, x*2);
//        b = 10.0*a;
//	}
//    // x is close to its lower bound
//    else if (x <= 2.0 * a) {
//        //ax = a;
//        //bx = Parameters_value(newps,0);
//        b = 5.0*x;
//	}
//    else {
//        a = 0.5*x;
//        //bx = xguess;
//        //cx = 2.0*Parameters_value(newps,0);
//        b = x*2;
//	}
//
//	/* ideally this range includes the true minimum, i.e.,
//     fb < fa and fb < fc
//     if not, we gradually expand the boundaries until it does,
//     or we near the boundary of the allowed range and use that
//     */
//    Parameters_set_value(newps, 0, a );
//	double fa = f(newps, NULL, data);
//    Parameters_set_value(newps, 0, x );
//	fx = f(newps, NULL, data);
//    Parameters_set_value(newps, 0, b );
//	double fb = f(newps, NULL, data);
//    Parameters_set_value(newps, 0, x );
//    
//    //printf("a %e (%f) b %e (%f) c %e (%f)\n", a, fa, x, fx, b, fb);
//    
//	while(fa < fx && a > xmin) {
//        a = (a+xmin)/2.0;
//        if (a < 2.0*xmin)	/* give up on shrinking the region */
//            a = xmin;
//        Parameters_set_value(newps, 0, a );
//        fa = f(newps, NULL, data);
//        //fa = (*f)(ax,data);
//	}
//	while(fb < fx && b < xmax) {
////        b = (b+xmax)/2.0;
////        if (b > xmax * 0.95)
////            b = xmax;
//        b *= 2.;
//        Parameters_set_value(newps, 0, b );
//        fb = f(newps, NULL, data);
//        //fc = (*f)(cx,data);
//	}
//    Parameters_set_value(newps, 0, x );
    
    //printf("a %e (%f) b %e (%f) c %e (%f)\n", a, fa, x, fx, b, fb);

	x = w = v = Parameters_value(ps,0);
	fw = fv = fx = f(ps, NULL, data);

    stop->f_eval_current = 1;
	
    stop->iter = 0;
    
	while ( *iter <= stop->iter_max ) {
        
		xm = 0.5*(a+b);
		tol2 = 2.0*(tol1 = tol*fabs(x) + ZEPS);
		
		if ( fabs(x - xm) <= (tol2 - 0.5*(b-a)) ) {
			Parameters_set_value(ps, 0, x );

			*fmin = fx;
            //printf("%s xmin %e (%f)\n\n", Parameters_name(ps, 0), x, fx);
			return OPT_SUCCESS;
		}
		
		// Construct a trial parabolic fit.
		if ( fabs(e) > tol1 ) {
			r = ( x - w )*( fx - fv );
			q = ( x - v )*( fx - fw );
			p = ( x - v )*q - ( x - w )*r;
			q = 2.0*( q - r );
			if (q > 0.0) p = -p;
			q = fabs(q);
			
			etemp = e;
			e = d;
			
			// Determine the acceptability of the parabolic fit (e.g. collinearity of the 3 points)
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)){
				d = CGOLDEN_RATIO * (e = (x >= xm ? a-x : b-x)); // Take the golden section step into the larger of the two segments.
			}
			//Take the parabolic step.
			else {
				d = p/q;
				u = x+d;
				if ( u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm-x);
			}
		}
		// Golden section search
		else {
			d = CGOLDEN_RATIO * (e = (x >= xm ? a-x : b-x));
		}
		
		u  = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		//This is the one function evaluation per iteration.
		Parameters_set_value(ps, 0, u );
		fu = f(ps, NULL, data);
        stop->f_eval_current++;
		
		if (fu <= fx) {
			if (u >= x) a = x;
			else b = x;
			shift4(&v,&w,&x,u);
			shift4(&fv,&fw,&fx,fu);
		}
		else {
			if (u < x) a = u;
			else b = u;
			if (fu <= fw || w == x) {
				//Now decide what to do with our function evaluation.
				//Housekeeping follows:
				v  = w;
				w  = u;
				fv = fw;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v  = u;
				fv = fu;
			}
		}
        (*iter)++;
	}
	fprintf(stderr,"Too many iterations in brent_one_d:%s\n",Parameters_name(ps,0));
	Parameters_set_value(ps, 0, x );
	*fmin = fx;
	return OPT_MAXITER;
}

