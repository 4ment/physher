/*
 *  bracket.c
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

#include "bracket.h"

#include <math.h>

#include "optimizer.h"
#include "parameters.h"
#include "utils.h"

#define TINY 1.0e-20
#define GOLD 1.618034  // GOLD is the default ratio by which successive intervals are magnified;
#define GLIMIT 100.0   // GLIMIT is the maximum magnification allowed for a parabolic-fit step.



inline static double sign( double a, double b){
	return b >= 0.0 ? fabs(a) : -fabs(a);
}



/*! Initially bracketting a minimum.
 * Given a function func, and given distinct initial points ax and bx, this routine searches in the downhill direction 
 * (defined by the function as evaluated at the initial points) and returns new points ax, bx, cx that bracket a minimum of the function.
 * Also returned are the function values at the three points, fa, fb, and fc.
 * a < b < c
 * fb < (fa,fb)
 */
//void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double)){
//void bracket(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, directionfunction *func ){
void bracket2( Parameters *ps, double *ax, double *bx, double *cx, opt_func f, void *data ){
	double fa=0, fb=0, fc=0;
	double ulim,u,r,q,fu;
	
	Parameter *p = Parameters_at(ps,0);
	
	Parameter_set_value(p, *ax);
	fa = f(ps, NULL, data);
	//fprintf(stderr, "fa -> a:%f (%f) b:%f (%f) c:%f (%f)\n",*ax,fa,*bx,fb,*cx,fc);
	
	Parameter_set_value(p, *bx);
	fb = f(ps, NULL, data);
	//fprintf(stderr, "fb -> a:%f (%f) b:%f (%f) c:%f (%f)\n",*ax,fa,*bx,fb,*cx,fc);
	
	//Switch roles of a and b so that we can go downhill in the direction from a to b.
	if (fb > fa) {
		dswap(ax, bx);
		dswap(&fb, &fa);
	} 
	
	
	*cx = (*bx)+GOLD*(*bx-*ax);
	Parameter_set_value(p, *cx);//First guess for c.
	fc = f(ps, NULL, data);
	//fprintf(stderr, "fc -> %f (%f) %f (%f) %f (%f)\n",*ax,fa,*bx,fb,*cx,fc);
	
	//Keep returning here until we bracket. Compute u by parabolic extrapolation from a, b, c. TINY is used to prevent any possible division by zero.
	while (fb > fc) {
		r = (*bx-*ax)*(fb-fc);
		q = (*bx-*cx)*(fb-fa);
		
		//if( q == r ) q+= DBL_EPSILON;
		
		u = (*bx)-((*bx-*cx)*q - (*bx-*ax)*r)/ (2.0*sign(dmax(fabs(q-r),TINY),q-r));
		//fprintf(stderr, "bracket.c: u=%f - %f %f %f %f %f fa=%f fb=%f fc=%f -%f %f\n",u, *cx,q,*bx,*ax,r, fa, fb,fc,   ((*bx)-((*bx-*cx)*q - (*bx-*ax)*r)), (2.0*sign(dmax(fabs(q-r),TINY),q-r)) );
		u = check_value(p->cnstr, u);
		
		
		ulim = (*bx)+GLIMIT*(*cx-*bx);
		//if (ulim < func->lower_bound) ulim = func->lower_bound;
		
		
		//fprintf(stdout, "%f %f\n", u,ulim);
		
		// Parabolic u is between b and c
		if ((*bx-u)*(u-*cx) > 0.0) {
			Parameter_set_value(p, u);
			fu = f(ps, NULL, data);
			
			// Got a minimum between b and c.
			if (fu < fc) {
				*ax = (*bx);
				*bx = u;
				fa = fb; 
				fb = fu;
				//fprintf(stderr, "1] %f (%f) %f (%f) %f (%f)\n",*ax,fa,*bx,fb,*cx,fc);
				return;
			}
			//Got a minimum between between a and u.
			else if (fu > fb) { 
				*cx = u;
				fc = fu;
				//fprintf(stderr, "2] %f (%f) %f (%f) %f (%f)\n",*ax,fa,*bx,fb,*cx,fc);
				return;
			} 
			//Parabolic failed, use default magnification.
			u = (*cx)+GOLD*(*cx-*bx);
			Parameter_set_value(p, u);
			fu = f(ps, NULL, data);
			
		} 
		//Parabolic fit is between c and its allowed limit.
		else if ((*cx-u)*(u-ulim) > 0.0) {
			Parameter_set_value(p, u);
			fu = f(ps, NULL, data);
			// u is not a minimum
			if (fu < fc) {
				//shift4(bx, cx, &u, check_value(p->cnstr, *cx+GOLD*(*cx-*bx) ) );
				shift4(bx, cx, &u, check_value(p->cnstr, u+GOLD*(u- *cx) ) );
				Parameter_set_value(p, u);
				shift4(&fb, &fc, &fu, f(ps, NULL, data) );
			} 
		}
		//Limit parabolic u to maximum allowed value.
		else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u = ulim;
			Parameter_set_value(p, u);
			fu = f(ps, NULL, data);
		}
		else { //Reject parabolic u, use default magnification.
			u = (*cx)+GOLD*(*cx-*bx) ;
			Parameter_set_value(p, u);
			fu = f(ps, NULL, data);
		} 
		shift4(ax, bx, cx, u); //Eliminate oldest point and continue.
		shift4(&fa, &fb, &fc, fu);
		//fprintf(stderr, " -- %f (%f) %f (%f) %f (%f)\n",*ax,fa,*bx,fb,*cx,fc);
	}
	
	//fprintf(stderr, "3] %f (%f) %f (%f) %f (%f)\n",*ax,fa,*bx,fb,*cx,fc);
}
