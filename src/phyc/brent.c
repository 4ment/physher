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

opt_result brent_optimize2( Parameter *parameter, size_t index, opt_func f, void *data, OptStopCriterion *stop, double *fminp );

/*
 * 
 */

opt_result serial_brent_optimize( Parameters *ps, opt_func f, void *data, OptStopCriterion *stop, double *fmin ){
	//for(int j = 0; j < stop->iter_min; j++)
	for(int i = 0; i < Parameters_count(ps); i++){
		Parameter* parameter = Parameters_at(ps, i);
		stop->iter = 0;
		stop->f_eval_current = 0;
		for(size_t j = 0; j < Parameter_size(parameter); j++){
			opt_result status = brent_optimize2(parameter, j, f, data, stop, fmin);
		}
	}
	return OPT_SUCCESS;
}

//opt_result brent_optimize( Parameters *ps, opt_func f, void *data, const int maxeval, const double tol, double *fmin ){
opt_result brent_optimize( Parameters *ps, opt_func f, void *data, OptStopCriterion *stop, double *fminp ){
	Parameter* parameter = Parameters_at(ps, 0);
	return brent_optimize2(parameter, 0, f, data, stop, fminp);
}


/**
 * Find a bracket [a,b,c] for a minimum of f(x)
 * given an initial guess x and bounds (lower, upper).
 *
 * Works for domains (0,1), [0,∞), (-∞,∞), etc.
 * Detects boundary minima and gracefully handles x at a boundary.
 */
void find_bracket(opt_func f, void* data, Parameter *xx, size_t index,
                  double *a, double *b, double *c,
                  double step_factor, int max_iter)
{
    if (step_factor <= 1.0){
        step_factor = 1.618034;  // golden ratio default
	}

	double lower = Parameter_lower(xx);
	double upper = Parameter_upper(xx);

    double eps = 1e-8;
    double width = (isinf(lower) || isinf(upper)) ? 1.0 : upper - lower;

	double x = Parameter_value_at(xx, index);
	if(isinf(x)){
		// GUARD: the optimizer handed us an infinite starting point. This happens when
		// an unconstrained coordinate runs off to +/-inf because the constrained
		// quantity it maps from has reached the edge of its domain -- e.g. a simplex
		// proportion collapsing to 0 (or 1) under the logit stick-breaking transform
		// drives its free coordinate to logit(0) = -inf (see simplex.c). Rather than
		// aborting, snap x back to a finite point on the same side so Brent has a valid
		// bracket origin and can pull the parameter back into the interior.
		//
		// Target preference: a finite hard bound, else a moderate finite magnitude.
		// RESET_MAGNITUDE is deliberately modest: large enough not to clobber a
		// genuinely large value, but small enough that for the logit-style transforms
		// that trigger this it still sits where the transform has usable gradient, so
		// the parameter can actually recover.
		const double RESET_MAGNITUDE = 20.0;
		if (x < 0) {
			x = isfinite(lower) ? lower : -RESET_MAGNITUDE;
		} else {
			x = isfinite(upper) ? upper : RESET_MAGNITUDE;
		}
		Parameter_set_value_at(xx, x, index);
		// fprintf(stderr, "find_bracket: initial x was infinite for %s index %zu, reset to %g\n",
		//         Parameter_name(xx), index, x);
	}

    // --- 1️⃣ Handle cases where x is at or near a boundary ---
	if (x <= lower + eps * width) {
		// printf("x at lower boundary %e %e %s\n", x, lower, Parameter_name(xx));
		double f0 = f(NULL, NULL, data);
		
        // double delta = 0.1 * (isinf(upper) ? fmax(1.0, fabs(x)) : width);
		double delta = 10 * fabs(x);
        double x1 = fmin(upper, x + delta);
		Parameter_set_value_at(xx, x1, index);
        double f1 = f(NULL, NULL, data);

        if (f0 <= f1) {
			// printf("lower f0 <= f1 %e < %e x0: %e x1: %e\n", f0, f1, x, x1);
            // Minimum at lower boundary
            // *a = *b = lower;
			*a = x;
			*b = x;
            *c = x1;
            return;
        } else {
            // Move into interior and expand upward
			// printf("%s.%lu lower f0 > f1 %e > %e x0: %e x1: %e\n", Parameter_name(xx), index, f0, f1, x, x1);
			while(f0 > f1){
				*a = x;
				f0 = f1;
				x = x1;
				x1 = fabs(x1)*10.0;
				Parameter_set_value_at(xx, x1, index);
        		f1 = f(NULL, NULL, data);
			}
			*b = x;
			*c = x1;
			// printf("%e %e %e == %e < %e\n\n", *a, *b, *c, f0, f1);
            // *a = lower;
            // *b = x1;
            // double delta2 = step_factor * delta;
            // *c = fmin(upper, *b + delta2);
            return;
        }
		assert(f0 < f1);
		assert(*a <= *b && *b < *c);
		if(f0 > f1 || *a > *b || *b >= *c){
			fprintf(stderr, "Warning: invalid bracket found: [%g, %g, %g] f0: %e  f1: %e\n", *a, *b, *c, f0, f1);
		}
		return;
    }

    if (x >= upper - eps * width) {
		// printf("x at upper boundary %e %e %s\n", x, upper, Parameter_name(xx));
        // x = upper;
		// Parameter_set_value_at(xx, x, index);
        double f0 = f(NULL, NULL, data);
        // double delta = 0.1 * (isinf(lower) ? fmax(1.0, fabs(x)) : width);
		double delta = 0.1 * fabs(x);
        double x1 = fmax(lower, x - delta);
		Parameter_set_value_at(xx, x1, index);
        double f1 = f(NULL, NULL, data);

        if (f0 <= f1) {
            // Minimum at upper boundary
			// printf("upper f0 <= f1 %e < %e x0: %e x1: %e\n", f0, f1, x, x1);
            // *a = x1;
            // *b = *c = upper;
			*a = x1;
			*b = x;
            *c = x;
        } else {
            // Move into interior and expand downward
			// printf("upper f0 > f1 %e > %e x0: %e x1: %e\n", f0, f1, x, x1);
            // *c = upper;
            // *b = x1;
            // double delta2 = step_factor * delta;
            // *a = fmax(lower, *b - delta2);
			while(f0 > f1){
				*c = x;
				f0 = f1;
				x = x1;
				x1 = fabs(x1)*0.1;
				Parameter_set_value_at(xx, x1, index);
        		f1 = f(NULL, NULL, data);
			}
			*b = x;
			*a = x1;
        }

		assert(f0 < f1);
		assert(*a < *b && *b <= *c);
		if(f0 > f1 || *a >= *b || *b > *c){
			fprintf(stderr, "Warning: invalid bracket found: [%g, %g, %g] f0: %e  f1: %e\n", *a, *b, *c, f0, f1);
		}
		return;
	}

    // --- 2️⃣ Interior starting point ---
    double delta;
    if (isinf(lower) && isinf(upper)) {
        delta = 0.1 * fmax(1.0, fabs(x));
	} else if (isfinite(lower) && isinf(upper)) {
		// Upper is unbounded, so the local scale is |x|, not the (possibly
		// far-away) distance to the lower bound. Capping at fmax(1, |x|) keeps
		// the bracket near x when lower is very negative, while leaving the
		// lower >= 0 case unchanged (there x - lower <= x <= fmax(1, |x|)).
		delta = 0.1 * fmin(fmax(1.0, fabs(x)), x - lower);
    } else if (isinf(lower) && isfinite(upper)) {
        delta = 0.1 * fmax(1.0, upper - x);
    } else if (lower == 0.0 && upper == 1.0) {
        delta = 0.05 * fmin(x, 1.0 - x);
    } else {
        delta = 0.1 * (upper - lower);
    }

    *a = fmax(lower, x - delta);
    *b = x;
    *c = fmin(upper, x + delta);
	// fprintf(stderr, "Initial bracket: a=%g b=%g c=%g\n", *a, *b, *c);

	Parameter_set_value_at(xx, *a, index);
    double fa = f(NULL, NULL, data);
	Parameter_set_value_at(xx, *b, index);
    double fb = f(NULL, NULL, data);
	Parameter_set_value_at(xx, *c, index);
    double fc = f(NULL, NULL, data);

    // --- 3️⃣ Check for boundary minima (both finite sides) ---
    // if (!isinf(lower)) {
	// 	Parameter_set_value_at(xx, lower, index);
    //     double fl = f(NULL, NULL, data);
	// 	Parameter_set_value_at(xx, x + eps * fmax(1.0, fabs(lower)), index);
    //     double fl_eps = f(NULL, NULL, data);
    //     if (fl <= fl_eps && fl <= fb) {
    //         *a = lower;
    //         *b = lower;
    //         *c = fmin(lower + delta, upper);
	// 		printf("lower %s\n", Parameter_name(xx));
    //         return;
    //     }
    // }
    // if (!isinf(upper)) {
	// 	Parameter_set_value_at(xx, upper, index);
    //     double fu = f(NULL, NULL, data);
	// 	Parameter_set_value_at(xx, upper - eps * fmax(1.0, fabs(upper)), index);
    //     double fu_eps = f(NULL, NULL, data);
    //     if (fu <= fu_eps && fu <= fb) {
	// 		printf("upper %s\n %e < %e < %e\n", Parameter_name(xx), *a, *b, *c, fu, fu_eps, fb);
    //         *a = fmax(lower, upper - delta);
    //         *b = upper;
    //         *c = upper;
    //         return;
    //     }
    // }

    // --- 4️⃣ Already bracketed? ---
    if (fb < fa && fb < fc)
        return;

    // Expanding towards an infinite bound, the objective may keep decreasing
    // yet flatten onto an asymptote (bounded below as x -> +/-inf). There is
    // then no interior minimum to bracket -- the infimum is at the unreachable
    // bound. Without a guard the expansion marches off to +/-1e40, emits a
    // "could not find a valid bracket" warning, and leaves the parameter pegged
    // at a value from which the outer optimiser cannot recover (the local slope
    // there is numerically zero). One way this arises is a free coordinate that
    // maps onto a constrained quantity collapsing to its boundary, but it is a
    // generic property of the objective and does not assume any transform.
    // We detect the flattening and pin a boundary-style bracket at the current
    // best finite point, which is at the asymptote to within flat_tol and stays
    // in a region where the parameter remains recoverable. MAX_REACH is a hard
    // backstop on how far past the start we are willing to travel.
    const double MAX_REACH = 60.0;
    const double flat_tol = 1e-9 * fmax(1.0, fabs(fb));

    // --- 5️⃣ Expand towards decreasing direction ---
    if (fc < fb) {
        for (int i = 0; i < max_iter; ++i) {
            *a = *b;
            fa = fb;
            *b = *c;
            fb = fc;
            delta *= step_factor;
            *c = fmin(upper, *b + delta);
			Parameter_set_value_at(xx, *c, index);
            fc = f(NULL, NULL, data);
            if (fb < fa && fb < fc) return;
            if (isinf(upper) && (fb - fc <= flat_tol || *b > x + MAX_REACH)) {
                // Asymptotic minimum against the +inf bound: pin at best point.
                *c = *b;  // [interior(higher), best, best], fa > fb
                return;
            }
            if (*c >= upper) return;
        }
    } else {
        for (int i = 0; i < max_iter; ++i) {
            *c = *b;
            fc = fb;
            *b = *a;
            fb = fa;
            delta *= step_factor;
            *a = fmax(lower, *b - delta);
			Parameter_set_value_at(xx, *a, index);
            fa = f(NULL, NULL, data);
            if (fb < fa && fb < fc) return;
            if (isinf(lower) && (fb - fa <= flat_tol || *b < x - MAX_REACH)) {
                // Asymptotic minimum against the -inf bound: pin at best point.
                *a = *b;  // [best, best, interior(higher)], fc > fb
                return;
            }
            if (*a <= lower) return;
        }
    }

    fprintf(stderr, "Warning: could not find a valid bracket around x=%g [%g, %g]\n",
            x, Parameter_lower(xx), Parameter_upper(xx));
}



opt_result brent_optimize2( Parameter *parameter, size_t index, opt_func f, void *data, OptStopCriterion *stop, double *fminp ){

	double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e = 0.0;
	double d = 0;
	
	double tol = stop->tolx;
    //int *numFun = &stop.f_eval_current;
    size_t *iter = &stop->iter;
	// double lbound = Parameter_lower(parameter);
	// double ubound = Parameter_upper(parameter);
	// double xmin = ( Parameter_flower(parameter) < Parameter_fupper(parameter) ? Parameter_flower(parameter) : Parameter_fupper(parameter));
	// double xmax = ( Parameter_flower(parameter) > Parameter_fupper(parameter) ? Parameter_flower(parameter) : Parameter_fupper(parameter));
	x = Parameter_value_at(parameter, index);

// // 	if(isinf(Parameter_lower(parameter)) && isinf(Parameter_upper(parameter))){
// // 		if(x == 0.0){
// // 			a = x - 1;
// // 			b = x + 1;
// // 		}else{
// // 			a = x - fabs(x)/2.0;
// // 			b = x + fabs(x)*2.0;
// // 		}
// // 		Parameter_set_value_at(parameter, a, index);
// // 		double fa = f(NULL, NULL, data);
// // 		Parameter_set_value_at(parameter, b, index);
// // 		double fb = f(NULL, NULL, data);
// // 		Parameter_set_value_at(parameter, x, index);
// // 		fx = f(NULL, NULL, data);
// // 		// printf("%zu Brent optimization: a = %f, b = %f, x = %f fa = %f, fb = %f, fx = %f\n", index, a, b, x, fa, fb, fx);
// // 		int i = 0;
// // 		while(fa < fx && a > lbound) {
// // 			// a = max(xmin, a - (c - a) * expansion)
// // 			a = fmax(lbound, a - (b - a) * 2.0);
// //             Parameter_set_value_at(parameter, a, index);
// // 			fa = f(NULL, NULL, data);
// // 			printf("%f %f\n", a, fa);
// // 			i++;
// // 			// if(i > 10){
// // 			// 	break;
// // 			// }
// // 		}
// // 		while(fb < fx && b < ubound) {
// // 			// b = min(xmax, b + (b - c) * expansion)
// // 			b = fmin(ubound, b + (b - a) * 2.0);
// // 			Parameter_set_value_at(parameter, b, index);
// // 			fb = f(NULL, NULL, data);
// // 		}
// // 		printf("Brent optimization: a = %f, b = %f, x = %f fa = %f, fb = %f, fx = %f\n\n", a, b, x, fa, fb, fx);

// // if(isinf(a)){
// // 	exit(3);
// // }
// // 	}
// // 	else{
// 	if(x == 0.0){
// 		a = x - 10;
// 		b = x + 10;
// 	}
// 	else{
// 		a = fmax(xmin, x - fabs(x)/2.0);
// 		b = fmin(xmax, x + fabs(x)*2.0);
// 	}
// 	x = a + (b - a)/2.0;
	
// 	Parameter_set_value_at(parameter, a, index);
// 	double fa = f(NULL, NULL, data);
//     Parameter_set_value_at(parameter, b, index);
// 	double fb = f(NULL, NULL, data);
//     Parameter_set_value_at(parameter, x, index);
// 	fx = f(NULL, NULL, data);

// 	while(fa < fx && a > xmin) {
// 		x = a;
// 		fx = fa;
//         a = (a+xmin)/2.0;
// 		if (a < 2.0*xmin){
//             a = xmin;
// 		}
//         Parameter_set_value_at(parameter, a, index);
//         fa = f(NULL, NULL, data);
// 		// printf("%f %f [%f %f]\n", a, fa, xmin, xmax);
// 	}

// 	while(fb < fx && b < xmax) {
// 		x = b;
// 		fx = fb;
//         b = (b+xmax)/2.0;
//         if (b > xmax * 0.95)
//             b = xmax;
//         Parameter_set_value_at(parameter, b, index);
//         fb = f(NULL, NULL, data);
// 	}
//     Parameter_set_value_at(parameter, x, index);
// // }

	find_bracket(f, data, parameter, index, &a, &x, &b, 1.618034, 100);
	Parameter_set_value_at(parameter, x, index);
	fx = f(NULL, NULL, data);
	// if(isinf(a) || isinf(x) || isinf(b)|| a >= b || x <= a || x >= b || isinf(fx)){
	// 	fprintf(stderr,"Brent optimization: invalid bracket fx: %e x: %e [%e, %e] for %s index %zu\n", fx, x, a, b, Parameter_name(parameter), index);
	// 	exit(2);
	// }

	x = w = v = x;
	fw = fv = fx;

    stop->f_eval_current = 1;
	
    stop->iter = 0;
    
	while ( *iter <= stop->iter_max ) {
        
		xm = 0.5*(a+b);
		tol2 = 2.0*(tol1 = tol*fabs(x) + ZEPS);
		
		if ( fabs(x - xm) <= (tol2 - 0.5*(b-a)) ) {
			Parameter_set_value_at(parameter, x, index );

			*fminp = fx;
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
		Parameter_set_value_at(parameter, u, index );
		fu = f(NULL, NULL, data);
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
	fprintf(stderr,"Too many iterations in brent_one_d: %s index: %zu x: %f fx old %f -> new %f\n",Parameter_name(parameter), index, x, *fminp, fx);
	Parameter_set_value_at(parameter, x, index);
	*fminp = fx;
	return OPT_MAXITER;
}

