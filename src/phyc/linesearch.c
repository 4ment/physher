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
#include <string.h>

#include "parameters.h"
#include "utils.h"
#include "optimizer.h"
#include "matrix.h"

#define ALF 1.0e-4	//Ensures sufficient decrease in function value.
#define TOLX 1.0e-7	//Convergence criterion on ∆x.



static double dot(const double *a, const double *b, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

static void vec_add(double *out, const double *x, const double *p, double alpha, int n) {
    for (int i = 0; i < n; i++) out[i] = x[i] + alpha * p[i];
}

// Safeguarded interpolation of a trial step inside the bracket (a_lo, a_hi).
// Uses a quadratic model from f_lo, the projected slope g_lo_p at a_lo, and
// f_hi; falls back to bisection when the model is degenerate or the minimizer
// lands outside the bracket. a_lo and a_hi may be given in either order.
static double safeguarded_interp(double a_lo, double f_lo, double g_lo_p,
                                 double a_hi, double f_hi) {
    double lo = a_lo < a_hi ? a_lo : a_hi;
    double hi = a_lo < a_hi ? a_hi : a_lo;
    double width = hi - lo;
    double bisect = 0.5 * (a_lo + a_hi);

    double denom = 2.0 * (f_hi - f_lo - g_lo_p * (a_hi - a_lo));
    if (denom == 0.0) return bisect;

    double a = a_lo - g_lo_p * (a_hi - a_lo) * (a_hi - a_lo) / denom;
    double margin = 0.1 * width;
    if (!isfinite(a) || a <= lo + margin || a >= hi - margin) return bisect;
    return a;
}

// Stage 2 of the strong-Wolfe search: refine a step inside a bracket known to
// contain a point satisfying the conditions (Nocedal & Wright, Alg. 3.6).
static LineSearchResult zoom(Parameters* parameters, opt_func fun, opt_grad_func grad_f,
                             void* data, const double *x, const double *p,
                             double f0, double g0p,
                             double a_lo, double f_lo, double g_lo_p,
                             double a_hi, double f_hi,
                             double c1, double c2) {
    size_t n = Parameters_size(parameters);
    LineSearchResult res = {0};
    double *x_trial = malloc(n * sizeof(double));
    double *g_j = malloc(n * sizeof(double));
    int max_iter = 30;
    double a_j = a_lo;
    double f_j = f_lo;

    for (int i = 0; i < max_iter; i++) {
        a_j = safeguarded_interp(a_lo, f_lo, g_lo_p, a_hi, f_hi);
        vec_add(x_trial, x, p, a_j, n);
        Parameters_restore_value(parameters, x_trial);
        f_j = fun(parameters, NULL, data);
        res.nfev++;

        if (!isfinite(f_j) || f_j > f0 + c1 * a_j * g0p || f_j >= f_lo) {
            a_hi = a_j;
            f_hi = f_j;
        } else {
            grad_f(parameters, g_j, data);
            double gtp = dot(g_j, p, n);
            if (fabs(gtp) <= -c2 * g0p) {
                res.fx = f_j;
                res.alpha = a_j;
                res.status = 0;
                free(x_trial);
                free(g_j);
                return res;
            }
            if (gtp * (a_hi - a_lo) >= 0.0) {
                a_hi = a_lo;
                f_hi = f_lo;
            }
            a_lo = a_j;
            f_lo = f_j;
            g_lo_p = gtp;
        }
        if (fabs(a_hi - a_lo) < 1e-12) break;
    }

    res.fx = f_j;
    res.alpha = a_j;
    res.status = 1;
    free(x_trial);
    free(g_j);
    return res;
}

// Strong-Wolfe line search (bracketing stage; Nocedal & Wright, Alg. 3.5).
// Function values are obtained from `fun` and gradients from `grad_f` so the
// gradient sign always matches the minimized objective. `g` is the gradient at
// x, `p` a descent direction, `alpha0` the initial trial step.
LineSearchResult strong_wolfe_line_search(Parameters* parameters, opt_func fun,
                                          opt_grad_func grad_f, void* data, const double *x,
                                          const double *p, const double *g, double f0,
                                          double c1, double c2, double alpha0, double amax) {
    size_t n = Parameters_size(parameters);
    LineSearchResult res = {0};
    double g0p = dot(g, p, n);

    double *x_trial = malloc(n * sizeof(double));
    double *g_new = malloc(n * sizeof(double));

    // Not a descent direction: signal the caller to restart.
    if (g0p >= 0.0) {
        res.alpha = 0.0;
        res.fx = f0;
        res.status = 2;
        free(x_trial);
        free(g_new);
        return res;
    }

    double a_prev = 0.0;
    double f_prev = f0;
    double gprev_p = g0p;
    double a_cur = fmin(alpha0 > 0.0 ? alpha0 : 1.0, amax);

    for (int i = 0; i < 30; i++) {
        vec_add(x_trial, x, p, a_cur, n);
        Parameters_restore_value(parameters, x_trial);
        double f_cur = fun(parameters, NULL, data);
        res.nfev++;

        // Armijo failure (or non-finite value): minimizer is bracketed below.
        if (!isfinite(f_cur) || f_cur > f0 + c1 * a_cur * g0p ||
            (i > 0 && f_cur >= f_prev)) {
            LineSearchResult r = zoom(parameters, fun, grad_f, data, x, p, f0, g0p,
                                      a_prev, f_prev, gprev_p, a_cur, f_cur, c1, c2);
            r.nfev += res.nfev;
            free(x_trial);
            free(g_new);
            return r;
        }

        grad_f(parameters, g_new, data);
        double gtp = dot(g_new, p, n);

        if (fabs(gtp) <= -c2 * g0p) {  // strong Wolfe satisfied
            res.alpha = a_cur;
            res.fx = f_cur;
            res.status = 0;
            free(x_trial);
            free(g_new);
            return res;
        }

        if (gtp >= 0.0) {  // minimizer bracketed between a_cur and a_prev
            LineSearchResult r = zoom(parameters, fun, grad_f, data, x, p, f0, g0p,
                                      a_cur, f_cur, gtp, a_prev, f_prev, c1, c2);
            r.nfev += res.nfev;
            free(x_trial);
            free(g_new);
            return r;
        }

        a_prev = a_cur;
        f_prev = f_cur;
        gprev_p = gtp;
        res.fx = f_cur;
        res.alpha = a_cur;
        if (a_cur >= amax) break;
        a_cur = fmin(a_cur * 2.0, amax);
    }

    res.status = 1;
    free(x_trial);
    free(g_new);
    return res;
}


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
opt_result lnsrch( Parameters* parameters, double *x,  opt_func func, void *data, double fold, double *g, double *p, double *fmin, double stpmax, double alpha){
	int i;
	double a,alam,alamin,b,disc,rhs1,rhs2,slope,temp,test,tmplam;
	double alam2 = 0;
	double f2 = 0;
	
	opt_result status = OPT_ERROR;	
	
	int n = Parameters_size(parameters);
	double* xold = dvector(n);
	Parameters_store_value(parameters, xold);
	
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
	
	// Roundoff problem: not a descent direction
	if ( slope >= 0.0 ){
		*fmin = fold;
		Parameters_restore_value(parameters, xold);
		free(xold);
		return status;
	}
	
	test = 0.0;	//Compute λmin.
	for ( i = 0; i < n; i++ ) {					
		temp = fabs(p[i])/fmax(fabs(xold[i]), 1.0); 
		if (temp > test) test = temp;
	} 
	alamin = TOLX/test;
	alam = alpha;  // first trial step (full Newton step when alpha == 1)

	for (;;) {
		size_t index = 0;
		for (size_t i = 0; i < Parameters_count(parameters); i++ ){
			Parameter* param = Parameters_at(parameters, i);
			size_t idx = index;
			double lower = Parameter_flower(param);
			double upper = Parameter_fupper(param);
			for (size_t j = 0; j < Parameter_size(param); j++ ){
				x[index] = xold[index] + alam * p[index]; //x_i = xold_i + alam*p_i;
				// Project onto the feasible box: an unbounded quasi-Newton step can
				// otherwise drive parameters past their bounds (e.g. negative branch
				// lengths) where the objective is meaningless but still finite.
				if (x[index] < lower) x[index] = lower;
				else if (x[index] > upper) x[index] = upper;
				index++;
			}
			Parameter_set_values(param, x + idx);
		}
		*fmin = func(parameters, NULL, data);

		// Convergence on ∆x. For zero finding, the calling program should verify the convergence
		if (alam < alamin) {
			memcpy(x, xold, n*sizeof(double));
			*fmin = fold;
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
			// First backtrack: only f(alam) is known, so use a quadratic model.
			if (alam == alpha){
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
		alam  = fmax(tmplam,0.1*alam); //λ ≥ 0.1λ1
	}
	Parameters_restore_value(parameters, xold);
	free(xold);
	return status;
}

