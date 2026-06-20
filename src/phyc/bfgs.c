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

#define STPMX 100.0  // Scaled maximum step length allowed in line searches
#define EPS 3.0e-8   // Machine accuracy
#define TOLX (4*EPS) // Convergence criterion on x values

/*
 * Bound-constrained BFGS minimization of func, using the gradient supplied by
 * grad_f. Starting from p, it returns the minimizer in p and the minimum value
 * in fmin. Bounds are enforced with an active set (coordinates pinned on a bound
 * whose gradient points outward are frozen) together with clamping in the line
 * search (lnsrch). Stops on a small relative step (TOLX), a small projected
 * gradient (stop.tolg), or after stop.iter_max iterations.
 */

opt_result dfpmin_optimize( Parameters *p, opt_func f, opt_grad_func grad_f, void *data, OptStopCriterion stop, double *fmin, double alpha){

	int i,its,j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

	opt_result status = OPT_SUCCESS;

	int maxeval = stop.iter_max;
	double gtol = stop.tolg;

	size_t n = Parameters_size(p);
	size_t paramCount = Parameters_count(p);

	double *dg      = dvector(n);
	double *g       = dvector(n);
	double *hdg     = dvector(n);
	double **hessin = dmatrix(n,n);
	double *xi      = dvector(n);
	bool *active    = bvector(n);
	double *lwr     = dvector(n);  // per-scalar lower bounds
	double *upr     = dvector(n);  // per-scalar upper bounds
	double *pnew    = dvector(n);

	// Evaluate the function and gradient at the starting point.
	fp = f(NULL, NULL, data);
	grad_f(p, g, data);

	// Initialize the inverse-Hessian approximation to the identity.
	for (size_t i = 0; i < n; i++) {
		memset(hessin[i], 0, n*sizeof(double));
		hessin[i][i] = 1.0;
	}

	// Initial direction is steepest descent; record per-scalar bounds.
	size_t index = 0;
	for (i = 0; i < paramCount; i++ ) {
		Parameter* param = Parameters_at(p, i);
		const double* values = Parameter_values(param);
		memcpy(pnew + index, values, Parameter_size(param)*sizeof(double));

		double lower = Parameter_flower(param);
		double upper = Parameter_fupper(param);

		for (size_t j = 0; j < Parameter_size(param); j++) {
			xi[index] = -g[index];
			sum  += values[j]*values[j];
			lwr[index] = lower;
			upr[index] = upper;
			active[index] = true;
			index++;
		}
	}

	stpmax = STPMX * dmax(sqrt(sum),(double)n);

	// Main loop over the iterations.
	for ( its = 0; its < maxeval; its++ ) {
		status = lnsrch( p, pnew, f, data, fp, g, xi, fmin, stpmax, alpha );

		// The new function evaluation occurs in lnsrch; save it for the next one.
		fp = *fmin;

		// Move to the new point and recover the step actually taken (xi = s);
		// lnsrch may have clamped it to the feasible box.
		size_t index = 0;
		for (size_t i = 0; i < paramCount; i++) {
			Parameter* param = Parameters_at(p, i);
			const double* values = Parameter_values(param);
			size_t idx = index;
			for (size_t j = 0; j < Parameter_size(param); j++) {
				xi[index] = pnew[index] - values[j];
				index++;
			}
			Parameter_set_values(param, pnew + idx);
		}

		// Test for convergence on the relative step size.
		test = 0.0;
		for ( i = 0; i < n; i++ ) {
			temp = fabs(xi[i])/dmax(fabs(pnew[i]), 1.0);
			if (temp > test) test = temp;
		}
		if (test < TOLX) {
			status = OPT_SUCCESS;
			break;
		}

		// Save the old gradient, then compute the new one at the new point.
		memcpy(dg, g, n*sizeof(double));
		grad_f(p, g, data);

		// Update the active set. Freeze a coordinate that sits on a bound while
		// its gradient points out of the feasible region (its optimum lies on the
		// boundary); release one that has moved off its bound. Frozen coordinates
		// are excluded from the search direction and the Hessian update via the
		// active[] guards below, so the bound-constrained problem is solved on the
		// free subspace. Without this, clamping in the line search pins variables
		// at a bound and the shared step length couples them, stalling the search.
		// Done before the convergence test so a coordinate that is optimal on its
		// bound does not keep the test from being satisfied.
		for ( i = 0; i < n; i++ ) {
			bool at_lower = pnew[i] <= lwr[i] + EPS && g[i] > 0.0;
			bool at_upper = pnew[i] >= upr[i] - EPS && g[i] < 0.0;
			bool now_active = !(at_lower || at_upper);
			if ( now_active != active[i] ) {
				// status changed: drop stale curvature for this coordinate
				for ( j = 0; j < n; j++ ) { hessin[i][j] = 0.0; hessin[j][i] = 0.0; }
				hessin[i][i] = 1.0;
			}
			active[i] = now_active;
		}

		// Test for convergence on the projected (zero) gradient. Only free
		// coordinates count; a frozen one already satisfies the KKT conditions.
		test = 0.0;
		den = dmax(*fmin,1.);
		for ( i = 0; i < n; i++ ) {
			if ( !active[i] ) continue;
			temp = fabs(g[i])*dmax(fabs(pnew[i]), 1.)/den;
			if (temp > test) test = temp;
		}
		if (test < gtol) {
			status = OPT_SUCCESS;
			break;
		}

		// Compute difference of gradients (dg = y).
		for ( i = 0; i < n; i++ ) dg[i] = g[i] - dg[i];

		// On the first update, rescale the identity inverse-Hessian to
		// H0 = (s^T y)/(y^T y) I (Nocedal & Wright eq. 6.20). An identity H0
		// combined with a gradient whose magnitude differs from the parameters'
		// by orders of magnitude gives a badly scaled first step and the method
		// stalls; this rescaling fixes the conditioning.
		if (its == 0) {
			double sy = 0.0, yy = 0.0;
			for (i = 0; i < n; i++) {
				if (active[i]) { sy += xi[i]*dg[i]; yy += dg[i]*dg[i]; }
			}
			if (sy > 0.0 && yy > 0.0) {
				double scale = sy/yy;
				for (i = 0; i < n; i++) {
					if (active[i]) hessin[i][i] = scale;
				}
			}
		}

		// hdg = H * y
		for ( i = 0; i < n; i++ ) {
			hdg[i] = 0.0;
			if ( !active[i] ) continue;

			for ( j = 0; j < n; j++ ){
				if ( active[j] ) {
					hdg[i] += hessin[i][j]*dg[j];
				}
			}
		}

		// Dot products for the denominators: fac = y^T s, fae = y^T H y.
		fac = fae = sumdg = sumxi = 0.0;
		for ( i = 0; i < n; i++ ) {
			if ( active[i] ) {
				fac += dg[i]*xi[i];
				fae += dg[i]*hdg[i];
				sumdg += SQR(dg[i]);
				sumxi += SQR(xi[i]);
			}
		}
		// Curvature condition: skip the update unless y^T s > sqrt(eps ||y||^2 ||s||^2),
		// which keeps the approximation positive definite and avoids tiny denominators.
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

		// Calculate the next direction to go.
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
	free(pnew);
	free_dmatrix(hessin,n);
	free(hdg);
	free(g);
	free(dg);
	free(active);
	free(lwr);
	free(upr);

	if ( its == maxeval ) {
		status = OPT_MAXEVAL;
	}
	stop.iter = its;

	return status;
}
