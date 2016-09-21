/*
 *  frpmrn.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/15/10.
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


#include "frpmrn.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "matrix.h"
#include "mathconstant.h"
#include "linefunction.h"


void conjugate_gradient( const int n, double *sdir, const double *gvec, const double *gold, const bool *active, opt_algorithm algorithm );
double gradientProjection( const int n, const double *sdir, const double *gvec );
void steepestDescentDirection( const int n, double *sdir, double *gvec, bool *active);

static double _findStep( LineFunction *lf, double f0, double s0, double lastStep, int *numFun );

static double _computeDerivative(LineFunction *lf, double lambda, int *numFun );

/*! Conjugate gradient in multidimension.
 * Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func, 
 * using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is input as ftol.
 * Returned quantities are p (the location of the minimum), iter (the number of iterations that were performed), 
 * and fret (the minimum value of the function). The routine linmin is called to perform line minimizations.
*/
//void frprmn(double p[], int n, double ftol, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double [])){


opt_result frprmn_optimize( Parameters *xvec, opt_func f, void *data, OptStopCriterion stop, double *fmin, opt_algorithm algorithm ){
    int prin = 0, numArgs;
    double *backup = NULL;
    //Parameters *xvec = x;
    
	numArgs = Parameters_count(xvec);
	
    backup = dvector(numArgs);
    Parameters_store_value(xvec, backup);

	
	opt_result status = OPT_KEEP_GOING;
	
	// line function
	LineFunction *lf = new_LineFunction( xvec, f, data );
	
	// Force points to be between bounds
	LineFunction_force_within_bounds( lf, xvec );
	
	
	// function value and gradient at current guess
	int *numFun = &stop.f_eval_current;
	double fx = 0;
	//double numGrad = 0;
	double *gvec = dvector(numArgs);
	
    // Evaluate function and calculate the gradient
	fx = f(xvec, gvec, data);
	(*numFun)++;
	*numFun += 2*numArgs;
    
    double fx_start = fx;
    
	double *gold = clone_dvector(gvec, numArgs);	

	// init stop condition
	opt_check_stop( &stop, xvec, fx );
	
	// currently active variables
	bool *active = bvector(numArgs);
	
    // Set variable as not active if its value are outside its bounds
	int numActive = LineFunction_set_active_parameters(lf, xvec, gvec, active);
	
	// if no variables are active return
	if ( numActive == 0 ){
		free(active);
		free(gvec);
		free(gold);
		free_LineFunction(lf);
		*fmin = fx;
        free(backup);
		return OPT_NEED_CHECK;
	}
	
	// initial search direction (steepest descent)
	double *sdir = dvector(numArgs);
    
    // Compute direction of steepest descent (active == true ? -gvec : 0)
	steepestDescentDirection(numArgs, sdir, gvec, active);

	LineFunction_update(lf, xvec, sdir);
	
	// Slope at start point in initial direction
	double slope = gradientProjection(numArgs, sdir, gvec);
	
	
	if (prin > 0){
		fprintf(stdout,"--- starting minimization ---\n");
		fprintf(stdout,"... current parameter settings ...\n");
		fprintf(stdout,"...   numArgs   ... %d\n", numArgs);
		fprintf(stdout,"...   tolx   ... %f\n", stop.tolx);
		fprintf(stdout,"...   tolfx   ... %f\n", stop.tolfx);
		fprintf(stdout,"... maxFun  ... %d\n", stop.f_eval_max);
		fprintf(stdout,"... numActive  ... %d\n", numActive);
		fprintf(stdout,"... slope  ... %f\n\n", slope);
		
		fprintf(stdout,"... start vector ...\n");
        for ( int i = 0; i < numArgs; i++ ) {
            fprintf(stdout,"%s value: %f active: %d direction: %f\n", Parameters_name(xvec, i), Parameters_value(xvec, i), active[i], sdir[i]);
        }
		fprintf(stderr, "\n\n");
	}
	
	int numLin = 0;
	
	double defaultStep = 1.;
	double lastStep = defaultStep;

	while( status == OPT_KEEP_GOING ){
		// determine an appropriate step
        if( prin > 0 ) printf("Find step lower: %e upper: %e\n", lf->lower, lf->upper);
		double step = _findStep( lf, fx, slope, lastStep, numFun );
		lastStep = step;
		numLin++;
        
		// update xvec using step
		LineFunction_set_parameters(lf, step, xvec);
        // force xvec to be within bounds (lf not needed)
		LineFunction_force_within_bounds( lf, xvec );
        
		if( prin > 0 ) printf("Step: %e LnL %f\n", step, fx);
		// Evaluate function using xvec
		fx = f(xvec, gvec, data);
        if( prin > 0 ) printf("LnL %f\n\n", fx);
		(*numFun)++;
		*numFun += 2*numArgs;
        
		// test for for convergence		
		if ( (status = opt_check_stop( &stop, xvec, fx )) != OPT_KEEP_GOING ){
			//fprintf(stderr, "opt_check_stop status %d iter = %d\n", status, stop.iter_current );
			break;
		}
		
		// Determine which parameters should be active based on their gradient
		numActive = LineFunction_set_active_parameters(lf, xvec, gvec, active);
		
		// if all variables are inactive return
		if (numActive == 0){
			break;
		}
		
		// Determine new search direction (sdir)
		conjugate_gradient(numArgs, sdir, gvec, gold, active, algorithm);
		
        // Directions (sdir) that point outside boundaries are set to 0 (lf not needed)
		LineFunction_constrain_direction(lf, xvec, sdir);
		
		// compute slope in new direction
		slope = gradientProjection(numArgs, sdir, gvec);
		
		if (slope >= 0){
			//reset to steepest descent direction
			steepestDescentDirection(numArgs, sdir, gvec, active);
			
			// compute slope in new direction
			slope = gradientProjection(numArgs, sdir, gvec);
			
			// reset to default step length
			lastStep = defaultStep;
		}
		
		
		// other updates
		LineFunction_update(lf, xvec, sdir);
	
		memcpy(gold, gvec, numArgs * sizeof(double) );
		
		if (prin > 1){
            double ff = f(xvec, NULL, data);
			fprintf(stderr, "\nFunction value: %f recalculated %f diff %e\n", fx, ff, (ff-fx) );
			fprintf(stderr, "... new vector ...\n");
            for ( int i = 0; i < numArgs; i++ ) {
                fprintf(stdout,"%s value: %f active: %d direction: %f\n", Parameters_name(xvec, i), Parameters_value(xvec, i), active[i], sdir[i]);
            }
			fprintf(stderr, "... numFun  ... %d\n", *numFun);
			
			fprintf(stderr, "... numLin  ... %d\n", numLin);
            fprintf(stdout,"... numActive  ... %d\n", numActive);
            fprintf(stdout,"------------------------------------------------------------------------------------------\n\n");
		}
		
	}
	
	if (prin > 0){
		fprintf(stderr, "\n");
		fprintf(stderr, "... final vector ...\n");
		fprintf(stderr,"... numFun  ... %d", *numFun);
		fprintf(stderr,"... numLin  ... %d", numLin);
		fprintf(stderr,"... LnL  ... %f (%f)", fx, fx_start);
        fprintf(stdout,"... numActive  ... %d\n", numActive);
		fprintf(stderr,"\n--- end of minimization ---\n\n");
	}
	
    if( fx > fx_start ){
        Parameters_restore_value(xvec, backup);
        status = OPT_FAIL;
        *fmin = fx_start;
    }
    else {
        *fmin = fx;
    }
    
	
	free(active);
	free(gvec);
	free(gold);
	free(sdir);
	free_LineFunction(lf);
    free(backup);
	return status;
}



void conjugate_gradient( const int n, double *sdir, const double *gvec, const double *gold, const bool *active, opt_algorithm algorithm ){
	double gg, dgg; 
	dgg = gg = 0.0;
	int i = 0;
	for ( i = 0; i < n; i++ ) {
		if ( active [i] ){
			switch ( algorithm ) {
				case OPT_CG_FR:{
					gg  += gvec[i] * gvec[i];
					dgg += gold[i] * gold[i];					
					break;
				}
				case OPT_CG_PR:{
					gg += gold[i] * gold[i];
					dgg += (gvec[i] - gold[i]) * gvec[i];					
					break;
				}
				default:
					assert(0);
			}
		}
	}
	
	double gam = dgg/gg;
	
	if( gam < 0 || gg == 0 ){
		gam = 0;
	}
	for ( i = 0; i < n; i++ ) {
		if ( active[i] ) {
			sdir[i] = -gvec[i] + gam * sdir[i];
		}
		else {
			sdir[i] = 0.0;
		}
		
	}
}

double gradientProjection( const int n, const double *sdir, const double *gvec ){
	double slope = 0;
	for (int i = 0; i < n; i++){
		slope += gvec[i] * sdir[i];
	}
	return slope;
}

// Direction of steepest descent
// r = -f'(x)
void steepestDescentDirection( const int n, double *sdir, double *gvec, bool *active){
	for (int i = 0; i < n; i++){
		if ( active[i] ){
			sdir[i] = -gvec[i];
		}
		else{
			sdir[i] = 0;
		}
	}
}

// Find alpha that minimize f(x_1) such that f(x_1) = f(x_0) + alpha * r_0
// where r_0 = -f'(x_0) the direction of steepest descent
double _findStep( LineFunction *lf, double f0, double s0, double lastStep, int *numFun ){
	// f0 function value at step = 0
	// s0 slope at step = 0
	
	double step;
	double maxStep = lf->upper;
    //printf("maxStep %e\n", maxStep);
	
	if (maxStep <= 0 || s0 == 0){
		return 0.0;
	}
	
	//step = Math.abs(lf.findMinimum());
	
	
	// growing/shrinking factors for bracketing
	
	double g1 = 2.0;
	double g2 = 1.25;
	double g3 = 0.5;
	
	// x1 and x2 try to bracket the minimum
	
	double x1 = 0;
	double s1 = s0;
	double x2 = lastStep*g2;
	if(x2 > maxStep){
		x2 = maxStep*g3;
	}
	double s2 = _computeDerivative( lf, x2, numFun );
	// we need to go further to bracket minimum
	bool boundReached = false;
	while (s2 <= 0 && !boundReached){
		x1 = x2;
		s1 = s2;
		x2 = x2*g1;
		if (x2 > maxStep){
			x2 = maxStep;
			boundReached = true;
		}
		s2 = _computeDerivative( lf, x2, numFun );
	}
	
	// determine step length by quadratic interpolation
	// for minimum in interval [x1,x2]
	
	if (s2 <= 0){
		// true local minimum could NOT be bracketed
		// instead we have a local minimum on a boundary
		
		step = x2;
	}
	else{
		// minimum is bracketed
		
		step = (x1*s2-x2*s1)/(s2-s1);
		// note that nominator is always positive
	}
	
	// just to be safe - should not be necessary
	if (step >= maxStep){
		step = maxStep;
	}
	if (step < 0){
		step = 0;
	}
	//printf("step %e [%d]\n", step,boundReached);
	return step;
}

// This function has changed the current parameters.
double _computeDerivative(LineFunction *lf, double x, int *numFun ){
	*numFun += 2;
	double h = SQRT_EPS*(fabs(x) + 1.0);
    
    double fxplus = LineFunction_evaluate( lf, dmin(x+h,lf->upper));
    double fxminus = LineFunction_evaluate( lf, x-h);
    
    //printf("_computeDerivative %e h %e d %e\n",x,h,(fxplus - fxminus)/(2.0*h));
	return (fxplus - fxminus)/(2.0*h);
}


void gradient( opt_func f, void *data, Parameters *x, double *grad ){	
	for (int i = 0; i < Parameters_count(x); i++){
		double h = SQRT_EPS*(fabs(Parameters_value(x,i)) + 1.0);
		
		double oldx = Parameters_value(x,i);
		Parameters_set_value(x,i, oldx + h);
		double fxplus = f(x, NULL, data);
		Parameters_set_value(x,i, oldx - h);
		double fxminus = f(x, NULL, data);
		Parameters_set_value(x,i, oldx);
		
		// Centered first derivative
		grad[i] = (fxplus-fxminus)/(2.0*h);
	}
}
