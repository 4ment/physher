/*
 *  linefunction.c
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

#include "linefunction.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "utils.h"
#include "matrix.h"
#include "mathconstant.h"

static void _LineFunction_compute_bounds( LineFunction *lf );

LineFunction *new_LineFunction( Parameters *x, opt_func f, void *data ){
	LineFunction *linfunc = (LineFunction*)malloc( sizeof(LineFunction) );
	assert(linfunc);
	
	linfunc->dim = Parameters_count(x);
	
	linfunc->parameters = x; // original parameters. I should clone x in linefunc->x in order to keep the boundaries
	
	linfunc->x = new_Parameters(Parameters_count(x));
	for (int i = 0; i < Parameters_count(x); i++) {
		Parameters_add(linfunc->x, new_Parameter(Parameters_name(x, i), 0, NULL) );
	}
	
	linfunc->s = dvector(Parameters_count(x));	
	linfunc->xi = dvector(Parameters_count(x));
	
	linfunc->f = f;
	linfunc->data = data;
	return linfunc;
}

void free_LineFunction( LineFunction *lf ){
	free_Parameters(lf->x);
	lf->parameters = NULL;
	free(lf->s);
	free(lf->xi);
	lf->f = NULL;
	lf->data = NULL;
	free(lf);
	lf = NULL;
}

//WTF
double LineFunction_minimize( LineFunction *lf ){
	Optimizer *brent = new_Optimizer(OPT_BRENT);
	double fret = 0;
	free_Optimizer(brent);
	return fret;
}

double LineFunction_evaluate( LineFunction *lf, double lambda ){
	LineFunction_set_parameters( lf, lambda, lf->x );
	return lf->f(lf->x, NULL, lf->data);
}


// p = s + lambda * xi
// where are s and xi are local variables
// Does not modify the LineFunction object
void LineFunction_set_parameters( const LineFunction *lf, const double lambda, Parameters *p ){
	for (int i = 0; i < lf->dim; i++){
		Parameters_set_value(p, i, lf->s[i] + lambda * lf->xi[i] );
	}
}


/**
 * check (and modify, if necessary) whether a point lies properly
 * within the predefined bounds
 *
 * @param p coordinates of point
 *
 * @return true if p was modified, false otherwise
 */
bool LineFunction_force_within_bounds( const LineFunction *lf, Parameters *p ){
	bool modified = false;
	//fprintf(stderr, "-------------------------\nLineFunction_check_point\n");
	for (int i = 0; i < lf->dim; i++){
		//fprintf(stderr, "%s %f\n",Parameters_name(p, i) , Parameters_value(p, i));
		if( Parameters_value(p, i) < Parameters_lower(lf->parameters,i) ){
			Parameters_set_value(p, i, Parameters_lower(lf->parameters,i) );
			modified = true;
		}
		if( Parameters_value(p, i) > Parameters_upper(lf->parameters,i) ){
			Parameters_set_value(p, i, Parameters_upper(lf->parameters,i) );
			modified = true;
		}
	}
	//fprintf(stderr, "%s\n\n", (modified ? "modified" : "not modified"));
	return modified;
}

// Update lf->s and lf->xi vectors with p and xi
void LineFunction_update( LineFunction *lf, Parameters *p, double *xi ){
	for (int i = 0; i < lf->dim; i++){
		lf->s[i] = Parameters_value(p, i);
	}
	memcpy(lf->xi, xi, lf->dim * sizeof(double) );
	
	_LineFunction_compute_bounds( lf );
}

int LineFunction_set_active_parameters( const LineFunction *lf, const Parameters *p, const double *grad, bool *active){
	int numActive = 0;
	//fprintf(stderr, "-------------------------\nLineFunction_check_variables\n");
	for (int i = 0; i < lf->dim; i++){
		active[i] = true;
		if ( Parameters_value(p, i) <= Parameters_lower(lf->parameters,i)+EPS ){
			// no search towards lower boundary
			if ( grad[i] > 0 ){
				active[i] = false;
			}
		}
		else if ( Parameters_value(p, i) >= Parameters_upper(lf->parameters,i)-EPS ){
			// no search towards upper boundary
			if ( grad[i] < 0 ){
				active[i] = false;
			}
		}
		else{
			numActive++;
		}
		//ftolfx(stderr, "%s %d\n", Parameters_name(p, i), active[i] );
	}
	//fprintf(stderr, "\n");
	return numActive;
}


int LineFunction_constrain_direction( const LineFunction *lf, const Parameters *p, double *dir ){
	int n = 0;
	for (int i = 0; i < lf->dim; i++){
		// no search towards lower boundary
		if (  Parameters_value(p, i) <= Parameters_lower(lf->parameters,i)+ EPS ){
			if( dir[i] < 0 ){
				dir[i] = 0;
				n++;
			}
		}
		// no search towards upper boundary
		else if ( Parameters_value(p, i) >= Parameters_upper(lf->parameters,i) - EPS ){
			if( dir[i] > 0 ){
				dir[i] = 0;
				n++;
			}
		}
	}
	
	return n;
}

void _LineFunction_compute_bounds( LineFunction *lf ){
	bool firstVisit = true;
	double upper, lower;
	for (int i = 0; i < lf->dim; i++){
		if ( lf->xi[i] != 0){
			upper = ( Parameters_upper(lf->parameters,i) - lf->s[i])/lf->xi[i];
			lower = ( Parameters_lower(lf->parameters,i) - lf->s[i])/lf->xi[i];
			if (lower > upper){
				dswap(&upper, &lower);
			}
			
			if (firstVisit){
				lf->lower = lower;
				lf->upper = upper;
				firstVisit = false;
			}
			else {
				if ( lower > lf->lower ){
					lf->lower = lower;
                    //printf("%s %f l\n", Parameters_name(lf->parameters, i), lower);
				}
				if ( upper < lf->upper ){
                    lf->upper = upper;
                    //printf("%s %f u %e %e %e\n", Parameters_name(lf->parameters, i), upper,Parameters_lower(lf->parameters,i), lf->s[i],lf->xi[i]);
				}
			}
		}
	}
}
