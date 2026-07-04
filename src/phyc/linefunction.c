// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
	linfunc->size = Parameters_size(x);
	
	linfunc->parameters = x; // original parameters. I should clone x in linefunc->x in order to keep the boundaries
	
	// linfunc->x = new_Parameters(Parameters_count(x));
	// for (int i = 0; i < Parameters_count(x); i++) {
	// 	Parameters_move(linfunc->x, new_Parameter(Parameters_name(x, i), 0, NULL) );
	// }
	
	linfunc->s = dvector(linfunc->size);
	linfunc->xi = dvector(linfunc->size);

	linfunc->f = f;
	linfunc->data = data;
	linfunc->lower = -INFINITY;
	linfunc->upper = INFINITY;
	return linfunc;
}

void free_LineFunction( LineFunction *lf ){
	// free_Parameters(lf->x);
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
	LineFunction_set_parameters( lf, lambda );
	return lf->f(NULL, NULL, lf->data);
}


// p = s + lambda * xi
// where are s and xi are local variables
// Does not modify the LineFunction object
void LineFunction_set_parameters( const LineFunction *lf, const double lambda ){
	for(size_t i = 0; i < lf->dim; i++){
		Parameter* param = Parameters_at(lf->parameters, i);
		double* temp = dvector( Parameter_size(param) );
		for(size_t j = 0; j < Parameter_size(param); j++){
			temp[j] = lf->s[j] + lambda * lf->xi[j];
		}
		Parameter_set_values(param, temp);
		free(temp);
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
bool LineFunction_force_within_bounds( const LineFunction *lf ){
	bool modified = false;
	//fprintf(stderr, "-------------------------\nLineFunction_check_point\n");
	for (int i = 0; i < lf->dim; i++){
		Parameter* param = Parameters_at(lf->parameters, i);
		const double* values = Parameter_values(param);
		double lower = Constraint_flower(param->cnstr);
		double upper = Constraint_fupper(param->cnstr);
	
		for(size_t j = 0; j < Parameter_size(param); j++){
			//fprintf(stderr, "%s %f\n",Parameters_name(p, i) , Parameters_value(p, i));
			if( values[j] < lower ){
				Parameter_set_value_at(param, lower, j);
				modified = true;
			}
			if( values[j] > upper ){
				Parameter_set_value_at(param, upper, j);
				modified = true;
			}
		}
	}
	//fprintf(stderr, "%s\n\n", (modified ? "modified" : "not modified"));
	return modified;
}

// Update lf->s and lf->xi vectors with p and xi
void LineFunction_update( LineFunction *lf, double *xi ){
	Parameters_store_value(lf->parameters, lf->s);
	memcpy(lf->xi, xi, lf->size * sizeof(double));
	
	_LineFunction_compute_bounds( lf );
}

int LineFunction_set_active_parameters( const LineFunction *lf, const double *grad, bool *active){
	int numActive = 0;
	size_t index = 0;
	//fprintf(stderr, "-------------------------\nLineFunction_check_variables\n");
	for (int i = 0; i < lf->dim; i++){
		Parameter* param = Parameters_at(lf->parameters, i);
		const double* values = Parameter_values(param);
		double lower = Constraint_flower(param->cnstr);
		double upper = Constraint_fupper(param->cnstr);

		for(size_t j = 0; j < Parameter_size(param); j++){
		active[index] = true;
		if ( values[j] <= lower + EPS ){
			// no search towards lower boundary
			if ( grad[index] > 0 ){
				active[index] = false;
			}
		}
		else if ( values[j] >= upper - EPS ){
			// no search towards upper boundary
			if ( grad[index] < 0 ){
				active[index] = false;
			}
		}
		else{
			numActive++;
		}
		index++;
		//ftolfx(stderr, "%s %d\n", Parameters_name(p, i), active[i] );
		}
	}
	//fprintf(stderr, "\n");
	return numActive;
}


int LineFunction_constrain_direction( const LineFunction *lf, double *dir ){
	int n = 0;
	size_t index = 0;
	for (int i = 0; i < lf->dim; i++){
		Parameter* param = Parameters_at(lf->parameters, i);
		const double* values = Parameter_values(param);
		double lower = Constraint_flower(param->cnstr);
		double upper = Constraint_fupper(param->cnstr);
		
		for(size_t j = 0; j < Parameter_size(param); j++){
		// no search towards lower boundary
		if (  values[j] <= lower + EPS ){
			if( dir[index] < 0 ){
				dir[index] = 0;
				n++;
			}
		}
		// no search towards upper boundary
		else if ( values[j] >= upper - EPS ){
			if( dir[index] > 0 ){
				dir[index] = 0;
				n++;
			}
		}
		index++;
		}
	}
	
	return n;
}

void _LineFunction_compute_bounds( LineFunction *lf ){
	bool firstVisit = true;
	double lower = -INFINITY;
	double upper = INFINITY;
	size_t index = 0;
	for (int i = 0; i < lf->dim; i++){
		Parameter* param = Parameters_at(lf->parameters, i);
		double flower = Constraint_flower(param->cnstr);
		double fupper = Constraint_fupper(param->cnstr);

		for(size_t j = 0; j < Parameter_size(param); j++){
		if ( lf->xi[index] != 0){
			upper = ( fupper - lf->s[i])/lf->xi[index];
			lower = ( flower - lf->s[i])/lf->xi[index];
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
		index++;
		}
	}
}
