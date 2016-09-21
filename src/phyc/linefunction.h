/*
 *  linefunction.h
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

#ifndef _LINE_FUNCTION_H_
#define _LINE_FUNCTION_H_

#include "parameters.h"
#include "optimizer.h"

typedef struct LineFunction{
	Parameters *parameters;
	
	Parameters *x;
	
	double *s;
	double *xi; // d
	
	int dim;
	
	double lower;
	double upper;
	
	//double (*evaluate)( struct LineFunction *, double );
	
	opt_func f; // the original objective function (e.g SingleLikelihood update and lk calculation)
	void *data; // and data (e.g. SingleLikelihood)
} LineFunction;

LineFunction *new_LineFunction( Parameters *x, opt_func f, void *data );

void free_LineFunction( LineFunction *lf );

double LineFunction_minimize( LineFunction *lf );

double LineFunction_evaluate( LineFunction *lf, double lambda );

void LineFunction_set_parameters( const LineFunction *lf, const double lambda, Parameters *p );

bool LineFunction_force_within_bounds( const LineFunction *lf, Parameters *p);

void LineFunction_update( LineFunction *lf, Parameters *p, double *xi );

int LineFunction_constrain_direction( const LineFunction *lf, const Parameters *p, double *xi );

int LineFunction_set_active_parameters( const LineFunction *lf, const Parameters *p, const double *grad, bool *active);

#endif
