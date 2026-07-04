// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _LINE_FUNCTION_H_
#define _LINE_FUNCTION_H_

#include "parameters.h"
#include "optimizer.h"

typedef struct LineFunction{
	Parameters *parameters;
	
	// Parameters *x;
	
	double *s;
	double *xi; // d
	
	size_t dim; // number of parameters
	size_t size; // total number of values
	
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

void LineFunction_set_parameters(const LineFunction *lf, const double lambda);

bool LineFunction_force_within_bounds(const LineFunction *lf);

void LineFunction_update( LineFunction *lf, double *xi );

int LineFunction_constrain_direction( const LineFunction *lf, double *xi );

int LineFunction_set_active_parameters( const LineFunction *lf, const double *grad, bool *active);

#endif
