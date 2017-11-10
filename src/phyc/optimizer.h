/*
 *  optimizer.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/4/11.
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


#ifndef _OPTIMIZER_H_
#define _OPTIMIZER_H_

#include "parameters.h"

#include "time.h"

#define OPT_XTOL 1.0e-3

struct _Optimizer;
typedef struct _Optimizer Optimizer;

typedef enum opt_result{
	OPT_MAXEVAL    = -1,
	OPT_MAXTIME    = -2,
    OPT_MAXITER    = -3,
	OPT_ERROR      = -4,
	OPT_NEED_CHECK = -5,
	OPT_FAIL       = -6, // when the likelihood is lower than then input
	
	OPT_KEEP_GOING =  0,
	
	OPT_SUCCESS     = 1
}opt_result;

typedef enum opt_algorithm {
	OPT_POWELL = 0,
	OPT_BRENT,
	OPT_SERIAL_BRENT,
	OPT_BFGS,
	OPT_CG_PR,
	OPT_CG_FR,
	OPT_META,
}opt_algorithm;

static const char *OPT_ALGORITHMS[5] = {"POWELL","BRENT","BFGS","CG_PR","CG_FR"};

typedef struct OptStopCriterion{
	time_t time_start;
	time_t time_current;
	time_t time_end;
	double time_max;
	
	// Number of iteration
	int iter_min;
	int iter_max;
	int iter;
	
	// Number of evaluation of objective function
	int f_eval_max;
	int f_eval_current;
	
	double tolfx;
	double tolx;
	double tolg; // gradient tolerance
	
	double oldfx;
	double *oldx;
	
	int count;
}OptStopCriterion;


typedef double (*opt_func)( Parameters *x, double *gradient, void *data );

typedef bool (*opt_update_data)( void *data, Parameters *p);

struct _OptimizerSchedule;
typedef struct _OptimizerSchedule OptimizerSchedule;

typedef bool(*OptimizerSchedule_post)(OptimizerSchedule*, double before, double after);

struct _OptimizerSchedule{
	Optimizer** optimizers;
	Parameters** parameters;
	OptimizerSchedule_post* post;
	int* rounds;
	int count;
	int capacity;
};


Optimizer * new_Optimizer( opt_algorithm algorithm );

void free_Optimizer( Optimizer *opt );

Optimizer* clone_Optimizer(Optimizer *opt, void* data, Parameters* parameters);


opt_result opt_optimize( Optimizer *opt, Parameters *ps, double *fmin );

opt_result opt_optimize_univariate( Optimizer *opt, Parameter *p, double *fmin );

opt_result opt_maximize( Optimizer *opt, Parameters *ps, double *fmin );

opt_result opt_maximize_univariate( Optimizer *opt, Parameter *p, double *fmin );

void opt_set_objective_function( Optimizer *opt, opt_func f );

void opt_set_data( Optimizer *opt, void *data );

void opt_set_max_evaluation( Optimizer *opt, const int maxeval );

void opt_set_max_iteration( Optimizer *opt, const int maxiter );

void opt_set_time_max( Optimizer *opt, const double maxtime );
void opt_set_time_max_minutes( Optimizer *opt, const double maxtime );
void opt_set_time_max_hours( Optimizer *opt, const double maxtime );
void opt_set_time_max_relative( Optimizer *opt, const time_t time, const double factor );

void opt_set_tolfx( Optimizer *opt, const double tolfx );

void opt_set_tolx( Optimizer *opt, const double tolx );

double opt_tolx( Optimizer *opt );

void opt_set_tolg( Optimizer *opt, const double tolg );

int opt_iterations( Optimizer *opt );

int opt_f_evaluations( Optimizer *opt );

void opt_set_update_data_function( Optimizer *opt, opt_update_data uf );

opt_result opt_check_stop( OptStopCriterion *stop, Parameters *x, double fx );

void opt_add_optimizer(Optimizer *opt_meta, Optimizer *opt, const Parameters* parameters);

OptimizerSchedule* opt_get_schedule(Optimizer *opt_meta);

#endif
