/*
 *  optimimize.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/12/10.
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

#ifndef _OPTIMIZE_H_
#define _OPTIMIZE_H_

#include "treelikelihood.h"



typedef double (*function_f)(void *);

//double optimize_treelikelihood( TreeLikelihood *tlk );

double optimize_singletreelikelihood( SingleTreeLikelihood *stlk );

double optimize_singletreelikelihood2( SingleTreeLikelihood *stlk );

typedef struct MultivariateData{
	SingleTreeLikelihood *tlk;
    function_f f;
	int *map_params;
    void *model; // extra pointer if we use coalescent or penalized likelihood
}MultivariateData;

typedef struct BrentData{
	SingleTreeLikelihood *tlk;
    function_f f;
	int index_param;
	double *backup;
    void *model; // extra pointer if we use coalescent or penalized likelihood
    double threshold;
} BrentData;



BrentData * new_BrentData( SingleTreeLikelihood *tlk);
void free_BrentData( BrentData *data );

BrentData * new_MultiBrentData( SingleTreeLikelihood *tlk, int *map );
void free_MultiBrentData( BrentData *data );

MultivariateData * new_MultivariateData( SingleTreeLikelihood *tlk, int *map );
void free_MultivariateData( MultivariateData *data );



double optimize_brent_branch_length_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, int iterations );

double optimize_brent_branch_length( Parameters *params, double *grad, void *data );


double standard_loglikelihood_brent( void *data );

double standard_loglikelihood_upper_brent( void *data );

// Time

double optimize_brent_height_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param );

double optimize_brent_rate_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param );

double optimize_scale_rate_height_strict( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk );
double optimize_scale_rate_height_strict2( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk );



void calculate_constrained(Node *node, BranchModel *bm);
#endif
