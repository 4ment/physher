/*
 *  demographicmodels.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/2/10.
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


#ifndef _DEMOGRAPHIC_MODELS_H_
#define _DEMOGRAPHIC_MODELS_H_

#include "parameters.h"
#include "tree.h"
#include "discreteparameter.h"
#include "utils.h"

typedef enum demography{
	COALESCENT_CONSTANT,
	COALESCENT_EXPONENTIAL,
	COALESCENT_SKYGRID,
	COALESCENT_SKYLINE,
	COALESCENT_SKYLINE_CLASSIC,
	COALESCENT_SKYRIDE
}demography;


typedef struct Coalescent{
    Tree *tree;
	demography type;
	Parameters *p;
	double logP;
	double stored_logP;
    int *lineages;
	int *stored_lineages;
	double *times;
	double *stored_times;
	double_int_pair_t** nodes; // indexes of nodes corresponding to interval
    bool *iscoalescent;
	bool *stored_iscoalescent;
	int n;
	double (*calculate)( struct Coalescent* );
	double (*dlogP)( struct Coalescent*, const Parameter* );
	double (*d2logP)( struct Coalescent*, const Parameter* );
	double (*ddlogP)( struct Coalescent*, const Parameter*, const Parameter* );
    bool need_update;
	bool need_update_intervals;
	double* grid;
	size_t gridCount;
	DiscreteParameter* groups;
}Coalescent;


#pragma mark Coalescent

Model* new_CoalescentModel_from_json(json_node* node, Hashtable* hash);

void free_Coalescent( Coalescent *coalescent );

#pragma mark -

Coalescent * new_ConstantCoalescent( Tree* tree, Parameter* theta );
Coalescent * new_ConstantCoalescent_with_data( Parameter* theta, double* times, bool* coalescent, int size );

Coalescent * new_ExponentialCoalescent( Tree *tree, Parameters* parameters );
Coalescent * new_ExponentialCoalescent_with_data( Parameters* parameters, double* times, bool* coalescent, int size );

Coalescent * new_ClassicalSkylineCoalescent_with_parameters( Tree *tree, Parameters* parameters);

Coalescent * new_SkylineCoalescent( Tree *tree, Parameters* parameters, DiscreteParameter* groups);
Coalescent * new_SkylineCoalescent_with_data(Parameters* parameters, double* times, bool* coalescent, int size, DiscreteParameter* groups);

Coalescent * new_SkyrideCoalescent( Tree *tree, Parameters* parameters);
Coalescent * new_SkyrideCoalescent_with_data(Parameters* parameters, double* times, bool* coalescent, int size);

Coalescent * new_GridCoalescent( Tree *tree, Parameters* parameters, int grid, double cutoff );
Coalescent * new_GridCoalescent_with_data(Parameters* parameters, double* times, bool* coalescent, int size, int grid, double cutoff);

#endif
