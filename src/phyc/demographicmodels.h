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

typedef enum demography{CONSTANT_DEMOGRAPHY, EXP_DEMOGRAPHY}demography;

typedef struct Coalescent{
    Tree *tree;
	demography type;
    Parameters *p;
	double logP;
    int *lineages;
	int *stored_lineages;
	double *times;
	double *stored_times;
	int* nodes; // indexes of nodes corresponding to interval
    bool *iscoalescent;
	bool *stored_iscoalescent;
	int n;
	double (*calculate)( struct Coalescent* );
	double (*dlogP)( struct Coalescent*, const Parameter* );
	double (*d2logP)( struct Coalescent*, const Parameter* );
	double (*ddlogP)( struct Coalescent*, const Parameter*, const Parameter* );
    bool need_update;
	bool need_update_intervals;
}Coalescent;


#pragma mark Coalescent

Model* new_CoalescentModel_from_json(json_node* node, Hashtable* hash);

void free_Coalescent( Coalescent *coalescent );

#pragma mark -
#pragma mark Constant coalescent

Coalescent * new_ConstantCoalescent_with_parameter( Tree *tree, Parameter* theta );


Coalescent * clone_ConstantCoalescent( const Coalescent *coal, Tree *tree );

#pragma mark -
#pragma mark Constant coalescent

Coalescent * new_ExponentialCoalescent_with_parameters( Tree *tree, Parameters* parameters );


#endif
