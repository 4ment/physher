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
    int *lineages;
    double *times;
    bool *iscoalescent;
    int n;
    double (*calculate)( struct Coalescent* );
    bool need_update;
}Coalescent;

typedef struct exponentialdemography{
	demography type;
	char *name;
	double n0;
	double lower;
	double upper;
	double r;
}exponentialdemography;

#pragma mark Coalescent

void free_Coalescent( Coalescent *coalescent );

Coalescent * clone_Coalescent( const Coalescent *coal, Tree *tree );

#pragma mark -
#pragma mark Constant coalescent

Coalescent * new_ConstantCoalescent( Tree *tree, double theta );

void free_ConstantCoalescent( Coalescent *coal );

Coalescent * clone_ConstantCoalescent( const Coalescent *coal, Tree *tree );


#endif
