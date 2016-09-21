/*
 *  mixedtreelikelihood.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/12/2015.
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

#ifndef PhyC_mixedtreelikelihood_h
#define PhyC_mixedtreelikelihood_h

#include "treelikelihood.h"
#include "parameters.h"

typedef struct MixedTreeLikelihood{
    SingleTreeLikelihood **tlks;
	int n;
	Parameters *params;
    
	double *categories;
	double *proportions;
    
    bool need_update;
    
    double   (*get_category)( struct MixedTreeLikelihood *, const int );
	double   (*get_proportion)( struct MixedTreeLikelihood *, const int );
    
	double   (*calculate)( struct MixedTreeLikelihood * );
    double lnl;
	
} MixedTreeLikelihood;

double optimize_mixedtreelikelihoods( MixedTreeLikelihood *mixed );

MixedTreeLikelihood * new_MixedStrictTreeLikelihood( SingleTreeLikelihood *tlk, int n );

void free_MixedTreeLikelihood( MixedTreeLikelihood *mixed );

#endif
