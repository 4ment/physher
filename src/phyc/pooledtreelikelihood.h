/*
 *  pooledtreelikelihood.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 19/10/12.
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


#ifndef PhyC_pooledtreelikelihood_h
#define PhyC_pooledtreelikelihood_h

#include "treelikelihood.h"
#include "parameters.h"

struct _PooledTreeLikelihood;

typedef struct _PooledTreeLikelihood PooledTreeLikelihood;

typedef void (*pFreePooledTreeLikelihood)(PooledTreeLikelihood*);

struct _PooledTreeLikelihood{
	SingleTreeLikelihood **tlks;
	
	unsigned count;
	unsigned capacity;
	Parameters **backup;
	unsigned backup_count;
	bool share_sitepattern;
	bool share_sitemodel;
	
	pFreePooledTreeLikelihood free;
};

PooledTreeLikelihood * new_PooledTreeLikelihood( SingleTreeLikelihood *tlk, const unsigned count, bool share_sitepattern, bool share_sitemodel );

void PooledTreeLikelihood_add_paramters( PooledTreeLikelihood *pool, Parameters *ps );

PooledTreeLikelihood * new_PooledTreeLikelihood_List( const unsigned count );

void PooledTreeLikelihood_add( PooledTreeLikelihood *pool, SingleTreeLikelihood *tlk );

#endif
