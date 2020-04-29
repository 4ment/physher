/*
 *  pooledtreelikelihood.c
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

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "pooledtreelikelihood.h"



double _calculate_pooled_likelihood( PooledTreeLikelihood *pool );

#pragma mark -
// MARK: PooledTreeLikelihood

static void _free_PooledTreeLikelihood( PooledTreeLikelihood *pool );

/** A collection of independent Likelihood objects where the site pattern can be shared.
 * The first Likelihood object is the original one (tlk) and the others are clones of tlk
 * Used by OPENMP ready algorithms
 */
PooledTreeLikelihood * new_PooledTreeLikelihood( SingleTreeLikelihood *tlk, const unsigned count, bool share_sitepattern, bool share_sitemodel ){
	assert( count > 0 );
	PooledTreeLikelihood *pool = (PooledTreeLikelihood*)malloc( sizeof(PooledTreeLikelihood) );
	assert(pool);
	
	pool->count = 0;
	pool->capacity = count;
	pool->tlks = (SingleTreeLikelihood**)malloc( pool->capacity * sizeof(SingleTreeLikelihood*) );
	assert(pool->tlks);
	
	// use the original tlk
	pool->tlks[0] = tlk;
	pool->count++;
	
	for (int i = 1; i < count; i++) {
		pool->tlks[i] = clone_SingleTreeLikelihood_share(tlk, share_sitepattern, share_sitemodel );
		pool->count++;
	}
	
	pool->share_sitepattern = share_sitepattern;
	pool->share_sitemodel = share_sitemodel;
	pool->backup = NULL;
	pool->backup_count = 0;
	pool->free = _free_PooledTreeLikelihood;
	
	return pool;
}

PooledTreeLikelihood * new_PooledTreeLikelihood_List( const unsigned count ){
	PooledTreeLikelihood *pool = (PooledTreeLikelihood*)malloc( sizeof(PooledTreeLikelihood) );
	assert(pool);
	pool->capacity = count;
	pool->count = 0;
	pool->tlks = (SingleTreeLikelihood**)malloc( count * sizeof(SingleTreeLikelihood*));
	assert(pool->tlks);
	for (int i = 0; i < count; i++) {
		pool->tlks[i] = NULL;
	}
	pool->share_sitepattern = false;
	pool->share_sitemodel = false;
	pool->backup = NULL;
	pool->backup_count = 0;
	return pool;
}



void PooledTreeLikelihood_add( PooledTreeLikelihood *pool, SingleTreeLikelihood *tlk ){
	if( pool->count == pool->capacity){
		pool->capacity++;
		pool->tlks = realloc(pool->tlks, pool->capacity * sizeof(SingleTreeLikelihood*) );	
	}
	pool->tlks[pool->count] = tlk;
	pool->count++;
}

void PooledTreeLikelihood_remove( PooledTreeLikelihood *pool, const unsigned index ){
	int i = index;
	while ( i < pool->count-1) {
		pool->tlks[i] = pool->tlks[i+1];
	}
	pool->tlks[i] = NULL;
	pool->count--;
}

// The original Likelihood is conserved
static void _free_PooledTreeLikelihood( PooledTreeLikelihood *pool ){
	for (int i = 1; i < pool->count; i++){
		free_SingleTreeLikelihood_share( pool->tlks[i], pool->share_sitepattern, pool->share_sitemodel );
	}
	
	//	if ( pool->share_sitepattern ) {
	//		for (int i = 1; i < pool->count; i++) {
	//			free_SingleTreeLikelihood_SitePattern_shared( pool->tlks[i] );
	//		}
	//	}
	//	else {
	//		for (int i = 1; i < pool->count; i++) {
	//			free_SingleTreeLikelihood( pool->tlks[i] );
	//		}
	//	}
	
	free(pool->tlks);
	for (int i = 0; i < pool->backup_count; i++) {
		free_Parameters( pool->backup[i] );
	}
	free(pool->backup);
	free(pool);
	pool = NULL;
}

void PooledTreeLikelihood_add_paramters( PooledTreeLikelihood *pool, Parameters *ps ){
	if ( pool->backup == NULL ) {
		pool->backup = (Parameters**)malloc( sizeof(Parameters *) );
		assert(pool->backup);
	}
	else {
		pool->backup = realloc( pool->backup, (pool->count+1) * sizeof(Parameters *) );
		assert(pool->backup);
	}
	pool->backup[pool->backup_count++] = ps;	
}
