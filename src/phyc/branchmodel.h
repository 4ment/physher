/*
 *  branchmodel.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/8/10.
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


#ifndef _BRANCH_MODEL_H_
#define _BRANCH_MODEL_H_

#include "parameters.h"
#include "utils.h"
#include "tree.h"
#include "node.h"

#define BRANCHMODEL_RATE_MIN 1e-15
#define BRANCHMODEL_RATE_MAX 0.1

typedef enum branchmodel { CLOCK_STRICT, CLOCK_LOCAL, CLOCK_DISCRETE, CLOCK_RELAXED } branchmodel;

typedef enum relaxed_clock{RELAXED_LOGNORMAL, RELAXED_EXPONENTIAL, RELAXED_DISCRETE} relaxed_clock;

typedef struct BranchModel{
	
	branchmodel name;
	int id;
	Tree *tree;
	Parameters *rates;
	
	double (*get)( struct BranchModel *, Node * );
	void (*set)( struct BranchModel *, const int, const double );
    
    
	void (*free)( struct BranchModel *, bool );

	bool need_update;
	
	// LOCAL clock
	bool *indicators; // location of local clcoks indexed by id
	double *unscaled_rates; // and relaxed
	double scalefactor;
	
	// RELAXED
	relaxed_clock type;
	
	// local and discrete
	// each element represents a node indexed by id
	// Each value contains the index of the corresponding rate
	// These indexes are set in PREORDER from indicators. This map is not necessarily ordered (e.g. 1203 instead of 0123)
	DiscreteParameter *map;
} BranchModel;

#pragma mark -
#pragma mark BranchModel

Model* new_BranchModel_from_json(json_node*node, Hashtable*hash);

Model * new_BranchModel2( const char* name, BranchModel *bm, Model* tree);

BranchModel * new_BranchModel( Tree *tree, branchmodel type );

BranchModel * clone_BranchModel(const BranchModel *bm, Tree *tree );

Parameters * BranchModel_save_rates( BranchModel *bm );

void BranchModel_restore_rates( BranchModel *bm, const Parameters *rates );

void BranchModel_rates_to_vector( BranchModel *bm, double *rates );

void BranchModel_vector_to_rates( BranchModel *bm, const double *rates );

void BranchModel_value_to_rates( const double rate, BranchModel *bm );

int BranchModel_n_rate( BranchModel *bm );

#pragma mark -
#pragma mark StrickClock

BranchModel * new_StrictClock( Tree *tree );

BranchModel * new_StrictClock_with_parameters( Tree *tree, const Parameters *rates );

BranchModel * new_StrictClock_with_parameter( Tree *tree, Parameter *rate );


#pragma mark -
#pragma mark LocalClock

BranchModel * new_LocalClock( Tree *tree, const int nLocalClocks );

BranchModel * new_LocalClock_with_parameters( Tree *tree, const Parameters *rates );

BranchModel * new_LocalClock_from_tree( Tree *tree );

void localclock_set_random_clock_indicators( BranchModel *bm, const int nLocalClocks );

void LocalClock_set_number_of_clocks( BranchModel *bm, const int nLocalClocks );


void localclock_rebuild_map( BranchModel *bm );

void LocalClock_indicator_to_map( const bool *indicators, unsigned int *map, Node *node, int *index );

void localclock_set_indicators( BranchModel *bm, const bool *indicators );

void localclock_set_indicators2( BranchModel *bm, const unsigned int *positions );


void LocalClock_get_indexes( const BranchModel *bm,  unsigned *indexes);

#pragma mark -
#pragma mark DiscreteClock

BranchModel * new_DiscreteClock2( Tree *tree, const int n );

BranchModel * new_DiscreteClock( Tree *tree, const int n );

BranchModel * new_DiscreteClock_with_parameters( Tree *tree, const Parameters *rates );

BranchModel * new_DiscreteClock_from_tree( Tree *tree );

BranchModel * new_DiscreteClock_from_LocalClock_tree( Tree *tree );

BranchModel * new_DiscreteClock_from_LocalClock( const BranchModel *localBm );

void DiscreteClock_set_number_of_rate_classes( BranchModel *bm, const int nClasses );

void DiscreteClock_set_classes( BranchModel *bm, const unsigned int *classes );

void DiscreteClock_set_random_branch_assigment( BranchModel *bm );


#pragma mark -
#pragma mark RelaxedClock

BranchModel * new_RelaxedClock( Tree *tree, const relaxed_clock type, const int n, ... );

BranchModel * new_RelaxedClock_with_parameters( Tree *tree, const Parameters *rates, const relaxed_clock type );

BranchModel * new_RelaxedClock_from_tree( Tree *tree, double center );

BranchModel * new_LocalClockFromTree( Tree *tree );

void RelaxedClock_set_classes( BranchModel *bm, const unsigned int *classes );

void RelaxedClock_set_random_branch_assigment( BranchModel *bm );

void RelaxedClock_set_random_branch_assigment2( BranchModel *bm, const int cat_count );


#pragma mark -
#pragma mark Misc

void infer_distance_from_rate_height( BranchModel *bm );

void print_rate_map( BranchModel *bm );

void print_compare_bms( const BranchModel *bm1, const BranchModel *bm2);

void compare_branchmodel( const BranchModel *bm1, const BranchModel *bm2 );

double BranchModel_mean_rate_scaled( BranchModel *bm );

double BranchModel_mean_rate_tips_scaled( BranchModel *bm );

double BranchModel_mean_rate_internal_scaled( BranchModel *bm );

double BranchModel_mean_rate( BranchModel *bm, double *min, double *max );

double BranchModel_correlation( BranchModel *bm );

double BranchModel_correlation_distance( BranchModel *bm );

void BranchModel_check_outliers( BranchModel * bm );

void BranchModel_to_distance( BranchModel *bm );

#endif
