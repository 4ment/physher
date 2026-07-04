// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later


#ifndef _BRANCH_MODEL_H_
#define _BRANCH_MODEL_H_

#include "parameters.h"
#include "utils.h"
#include "tree.h"
#include "node.h"
#include "discreteparameter.h"

#define BRANCHMODEL_RATE_MIN 1e-15
#define BRANCHMODEL_RATE_MAX 0.1

typedef enum branchmodel { NO_CLOCK, CLOCK_STRICT, CLOCK_LOCAL, CLOCK_DISCRETE, CLOCK_RELAXED, CLOCK_ARBITRARY } branchmodel;

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
	
	// SSVS
	DiscreteParameter* ssvs;
	unsigned* ssvs_map; // discrete index to node index
	unsigned* ssvs_map2; // node index to discrete index
} BranchModel;

#pragma mark -
#pragma mark BranchModel

Model* new_BranchModel_from_json(json_node*node, Hashtable*hash);

Model * new_BranchModel2( const char* name, BranchModel *bm, Model* tree, Model* ssvs);

BranchModel * clone_BranchModel(const BranchModel *bm, Tree *tree, DiscreteParameter* dp );

void BranchModel_vector_to_rates( BranchModel *bm, const double *rates );

void BranchModel_backward(BranchModel *bm, Parameters* parameters, const double* ingrad);

#pragma mark -
#pragma mark NoClock

BranchModel * new_NoClock( Tree *tree );

#pragma mark -
#pragma mark StrickClock

BranchModel * new_StrictClock( Tree *tree );

BranchModel * new_StrictClock_with_parameter( Tree *tree, Parameter *rate );

#pragma mark -
#pragma mark ArbitraryClock

BranchModel * new_ArbitraryClock_with_parameters( Tree *tree, Parameter *rates, Parameter *location, Parameter *scale );

#pragma mark -
#pragma mark LocalClock

BranchModel * new_LocalClock( Tree *tree, const int nLocalClocks );

BranchModel * new_LocalClock_with_parameter( Tree *tree, Parameter *rates );

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

BranchModel *new_DiscreteClock_with_parameters(Tree *tree, Parameter *rates,
                                               DiscreteParameter *map);

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

double BranchModel_mean_rate_scaled( BranchModel *bm );

double BranchModel_mean_rate_tips_scaled( BranchModel *bm );

double BranchModel_mean_rate_internal_scaled( BranchModel *bm );

double BranchModel_mean_rate( BranchModel *bm, double *min, double *max );

double BranchModel_correlation( BranchModel *bm );

double BranchModel_correlation_distance( BranchModel *bm );

void BranchModel_check_outliers( BranchModel * bm );

void BranchModel_to_distance( BranchModel *bm );

#endif
