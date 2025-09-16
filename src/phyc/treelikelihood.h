/*
 *  treelikelihood.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/2/10.
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


#ifndef _TREELIKELIHOOD_H_
#define _TREELIKELIHOOD_H_

#include "tree.h"
#include "node.h"
#include "sitepattern.h"
#include "substmodel.h"
#include "sitemodel.h"
#include "branchmodel.h"
#include "mjson.h"

static double TWENTY_DOUBLE_ONES[20] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

#define TREELIKELIHOOD_FLAG_TREE_MODEL         1 << 0
#define TREELIKELIHOOD_FLAG_SITE_MODEL         1 << 1
#define TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL 1 << 2
#define TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_RATES 1 << 3
#define TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_FREQUENCIES 1 << 4
#define TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_UNCONSTRAINED 1 << 5
#define TREELIKELIHOOD_FLAG_BRANCH_MODEL       1 << 6

struct _SingleTreeLikelihood;

typedef struct _SingleTreeLikelihood SingleTreeLikelihood;

typedef double (*calculate_upper_t)( SingleTreeLikelihood *, Node *);

struct _SingleTreeLikelihood{
	
	Tree        *tree;
	SubstitutionModel *m;
	SiteModel   *sm;
	SitePattern *sp;
	BranchModel *bm;
	
	int cat_count;
	int pattern_count;
	bool use_tip_states;
	
	int *mapping; // index of node to index of sequence (-1 if no sequence)
	
	int partials_dim; // first dimension of partials
	int partials_size; // second dimension of partials
	double ***partials;
	unsigned *current_partials_indexes;
	unsigned *stored_partials_indexes;
	
	double ***scaling_factors;

	bool scale;
	double scaling_threshold;
	
	int matrix_dim;
	int matrix_size;
	double ***matrices;
	unsigned *current_matrices_indexes;
	unsigned *stored_matrices_indexes;

	
	
	bool *update_nodes;
	bool update;
	
	int root_partials_size;
	int pattern_lk_size;
	double *root_partials;
	double *pattern_lk;
	double lk;
	double stored_lk;
	
	
	double (*calculate)( SingleTreeLikelihood *);
	void (*update_partials)( SingleTreeLikelihood *, int, int, int, int, int );
	void (*update_partials_flexible)( SingleTreeLikelihood *, double*, int, double*, double*, int, double*, double* );
	void (*integrate_partials)( const SingleTreeLikelihood *, const double *, const double *, double * );
	void (*node_log_likelihoods)( const SingleTreeLikelihood *, const double *, const double *, double * );

#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	bool use_SIMD;
#endif

	int node_id; // use to optimize subtrees
    
    int nthreads;

    // Calculation using lower and upper
    bool use_upper;
    Node *node_upper;
	bool tripod;
	
    calculate_upper_t calculate_upper;
	
	void (*calculate_per_cat_partials)(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);
	
	int* upper_partial_indexes;
	bool update_upper; // for derivatives
	
	const double* (*get_root_frequencies)(SingleTreeLikelihood*);
	double* root_frequencies;
	bool include_jacobian;
	
	int prepared_gradient;
	double* gradient;
	size_t gradient_length;
	bool include_root_freqs;
};

Model * new_TreeLikelihoodModel( const char* name, SingleTreeLikelihood *tlk, Model *tree, Model *m, Model *sm, Model *bm );

Model * new_TreeLikelihoodModel_from_json(json_node*node, Hashtable*hash);

size_t TreeLikelihood_initialize_gradient(Model *self, int flags);
void TreeLikelihood_gradient(Model *self, int flags, double* grads);


#pragma mark -
#pragma mark SingleTreeLikelihood


SingleTreeLikelihood * new_SingleTreeLikelihood( Tree *tree, SubstitutionModel *m, SiteModel *sm, SitePattern *sp, BranchModel *bm, bool use_tip_states );

// does not free tree, branchmodel, sitemodel
void free_SingleTreeLikelihood_internals( SingleTreeLikelihood *tlk );

void free_SingleTreeLikelihood( SingleTreeLikelihood *tlk );


SingleTreeLikelihood * clone_SingleTreeLikelihood( SingleTreeLikelihood *tlk );

SingleTreeLikelihood * clone_SingleTreeLikelihood_with( SingleTreeLikelihood *tlk, Tree *tree, SubstitutionModel *m, SiteModel *sm, SitePattern *sp, BranchModel *bm);



double SingleTreeLikelihood_calculate_at_node( SingleTreeLikelihood *tlk, const Node *node );


void SingleTreeLikelihood_update_all_nodes( SingleTreeLikelihood *tlk );

void SingleTreeLikelihood_update_one_node( SingleTreeLikelihood *tlk, const Node *node );

void SingleTreeLikelihood_update_three_nodes( SingleTreeLikelihood *tlk, const Node *node );

void SingleTreeLikelihood_update_Q(SingleTreeLikelihood* tlk, Node* n);

bool SingleTreeLikelihood_rescaling( SingleTreeLikelihood *tlk );
void SingleTreeLikelihood_use_rescaling( SingleTreeLikelihood *tlk, bool use );

double getLogScalingFactor( const SingleTreeLikelihood *tlk, int pattern );

void SingleTreeLikelihood_scalePartials( SingleTreeLikelihood *tlk, int nodeIndex, int childIndex1, int childIndex2 );

void SingleTreeLikelihood_enable_SSE( SingleTreeLikelihood *tlk, bool value );
bool SingleTreeLikelihood_SSE( SingleTreeLikelihood *tlk );

const double* get_root_frequencies(SingleTreeLikelihood* tlk);

const double* get_root_frequencies_fixed(SingleTreeLikelihood* tlk);


void SingleTreeLikelihood_set_nthreads( SingleTreeLikelihood *tlk, int nthreads );

double calculate_dlnl_dQ( SingleTreeLikelihood *tlk, int index, const double* pattern_likelihoods );

double calculate_dlnl_dWeibull( SingleTreeLikelihood *tlk, const double* pattern_likelihoods );


#pragma mark -
#pragma mark Upper Likelihood

void update_upper_partials(SingleTreeLikelihood *tlk, Node* node, bool update_all_partials);
void update_upper_partials2(SingleTreeLikelihood *tlk, Node* node);

void calculate_dldt_uppper( SingleTreeLikelihood *tlk, Node *node, double* pattern_dlikelihoods );

void calculate_dldh_uppper( SingleTreeLikelihood *tlk, Node *node, double* pattern_dlikelihoods );

double dlnldt_uppper( SingleTreeLikelihood *tlk, const double* pattern_likelihoods, const double* pattern_dlikelihoods );

double d2lnldt2_uppper( SingleTreeLikelihood *tlk, Node *node, const double* pattern_likelihoods, const double* pattern_dlikelihoods);

void SingleTreeLikelihood_update_uppers(SingleTreeLikelihood *tlk);
void SingleTreeLikelihood_update_uppers2(SingleTreeLikelihood *tlk);

#endif
