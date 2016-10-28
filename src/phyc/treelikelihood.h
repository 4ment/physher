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

#include "utils.h"
#include "tree.h"
#include "node.h"
#include "sitepattern.h"
#include "substmodel.h"
#include "sitemodel.h"
#include "branchmodel.h"
#include "optimizer.h"
#include "treesearch.h"

//#define SCALING_THRESHOLD 1.0e-100


//#define SSE3_ENABLED 1

typedef enum treelikelihood_approximation {
    TREELIKELIHOOD_APPROXIMATION_NONE,
    TREELIKELIHOOD_APPROXIMATION_HESSIAN_DIAGONAL,
    TREELIKELIHOOD_APPROXIMATION_HESSIAN,
    TREELIKELIHOOD_APPROXIMATION_DERIVATIVE_DIAGONAL
}treelikelihood_approximation;

// 32 bytes
typedef struct Optimizable{
	double tolx;
	double tolfx;
	opt_algorithm method;
	int max_iteration;
	bool optimize;

	
} Optimizable;

// 244 bytes
typedef struct OptConfig{
	Optimizable bl;
	Optimizable freqs;
	Optimizable relative_rates;
	Optimizable gamma;
	Optimizable pinv;
	Optimizable rates;
	Optimizable heights;
    tree_search_algorithm topology_alogrithm;
    int topology_threads;
    bool topology_optimize;
	double precision;
	int max_rounds;
	int verbosity;
	bool interruptible;
}OptConfig;

struct _SingleTreeLikelihood;

typedef struct _SingleTreeLikelihood SingleTreeLikelihood;

typedef double (*calculate_upper_t)( SingleTreeLikelihood *, Node *);

struct _SingleTreeLikelihood{
	OptConfig opt;
	
	Tree        *tree;
	SiteModel   *sm;
	SitePattern *sp;
	BranchModel *bm;
	
	int id;
	int *mapping; // index of node to index of sequence (-1 if no sequence)
	
	int partials_size;
	double **partials;
	
	double **scaling_factors;
	bool scale;
	double scaling_threshold;
	
	int matrix_size;
	double **matrices;
	
	
	bool *update_nodes;
	bool update;
	
	double *root_partials;
	double *pattern_lk;
	double lk;
	
	
	double (*calculate)( SingleTreeLikelihood *);
	void (*update_partials)( SingleTreeLikelihood *, int, int, int );
	void (*integrate_partials)( const SingleTreeLikelihood *, const double *, const double *, double * );
	void (*node_log_likelihoods)( const SingleTreeLikelihood *, const double *, const double *, double * );

    treelikelihood_approximation approx;
	double *hessian; // used for the Taylor approximation [nNodes x nNodes]
    int hessian_length;
    double lnl_bl;

#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	bool use_SIMD;
#endif

	int node_id; // use to optimize subtrees
    
    int nthreads;

    // Calculation using lower and upper
    bool use_upper;
    double **partials_upper;
	double *node_pattern_lk;
    //int index_upper;
    Node *node_upper;
    
	void   (*update_partials_upper)( SingleTreeLikelihood *, Node * );
    calculate_upper_t calculate_upper;
	void   (*node_log_likelihoods_upper)( const SingleTreeLikelihood *, Node * );
    bool use_ambiguity;
};


#pragma mark -
#pragma mark Optimizable

void Optimizable_init( Optimizable *opt, bool optimize, opt_algorithm method, const int max_iteration, const double tolx, const double tolfx );

void Optimizable_copy( const Optimizable *src, Optimizable *dst );

void OptConfig_init( SingleTreeLikelihood *tlk );

void OptConfig_copy( const OptConfig *src, OptConfig *dst );


#pragma mark -
#pragma mark SingleTreeLikelihood


SingleTreeLikelihood * new_SingleTreeLikelihood( Tree *tree, SiteModel *sm, SitePattern *sp, BranchModel *bm );

void free_SingleTreeLikelihood( SingleTreeLikelihood *tlk );

void free_SingleTreeLikelihood_share2( SingleTreeLikelihood *tlk, bool shared_tree, bool shared_sitemodel, bool shared_sitepattern, bool shared_branchmodel );

void free_SingleTreeLikelihood_share( SingleTreeLikelihood *tlk, bool shared_sitepattern, bool shared_sitemodel );

void SingleTreeLikelihood_free_upper( SingleTreeLikelihood *tlk );

SingleTreeLikelihood * clone_SingleTreeLikelihood( SingleTreeLikelihood *tlk );

SingleTreeLikelihood * clone_SingleTreeLikelihood_share( SingleTreeLikelihood *tlk, bool share_sitepattern, bool share_sitemodel );

SingleTreeLikelihood * clone_SingleTreeLikelihood_with( SingleTreeLikelihood *tlk, SitePattern *sp, SiteModel *sm );


void SingleTreeLikelihood_use_ambiguity(SingleTreeLikelihood *tlk);

void SingleTreeLikelihood_set_BranchModel( SingleTreeLikelihood *stlk, BranchModel *bm, bool removeTree );

void SingleTreeLikelihood_remove_BranchModel( SingleTreeLikelihood *tlk, bool removeTree );

void SingleTreeLikelihood_set_Tree( SingleTreeLikelihood *stlk, Tree *tree );

void SingleTreeLikelihood_remove_Tree( SingleTreeLikelihood *stlk );


double SingleTreeLikelihood_calculate_at_node( SingleTreeLikelihood *tlk, const Node *node );


void SingleTreeLikelihood_copy_partials( SingleTreeLikelihood *src, SingleTreeLikelihood *dst );

void SingleTreeLikelihood_rearrange_partials( SingleTreeLikelihood *tlk );
void SingleTreeLikelihood_rearrange( SingleTreeLikelihood *tlk, Node *node1, Node *node2 );

// using upper partial likelihood
// should notbe called at the root
double SingleTreeLikelilhood_calculate_uppper( SingleTreeLikelihood *tlk, int node_index );
//void calculate_uppper_partials( SingleTreeLikelihood *tlk );
double SingleTreeLikelilhood_calculate_uppper2( SingleTreeLikelihood *tlk, Node *node );


void TreeLikelihood_set_calculate( SingleTreeLikelihood *tlk, Node *node );

void update_branch_length( SingleTreeLikelihood *tlk, Parameters *ps, int index );

void update_branches( SingleTreeLikelihood *tlk, Parameters *ps );

void SingleTreeLikelihood_update_all_nodes( SingleTreeLikelihood *tlk );

//void SingleTreeLikelihood_update_one_node( SingleTreeLikelihood *tlk, const int index );
void SingleTreeLikelihood_update_one_node( SingleTreeLikelihood *tlk, const Node *node );

void SingleTreeLikelihood_update_three_nodes( SingleTreeLikelihood *tlk, const Node *node );


void compare_singlelikelihood( const SingleTreeLikelihood *tlk1, const SingleTreeLikelihood *tlk2 );

void SingleTreeLikelihood_use_rescaling( SingleTreeLikelihood *tlk, bool use );

int SingleTreeLikelihood_df_count( const SingleTreeLikelihood *stlk );

double getLogScalingFactor( const SingleTreeLikelihood *tlk, int pattern );

void SingleTreeLikelihood_scalePartials( SingleTreeLikelihood *tlk, int nodeIndex );

double ** SingleTreeLikelihood_posterior_sites( SingleTreeLikelihood *tlk );

void SingleTreeLikelihood_enable_SSE( SingleTreeLikelihood *tlk, bool value );
bool SingleTreeLikelihood_SSE( SingleTreeLikelihood *tlk );


void SingleTreeLikelihood_add_height( SingleTreeLikelihood *tlk, Node *node, double value );
void SingleTreeLikelihood_scaler( SingleTreeLikelihood *tlk, Node *node, const double scaler );
void SingleTreeLikelihood_scale_root( SingleTreeLikelihood *tlk, const double scaler );

void SingleTreeLikelihood_set_nthreads( SingleTreeLikelihood *tlk, int nthreads );

void calculate_all_dlnl_dt( SingleTreeLikelihood *tlk, double *dlnl );

void calculate_all_d2lnl_d2t( SingleTreeLikelihood *tlk, double *d2lnl );

void calculate_hessian_branches( SingleTreeLikelihood *tlk, double *d2lnl );

#pragma mark -
#pragma mark Second order Taylor series approximation

void SingleTreeLikelihood_init_approximation( SingleTreeLikelihood *tlk, const double *a );

void SingleTreeLikelihood_fisher_information( SingleTreeLikelihood *tlk, double *hessian );

void SingleTreeLikelihood_covariance( SingleTreeLikelihood *tlk, double *hessian );

bool SingleTreeLikelihood_Hessian( SingleTreeLikelihood *tlk, double *hessian, double *gradient );

bool SingleTreeLikelihood_Hessian_diag( SingleTreeLikelihood *tlk, double *hessian, int *len );

double SingleTreeLikelihood_d2_distance( SingleTreeLikelihood *tlk, Node *node );

void SingleTreeLikelihood_dt( SingleTreeLikelihood *tlk, double *dlnls );

#pragma mark -
#pragma mark Upper Likelihood

void SingleTreeLikelihood_set_upper_function( SingleTreeLikelihood *tlk, calculate_upper_t function );

void SingleTreeLikelihood_use_upper( SingleTreeLikelihood *tlk, bool use_upper );

double calculate_uppper_2nodes( SingleTreeLikelihood *tlk, Node *node );


/*char * SingleTreeLikelihood_stringify( SingleTreeLikelihood *stlk );
 
 StringBuffer * SingleTreeLikelihood_bufferize( StringBuffer *buffer, SingleTreeLikelihood *stlk );
 
 void * SingleTreeLikelihood_SML_to_object( ObjectStore *store, SMLNode node );*/

#endif
