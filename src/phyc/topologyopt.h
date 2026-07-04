// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_topologyopt_h
#define PhyC_topologyopt_h

#include "treelikelihood.h"
#include "treesearch.h"



typedef struct TopologyOptimizer{
	Model* model;
	Model *tlk;
	Tree* tree;
    double best_lnl;
    double *lnls;
    double *branches;
    int    *positions;
    tree_search_algorithm algorithm;
    
    double (*optimize)( struct TopologyOptimizer * );
	int max_distance; // for SPR
	int max_failures;
	int failures;
    int moves;
    double K;
    int threads;
	int verbosity;
} TopologyOptimizer;

TopologyOptimizer * new_TopologyOptimizer( Model *tlk );

void free_TopologyOptimizer( TopologyOptimizer *opt );

void TopologyOptimizer_set_algorithm( TopologyOptimizer *opt, tree_search_algorithm algorithm );

void TopologyOptimizer_set_nthreads( TopologyOptimizer *opt, int nthreads );

TopologyOptimizer* new_TopologyOptimizer_from_json(json_node* node, Hashtable* hash);

#endif
