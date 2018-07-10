/*
 *  topologyopt.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/08/13.
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
