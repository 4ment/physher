/*
 *  topologyopt.c
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

#include "topologyopt.h"

#include <stdio.h>
#include <assert.h>

#include "matrix.h"
#include "tree.h"
#include "optimize.h"
#include "treeio.h"
#include "parsimony.h"
#include "spropt.h"
#include "nniopt.h"


//#define DEBUG_TOPOLOGY2 1
//#define DEBUG_TOPOLOGY 1
//#define DEBUG_TOPOLOGY_NNNI 1


TopologyOptimizer * new_TopologyOptimizer( SingleTreeLikelihood *tlk, tree_search_algorithm algorithm ){
    TopologyOptimizer *opt = (TopologyOptimizer*)malloc(sizeof(TopologyOptimizer));
    assert(opt);
    opt->lnls      = NULL;
    opt->branches  = NULL;
    opt->positions = NULL;
    opt->best_lnl  = -INFINITY;
    opt->algorithm = algorithm;
    opt->tlk = tlk;
    opt->moves = 0;
    opt->K = 0.75;
    opt->threads = 1;

    switch ( algorithm ) {
        case TREE_SEARCH_NNI:{
            if ( tlk->bm == NULL ){
                opt->optimize = nni_optimize_bl;
            }
            else {
                opt->optimize = nni_optimize_heights;
            }
            break;
        }
        case TREE_SEARCH_NNNI:{
            if ( tlk->bm == NULL ){
                opt->optimize = nnni_optimize_bl;
            }
            break;
        }
        case TREE_SEARCH_PARSIMONY_NNI:{
            opt->optimize = nni_optimize_bl_parsimony;
            
            break;
        }
        case TREE_SEARCH_PARSIMONY_SPR:{
            opt->optimize = spr_optimize_bl_parsimony_only;
            
            break;
        }
        case TREE_SEARCH_SPR:{
            opt->optimize = spr_optimize_bl_openmp;
            
            break;
        }
        default:{
            free(opt);
            return NULL;
            break;
        }
    }
    
    TopologyOptimizer_set_algorithm(opt, algorithm);
    
    return opt;
}

void free_TopologyOptimizer( TopologyOptimizer *opt ){
    
    if( opt->lnls != NULL )free(opt->lnls);
    if( opt->branches != NULL )free(opt->branches);
    if( opt->positions != NULL ) free(opt->positions);
    free(opt);
}

void TopologyOptimizer_set_algorithm( TopologyOptimizer *opt, tree_search_algorithm algorithm ){
    
    int size = Tree_tip_count(opt->tlk->tree);
    
    switch ( algorithm ) {
        case TREE_SEARCH_NNI:{
            opt->K = 0.75;
            if ( opt->tlk->bm == NULL ){
                opt->optimize = nni_optimize_bl;
            }
            else {
                opt->optimize = nni_optimize_heights;
            }
            break;
        }
        case TREE_SEARCH_NNNI:{
            opt->K = 0.75;
            if ( opt->tlk->bm == NULL ){
                opt->optimize = nnni_optimize_bl;
            }
            break;
        }
        case TREE_SEARCH_PARSIMONY_NNI:{
            opt->K = 0.75;
            opt->optimize = nni_optimize_bl_parsimony;
            
            break;
        }
        case TREE_SEARCH_PARSIMONY_SPR:{
            size = Tree_node_count(opt->tlk->tree);
            opt->optimize = spr_optimize_bl_parsimony_only;
            break;
        }
        case TREE_SEARCH_SPR:{
            size = Tree_node_count(opt->tlk->tree);
            opt->optimize = spr_optimize_bl_openmp;
            break;
        }
        default:{
            opt->K = 0;
            break;
        }
    }
    
    if( opt->lnls == NULL ){
        opt->lnls      = dvector(size);
        opt->branches  = dvector(size);
        opt->positions = ivector(size);
    }
    else if( opt->algorithm != algorithm ){
        opt->lnls = realloc(opt->lnls, size*sizeof(double));
        opt->branches = realloc(opt->branches, size*sizeof(double));
        opt->positions = realloc(opt->positions, size*sizeof(int));
    }
    
    opt->algorithm = algorithm;
}

void TopologyOptimizer_set_nthreads( TopologyOptimizer *opt, int nthreads ){
#if defined (_OPENMP)
    opt->threads = nthreads;
#endif
}


