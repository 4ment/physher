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
#include <strings.h>

#include "matrix.h"
#include "tree.h"
#include "optimize.h"
#include "treeio.h"
#include "parsimony.h"
#include "spropt.h"
#include "nniopt.h"
#include "compoundmodel.h"


//#define DEBUG_TOPOLOGY2 1
//#define DEBUG_TOPOLOGY 1
//#define DEBUG_TOPOLOGY_NNNI 1


TopologyOptimizer * new_TopologyOptimizer( Model *mmodel ){
    TopologyOptimizer *opt = (TopologyOptimizer*)malloc(sizeof(TopologyOptimizer));
    assert(opt);
    opt->lnls      = NULL;
    opt->branches  = NULL;
    opt->positions = NULL;
    opt->best_lnl  = -INFINITY;
    opt->model = mmodel;
	opt->tlk = NULL;
    opt->moves = 0;
    opt->K = 0.75;
    opt->threads = 1;
    
    return opt;
}

void free_TopologyOptimizer( TopologyOptimizer *opt ){
    
    if( opt->lnls != NULL )free(opt->lnls);
    if( opt->branches != NULL )free(opt->branches);
    if( opt->positions != NULL ) free(opt->positions);
	if(opt->tlk != NULL)opt->tlk->free(opt->tlk);
	opt->model->free(opt->model);
    free(opt);
}

void TopologyOptimizer_set_algorithm( TopologyOptimizer *opt, tree_search_algorithm algorithm ){
    int size = Tree_tip_count(opt->tree);
    
    switch ( algorithm ) {
        case TREE_SEARCH_NNI:{
            opt->K = 0.75;
			SingleTreeLikelihood* tlk = opt->tlk->obj;
            if ( tlk->bm == NULL ){
                opt->optimize = nni_optimize_bl;
            }
            else {
                opt->optimize = nni_optimize_heights;
            }
            break;
        }
        case TREE_SEARCH_NNNI:{
            opt->K = 0.75;
			SingleTreeLikelihood* tlk = opt->tlk->obj;
            if ( tlk->bm == NULL ){
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
            opt->optimize = spr_optimize_bl_parsimony_only;
            break;
        }
        case TREE_SEARCH_SPR:{
			SingleTreeLikelihood* tlk = opt->tlk->obj;
            size = Tree_node_count(tlk->tree);
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


TopologyOptimizer* new_TopologyOptimizer_from_json(json_node* node, Hashtable* hash){
	char* algorithm_string = get_json_node_value_string(node, "move");
	char* criterion = get_json_node_value_string(node, "criterion");
	int nthreads = get_json_node_value_int(node, "threads", 1);
	int verbosity = get_json_node_value_int(node, "verbosity", 1);
	//	int maxiter = get_json_node_value_int(node, "maxiter", 100);
	json_node* tlk_node = get_json_node(node, "treelikelihood"); // treelikelihood with the tree
	json_node* model_node = get_json_node(node, "model"); // model to optimize
	
	tree_search_algorithm algorithm = TREE_SEARCH_NNI;
	if (strcasecmp(algorithm_string, "nnni") == 0) {
		algorithm = TREE_SEARCH_NNNI;
	}
	else if(strcasecmp(algorithm_string, "spr") == 0 && criterion != NULL && strcasecmp(criterion, "parsimony") == 0){
		algorithm = TREE_SEARCH_PARSIMONY_SPR;
	}
	else if(strcasecmp(algorithm_string, "spr") == 0){
		algorithm = TREE_SEARCH_SPR;
	}
	
	Model* likelihood = NULL;
	if(tlk_node != NULL){
		if (tlk_node->node_type == MJSON_OBJECT) {
			likelihood = new_TreeLikelihoodModel_from_json(tlk_node, hash);
			char* id = get_json_node_value_string(tlk_node, "id");
			Hashtable_add(hash, id, likelihood);
		}
		else if(tlk_node->node_type == MJSON_STRING){
			char* ref = (char*)tlk_node->value;
			likelihood = Hashtable_get(hash, ref+1);
			likelihood->ref_count++;
		}
		else{
			exit(10);
		}
	}
	
	Model* model = NULL; // full model (could be treelikelihood)
	if (model_node->node_type == MJSON_OBJECT) {
		json_node* type_node = get_json_node(model_node, "type");
		char* id = get_json_node_value_string(model_node, "id");
		
		if (strcasecmp((char*)type_node->value, "compound") == 0) {
			model = new_CompoundModel_from_json(model_node, hash);
		}
		else if(strcasecmp((char*)type_node->value, "treelikelihood") == 0){
			model = new_TreeLikelihoodModel_from_json(model_node, hash);
		}
		else if(strcasecmp((char*)type_node->value, "parsimony") == 0){
			model = new_ParsimonyModel_from_json(model_node, hash);
		}
		else{
			exit(10);
		}
		Hashtable_add(hash, id, model);
	}
	else if(model_node->node_type == MJSON_STRING){
		char* ref = (char*)model_node->value;
		model = Hashtable_get(hash, ref+1);
		model->ref_count++;
	}
	else{
		exit(10);
	}
	
	TopologyOptimizer *opt = new_TopologyOptimizer(model);
	opt->tlk = likelihood;
	opt->verbosity = verbosity;
	opt->max_distance = get_json_node_value_int(node, "radius", 20);
	
	if (likelihood == NULL) {
		opt->tree = ((Parsimony*)(model->obj))->tree;
	}
	else{
		opt->tree = ((SingleTreeLikelihood*)(likelihood->obj))->tree;
	}
	
	TopologyOptimizer_set_algorithm(opt, algorithm);
	
	return opt;
}
