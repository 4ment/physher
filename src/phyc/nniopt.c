/*
 *  nniopt.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/04/2014.
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

#include "nniopt.h"

#include "topologyopt.h"
#include "optimize.h"
#include "matrix.h"
#include "parsimony.h"

#include "pooledtreelikelihood.h"

#include "treeio.h"


static double _brent_optimize_height( Parameters *params, double *grad, void *data ){
	BrentData *opt = (BrentData*)data;
    Node *node = Tree_get_node(opt->tlk->tree, POSTORDER, opt->index_param);
    Node_set_height(node, Parameters_value(params, 0));
    
	SingleTreeLikelihood_update_three_nodes(opt->tlk, node );
#ifdef USE_UPPER_PARTIALS
    double lnl = fabs(opt->tlk->calculate_upper(opt->tlk, Tree_get_node(opt->tlk->tree, POSTORDER, opt->index_param)));
    //printf("%d] %f\n",opt->index_param,lnl);
    return lnl;
#else
	return fabs(opt->tlk->calculate(opt->tlk));
#endif
}

static void _sort_decreasing_double( double *values, int *order, int size ){
    
    for ( int i = 0; i < size; i++ ) {
        order[i] = i;
    }
    
    bool done = false;
    while ( !done ) {
        done = true;
        for ( int i = 0 ; i < size-1 ; i++ ) {
            if ( values[ order[i] ] < values[ order[i+1] ] ) {
                done = false;
                swap_int(&order[i], &order[i+1]);
                dswap(&values[ order[i] ], &values[ order[i+1] ]);
            }
        }
        size--;
    }
}

static void _sort_decreasing_lnls( /*int *v,*/ TopologyOptimizer *opt, int *map, int size ){
    bool done = false;
    while ( !done ) {
        done = true;
        for ( int i = 0 ; i < size-1 ; i++ ) {
            if ( opt->lnls[i] < opt->lnls[i+1] ) {
                done = false;
                //swap_int(&v[i], &v[i+1]);
                swap_int(&(opt->positions[i]), &(opt->positions[i+1]));
                dswap(&(opt->lnls[i]), &(opt->lnls[i+1]));
                dswap(&(opt->branches[i]), &(opt->branches[i+1]));
                swap_int(&map[i], &map[i+1]);
            }
        }
        size--;
    }
}

// descreasing nni
static void _sort_by_nni( /*int *v,*/ TopologyOptimizer *opt, int *map, int size ){
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( opt->positions[i] < opt->positions[i+1] ) {
				done = false;
                //swap_int(&v[i], &v[i+1]);
                swap_int(&(opt->positions[i]), &(opt->positions[i+1]));
                dswap(&(opt->lnls[i]), &(opt->lnls[i+1]));
                dswap(&(opt->branches[i]), &(opt->branches[i+1]));
                swap_int(&map[i], &map[i+1]);
			}
		}
		size--;
	}
}


// remove adjacent good NNIs
void remove_adjacent_nnis( TopologyOptimizer *opt, Node *node, int *map ){
    
    if( Node_isleaf(node)){
        return;
    }
    
    remove_adjacent_nnis(opt, node->left, map);
    remove_adjacent_nnis(opt, node->right, map);
    
    if( map[ Node_id(node) ] == -1 ) return;
    
    if( opt->positions[ map[ Node_id(node) ] ] != 0 ){
        // consider at least left if it is not a tip
        if ( !Node_isleaf(Node_left(node)) && opt->positions[ map[  Node_id(Node_left(node)) ] ] != 0 ) {
            
            // node better than left
            if( opt->lnls[ map[ Node_id(node) ] ] >= opt->lnls[ map[  Node_id(Node_left(node)) ] ] ){
                opt->positions[ map[  Node_id(Node_left(node)) ] ] = 0;
                
                // no tips under node
                if ( !Node_isleaf(Node_right(node)) && opt->positions[ map[  Node_id(Node_right(node)) ] ] != 0 ){
                    // node better than right
                    if( opt->lnls[ map[ Node_id(node) ] ] >= opt->lnls[ map[ Node_id(Node_right(node)) ] ] ){
                        opt->positions[ map[ Node_id(Node_right(node)) ] ] = 0;
                    }
                    // right better than node
                    else {
                        opt->positions[ map[  Node_id(node) ] ] = 0;
                    }
                }
            }
            // left better than node
            else {
                opt->positions[ map[  Node_id(node) ] ] = 0;
                // consider right if is not a tip
                if ( !Node_isleaf(Node_right(node)) && opt->positions[ map[  Node_id(Node_right(node)) ] ] != 0 ){
                    // left better than right
                    if( opt->lnls[  Node_id(Node_left(node)) ] >= opt->lnls[ map[  Node_id(Node_right(node)) ] ] ){
                        opt->positions[ map[  Node_id(Node_right(node)) ] ] = 0;
                    }
                    // right better than left
                    else {
                        opt->positions[ map[  Node_id(Node_left(node)) ] ] = 0;
                    }
                }
            }
        }
        // consider only right
        else if ( !Node_isleaf(Node_right(node)) && opt->positions[ map[  Node_id(Node_right(node)) ] ] != 0 ) {
            
            // node better than right
            if( opt->lnls[ map[ Node_id(node) ] ] >= opt->lnls[ map[ Node_id(Node_right(node)) ] ] ){
                opt->positions[ map[  Node_id(Node_right(node)) ] ] = 0;
            }
            // right better than node
            else {
                opt->positions[ map[ Node_id(node) ] ] = 0;
            }
            
        }
    }
}

//#define DEBUG_TOPOLOGY 2

// there is 2n-6 NNIs
double nni_optimize_bl( struct TopologyOptimizer * opt ){
	double lnl = opt->model->logP(opt->model);
	
    if(opt->failures >= opt->max_failures) return lnl;
	
    unsigned nthreads = opt->threads;
    assert(nthreads > 0);
	
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	Model** pool = malloc(sizeof(Model*)*nthreads);
	Model** tlks = malloc(sizeof(Model*)*nthreads);
	pool[0] = opt->model;
	tlks[0] = opt->tlk;
	for (size_t i = 1; i < nthreads; i++) {
		pool[i] = opt->model->clone(opt->model, hash);
		tlks[i] = Hashtable_get(hash, opt->tlk->name); // not ref++
		Hashtable_empty(hash);
	}
	free_Hashtable(hash);
	
    Optimizer **opt_bls = (Optimizer**)malloc(nthreads*sizeof(Optimizer*));
    assert(opt_bls);
	Parameters **oneparameters = (Parameters**)malloc(nthreads*sizeof(Parameters*));
    assert(oneparameters);
	
	for ( int i = 0; i < nthreads; i++ ) {
        opt_bls[i] = new_Optimizer(OPT_BRENT);
		
		oneparameters[i] = new_Parameters(1);
		opt_set_data(opt_bls[i], pool[i]);
		opt_set_objective_function(opt_bls[i], model_negative_logP);
		opt_set_treelikelihood(opt_bls[i], tlks[i]->obj);
//		opt_set_max_iteration(opt, max);
//		opt_set_tolx(opt, precision);
    }
	
	Optimizer* full_opt = new_Optimizer(OPT_SERIAL_BRENT);
	opt_set_treelikelihood(full_opt, opt->tlk);
	opt_set_data(full_opt, opt->model);
	opt_set_objective_function(full_opt, model_negative_logP);
	
    int nNodes = Tree_node_count(((SingleTreeLikelihood*)(opt->tlk->obj))->tree);
    
    int *map = ivector(nNodes);
    for( int i = 0; i < nNodes; i++ ){
        map[i] = -1;
    }
    
    int *map_index_to_post = ivector(nNodes);
    int index = 0;
    bool failed = false;
    opt_result status;
    
    double *branches = dvector(nNodes); // backup branch length in case we need to backtrack
    Tree_branch_length_to_vector(((SingleTreeLikelihood*)(opt->tlk->obj))->tree, branches);
    
    time_t start_time, end_time;
    double diff_time;
    
	if (opt->verbosity > 0) {
    	printf("\nStart NNI optimization LnL: %f\n\n", lnl);
	}
	
    opt->moves = 0;
    
    while ( !failed ) {

        index = 0;
        int count = 0;
        
        time(&start_time);
        
        // do not do NNIs on the root and one of its children
        // choose NNIs than do not change the 'fake' rooting, except for the first child of the root
        #pragma omp parallel for num_threads(nthreads)
        for ( int i = 0; i < nNodes; i++ ) {
            
            int tid = 0;
#if defined (_OPENMP)
            tid = omp_get_thread_num();
#endif
            SingleTreeLikelihood *tlk = tlks[tid]->obj;
            Node *node = Tree_node(tlk->tree, i);
            
            if ( Node_isroot(node) || ( Node_isroot(Node_parent(node)) && node == Node_right(Tree_root(tlk->tree))) ) continue;
            
            
            Optimizer *opt_bl = opt_bls[tid];
            Parameters *oneparameter = oneparameters[tid];
            
            if ( !Node_isleaf(node) ) {
                
                // Do not do NNI on a topology where the sibling is a leaf (d) and its parent is the root
                //        - a
                //       |
                //      -
                //   * | |
                //   --|  - b
                //  |  |
                //  |   - c
                //   - d
                if( Node_isroot(Node_parent(node)) && Node_isleaf(Node_sibling(node)) ){
                    continue;
                }
                
                Node *left    = Node_left(node);
                Node *right   = Node_right(node);
                Node *sibling = Node_sibling(node);
                
                // Use the left node of the sibling (d)
                // NNI 1: Exchange (a,b) with d
                // NNI 2: Exchange c with d
                //        - a
                //       |
                //      -
                //   * | |
                //   --|  - b
                //  |  |
                //  |   - c
                //  |  - d
                //  | |
                //   -
                //    |
                //     - e
                if( Node_isroot(Node_parent(node)) ){
                    sibling = Node_left(sibling);
                }
                
#ifdef DEBUG_TOPOLOGY2
                printf("\nNNI around %s (%d)\n", Node_name(node), node->id);
                printf("  %f\n", lnl);
                
                printf("  rearragement 1 %s %s\n", left->name, sibling->name);
                printf("  rearragement 2 %s %s\n", right->name, sibling->name);
#endif
                double original_bl = Node_distance(node);
                Parameters_add(oneparameter, node->distance);

                
                // NNI 1
				Tree_save(tlk->tree);
				NNI_move(tlk->tree, sibling, left);
				
				tlk->use_upper = false;
				pool[0]->logP(pool[0]);
				
				Node* n = node;
				Node* root = Tree_root(tlk->tree);
				while (n != root) {
					tlk->update_nodes[Node_id(n)] = true;
					n = n->parent;
				}
				SingleTreeLikelihood_update_uppers2(tlk);
				
                tlk->node_upper = node;
				
                double nni_1 = 0;
                
                status = opt_maximize( opt_bl, oneparameter, &nni_1);
                if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
                
                double bl_1 = Node_distance(node);
                Tree_revert(tlk->tree);
				
                // NNI 2
				Tree_save(tlk->tree);
                NNI_move(tlk->tree, sibling, right);
				
				tlk->use_upper = false;
				pool[0]->logP(pool[0]);
				
				n = node;
				while (n != root) {
					tlk->update_nodes[Node_id(n)] = true;
					n = n->parent;
				}
				SingleTreeLikelihood_update_uppers2(tlk);
				
				tlk->node_upper = node;
                
                
                double nni_2 = 0;
                
                status = opt_maximize( opt_bl, oneparameter, &nni_2);
                if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");

				Parameters_pop(oneparameter);
                
                double bl_2 = Node_distance(node);
				Tree_revert(tlk->tree);
				
                double local_lnl = lnl;
                int local_position = 0;
                double local_bl = original_bl;

                if ( nni_1-0.1 > lnl ) {
                    // NNI2 is the best
                    if ( nni_2 > nni_1) {
                        local_lnl      = nni_2;
                        local_position = 2;
                        local_bl = bl_2;
                    }
                    // NNI1 is the best
                    else {
                        local_lnl      = nni_1;
                        local_position = 1;
                        local_bl = bl_1;
                    }
                }
                else {
                    // NNI2 is the best
                    if ( nni_2-0.1 > lnl ) {
                        local_lnl      = nni_2;
                        local_position = 2;
                        local_bl = bl_2;
                    }
                }

                int local_index;
#pragma omp critical
{
                local_index = index++;
    
                if( local_position != 0 ){
                    count++;
                }
}
                opt->lnls[local_index]      = local_lnl;
                opt->positions[local_index] = local_position;
                opt->branches[local_index] = local_bl;
                
                map_index_to_post[local_index] = i;
                map[i] = local_index;
            }
        }
        
        time(&end_time);
        diff_time = difftime(end_time, start_time);
        
		if (opt->verbosity > 0) {
        	printf("\n%d potentials NNIs [%f]\n", count, diff_time);
		}
		
		SingleTreeLikelihood* tlk = opt->tlk->obj;
		Tree* tree = tlk->tree;

        if( count != 0 ){
            
            if(count > 1){
                remove_adjacent_nnis(opt, Tree_root(tlk->tree), map);
                
                count = 0;
                for ( int j = 0; j < Tree_tip_count(tlk->tree)-3; j++ ) {
                    if ( opt->positions[j] != 0 ) {
                        count++;
                    }
                }
            }
            
            // 0..count-1 are valid nnis but not ordererd by LnL (should be done even if there is only one nni)
            _sort_by_nni(/*order,*/ opt, map_index_to_post, Tree_tip_count(tree)-3);
            
            // 0..count-1 are valid nnis ordererd by LnL
            _sort_decreasing_lnls(/*order,*/ opt, map_index_to_post, count);
            
        
			if (opt->verbosity > 0) {
            	printf("%d potentials non adjacent NNIs\n", count);
			}
			
            
            Node **nodes = Tree_nodes(tlk->tree);
			
//            count *= opt->K;
//            count = imax(1, count); // in case we have count == 1, we would get (count * 0.75) == 0
            double nni_lnl = lnl;

#ifdef DEBUG_TOPOLOGY
            fprintf(stderr, "%d actual non adjacent NNIs\n", count);
#endif
            
            while ( count > 0 ) {
				Tree_save(tlk->tree);
                // apply NNIs
                for ( int j = 0; j < count; j++ ) {
                    
                    int node_index = map_index_to_post[j]; // NNI node
                    
#ifdef DEBUG_TOPOLOGY
                    printf("%s %f %d d: %f\n", nodes[node_index]->name, opt->lnls[j], opt->positions[j], opt->branches[j]);
#endif
                    
                    Node *sibling = Node_sibling(nodes[node_index]);
                    
                    if( Node_isroot(Node_parent(nodes[node_index])) ){
                        sibling = Node_left(sibling);
                    }
                    
                    Node *child   = NULL;
                    
                    // sibling and left
                    if( opt->positions[j] == 1 ){
                        child = Node_left(nodes[node_index]);
                    }
                    // sibling and right
                    else {
                        child = Node_right(nodes[node_index]);
                    }
                    
                    Node_swap_parents(sibling, child);
                    
                    //Node_set_distance(nodes[map_index_to_post[j]], opt->branches[map_index_to_post[j]]);
                    Node_set_distance(nodes[map_index_to_post[j]], opt->branches[j]);
                }
#ifdef DEBUG_TOPOLOGY
                for ( int j = count; j < count2; j++ ) {
                    int node_index = map_index_to_post[j]; // NNI node
                    printf("%s %f %d *\n", nodes[node_index]->name, opt->lnls[j], opt->positions[j]);
                }
#endif
                SingleTreeLikelihood_update_all_nodes(tlk);
                Tree_set_topology_changed(tlk->tree);

//                if( count == 1 ){
//                    data_brent->index_param = Node_id(nodes[map_index_to_post[0]]);
//                    Parameters_set(oneparameter, 0, nodes[map_index_to_post[0]]->distance);
//                    SingleTreeLikelihood_update_all_nodes(tlk);
//                    status = opt_maximize( opt_bl, oneparameter, &nni_lnl);
//                    if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
//                }
                
                // Optimize all branches
#ifdef DEBUG_TOPOLOGY
                time(&start_time);
#endif
				
				status = opt_maximize( full_opt, NULL, &nni_lnl);
				if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
                
				if (opt->verbosity > 0) {
					time(&end_time);
					diff_time = difftime(end_time, start_time);
		
					printf("nni LnL: %f LnL: %f (count: %d) [%f]\n", nni_lnl, lnl,count,diff_time);
				}

                if( nni_lnl > lnl || count == 1 ){
                    opt->moves += count;
                    break;
                }
                
				if (opt->verbosity > 0) {
                	printf("Backtrack...\n");
				}
				Tree_revert(tlk->tree);
				
                //count *= 0.5;
                count *= opt->K;
                count = imax(1, count); // in case we have count == 1, we would get (count * 0.75) == 0
            }
            
            Tree_branch_length_to_vector(tlk->tree, branches);
            
            // Apply the NNIs to the others
            for ( int i = 1; i < opt->threads; i++ ) {
                SingleTreeLikelihood *tlk = tlks[i]->obj;
                Node **nodes = Tree_nodes(tlk->tree);
                
                for ( int j = 0; j < count; j++ ) {
                    int node_index = map_index_to_post[j]; // NNI node
                    
                    Node *sibling = Node_sibling(nodes[node_index]);
                    
                    if( Node_isroot(Node_parent(nodes[node_index])) ){
                        sibling = Node_left(sibling);
                    }
                    
                    Node *child   = NULL;
                    
                    // sibling and left
                    if( opt->positions[j] == 1 ){
                        child = Node_left(nodes[node_index]);
                    }
                    // sibling and right
                    else {
                        child = Node_right(nodes[node_index]);
                    }
                    
                    Node_swap_parents(sibling, child);
                }
                Tree_vector_to_branch_length(tlk->tree, branches);
                Tree_set_topology_changed(tlk->tree);
            }
            
            lnl = nni_lnl;
        }
        else {
            failed = true;
		}
    }
	
	if (opt->verbosity > 0) {
    	printf("NNI round LnL: %f\n", lnl);
	}
	
    for ( int i = 0; i < opt->threads; i++ ) {
        free_Optimizer(opt_bls[i]);
        free_Parameters(oneparameters[i]);
    }
	free_Optimizer(full_opt);
	
    free(opt_bls);
    free(oneparameters);
    
    free(map);
    free(map_index_to_post);
    free(branches);
    //free(order);
	for (int i = 1; i < nthreads; i++) {
		pool[i]->free(pool[i]);
	}
	free(pool);
	free(tlks);
	
    opt->best_lnl = lnl;
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)(opt->tlk->obj);
	tlk->use_upper = false;
	if (opt->moves == 0) {
		opt->failures++;
	}
	else{
		opt->failures = 0;
	}
    return lnl;
}


//TODO: we don't need to update all lnl
double optimize_brent_branch_length_2( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, Node *i, Node *j ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = stlk->calculate(stlk);
    double lnl2 = lnl+10;
    data->index_param = -1;
    
    SingleTreeLikelihood_update_all_nodes(stlk);
    
    while ( 1 ) {
        
        //SingleTreeLikelihood_update_all_nodes(stlk);
        SingleTreeLikelihood_update_one_node(stlk, i);
        data->index_param = Node_id(i);
        Parameters_add(param, i->distance);
        
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        Node_set_distance(i, Parameters_value(param, 0));
		Parameters_pop(param);
        
        
        //SingleTreeLikelihood_update_all_nodes(stlk);
        SingleTreeLikelihood_update_one_node(stlk, i);
        SingleTreeLikelihood_update_one_node(stlk, j);
        data->index_param = Node_id(j);
        Parameters_add(param, j->distance);
        
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        Node_set_distance(j, Parameters_value(param, 0));
		Parameters_pop(param);
        
        //SingleTreeLikelihood_update_all_nodes(stlk);
        SingleTreeLikelihood_update_one_node(stlk, j);
        
        if ( -lnl2 - lnl < 0.01 ){
            lnl = -lnl2;
            break;
        }
        lnl = -lnl2;
    }
    
    //SingleTreeLikelihood_update_all_nodes(stlk);
	//return stlk->calculate(stlk);
    return lnl;
}

//#define DEBUG_TOPOLOGY_NNNI 1

// there is LESS than 2*(N-3) NNIs
double nnni_optimize_bl( struct TopologyOptimizer * opt ){
	double lnl = opt->model->logP(opt->model);
	
	unsigned nthreads = opt->threads;
	assert(nthreads > 0);
	
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	Model** pool = malloc(sizeof(Model*)*nthreads);
	Model** tlks = malloc(sizeof(Model*)*nthreads);
	pool[0] = opt->model;
	tlks[0] = opt->tlk;
	for (size_t i = 1; i < nthreads; i++) {
		pool[i] = opt->model->clone(opt->model, hash);
		tlks[i] = Hashtable_get(hash, opt->tlk->name); // not ref++
		Hashtable_empty(hash);
	}
	free_Hashtable(hash);
	
	Optimizer **opt_bls = (Optimizer**)malloc(nthreads*sizeof(Optimizer*));
	assert(opt_bls);
	Parameters **oneparameters = (Parameters**)malloc(nthreads*sizeof(Parameters*));
	assert(oneparameters);
	
	for ( int i = 0; i < nthreads; i++ ) {
		opt_bls[i] = new_Optimizer(OPT_BRENT);
		
		oneparameters[i] = new_Parameters(1);
		opt_set_data(opt_bls[i], pool[i]);
		opt_set_objective_function(opt_bls[i], model_negative_logP);
		opt_set_treelikelihood(opt_bls[i], tlks[i]->obj);
		//		opt_set_max_iteration(opt, max);
		//		opt_set_tolx(opt, precision);
	}
	
	Optimizer* full_opt = new_Optimizer(OPT_SERIAL_BRENT);
	opt_set_treelikelihood(full_opt, opt->tlk);
	opt_set_data(full_opt, opt->model);
	opt_set_objective_function(full_opt, model_negative_logP);
	
    int nNodes = Tree_node_count(((SingleTreeLikelihood*)(opt->tlk->obj))->tree);
    
    Node **temp_node_list[2];
    temp_node_list[0] = (Node**)malloc(sizeof(Node*) * nNodes );
    assert(temp_node_list[0]);
    temp_node_list[1] = (Node**)malloc(sizeof(Node*) * nNodes );
    assert(temp_node_list[1]);
    
    for ( int i = 0; i < nNodes; i++ ) {
        temp_node_list[0][i] = NULL;
        temp_node_list[1][i] = NULL;
    }
    
    
    int *map = ivector(nNodes);
    for( int i = 0; i < nNodes; i++ ){
        map[i] = -1;
    }
    
    int *map_index_to_post = ivector(nNodes);
    //int *order       = ivector(Tree_tip_count(opt->tlk->tree)-3);
    int index = 0;
    bool failed = false;
    opt_result status;
    
    double *branches = dvector(nNodes); // backup branch length in case we need to backtrack
    Tree_branch_length_to_vector(((SingleTreeLikelihood*)(opt->tlk->obj))->tree, branches);
    
    
    
#ifdef DEBUG_TOPOLOGY_NNNI
    printf("\nStart NNNI optimization LnL: %f\n\n", lnl);
#endif
    
    opt->moves = 0;
    
    while ( !failed ) {
        
        index = 0;
        int count = 0;
        
        // do not do NNIs on the root and one of its children
        // choose NNIs than do not change the 'fake' rooting, except for the first child of the root
        #pragma omp parallel for num_threads(nthreads)
        for ( int i = 0; i < nNodes; i++ ) {
            
            int tid = 0;
#if defined (_OPENMP)
            tid = omp_get_thread_num();
#endif
			SingleTreeLikelihood *tlk = tlks[tid]->obj;
            Node *node = Tree_node(tlk->tree, i);
            
            if ( Node_isroot(node) || ( Node_isroot(Node_parent(node)) && node == Node_right(Tree_root(tlk->tree))) ) continue;
            
            if ( !Node_isleaf(node) ) {
                
                // Do not do NNNI on a topology where the sibling is a leaf (d) and its parent is the root
                //        - a
                //       |
                //      -
                //   * | |
                //   --|  - b
                //  |  |
                //  |   - c
                //   - d
                //
                // Do not do NNNI on a topology where the children of sibling are both leaves (d) and its parent is the root (that is a regular NNI)
                //        - a
                //       |
                //      -
                //   * | |
                //   --|  - b
                //  |  |
                //  |   - c
                //  |   - d
                //   --|
                //      - e
                if( Node_isroot(Node_parent(node)) && (Node_isleaf(Node_sibling(node)) || (Node_isleaf(Node_sibling(node)->left) &&  Node_isleaf(Node_sibling(node)->right) )  ) ){
                    map[i] = -1;
                    continue;
                }
                
                Node *left    = Node_left(node);
                Node *right   = Node_right(node);
                Node *uncle   = Node_sibling(Node_parent(node));
                
                // Use the uncle (D)
                // NNI 1: Exchange A with D
                // NNI 2: Exchange B with D
                //          - A
                //       * |
                //       --
                //      |  |
                //    --|   - B
                //   |  |
                // --|   - c
                //   |
                //    -- D
                //
                //
                
                // Use the uncle (D)
                // NNI 1: Exchange A with D
                // NNI 2: Exchange B with D
                //          - A
                //         |
                //       --
                //     *|  |
                //    --|   - B
                //   |  |
                // --|   - c
                //   |   - d
                //   |  |
                //    --
                //      |  - e
                //       -|
                //         - f
                if( Node_isroot(Node_parent(node)) ){
                    map[i] = -1;
                    continue;
                    Node *sibling = Node_sibling(Node_parent(node));
                    if( Node_isleaf( Node_left(sibling) ) ){
                        uncle = Node_right(sibling);
                    }
                    else {
                        uncle = Node_left(sibling);
                    }
                }
                
                double lnl = tlk->calculate(tlk);
                
#ifdef DEBUG_TOPOLOGY_NNNI
                printf("\nNNI around %s (%d)\n", Node_name(node), node->postorder_idx);
                printf("  %f\n", lnl);
                
                printf("  rearragement 1 %s %s\n", left->name, uncle->name);
                printf("  rearragement 2 %s %s\n", right->name, uncle->name);
#endif
                double original_bl1 = Node_distance(node);
                double original_bl2 = Node_distance(Node_parent(node));
                
                // NNI 1
				NNI_move(tlk->tree, uncle, left);
				SingleTreeLikelihood_update_one_node(tlk, uncle);
				SingleTreeLikelihood_update_one_node(tlk, left);
				SingleTreeLikelihood_update_uppers(tlk);
				tlk->node_upper = node;
				
				// should not be commented
				double nni_1 = 0;//optimize_brent_branch_length_2(tlk, opt_bls[tid], data_brents[tid], oneparameters[tid], node, Node_parent(node));
                
                double bl_1 = 0;//Parameters_value(oneparameter, 0);
				NNI_move(tlk->tree, uncle, left);
				SingleTreeLikelihood_update_one_node(tlk, uncle);
				SingleTreeLikelihood_update_one_node(tlk, left);
				
                Tree_vector_to_branch_length(tlk->tree, branches);
                
                // NNI 2
				NNI_move(tlk->tree, uncle, right);
				SingleTreeLikelihood_update_one_node(tlk, uncle);
				SingleTreeLikelihood_update_one_node(tlk, right);
				SingleTreeLikelihood_update_uppers(tlk);
				tlk->node_upper = node;
				
				// should not be commented
				double nni_2 = 0;//optimize_brent_branch_length_2(tlk, opt_bls[tid], data_brents[tid], oneparameters[tid], node, Node_parent(node));
                
                double bl_2 = 0;//Parameters_value(oneparameter, 0);
                
                
				NNI_move(tlk->tree, uncle, right);
				SingleTreeLikelihood_update_one_node(tlk, uncle);
				SingleTreeLikelihood_update_one_node(tlk, right);
                
                Tree_vector_to_branch_length(tlk->tree, branches);
                
                double local_lnl = lnl;
                int local_position = 0;
                double local_bl = original_bl1;
                
                if ( nni_1-0.1 > lnl ) {
                    // NNI2 is the best
                    if ( nni_2 > nni_1) {
                        local_lnl      = nni_2;
                        local_position = 2;
                        local_bl = bl_2;
                    }
                    // NNI1 is the best
                    else {
                        local_lnl      = nni_1;
                        local_position = 1;
                        local_bl = bl_1;
                    }
                }
                else {
                    // NNI2 is the best
                    if ( nni_2-0.1 > lnl ) {
                        local_lnl      = nni_2;
                        local_position = 2;
                        local_bl = bl_2;
                    }
                }
#ifdef DEBUG_TOPOLOGY_NNNI
                printf("nn1: %f nni2: %f LnL: %f conf: %d\n", nni_1, nni_2, lnl, opt->positions[index]);
#endif
                int local_index;
#pragma omp critical
                {
                    local_index = index++;
                    
                    if( local_position != 0 ){
                        count++;
                    }
                }
                opt->lnls[local_index]      = local_lnl;
                opt->positions[local_index] = local_position;
                opt->branches[local_index] = local_bl;
                
                map_index_to_post[local_index] = i;
                map[i] = local_index;
            }
        }
        
#ifdef DEBUG_TOPOLOGY_NNNI
        fprintf(stderr, "%d potentials NNIs\n", count);
#endif
		SingleTreeLikelihood* tlk = opt->tlk->obj;
		Tree* tree = tlk->tree;
		//tlk->use_upper = false;
		
		if( count != 0 ){
			Optimizer *opt_bl = opt_bls[0];
			Parameters *oneparameter = oneparameters[0];
            
#ifdef DEBUG_TOPOLOGY_NNNI
            for ( int i = 0; i < Tree_tip_count(tlk->tree)-3; i++ ) {
                if( opt->positions[i] != 0 ){
                    printf("%s %f %d\n", Tree_node(tlk->tree,map_index_to_post[i])->name, opt->lnls[i], opt->positions[i]);
                }
            }
#endif
            if(count > 1){
                remove_adjacent_nnis(opt, Tree_root(tlk->tree), map);
            }
            
            count = 0;
            for ( int j = 0; j < Tree_tip_count(tlk->tree)-3; j++ ) {
                if ( opt->positions[j] != 0 ) {
                    count++;
                }
                //order[j] = j;
            }
#ifdef DEBUG_TOPOLOGY_NNNI
            fprintf(stderr, "\n%d potentials non adjacent NNIs\n", count);
            for ( int i = 0; i < Tree_tip_count(tlk->tree)-3; i++ ) {
                if( opt->positions[i] != 0 ){
                    printf("%s %f %d\n", Tree_node(tlk->tree,map_index_to_post[i])->name, opt->lnls[i], opt->positions[i]);
                }
            }
#endif
            
            _sort_by_nni(/*order,*/ opt, map_index_to_post, Tree_tip_count(tlk->tree)-3);
            
            _sort_decreasing_lnls(/*order,*/ opt, map_index_to_post, count);
            
#ifdef DEBUG_TOPOLOGY_NNNI
            fprintf(stderr, "\nBest NNI %f\n", opt->lnls[0]);
#endif
            Node **nodes = Tree_nodes(tlk->tree);
            
            // used to backtrack
            for ( int i = 0; i < (Tree_node_count(tlk->tree)); i++ ) {
                temp_node_list[0][i] = NULL;
                temp_node_list[1][i] = NULL;
            }
            
            count *= opt->K;
            count = imax(1, count); // in case we have count == 1, we would get (count * 0.75) == 0
            double nni_lnl = lnl;
            
            while ( count > 0 ) {
                // apply NNIs
                for ( int j = 0; j < count; j++ ) {
                    int node_index = map_index_to_post[j]; // NNI node
                    
                    Node *uncle   = Node_sibling(Node_parent(nodes[node_index]));
                    
                    //                    if( Node_isroot(Node_parent(nodes2[node_index])) ){
                    //                        sibling = Node_left(sibling);
                    //                    }
                    
                    Node *child   = NULL;
                    
                    // uncle and left
                    if( opt->positions[j] == 1 ){
                        child = Node_left(nodes[node_index]);
                    }
                    // uncle and right
                    else {
                        child = Node_right(nodes[node_index]);
                    }
                    
                    temp_node_list[0][j] = uncle;
                    temp_node_list[1][j] = child;
                    
                    Node_swap_parents(uncle, child);
                }
                
                SingleTreeLikelihood_update_all_nodes(tlk);
                Tree_set_topology_changed(tlk->tree);
				
				status = opt_maximize( full_opt, NULL, &nni_lnl);
				if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
                
#ifdef DEBUG_TOPOLOGY_NNNI
                printf("nni LnL: %f (%f)\n", nni_lnl, lnl);
#endif
                
                
                // nni_lnl can be > lnl but < best single nni_lnl
                //if( nni_lnl > opt->lnls[0] || count == 1 ){
                if( nni_lnl > lnl || count == 1 ){
#ifdef DEBUG_TOPOLOGY_NNNI
                    printf("LnL %f count: %d\n", nni_lnl, count);
#endif
                    opt->moves += count;
                    break;
                }
                
#ifdef DEBUG_TOPOLOGY_NNNI
                printf("backtrack\n");
#endif
                // backtrack to the previous tree topology using tlk
                for ( int j = count-1; j >= 0; j-- ) {
                    Node_swap_parents(temp_node_list[0][j], temp_node_list[1][j]);
                }
                
                Tree_set_topology_changed(tlk->tree);
                Tree_vector_to_branch_length(tlk->tree, branches);
                
                count *= 0.5;
            }
            
            Tree_branch_length_to_vector(tlk->tree, branches);
            
            // Apply the NNIs to the others
			for ( int i = 1; i < opt->threads; i++ ) {
				SingleTreeLikelihood *tlk = tlks[i]->obj;
                Node **nodes = Tree_nodes(tlk->tree);
                
                for ( int j = 0; j < count; j++ ) {
                    int node_index = map_index_to_post[j]; // NNI node
                    
                    Node *uncle   = Node_sibling(Node_parent(nodes[node_index]));
                    
                    //                    if( Node_isroot(Node_parent(nodes[node_index])) ){
                    //                        sibling = Node_left(sibling);
                    //                    }
                    
                    Node *child   = NULL;
                    
                    // uncle and left
                    if( opt->positions[j] == 1 ){
                        child = Node_left(nodes[node_index]);
                    }
                    // uncle and right
                    else {
                        child = Node_right(nodes[node_index]);
                    }
                    
                    Node_swap_parents(uncle, child);
				}
				Tree_vector_to_branch_length(tlk->tree, branches);
				Tree_set_topology_changed(tlk->tree);
            }
			
            lnl = nni_lnl;
            
            //break;// only one round
        }
        else {
            failed = true;
        }
    }
    
#ifdef DEBUG_TOPOLOGY_NNNI
    fprintf(stderr, "NNNI round LnL: %f\n", lnl);
#endif
    
    for ( int i = 0; i < opt->threads; i++ ) {
        free_Optimizer(opt_bls[i]);
        free_Parameters(oneparameters[i]);
    }
    free_Optimizer(full_opt);
    free(opt_bls);
    free(oneparameters);
    
    free(map);
    free(map_index_to_post);
    free(temp_node_list[0]);
    free(temp_node_list[1]);
    free(branches);
    //free(order);
    for (int i = 1; i < nthreads; i++) {
        pool[i]->free(pool[i]);
    }
    free(pool);
    free(tlks);
    
    opt->best_lnl = lnl;
    SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)(opt->tlk->obj);
    tlk->use_upper = false;
    
    return lnl;
}

//TODO: do not allow NNIs on constrained nodes
// there is 2n-6 NNIs in a rooted tree
double nni_optimize_heights( struct TopologyOptimizer * opt ){
//    SingleTreeLikelihood *tlk = opt->tlk;
//    double lnl = tlk->calculate(tlk);
//    Node **nodes = NULL;
//
//    Node **temp_node_list[2];
//    temp_node_list[0] = (Node**)malloc(sizeof(Node*) * (Tree_node_count(tlk->tree)-2) );
//    assert(temp_node_list[0]);
//    temp_node_list[1] = (Node**)malloc(sizeof(Node*) * (Tree_node_count(tlk->tree)-2) );
//    assert(temp_node_list[1]);
//
//    for ( int i = 0; i < (Tree_node_count(tlk->tree)-2); i++ ) {
//        temp_node_list[0][i] = NULL;
//        temp_node_list[1][i] = NULL;
//    }
//
//
//	Optimizer *opt_height = new_Optimizer(OPT_BRENT);
//    BrentData *data_brent = new_BrentData( tlk );
//    opt_set_data(opt_height, data_brent);
//    opt_set_objective_function(opt_height, _brent_optimize_height);
//    opt_set_max_iteration(opt_height, tlk->opt.heights.max_iteration);
//    opt_set_tolx( opt_height, tlk->opt.heights.tolx);
//
//	Parameters *oneparameter = new_Parameters(1);
//
//    int nNodes = Tree_node_count(tlk->tree);
//    int *map = ivector(nNodes);
//    map[nNodes-1] = -1;
//    int *map_index_to_post = ivector(Tree_node_count(tlk->tree)-1);
//    //int *order       = ivector(Tree_tip_count(tlk->tree)-2);
//    int index = 0;
//    bool failed = false;
//
//#ifdef DEBUG_TOPOLOGY
//    printf("\nStart NNI optimization LnL: %f\n\n", lnl);
//#endif
//
//    opt->moves = 0;
//
//    double *heights = dvector(Tree_node_count(tlk->tree)); // backup
//    bool done = false;
//
//    while ( !failed ) {
//
//        nodes = Tree_get_nodes(tlk->tree, POSTORDER);
//
//        Tree_heights_to_vector(tlk->tree, heights);
//
//        index = 0;
//        int count = 0;
//        done = false;
//
//        // do not do NNIs on the root
//        for ( int i = 0; i < Tree_node_count(tlk->tree)-1; i++ ) {
//
//            if ( !Node_isleaf(nodes[i]) ) {
//                // The node with the branch on which we do NNIs
//                Node *node = Tree_get_node(tlk->tree, POSTORDER, i);
//
//                // We don't want to do NNIs on both root's children
//                // If one of the children is a leaf the other one is not and this condition does not apply
//                if ( Node_isroot(Node_parent(node)) && !done) {
//                    done = true;
//                }
//                else if( Node_isroot(Node_parent(node)) && done ){
//                    continue;
//                }
//
//                Node *left    = Node_left(node);
//                Node *right   = Node_right(node);
//                Node *sibling = Node_sibling(node);
//
//                double lnl = tlk->calculate(tlk);
//
//#ifdef DEBUG_TOPOLOGY2
//                printf("\nNNI around %s (%d)\n", Node_name(node), node->postorder_idx);
//                printf("  %f\n", lnl);
//
//                printf("  rearragement 1 %s %s\n", left->name, sibling->name);
//                printf("  rearragement 2 %s %s\n", right->name, sibling->name);
//#endif
//
//                // NNI 1
//
//                if( Node_height(node) < Node_height(sibling) ){
//                    Node_set_height(node, Node_height(sibling)+ Node_time_elapsed(sibling)*0.5);
//                }
//
//                Node_swap_parents(sibling, left);
//
//                Tree_update_topology(tlk->tree);
//                Tree_constraint_heights(tlk->tree);
//                SingleTreeLikelihood_rearrange_partials(tlk);
//
//                double nni_1 = optimize_brent_height_all(tlk, opt_height, data_brent, oneparameter);
//#ifdef DEBUG_TOPOLOGY2
//                printf("NNI 1: %f\n", nni_1);
//#endif
//                SingleTreeLikelihood_rearrange( tlk, sibling, left );
//                //SingleTreeLikelihood_copy_partials(tlk, tlk2);
//
//                Tree_vector_to_heights(heights, tlk->tree);
//                Tree_constraint_heights(tlk->tree);
//
//
//                // NNI 2
//
//                if( Node_height(node) < Node_height(sibling) ){
//                    Node_set_height(node, Node_height(sibling)+ Node_time_elapsed(sibling)*0.5);
//                }
//
//                Node_swap_parents(sibling, right);
//
//                Tree_update_topology(tlk->tree);
//                Tree_constraint_heights(tlk->tree);
//                SingleTreeLikelihood_rearrange_partials(tlk);
//
//
//                double nni_2 = optimize_brent_height_all(tlk, opt_height, data_brent, oneparameter);
//#ifdef DEBUG_TOPOLOGY2
//                printf("NNI 2: %f\n", nni_2);
//#endif
//                //SingleTreeLikelihood_rearrange( tlk2, sibling, right );
//                Node_swap_parents(sibling, right);
//
//                Tree_update_topology(tlk->tree);
//                Tree_vector_to_heights(heights, tlk->tree);
//                Tree_constraint_heights(tlk->tree);
//                SingleTreeLikelihood_rearrange_partials(tlk);
//
//
//                if ( nni_1-0.001 > lnl ) {
//                    // NNI2 is the best
//                    if ( nni_2 > nni_1) {
//                        opt->lnls[index]      = nni_2;
//                        opt->positions[index] = 2;
//                    }
//                    // NNI1 is the best
//                    else {
//                        opt->lnls[index]      = nni_1;
//                        opt->positions[index] = 1;
//
//                    }
//                    count++;
//                }
//                else {
//                    // NNI2 is the best
//                    if ( nni_2-0.001 > lnl ) {
//                        opt->lnls[index]      = nni_2;
//                        opt->positions[index] = 2;
//                        count++;
//                    }
//                    // no good NNI
//                    else {
//                        opt->lnls[index]      = lnl;
//                        opt->positions[index] = 0;
//
//                    }
//                }
//#ifdef DEBUG_TOPOLOGY
//                printf("nn1: %f nni2: %fLnL: %f\n", nni_1, nni_2, lnl);
//#endif
//
//                map_index_to_post[index] = i;
//                map[i] = index++;
//            }
//        }
//        // At this stage tlk and tlk2 should have the same topology and branch length
//
//#ifdef DEBUG_TOPOLOGY
//        fprintf(stderr, "%d potentials NNIs\n", count);
//#endif
//
//        if( count != 0 ){
//#ifdef DEBUG_TOPOLOGY
//            for ( int i = 0; i < Tree_tip_count(tlk->tree)-2; i++ ) {
//                if( opt->positions[i] != 0 ){
//                    printf("%s %f %d\n", nodes[map_index_to_post[i]]->name, opt->lnls[i], opt->positions[i]);
//                }
//            }
//#endif
//
//            remove_adjacent_nnis(opt, Tree_root(tlk->tree), map);
//
//            count = 0;
//            for ( int j = 0; j < Tree_tip_count(tlk->tree)-2; j++ ) {
//                if ( opt->positions[j] != 0 ) {
//                    count++;
//                }
//                //order[j] = j;
//            }
//#ifdef DEBUG_TOPOLOGY
//            fprintf(stderr, "%d potentials non adjacent NNIs\n", count);
//#endif
//
//            _sort_by_nni(/*order,*/ opt, map_index_to_post, Tree_tip_count(tlk->tree)-3);
//
//            _sort_decreasing_lnls(/*order,*/ opt, map_index_to_post, count);
//
//            nodes = Tree_get_nodes(tlk->tree, POSTORDER);
//
//            // used to backtrack
//            for ( int i = 0; i < (Tree_node_count(tlk->tree)-2); i++ ) {
//                temp_node_list[0][i] = NULL;
//                temp_node_list[1][i] = NULL;
//            }
//
//            count *= opt->K;
//            count = imax(1, count); // in case we have count == 1, we would get (count * 0.75) == 0
//            double nni_lnl = lnl;
//
//            while ( count > 0 ) {
//                // apply NNIs
//                for ( int j = 0; j < count; j++ ) {
//                    int node_index = map_index_to_post[j]; // NNI node
//
//                    Node *sibling = Node_sibling(nodes[node_index]);
//                    Node *child   = NULL;
//
//                    // sibling and left
//                    if( opt->positions[j] == 1 ){
//                        child = Node_left(nodes[node_index]);
//                    }
//                    // sibling and right
//                    else {
//                        child = Node_right(nodes[node_index]);
//                    }
//
//                    if( Node_height(nodes[node_index]) < Node_height(sibling) ){
//                        Node_set_height(nodes[node_index], Node_height(sibling)+ Node_time_elapsed(sibling)*0.5);
//                    }
//
//                    temp_node_list[0][j] = sibling;
//                    temp_node_list[1][j] = child;
//
//                    Node_swap_parents(sibling, child);
//                }
//
//                Tree_update_topology(tlk->tree);
//                SingleTreeLikelihood_rearrange_partials(tlk);
//
//                // Optimize all branches
//                nni_lnl = optimize_brent_height_all(tlk, opt_height, data_brent, oneparameter);
//
//#ifdef DEBUG_TOPOLOGY
//                printf("nni LnL: %f (%f)\n", nni_lnl, lnl);
//#endif
//
//                if( nni_lnl > lnl || count == 1 ){
//                    opt->moves += count;
//                    break;
//                }
//
//#ifdef DEBUG_TOPOLOGY
//                printf("backtrack\n");
//#endif
//                // backtrack to the previous tree topology
//                for ( int j = count-1; j >= 0; j-- ) {
//                    Node_swap_parents(temp_node_list[0][j], temp_node_list[1][j]);
//                }
//
//                Tree_update_topology(tlk->tree);
//                Tree_vector_to_heights(heights, tlk->tree);
//                Tree_constraint_heights(tlk->tree);
//                SingleTreeLikelihood_rearrange_partials(tlk);
//
//                count *= 0.5;
//            }
//
//            lnl = nni_lnl;
//        }
//        else {
//            failed = true;
//        }
//    }
//
//#ifdef DEBUG_TOPOLOGY2
//    fprintf(stderr, "NNI round LnL: %f\n", lnl);
//#endif
//
//    free_BrentData(data_brent);
//    free_Optimizer(opt_height);
//    free_Parameters(oneparameter);
//    free(map);
//    free(map_index_to_post);
//    free(temp_node_list[0]);
//    free(temp_node_list[1]);
//    //free(order);
//    free(heights);
//
//    opt->best_lnl = lnl;
//
//    return lnl;
	return 0;
}


// Use the sibling (c)
// NNI 1: Exchange a with c
// NNI 2: Exchange b with c
//        - a
//      *|
//      -
//     | |
//   --|  - b
//  |  |
//  |   - c
//  |  - d
//  | |
//   -
//    |
//     - e
double nni_optimize_bl_parsimony( struct TopologyOptimizer * opt ){
    SingleTreeLikelihood *tlk = opt->tlk;
    
    
    double lnl; //= tlk->calculate(tlk);
    
    Node **nodes = NULL;
    SingleTreeLikelihood *tlk2 = clone_SingleTreeLikelihood_share(tlk, true, false);
    Parsimony *parsimony = new_Parsimony(tlk2->sp, tlk2->tree);
    lnl = -parsimony->calculate(parsimony);
    
    Node **temp_node_list[2];
    temp_node_list[0] = (Node**)malloc(sizeof(Node*) * (Tree_node_count(tlk->tree)-2) );
    assert(temp_node_list[0]);
    temp_node_list[1] = (Node**)malloc(sizeof(Node*) * (Tree_node_count(tlk->tree)-2) );
    assert(temp_node_list[1]);
    
    for ( int i = 0; i < (Tree_node_count(tlk->tree)-2); i++ ) {
        temp_node_list[0][i] = NULL;
        temp_node_list[1][i] = NULL;
    }
    
    
    // fix the right node, set it to 0 and add its distance to the left node
    // we can't do that if the substitution model is not reversible
    Node *right = Node_right( Tree_root(tlk->tree) );
    Node *left  = Node_left( Tree_root(tlk->tree) );
    Node_set_distance(left, Node_distance(right)+Node_distance(left) );
    Node_set_distance(right, 0);
    Parameter_set_estimate(right->distance, false);
    
	Optimizer *opt_bl = new_Optimizer(OPT_BRENT);
    BrentData *data_brent = new_BrentData( tlk2 );
    data_brent->f = standard_loglikelihood_brent;
    if( tlk->use_upper ){
        data_brent->f = standard_loglikelihood_upper_brent;
    }
    
    opt_set_data(opt_bl, data_brent);
    opt_set_objective_function(opt_bl, optimize_brent_branch_length);
    opt_set_max_iteration(opt_bl, tlk->opt.bl.max_iteration);
    opt_set_tolx( opt_bl, tlk->opt.bl.tolx);
    
	Parameters *oneparameter = new_Parameters(1);
    
    
    int nNodes = Tree_node_count(tlk->tree);
    int *map = ivector(nNodes);
    map[nNodes-1] = map[nNodes-2] = -1;
    int *map_index_to_post = ivector(Tree_node_count(tlk->tree)-2);
    //int *order       = ivector(Tree_tip_count(tlk->tree)-3);
    int index = 0;
    bool failed = false;
    
#ifdef DEBUG_TOPOLOGY
    printf("\nStart NNI optimization LnL: %f\n\n", lnl);
#endif
    
    opt->moves = 0;
    
    while ( !failed ) {
        
        nodes = Tree_get_nodes(tlk->tree, POSTORDER);
        
        index = 0;
        int count = 0;
        
        // do not do NNIs on the root and one of its children
        // choose NNIs than do not change the 'fake' rooting, except for the first child of the root
        for ( int i = 0; i < Tree_node_count(tlk->tree)-2; i++ ) {
            
            if ( !Node_isleaf(nodes[i]) ) {
                // The node with the branch on which we do NNIs
                Node *node = Tree_get_node(tlk2->tree, POSTORDER, i);
                
                // Do not do NNI on a topology where the sibling is a leaf (d) and its parent is the root
                //        - a
                //       |
                //      -
                //   * | |
                //   --|  - b
                //  |  |
                //  |   - c
                //   - d
                if( Node_isroot(Node_parent(node)) && Node_isleaf(Node_sibling(node)) ){
                    continue;
                }
                
                Node *left    = Node_left(node);
                Node *right   = Node_right(node);
                Node *sibling = Node_sibling(node);
                
                // Use the left node of the sibling (d)
                // NNI 1: Exchange (a,b) with d
                // NNI 2: Exchange c with d
                //        - a
                //       |
                //      -
                //   * | |
                //   --|  - b
                //  |  |
                //  |   - c
                //  |  - d
                //  | |
                //   -
                //    |
                //     - e
                if( Node_isroot(Node_parent(node)) ){
                    sibling = Node_left(sibling);
                }
                
                double lnl;// = tlk2->calculate(tlk2);
                
                lnl = -parsimony->calculate(parsimony); // minus because we want to minimize this score
                
#ifdef DEBUG_TOPOLOGY
                printf("\nNNI around %s (%d)\n", Node_name(node), node->postorder_idx);
                printf("  %f\n", lnl);
                
                printf("  rearragement 1 %s %s\n", left->name, sibling->name);
                printf("  rearragement 2 %s %s\n", right->name, sibling->name);
#endif
                double original_bl = Node_distance(node);
                Parameters_add(oneparameter, node->distance);
                
                
                // NNI 1
                Tree_rearrange( tlk->tree, sibling, left );
                
                //SingleTreeLikelihood_rearrange_partials(tlk);
                
                double nni_1 = 0;
                
                //                data_brent->index_param = node->postorder_idx; // that's the updated postorder index in the new topology
                //                status = opt_maximize( opt_bl, oneparameter, &nni_1);
                //                if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
				
                nni_1 = -parsimony->calculate(parsimony);
                
                Tree_rearrange( tlk->tree, sibling, left );
                
                //SingleTreeLikelihood_rearrange_partials(tlk);
                
                Node_set_distance(node, original_bl);
                
                
                // NNI 2
                Tree_rearrange( tlk->tree, sibling, right );
                
                //SingleTreeLikelihood_rearrange_partials(tlk);
                
                double nni_2 = 0;
                
                //                data_brent->index_param = node->postorder_idx;
                //                status = opt_maximize( opt_bl, oneparameter, &nni_2);
                //                if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
				//				Parameters_pop(oneparameter);
                
                //parsimony->update(parsimony);
                nni_2 = -parsimony->calculate(parsimony);
                
                Tree_rearrange( tlk->tree, sibling, right );
                
                //SingleTreeLikelihood_rearrange_partials(tlk);
                
                Node_set_distance(node, original_bl);
                
                
                if ( nni_1 > lnl ) {
                    // NNI2 is the best
                    if ( nni_2 > nni_1) {
                        opt->lnls[index]      = nni_2;
                        opt->positions[index] = 2;
                    }
                    // NNI1 is the best
                    else {
                        opt->lnls[index]      = nni_1;
                        opt->positions[index] = 1;
                        
                    }
                    count++;
                }
                else {
                    // NNI2 is the best
                    if ( nni_2 > lnl ) {
                        opt->lnls[index]      = nni_2;
                        opt->positions[index] = 2;
                        count++;
                    }
                    // no good NNI
                    else {
                        opt->lnls[index]      = lnl;
                        opt->positions[index] = 0;
                        
                    }
                }
#ifdef DEBUG_TOPOLOGY
                printf("nn1: %f nni2: %fLnL: %f\n", nni_1, nni_2, lnl);
#endif
                
                map_index_to_post[index] = i;
                map[i] = index++;
            }
        }
        
        // At this stage tlk and tlk2 should have the same topology and branch length
        
#ifdef DEBUG_TOPOLOGY
        fprintf(stderr, "%d potentials NNIs\n", count);
#endif
        
        if( count != 0 ){
#ifdef DEBUG_TOPOLOGY
            for ( int i = 0; i < Tree_tip_count(tlk->tree)-3; i++ ) {
                if( opt->positions[i] != 0 ){
                    printf("%s (%d) %f %d\n", nodes[map_index_to_post[i]]->name, opt->lnls[i], opt->positions[i]);
                }
            }
#endif
            if(count > 1){
                remove_adjacent_nnis(opt, Tree_root(tlk->tree), map);
            }
            
            count = 0;
            for ( int j = 0; j < Tree_tip_count(tlk->tree)-3; j++ ) {
                if ( opt->positions[j] != 0 ) {
                    count++;
                }
                //order[j] = j;
            }
#ifdef DEBUG_TOPOLOGY
            fprintf(stderr, "%d potentials non adjacent NNIs\n", count);
#endif
            
            _sort_by_nni(/*order,*/ opt, map_index_to_post, Tree_tip_count(tlk->tree)-3);
            
            _sort_decreasing_lnls(/*order,*/ opt, map_index_to_post, count);
            
            Node **nodes2 = Tree_get_nodes(tlk2->tree, POSTORDER);
            
            // used to backtrack
            for ( int i = 0; i < (Tree_node_count(tlk->tree)-2); i++ ) {
                temp_node_list[0][i] = NULL;
                temp_node_list[1][i] = NULL;
            }
            
            count *= opt->K;
            count = imax(1, count); // in case we have count == 1, we would get (count * 0.75) == 0
            double nni_lnl = lnl;
            
            while ( count > 0 ) {
                // apply NNIs
                for ( int j = 0; j < count; j++ ) {
                    int node_index = map_index_to_post[j]; // NNI node
                    
                    Node *sibling = Node_sibling(nodes2[node_index]);
                    
                    if( Node_isroot(Node_parent(nodes2[node_index])) ){
                        sibling = Node_left(sibling);
                    }
                    
                    Node *child   = NULL;
                    
                    // sibling and left
                    if( opt->positions[j] == 1 ){
                        child = Node_left(nodes2[node_index]);
                    }
                    // sibling and right
                    else {
                        child = Node_right(nodes2[node_index]);
                    }
                    
                    //Node_set_distance(nodes2[node_index], branches[ opt->order[j] ]);
                    
                    temp_node_list[0][j] = sibling;
                    temp_node_list[1][j] = child;
                    
                    Node_swap_parents(sibling, child);
                }
                
                Tree_update_topology(tlk2->tree);
                //SingleTreeLikelihood_rearrange_partials(tlk2);
                
                // Optimize all branches
                //                nni_lnl = optimize_brent_branch_length_all(tlk2, opt_bl, data_brent, oneparameter, 1);
                
                nni_lnl = -parsimony->calculate(parsimony);
                
#ifdef DEBUG_TOPOLOGY
                printf("nni LnL: %f (%f)\n", nni_lnl, lnl);
#endif
                
                if( nni_lnl > lnl || count == 1 ){
                    opt->moves += count;
                    break;
                }
                
#ifdef DEBUG_TOPOLOGY
                printf("backtrack\n");
#endif
                // backtrack to the previous tree topology using tlk
                for ( int j = count-1; j >= 0; j-- ) {
                    Node_swap_parents(temp_node_list[0][j], temp_node_list[1][j]);
                    
                }
                Tree_update_topology(tlk2->tree);
                
                Tree_copy_distances(tlk->tree, tlk2->tree);
                
                //SingleTreeLikelihood_rearrange_partials(tlk2);
                
                //SingleTreeLikelihood_copy_partials(tlk, tlk2);
                
                count *= 0.5;
            }
            
            // tlk now has the old topology with lower log-likelihood
            // update the tree topology, partials... of tlk using tlk2
            for ( int j = 0; j <count; j++ ) {
                Tree_swap_parents_by_name( tlk->tree, Node_name(temp_node_list[0][j]), Node_name(temp_node_list[1][j]) );
            }
            
            Tree_update_topology(tlk->tree);
            
            Tree_copy_distances(tlk2->tree, tlk->tree);
            
            //SingleTreeLikelihood_rearrange_partials(tlk);
            
            //SingleTreeLikelihood_copy_partials(tlk2, tlk);
            
            lnl = nni_lnl;
        }
        else {
            failed = true;
        }
    }
    
#ifdef DEBUG_TOPOLOGY
    fprintf(stderr, "NNI round LnL: %f\n", lnl);
#endif
    
    free_Parsimony(parsimony);
    free_BrentData(data_brent);
    free_Optimizer(opt_bl);
    free_Parameters(oneparameter);
    free_SingleTreeLikelihood_share(tlk2, true, false);
    free(map);
    free(map_index_to_post);
    free(temp_node_list[0]);
    free(temp_node_list[1]);
    //free(order);
    
    opt->best_lnl = lnl;
    
    SingleTreeLikelihood_rearrange_partials(tlk);
    
    return lnl;
}
