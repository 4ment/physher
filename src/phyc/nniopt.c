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
#include "matrix.h"
#include "treeio.h"



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

