/*
 *  spr.c
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

#include "spropt.h"

#include <stdio.h>
#include <assert.h>

#include "treesearch.h"
#include "parsimony.h"
#include "optimize.h"
#include "matrix.h"

#include <time.h>

#include "treeio.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//#define DEBUG_TOPOLOGY_SPR 1


void compare_trees(Node *node1, Node *node2 ){
    if( node1 == NULL && node2 == NULL ){
        return;
    }
    else if( node1 == NULL && node2 == NULL ){
        
    }
    else if( node1 != NULL && node2 != NULL ){
        if( strcmp(node1->name, node2->name) != 0 ){
            printf("Error %s %s\n", node1->name, node2->name);
            exit(1);
        }
        if( !Node_isroot(node1) &&  !Node_isroot(node2) && strcmp(node1->parent->name, node2->parent->name) != 0 ){
            printf("Error parents %s (%s) %s (%s)\n", node1->name, node1->parent->name, node2->name, node2->parent->name);
            exit(1);
            
        }
        compare_trees(node1->left, node2->left);
        compare_trees(node1->right, node2->right);
    }
    else {
        printf("crap\n");
        exit(1);
    }
}

static void sort_asc_spr_parsimony( int *prunes, int *grafts, double *scores, int size ){
    bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( scores[i] > scores[i+1] ) {
				done = false;
                swap_int(&prunes[i], &prunes[i+1]);
                swap_int(&grafts[i], &grafts[i+1]);
                dswap(&scores[i], &scores[i+1]);
			}
		}
		size--;
	}
}

static void sort_asc_spr_parsimony2( int *prunes, int *grafts, double *scores, int *ds, int size ){
    bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( scores[i] > scores[i+1] ) {
				done = false;
                swap_int(&prunes[i], &prunes[i+1]);
                swap_int(&grafts[i], &grafts[i+1]);
                swap_int(&ds[i], &ds[i+1]);
                dswap(&scores[i], &scores[i+1]);
			}
		}
		size--;
	}
}

static void sort_desc_spr_parsimony( int *prunes, int *grafts, double *scores, int *ds, int size ){
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( scores[i] < scores[i+1] ) {
				done = false;
				swap_int(&prunes[i], &prunes[i+1]);
				swap_int(&grafts[i], &grafts[i+1]);
				swap_int(&ds[i], &ds[i+1]);
				dswap(&scores[i], &scores[i+1]);
			}
		}
		size--;
	}
}


// backup contains the two branches together
// params should be contrained between 0 and backup
// should not be called on the children of root
static double _optimize_brent_branch_length_constrained( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
	
    Node *node = Tree_get_node(mydata->tlk->tree, POSTORDER, mydata->index_param);
    
    Node_set_distance( node, Parameters_value(params, 0) );
    Node_set_distance( Node_parent(node), mydata->backup[0]-Parameters_value(params, 0) );
    
    SingleTreeLikelihood_update_one_node( mydata->tlk, node);
    double lnl = fabs(mydata->f(mydata));
    //printf("sf %f - %f %f %f\n",lnl, mydata->backup[0], Parameters_value(params, 0),mydata->backup[0]-Parameters_value(params, 0));
    
    return lnl;
}


// node node1 is the node that also changes its parent
// node node2 is the grafted node with meaningless branch length
double optimize_brent_branch_length_3nodes( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, Node *node1, Node *node2, Node *node3 ){
    opt_result status = OPT_NEED_CHECK;
    double lnl = stlk->calculate(stlk);
    double lnl2 = lnl+10;
    data->index_param = -1;
    //printf("-- lnl %f\n",lnl);
    
    data->backup[0] = Node_distance(node1) + Node_distance(Node_parent(node1));
    
    if( Node_isroot(node1) ){
        printf("Cannot be called on the root\n%s\n%s: %d\n", __func__, __FILE__, __LINE__);
        exit(1);
    }
    double upper = Parameter_upper(node1->distance);
    Parameter_set_upper(node1->distance, data->backup[0]);

    int i = 0;
    while ( i < 3 ) {
        
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node1->postorder_idx;
        Parameters_add(param, node1->distance);
        
        opt_set_objective_function(opt, _optimize_brent_branch_length_constrained);
//        SingleTreeLikelihood_set_upper_function(stlk, calculate_uppper_2nodes);
		
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_distance( node1, Parameters_value(param, 0) );
        Node_set_distance( Node_parent(node1), data->backup[0]-Parameters_value(param, 0) );
		Parameters_pop(param);
        
        //printf("== lnl %f (%f) %f  %s %f %s %f\n", -lnl2, lnl, data->backup[0], node1->name, Parameters_value(param, 0), node1->parent->name, (data->backup[0]-Parameters_value(param, 0)) );
        
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node2->postorder_idx;
        Parameters_add(param, node2->distance);
        
        opt_set_objective_function(opt, optimize_brent_branch_length);
//        SingleTreeLikelihood_set_upper_function(stlk, NULL);
		
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        Node_set_distance(node2, Parameters_value(param, 0));
        Parameters_pop(param);
        //printf("++ lnl %f (%f) %s %f\n", -lnl2, lnl,node2->name, Parameters_value(param, 0) );
        
        
        if ( -lnl2 - lnl < 0.001 ){
            //printf("  %f %f\n", lnl2,lnl);
            lnl = -lnl2;
            break;
        }
        lnl = -lnl2;
    }
    
    Parameter_set_upper(node1->distance, upper);
    SingleTreeLikelihood_update_all_nodes(stlk);
    
    //return stlk->calculate(stlk);
    return lnl;
}

// node node1 is the node that also changes its parent
// node node2 is the pruned node with meaningless branch length
double optimize_brent_branch_length_2b( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, Node *node1, Node *node2 ){
    opt_result status = OPT_NEED_CHECK;
    double lnl = stlk->calculate(stlk);
    double lnl2 = lnl+10;
    data->index_param = -1;
    //printf("-- lnl %f\n",lnl);
    
    data->backup[0] = Node_distance(node1) + Node_distance(Node_parent(node1));
    
    if( Node_isroot(node1) ){
        printf("Cannot be called on the root\n%s\n%s: %d\n", __func__, __FILE__, __LINE__);
        exit(1);
    }
    double upper = Parameter_upper(node1->distance);
    Parameter_set_upper(node1->distance, data->backup[0]);
    
    while ( 1 ) {
        
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node1->postorder_idx;
        Parameters_add(param, node1->distance);
        
        opt_set_objective_function(opt, _optimize_brent_branch_length_constrained);
//        SingleTreeLikelihood_set_upper_function(stlk, calculate_uppper_2nodes);
		
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_distance( node1, Parameters_value(param, 0) );
        Node_set_distance( Node_parent(node1), data->backup[0]-Parameters_value(param, 0) );
		Parameters_pop(param);
        
        //printf("== lnl %f (%f) %f  %s %f %s %f\n", -lnl2, lnl, data->backup[0], node1->name, Parameters_value(param, 0), node1->parent->name, (data->backup[0]-Parameters_value(param, 0)) );
        
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node2->postorder_idx;
        Parameters_add(param, node2->distance);
        
        opt_set_objective_function(opt, optimize_brent_branch_length);
//        SingleTreeLikelihood_set_upper_function(stlk, NULL);
		
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        Node_set_distance(node2, Parameters_value(param, 0));
        Parameters_pop(param);
        //printf("++ lnl %f (%f) %s %f\n", -lnl2, lnl,node2->name, Parameters_value(param, 0) );
        
        
        if ( -lnl2 - lnl < 0.001 ){
            //printf("  %f %f\n", lnl2,lnl);
            lnl = -lnl2;
            break;
        }
        lnl = -lnl2;
    }
    
    Parameter_set_upper(node1->distance, upper);
    SingleTreeLikelihood_update_all_nodes(stlk);
    
    //return stlk->calculate(stlk);
    return lnl;
}

// node node1 is the node that also changes its parent
// node node2 is the pruned node with meaningless branch length
double optimize_brent_branch_length_3( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, Node *node1, Node *node2 ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = stlk->calculate(stlk);
    double lnl2 = lnl+10;
    data->index_param = -1;
#ifdef DEBUG_TOPOLOGY_SPR
    printf("-- lnl %f\n",lnl);
#endif
    
    if( Node_isroot(node1) ){
        printf("Cannot be called on the root\n%s\n%s: %d\n", __func__, __FILE__, __LINE__);
        exit(1);
    }
    
#ifdef DEBUG_TOPOLOGY_SPR
    printf("%f %f %f %f\n", Node_distance(node1), Node_distance(node1->left), Node_distance(node1->right), Node_distance(node2));
#endif
    while ( 1 ) {
        
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node1->left->postorder_idx;
        Parameters_add(param, node1->left->distance);
        
        
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_distance( node1->left, Parameters_value(param, 0) );
#ifdef DEBUG_TOPOLOGY_SPR
        printf("** lnl %f (%f) %s %f\n", -lnl2, lnl,node1->left->name, Parameters_value(param, 0) );
#endif
        Parameters_pop(param);
        
        
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node1->right->postorder_idx;
        Parameters_add(param, node1->right->distance);
        
        
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_distance( node1->right, Parameters_value(param, 0) );
#ifdef DEBUG_TOPOLOGY_SPR
        printf("** lnl %f (%f) %s %f\n", -lnl2, lnl,node1->right->name, Parameters_value(param, 0) );
#endif
        Parameters_pop(param);
  
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node1->postorder_idx;
        Parameters_add(param, node1->distance);
        
        
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_distance( node1, Parameters_value(param, 0) );
#ifdef DEBUG_TOPOLOGY_SPR
        printf("** lnl %f (%f) %s %f\n", -lnl2, lnl,node1->name, Parameters_value(param, 0) );
#endif
        Parameters_pop(param);
        
        SingleTreeLikelihood_update_all_nodes(stlk);
        data->index_param = node2->postorder_idx;
        Parameters_add(param, node2->distance);
        
        
        status = opt_optimize( opt, param, &lnl2);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_distance( node2, Parameters_value(param, 0) );
#ifdef DEBUG_TOPOLOGY_SPR
        printf("** lnl %f (%f) %s %f\n", -lnl2, lnl,node2->name, Parameters_value(param, 0) );
#endif
		Parameters_pop(param);
		
        if ( -lnl2 - lnl < 0.01 ){
            lnl = -lnl2;
            break;
        }
        lnl = -lnl2;
    }
    
    SingleTreeLikelihood_update_all_nodes(stlk);
#ifdef DEBUG_TOPOLOGY_SPR
    printf("%f %f %f %f\n", Node_distance(node1), Node_distance(node1->left), Node_distance(node1->right), Node_distance(node2));
#endif
    //return stlk->calculate(stlk);
    return lnl;
}

bool is_SPR_valid( Node *prune, Node *graft ){
    if( prune == graft ) return false;
    
    Node *prune_parent = Node_parent(prune);
    Node *graft_parent = Node_parent(graft);
    
    if( Node_isroot( prune_parent ) ){
        
        // Adjacent edge
        // parent-child: no change in topology
        if( graft_parent == prune ) return false;
        
        // The sibling of prune is never evaluated (see for loop).
        // If we get there it means that the sibling of prune has a parent other than root
        
        // Adjacent edge (because of the fake rooting)
        // testing root == root
        if( Node_isroot( Node_parent( graft_parent ) ) ) return false;
        
        // If we get there it means that the sibling of prune has a grand parent other than root
        
        // NNI case
        // Skip if graft is the grand-nephew of prune
        // The equivalent NNIs will prune the grand newphews and regraft them on the reachable son of root
        if( Node_isroot( Node_parent(Node_parent( graft_parent ) ) ) ) return false;
        
    }
    else if( Node_isroot( graft_parent ) ){
        
        // Adjacent edge
        // parent-child: no change in topology
        if( prune_parent == graft ) return false;
        
        // Adjacent edge (because of the fake rooting)
        // testing root == root
        if( Node_isroot( Node_parent( prune_parent ) ) ) return false;
        
        // Adjacent edge
        // Child j of the root and its nephiew i
        if( Node_parent( prune_parent ) != NULL ){
            if( Node_isroot( Node_parent(Node_parent( prune_parent ) ) ) ) return false;
        }
        
    }
    else {
        
        // Adjacent edge
        // siblings: no change in topology
        if( Node_sibling(prune) == graft ) return false;
        
        // Adjacent edge
        // parent-child: no change in topology
        if( prune_parent == graft || graft_parent == prune ) return false;
        
//        // NNI case I
//        // Skip when prune is the grand child of graft (seprated by one edge)
//        // NNI #1: prune with its uncle graft
//        // NNI #2: sibling of prune with its uncle graft
//        if( Node_parent( prune_parent ) == graft ) return false;
//        
//        // NNI case II
//        // Skip when prune is the grand father of graft (seprated by one edge)
//        // The equivalent NNIs are taken care of in the NNIs described in case I
//        if( Node_parent( graft_parent ) == prune ) return false;
//        
//        // NNI case III
//        // Avoid repeating the same NNI
//        // We always prune a newphew and graft it on its uncle, not the other way around
//        if( Node_sibling(prune) == graft_parent ) return false;
//        
//        
//        // NNI case IV
//        // For cousins with root grand parent: only 1 node (out of 4: its sibling and 2 cousins) is pruned and regrafted
//        if( Node_isroot( Node_parent( prune_parent ) ) && Node_isroot( Node_parent( graft_parent )) ){
//            if( i > j || Node_left(prune_parent) == prune )return false;
//        }
    }
    return true;
}

// neighborhood size 2(n − 3)(2n − 7)
double spr_optimize_bl_parsimony( struct TopologyOptimizer * opt ){
    SingleTreeLikelihood *tlk = opt->tlk;
    
    
    double lnl = tlk->calculate(tlk);
    //double max = lnl;
    
    SingleTreeLikelihood *tlk2 = clone_SingleTreeLikelihood_share(tlk, true, false);
    
    Parsimony *parsimony = new_Parsimony(tlk2->sp, tlk2->tree);
    double score = parsimony->calculate(parsimony);
    if(tlk->opt.verbosity > 0){
        printf("\nParsimony score: %f\n\n", score);
    }
    
	Optimizer *opt_bl = new_Optimizer(OPT_BRENT);
    BrentData *data_brent = new_BrentData( tlk2 );
    data_brent->backup = dvector(1);
    data_brent->f = standard_loglikelihood_brent;
    if( tlk->use_upper ){
        data_brent->f = standard_loglikelihood_upper_brent;
    }
    
    opt_set_data(opt_bl, data_brent);
    opt_set_objective_function(opt_bl, optimize_brent_branch_length);
    opt_set_max_iteration(opt_bl, tlk->opt.bl.max_iteration);
    opt_set_tolx( opt_bl, tlk->opt.bl.tolx);
    
	Parameters *oneparameter = new_Parameters(1);
    
    double *branches = dvector(Tree_node_count(tlk->tree));
    
    int nNodes = Tree_node_count(tlk->tree);
    bool failed = false;
    
    double *scores = dvector(nNodes);
    int *prunes    = ivector(nNodes);
    int *grafts    = ivector(nNodes);
    int *ds        = ivector(nNodes);
    
    
    opt->moves = 0;
    
    int count = 0;
    
    int rounds = 0;
    
    while ( !failed ) {
        
        //nodes = Tree_get_nodes(tlk->tree, POSTORDER);
        Node *root = Tree_root(tlk2->tree);
        
        Tree_init_depth(tlk2->tree);
        
        
        count = 0;
        
        for ( int i = 0; i < nNodes; i++ ) {
            
            scores[count] = -1;
            
            // The node we want to prune
            //Node *prune = Tree_get_node(tlk2->tree, POSTORDER, i);
            Node *prune = Tree_node(tlk2->tree, i);
            
            // we don't prune the root node and its right child
            if ( Node_isroot(prune) || ( Node_isroot(Node_parent(prune)) && prune == Node_right(Tree_root(tlk2->tree))) ) continue;
            
            for ( int j = 0; j < nNodes; j++ ) {
                
                // same node
                if( i == j ) continue;
                
                //Node *graft = Tree_get_node(tlk2->tree, POSTORDER, j);
                Node *graft = Tree_node(tlk2->tree, j);
                
                // we don't graft on the root node and its right child
                if ( Node_isroot(graft) || ( Node_isroot(Node_parent(graft)) && graft == Node_right(Tree_root(tlk2->tree))) ) continue;
                
                if( Node_isroot( Node_parent(prune) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(graft) == prune ) continue;
                    
                    // The sibling of prune is never evaluated (see for loop).
                    // If we get there it means that the sibling of prune has a parent other than root
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(graft) ) == root ) continue;
                    
                    // If we get there it means that the sibling of prune has a grand parent other than root
                    
                    // NNI case
                    // Skip if graft is the grand-nephew of prune
                    // The equivalent NNIs will prune the grand newphews and regraft them on the reachable son of root
                    if( Node_parent(Node_parent( Node_parent(graft)) ) == root ) continue;
                    
                }
                else if( Node_isroot( Node_parent(graft) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft ) continue;
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(prune) ) == root ) continue;
                    
                    // Adjacent edge
                    // Child j of the root and its nephiew i
                    if( Node_parent( Node_parent(prune)) != NULL ){
                        if( Node_parent(Node_parent( Node_parent(prune)) ) == root ) continue;
                    }
                    
                }
                else {
                    
                    // Adjacent edge
                    // siblings: no change in topology
                    if( Node_sibling(prune) == graft ) continue;
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft || Node_parent(graft) == prune ) continue;
                    
                    // NNI case I
                    // Skip when prune is the grand child of graft (seprated by one edge)
                    // NNI #1: prune with its uncle graft
                    // NNI #2: sibling of prune with its uncle graft
                    if( Node_parent( Node_parent(prune) ) == graft ) continue;
                    
                    // NNI case II
                    // Skip when prune is the grand father of graft (seprated by one edge)
                    // The equivalent NNIs are taken care of in the NNIs described in case I
                    if( Node_parent( Node_parent(graft) ) == prune ) continue;
                    
                    // NNI case III
                    // Avoid repeating the same NNI
                    // We always prune a newphew and graft it on its uncle, not the other way around
                    if( Node_sibling(prune) == Node_parent(graft) ) continue;
                    
                    
                    // NNI case IV
                    // For cousins with root grand parent: only 1 node (out of 4: its sibling and 2 cousins) is pruned and regrafted
                    if( Node_isroot( Node_parent( Node_parent(prune)) ) && Node_isroot( Node_parent( Node_parent(graft) )) ){
                        if( i > j || prune->parent->left == prune )continue;
                    }
                }
                
                int d = Node_graph_distance(prune, graft);
                
                if( d > 10 ) continue;
                
                
                Node *regraft = Node_sibling(prune);
                bool rerooted = false;
                Node *n = prune;
                
                
                if( Node_isancestor(graft, prune) ){
                    if( !Node_isroot(Node_parent(prune)) ){
                        regraft = Node_parent(prune);
                    }
                    rerooted = true;
                    
                    // Get the location of the root (sister lineage at the root)
                    while (!Node_isroot(Node_parent(n)) ) n = Node_parent(n);
                    n = Node_sibling(n);
                }
                
                
                // SPR
                SPR_move(tlk2->tree, prune, graft);
                
                
                double spr_score = 0;
                Parsimony_update_all_nodes(parsimony);
                spr_score = parsimony->calculate(parsimony);
                
                /*
                 //if ( spr_score < score ) {
                 
                 double lnl2;
                 
                 if( rerooted ){
                 lnl2 = optimize_brent_branch_length_2b(tlk2, opt_bl, data_brent, oneparameter, graft, parent );
                 }
                 else {
                 //lnl2 = optimize_brent_branch_length_all(tlk2, opt_bl, data_brent, oneparameter, 1);
                 lnl2 = optimize_brent_branch_length_2b(tlk2, opt_bl, data_brent, oneparameter, graft, prune);
                 //lnl2 = optimize_brent_branch_length_3(tlk2, opt_bl, data_brent, oneparameter, parent, regraft);
                 }
                 
                 if(lnl2 > max ){
                 max = lnl2;
                 }
                 
                 if(lnl2 > lnl){
                 printf("%s %s lnl %f (%f) parsimony %f  %f %s\n\n",graft->name, prune->name, lnl2, lnl, spr_score, max, (lnl2 > lnl ? "*" : ""));
                 
                 }
                 
                 //                if( strcmp(Node_name(prune), "node73") == 0 && strcmp(Node_name(graft), "Avian|A/mallard/Sweden/57/2003|Sweden|H1N1_2003") == 0 ){
                 //                    printf("%s %s lnl %f (%f) %s\n\n",graft->name, parent->name, lnl2, lnl, (lnl2 > lnl ? "*" : ""));
                 //                }
                 
                 //}
                 */
                
                
                
                if(rerooted){
                    SPR_move(tlk2->tree, regraft, prune);
                    Tree_reroot(tlk2->tree, n); // this function uses Tree_update_topology. should use Tree_set_topology_changed instead
                }
                else {
                    SPR_move(tlk2->tree, prune, regraft);
                }
                
                Node *root  = Tree_root(tlk->tree);
                Node *root2 = Tree_root(tlk2->tree);
                
                if( strcmp(Node_name(root->left), Node_name(root2->left) ) != 0  ){
                    Node_rotate(root2);
                    Tree_set_topology_changed(tlk2->tree);
                }
                
                if ( spr_score < score  ) {
                    scores[count] = spr_score;
                    prunes[count] = Node_id(prune);
                    grafts[count] = Node_id(graft);
                    ds[count]     = d;
                }
            }
            
            if ( scores[count] > 0 ) {
                count++;
            }
            
        }
        
        
#ifdef DEBUG_TOPOLOGY_SPR
        fprintf(stderr, "%d potentials SPRs\n", count);
#endif
        
        if( count != 0 ){
            
            if(count > 1){
                sort_asc_spr_parsimony2(prunes, grafts, scores, ds, count);
            }
            
            
            // apply one by one until the score increases
            for ( int i = 0; i < count; i++ ) {
                Node *prune = Tree_node(tlk2->tree, prunes[i]);
                Node *graft = Tree_node(tlk2->tree, grafts[i]);
                
                if( prune == graft ) continue;
                
                if( Node_isroot(prune) || Node_isroot(graft)  ) continue;
                
                // Adjacent edge
                // siblings: no change in topology
                if( Node_sibling(prune) == graft ) continue;
                
                
                Node *prune_parent = Node_parent(prune);
                Node *graft_parent = Node_parent(graft);
                
                // Adjacent edge
                // parent-child: no change in topology
                if( prune_parent == graft || graft_parent == prune ) continue;
                
                bool rerooted = false;
                Node *n = prune;
                Node *regraft = Node_sibling(prune);
                
                
                if( Node_isancestor(graft, prune) ){
                    if( !Node_isroot(Node_parent(prune)) ){
                        regraft = Node_parent(prune);
                    }
                    rerooted = true;
                    
                    // Get the location of the root (sister lineage at the root)
                    while (!Node_isroot(Node_parent(n)) ) n = Node_parent(n);
                    n = Node_sibling(n);
                }
                else {
                    //Parsimony_update_node(parsimony, prune);
                    //Parsimony_update_node(parsimony, regraft);
                }
                
                
                SPR_move(tlk2->tree, Tree_node(tlk2->tree, prunes[i]), Tree_node(tlk2->tree, grafts[i]));
                
                double spr_score = 0;
                
                Parsimony_update_all_nodes(parsimony);
                spr_score = parsimony->calculate(parsimony);
                
#ifdef DEBUG_TOPOLOGY_SPR
                printf("Parsimony %f (%f) d(prune,graft)=%d\n", spr_score, scores[i], ds[i]);
#endif
                if ( spr_score > score ) {
                    
                    if(rerooted){
                        SPR_move(tlk2->tree, regraft, prune);
                        Tree_reroot(tlk2->tree, n);
                    }
                    else {
                        SPR_move(tlk2->tree, prune, regraft);
                    }
                    
                    Node *root  = Tree_root(tlk->tree);
                    Node *root2 = Tree_root(tlk2->tree);
                    
                    if( strcmp(Node_name(root->left), Node_name(root2->left) ) != 0  ){
                        Node_rotate(root2);
                        Tree_set_topology_changed(tlk2->tree);
                    }
                    break;
                }
                
                
                //                double lnl2 = 0;
                //                if( rerooted ){
                //                    lnl2 = optimize_brent_branch_length_2b(tlk2, opt_bl, data_brent, oneparameter, graft, prune_parent );
                //                }
                //                else {
                //                    lnl2 = optimize_brent_branch_length_2b(tlk2, opt_bl, data_brent, oneparameter, graft, prune);
                //                }
                //                printf(" LnL %f (%f)\n", lnl2, lnl);
                
                /*if( lnl2 < lnl ){
                 
                 if(rerooted){
                 SPR_move(tlk2->tree, regraft, prune);
                 Tree_reroot(tlk2->tree, n);
                 }
                 else {
                 SPR_move(tlk2->tree, prune, regraft);
                 }
                 
                 Node *root  = Tree_root(tlk->tree);
                 Node *root2 = Tree_root(tlk2->tree);
                 
                 if( strcmp(Node_name(root->left), Node_name(root2->left) ) != 0  ){
                 Node_rotate(root2);
                 }
                 
                 break;
                 }*/
                
                score = spr_score;
                //lnl   = lnl2;
                opt->moves++;
                
                SPR_move(tlk->tree, Tree_node(tlk->tree, prunes[i]), Tree_node(tlk->tree, grafts[i]));
            }
            
            
            Tree_copy_distances(tlk2->tree, tlk->tree);
            
            //Tree_print_newick(stdout, tlk2->tree, true);
            
            //lnl = optimize_brent_branch_length_all(tlk2, opt_bl, data_brent, oneparameter, 1);
            //printf(" LnL %f\n", lnl);
            
        }
        else {
            failed = true;
        }
        rounds++;
    }
    
    free_Parsimony(parsimony);
    free_BrentData(data_brent);
    free_Optimizer(opt_bl);
    free_Parameters(oneparameter);
    free_SingleTreeLikelihood_share(tlk2, true, false);
    
    free(branches);
    free(prunes);
    free(grafts);
    free(scores);
    free(ds);
    
    opt->best_lnl = lnl;
    
    return lnl;
}
//#define SPR_THREAD 1
//#define _OPENMP 1


// neighborhood size 2(n − 3)(2n − 7)
double spr_optimize_bl_parsimony_only( struct TopologyOptimizer * opt ){
	Parsimony* parsimony = opt->model->obj;
	Model* mtlk = opt->tlk;
	
	Tree *tree = parsimony->tree;
	Node *root  = NULL;
	
    double score = parsimony->calculate(parsimony);
    
    int total = (Tree_tip_count(tree)-3)*2*(2*Tree_tip_count(tree)-7);
    
	int max_distance = opt->max_distance;

    if(opt->verbosity > 0){
        printf("\nParsimony score: %f (%i)\n\n", score, total);
        printf("Maximum radius %d\n\n", max_distance);
    }
	
    int nNodes = Tree_node_count(tree);
    bool failed = false;
    
    double *scores = dvector(nNodes);
    int *prunes    = ivector(nNodes);
    int *grafts    = ivector(nNodes);
    int *ds        = ivector(nNodes);
    
    for ( int i = 0; i < nNodes; i++ ) {
        scores[i] = -1;
    }
    
    
    opt->moves = 0;
    
    int count = 0;
    
    int rounds = 0;
    
    while ( !failed && rounds < 5 ) {
		
        root = Tree_root(tree);
        
        Tree_init_depth(tree);
        
        
        count = 0;
        scores[count] = -1;
        int tot = 0;
		
//		Tree_print_newick(stderr, tree, true);
//		printf("\n");

        for ( int i = 0; i < nNodes; i++ ) {
      
            // The node we want to prune
            //Node *prune = Tree_get_node(tlk2->tree, POSTORDER, i);
            Node *prune = Tree_node(tree, i);
            
            // we don't prune the root node and its right child
            if ( Node_isroot(prune) || ( Node_isroot(Node_parent(prune)) && prune == Node_right(Tree_root(tree))) ) continue;
            
            for ( int j = 0; j < nNodes; j++ ) {
                
                // same node
                if( i == j ) continue;
                
                //Node *graft = Tree_get_node(tlk2->tree, POSTORDER, j);
                Node *graft = Tree_node(tree, j);
                
                
                //if( abs(Node_depth(prune)-Node_depth(graft)) > max_distance )continue;
                
                // we don't graft on the root node and its right child
                if ( Node_isroot(graft) || ( Node_isroot(Node_parent(graft)) && graft == Node_right(Tree_root(tree))) ) continue;
                
                if( Node_isroot( Node_parent(prune) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(graft) == prune ) continue;
                    
                    // The sibling of prune is never evaluated (see for loop).
                    // If we get there it means that the sibling of prune has a parent other than root
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(graft) ) == root ) continue;
                    
                    // If we get there it means that the sibling of prune has a grand parent other than root
                    
                    // NNI case
                    // Skip if graft is the grand-nephew of prune
                    // The equivalent NNIs will prune the grand newphews and regraft them on the reachable son of root
                    if( Node_parent(Node_parent( Node_parent(graft)) ) == root ) continue;
                    
                }
                else if( Node_isroot( Node_parent(graft) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft ) continue;
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(prune) ) == root ) continue;
                    
                    // Adjacent edge
                    // Child j of the root and its nephiew i
                    if( Node_parent( Node_parent(prune)) != NULL ){
                        if( Node_parent(Node_parent( Node_parent(prune)) ) == root ) continue;
                    }
                    
                }
                else {
                    
                    // Adjacent edge
                    // siblings: no change in topology
                    if( Node_sibling(prune) == graft ) continue;
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft || Node_parent(graft) == prune ) continue;
                    
                    // NNI case I
                    // Skip when prune is the grand child of graft (seprated by one edge)
                    // NNI #1: prune with its uncle graft
                    // NNI #2: sibling of prune with its uncle graft
                    if( Node_parent( Node_parent(prune) ) == graft ) continue;
                    
                    // NNI case II
                    // Skip when prune is the grand father of graft (seprated by one edge)
                    // The equivalent NNIs are taken care of in the NNIs described in case I
                    if( Node_parent( Node_parent(graft) ) == prune ) continue;
                    
                    // NNI case III
                    // Avoid repeating the same NNI
                    // We always prune a newphew and graft it on its uncle, not the other way around
                    if( Node_sibling(prune) == Node_parent(graft) ) continue;
                    
                    
                    // NNI case IV
                    // For cousins with root grand parent: only 1 node (out of 4: its sibling and 2 cousins) is pruned and regrafted
                    if( Node_isroot( Node_parent( Node_parent(prune)) ) && Node_isroot( Node_parent( Node_parent(graft) )) ){
                        if( i > j || prune->parent->left == prune )continue;
                    }
                }
                
                int d = Node_graph_distance(prune, graft);
                
                if( abs(prune->depth - graft->depth) > max_distance || d > max_distance ) continue;
                //printf("%s %s\n", prune->name, graft->name);
				
				Tree_save(tree);
                
                // SPR
                SPR_move(tree, prune, graft);
//				bool ok = check_tree_topology(tree);
//				if(!ok){
//					Tree_print_newick(stderr, tree, true);
//					printf("\n");
//					exit(123);
//				}
                double spr_score = parsimony->calculate(parsimony);
				
				Tree_revert(tree);
                
                if ( spr_score < score  ) {
                    scores[count] = spr_score;
                    prunes[count] = Node_id(prune);
                    grafts[count] = Node_id(graft);
                    ds[count]     = d;
                }
            }
            
            if ( scores[count] > 0 ) {
                count++;
                scores[count] = -1;
            }
            tot++;
            //max_distance = ( max_distance == 5 ? 5 : max_distance-1);
        }
		
		if(opt->verbosity > 0){
        	printf("%d potentials SPRs out of %d\n", count, tot);
		}
        
        if( count != 0 ){
            
            int accepted = 0;
            
            if(count > 1){
                sort_asc_spr_parsimony2(prunes, grafts, scores, ds, count);
            }
            
            
            // apply one by one until the score increases
            for ( int i = 0; i < count; i++ ) {
                Node *prune = Tree_node(tree, prunes[i]);
                Node *graft = Tree_node(tree, grafts[i]);
                
                if( prune == graft ) continue;
                
                if( Node_isroot(prune) || Node_isroot(graft)  ) continue;
                
                // Adjacent edge
                // siblings: no change in topology
                if( Node_sibling(prune) == graft ) continue;
                
                
                Node *prune_parent = Node_parent(prune);
                Node *graft_parent = Node_parent(graft);
                
                // Adjacent edge
                // parent-child: no change in topology
                if( prune_parent == graft || graft_parent == prune ) continue;
				
				// Adjacent edge
				// parent-child: no change in topology
				if( Node_isroot( prune_parent ) && (graft_parent == prune || Node_parent(graft_parent) == root) ) continue;
				if( Node_isroot( graft_parent ) && (prune_parent == graft || Node_parent(prune_parent) == root) ) continue;
				
				Tree_save(tree);
                SPR_move(tree, Tree_node(tree, prunes[i]), Tree_node(tree, grafts[i]));
                
                double spr_score = 0;
				
                spr_score = parsimony->calculate(parsimony);
                if(opt->verbosity > 1){
                    printf("Parsimony %f (%f) d(prune,graft)=%d max %f\n", spr_score, scores[i], ds[i], score);
                }
                // skip moves that increase the score
                if ( spr_score >= score ) {
					Tree_revert(tree);
                }
                else {
                    score = spr_score;
                    accepted++;
                    
                    opt->moves++;
                }
            }
            if(opt->verbosity > 0){
                printf("Parsimony score: %f #SPR: %d\n\n", score, accepted);
            }
        }
        else {
            failed = true;
        }
        rounds++;
    }
    
    free(prunes);
    free(grafts);
    free(scores);
    free(ds);
	
    opt->best_lnl = score;
	
	if (mtlk != NULL) {
		SingleTreeLikelihood* tlk = mtlk->obj;
		SingleTreeLikelihood_update_all_nodes(tlk);
		Tree_topology_changed(tlk->tree);
	}
    
    return score;
}

/* With threads, the ordering of scores can be different as moves with the same score will be added to the scores array
   in a different order. This can result in different parsimony scores with different number of threads. */

double spr_optimize_bl_parsimony_only_openmp( struct TopologyOptimizer * opt ){
    Parsimony* parsimony = opt->model->obj;
	Tree* tree = parsimony->tree;

	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
    //opt->threads = 1;
    Model **pool = malloc(opt->threads*sizeof(Model*));
    assert(pool);
    pool[0] = opt->model;

	for ( int i = 1; i < opt->threads; i++ ) {
		pool[i] = opt->model->clone(opt->model, hash);
		Hashtable_empty(hash);
	}
	free_Hashtable(hash);

    double score = pool[0]->logP(pool[0]);
    
    
    
    int total = (Tree_tip_count(tree)-3)*2*(2*Tree_tip_count(tree)-7);
    
    if(opt->verbosity > 0){
        printf("\nParsimony score: %f (%i)\n\n", score, total);
    }
    
    int nNodes = Tree_node_count(tree);
    bool failed = false;
    
    double *scores = dvector(nNodes);
    int *prunes    = ivector(nNodes);
    int *grafts    = ivector(nNodes);
    int *ds        = ivector(nNodes);
    
    for ( int i = 0; i < nNodes; i++ ) {
        scores[i] = -1;
    }
    
    
    opt->moves = 0;
    
    int count = 0;
    
    int rounds = 0;
    
	int max_distance = opt->max_distance;
	
    while ( !failed && rounds < 5 ) {
        
        
        
        count = 0;
        scores[count] = -1;
        int tot = 0;
        
#pragma omp parallel for num_threads(nthreads)
        for ( int i = 0; i < nNodes; i++ ) {
            
            int tid = 0;
#if defined (_OPENMP)
            tid = omp_get_thread_num();
#endif
            Parsimony *parsimony = pool[tid]->obj;
            Tree *tree = parsimony->tree;
            Node *root = Tree_root(tree);
            Tree_init_depth(tree);
            
            int local_score = -1;
            int local_prune;
            int local_graft;
            int local_d;
            
            // The node we want to prune
            Node *prune = Tree_node(tree, i);
            
            // we don't prune the root node and its right child
            if ( Node_isroot(prune) || ( Node_isroot(Node_parent(prune)) && prune == Node_right(Tree_root(tree))) ) continue;
            
            for ( int j = 0; j < nNodes; j++ ) {
                
                // same node
                if( i == j ) continue;
                
                //Node *graft = Tree_get_node(tlk2->tree, POSTORDER, j);
                Node *graft = Tree_node(tree, j);
                
                
                //if( abs(Node_depth(prune)-Node_depth(graft)) > max_distance )continue;
                
                // we don't graft on the root node and its right child
                if ( Node_isroot(graft) || ( Node_isroot(Node_parent(graft)) && graft == Node_right(Tree_root(tree))) ) continue;
                
                if( Node_isroot( Node_parent(prune) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(graft) == prune ) continue;
                    
                    // The sibling of prune is never evaluated (see for loop).
                    // If we get there it means that the sibling of prune has a parent other than root
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(graft) ) == root ) continue;
                    
                    // If we get there it means that the sibling of prune has a grand parent other than root
                    
                    // NNI case
                    // Skip if graft is the grand-nephew of prune
                    // The equivalent NNIs will prune the grand newphews and regraft them on the reachable son of root
                    if( Node_parent(Node_parent( Node_parent(graft)) ) == root ) continue;
                    
                }
                else if( Node_isroot( Node_parent(graft) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft ) continue;
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(prune) ) == root ) continue;
                    
                    // Adjacent edge
                    // Child j of the root and its nephiew i
                    if( Node_parent( Node_parent(prune)) != NULL ){
                        if( Node_parent(Node_parent( Node_parent(prune)) ) == root ) continue;
                    }
                    
                }
                else {
                    
                    // Adjacent edge
                    // siblings: no change in topology
                    if( Node_sibling(prune) == graft ) continue;
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft || Node_parent(graft) == prune ) continue;
                    
                    // NNI case I
                    // Skip when prune is the grand child of graft (seprated by one edge)
                    // NNI #1: prune with its uncle graft
                    // NNI #2: sibling of prune with its uncle graft
                    if( Node_parent( Node_parent(prune) ) == graft ) continue;
                    
                    // NNI case II
                    // Skip when prune is the grand father of graft (seprated by one edge)
                    // The equivalent NNIs are taken care of in the NNIs described in case I
                    if( Node_parent( Node_parent(graft) ) == prune ) continue;
                    
                    // NNI case III
                    // Avoid repeating the same NNI
                    // We always prune a newphew and graft it on its uncle, not the other way around
                    if( Node_sibling(prune) == Node_parent(graft) ) continue;
                    
                    
                    // NNI case IV
                    // For cousins with root grand parent: only 1 node (out of 4: its sibling and 2 cousins) is pruned and regrafted
                    if( Node_isroot( Node_parent( Node_parent(prune)) ) && Node_isroot( Node_parent( Node_parent(graft) )) ){
                        if( i > j || prune->parent->left == prune )continue;
                    }
                }
                
                int d = Node_graph_distance(prune, graft);
                
                if( d > max_distance ) continue;
                //printf("%s %s\n", prune->name, graft->name);
                
                Node *regraft = Node_sibling(prune);
                bool rerooted = false;
                Node *n = prune;
                
                
                if( Node_isancestor(graft, prune) ){
                    if( !Node_isroot(Node_parent(prune)) ){
                        regraft = Node_parent(prune);
                    }
                    rerooted = true;
                    
                    // Get the location of the root (sister lineage at the root)
                    while (!Node_isroot(Node_parent(n)) ) n = Node_parent(n);
                    n = Node_sibling(n);
                }
                
                
                // SPR
                bool simple = SPR_move(tree, prune, graft);
                
                
                double spr_score = 0;
                
                
                if( rerooted ){
                    Parsimony_update_node(parsimony, n);
                }
                else if( !simple ){
                    Parsimony_update_all_nodes(parsimony);
                }
                else {
                    Parsimony_update_node(parsimony, prune);
                    Parsimony_update_node(parsimony, regraft);
                }
                
                
                spr_score = parsimony->calculate(parsimony);
                
                
                if(rerooted){
                    simple = SPR_move(tree, regraft, prune);
                    Tree_reroot(tree, n);
                }
                else {
                    simple = SPR_move(tree, prune, regraft);
                }
                
                if( rerooted ){
                    Parsimony_update_node(parsimony, graft);
                    Parsimony_update_node(parsimony, n);
                }
                else if( !simple ){
                    Parsimony_update_all_nodes(parsimony);
                }
                else {
                    Parsimony_update_node(parsimony, prune);
                    Parsimony_update_node(parsimony, graft);
                }
                
                Node *root2 = Tree_root(tree);
                
                if( strcmp(Node_name(root->left), Node_name(root2->left) ) != 0  ){
                    Node_rotate(root2);
                    Tree_set_topology_changed(tree);
                }
                
                if ( spr_score < score  ) {
                    local_score = spr_score;
                    local_prune = Node_id(prune);
                    local_graft = Node_id(graft);
                    local_d     = d;
                }
            }
            
            
            if ( local_score > 0 ) {
#pragma omp critical
                {
                    scores[count] = local_score;
                    prunes[count] = local_prune;
                    grafts[count] = local_graft;
                    ds[count]     = local_d;
                    
                    count++;
                    //scores[count] = -1;
                }
            }
            tot++;
            //max_distance = ( max_distance == 5 ? 5 : max_distance-1);
        }
        //exit(1);
        
#ifdef DEBUG_TOPOLOGY_SPR
        fprintf(stderr, "%d potentials SPRs out of %d\n", count, tot);
#endif
        
        if( count != 0 ){
            
            int accepted = 0;
            
            if(count > 1){
                sort_asc_spr_parsimony2(prunes, grafts, scores, ds, count);
            }
            
            Parsimony *parsimony = pool[0]->obj;
            Tree *tree = parsimony->tree;
            Node *root = Tree_root(tree);
            
            // apply one by one until the score increases
            for ( int i = 0; i < count; i++ ) {
                Node *prune = Tree_node(tree, prunes[i]);
                Node *graft = Tree_node(tree, grafts[i]);
                
                if( prune == graft ) continue;
                
                if( Node_isroot(prune) || Node_isroot(graft)  ) continue;
                
                // Adjacent edge
                // siblings: no change in topology
                if( Node_sibling(prune) == graft ) continue;
                
                
                Node *prune_parent = Node_parent(prune);
                Node *graft_parent = Node_parent(graft);
                
                // Adjacent edge
                // parent-child: no change in topology
                if( prune_parent == graft || graft_parent == prune ) continue;
                
                bool rerooted = false;
                Node *n = prune;
                Node *regraft = Node_sibling(prune);
                
                
                if( Node_isancestor(graft, prune) ){
                    if( !Node_isroot(Node_parent(prune)) ){
                        regraft = Node_parent(prune);
                    }
                    rerooted = true;
                    
                    // Get the location of the root (sister lineage at the root)
                    while (!Node_isroot(Node_parent(n)) ) n = Node_parent(n);
                    n = Node_sibling(n);
                }
                else {
                    //Parsimony_update_node(parsimony, prune);
                    //Parsimony_update_node(parsimony, regraft);
                }
                
                
                SPR_move(tree, Tree_node(tree, prunes[i]), Tree_node(tree, grafts[i]));
                
                double spr_score = 0;
                
                
                Parsimony_update_all_nodes(parsimony);
                
                spr_score = parsimony->calculate(parsimony);
                
                printf("Parsimony %f (%f) d(prune,graft)=%d max %f\n", spr_score, scores[i], ds[i], score);
                
                // skip moves that increase the score
                if ( spr_score >= score ) {
                    
                    if(rerooted){
                        SPR_move(tree, regraft, prune);
                        Tree_reroot(tree, n);
                    }
                    else {
                        SPR_move(tree, prune, regraft);
                    }
                    
                    Node *root2 = Tree_root(tree);
                    
                    if( strcmp(Node_name(root->left), Node_name(root2->left) ) != 0  ){
                        Node_rotate(root2);
                        Tree_set_topology_changed(tree);
                    }
                    //break;
                }
                else {
#ifdef DEBUG_TOPOLOGY_SPR
                    printf("Parsimony %f (%f) d(prune,graft)=%d max %f\n", spr_score, scores[i], ds[i], score);
#endif
                    
                    Parsimony_update_all_nodes(pool[0]->obj);
//                    for( int j = 1; j < opt->threads; j++ ){
//                        SPR_move(pool[j]->tree, Tree_node(pool[j]->tree, prunes[i]), Tree_node(pool[j]->tree, grafts[i]));
//                        Parsimony_update_all_nodes(pool[j]->obj);
//                    }
					
                    score = spr_score;
                    accepted++;

                    opt->moves++;
                }
            }
            if(opt->verbosity > 0){
                printf("Parsimony score: %f #SPR: %d\n\n", score, accepted);
            }
        }
        else {
            failed = true;
        }
        rounds++;
    }
    
    free(prunes);
    free(grafts);
    free(scores);
    free(ds);
	
    
    if(opt->threads > 1 ){
        for ( int i = 1; i < opt->threads; i++ ) {
			pool[i]->free(pool[i]);
        }
    }
    free(pool);

    opt->best_lnl = score;
    
    return score;
}

void reorder_topology(Node* node, Node* backup){
	if(!Node_isleaf(backup)){
		if (strcmp(node->left->name, backup->left->name) != 0) {
			Node_rotate(node);
		}
		reorder_topology(node->left, backup->left);
		reorder_topology(node->right, backup->right);
	}
}

double optimize_tripod(Model* model, Node* centralNode){
	SingleTreeLikelihood* tlk = model->obj;
	
	double spr_score = 0;
	
	Optimizer* opt_3bl = new_Optimizer(OPT_BRENT);
	opt_set_data(opt_3bl, model);
	opt_set_objective_function(opt_3bl, model_negative_logP);
	Parameters *parameters = new_Parameters(1);
	
	tlk->use_upper = false;
//	SingleTreeLikelihood_update_all_nodes(tlk);
	double lnl2 = model->logP(model);
	
	Node* n = centralNode;
	while (n != NULL) {
		tlk->update_nodes[Node_id(n)] = true;
		n = n->parent;
	}
	tlk->update_nodes[Node_id(centralNode->left)] = true;
	tlk->update_nodes[Node_id(centralNode->right)] = true;
	SingleTreeLikelihood_update_uppers2(tlk);
	
	tlk->tripod = true;
	
	int idLeft = Node_id(centralNode->left);
	int idRight = Node_id(centralNode->right);
	int idCentral = Node_id(centralNode);
	
	do {
		// optimize top node branch
		Parameters_pop(parameters);
		Parameters_add(parameters, centralNode->distance);
		tlk->node_upper = centralNode;
		opt_result status = opt_maximize( opt_3bl, parameters, &spr_score);
		if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
//		if(logit)printf("%s %f %e %e %e\n", centralNode->name, spr_score, centralNode->distance->value, centralNode->left->distance->value, centralNode->right->distance->value);
		// make sure its matrices are updated
		SingleTreeLikelihood_update_Q(tlk, centralNode);
		// update upper of left node
		tlk->update_partials(tlk, tlk->upper_partial_indexes[idLeft],
							 tlk->upper_partial_indexes[idCentral], idCentral,
							 idRight, idRight);
		
		// optimize left node
		Parameters_pop(parameters);
		Parameters_add(parameters, centralNode->left->distance);
		tlk->node_upper = centralNode->left;
		status = opt_maximize( opt_3bl, parameters, &spr_score);
		if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
//		if(logit)printf("%s %f %e %e %e\n", centralNode->left->name, spr_score, centralNode->distance->value, centralNode->left->distance->value, centralNode->right->distance->value);
		
		SingleTreeLikelihood_update_Q(tlk, centralNode->left);
		// update upper of right node
		tlk->update_partials(tlk, tlk->upper_partial_indexes[idRight],
							 tlk->upper_partial_indexes[idCentral], idCentral,
							 idLeft, idLeft);
		
		// optimize right node
		Parameters_pop(parameters);
		Parameters_add(parameters, centralNode->right->distance);
		tlk->node_upper = centralNode->right;
		status = opt_maximize( opt_3bl, parameters, &spr_score);
		if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");

		if (fabs(spr_score-lnl2) < 0.1) {
			break;
		}
		SingleTreeLikelihood_update_Q(tlk, centralNode->right);
		// update lower of central node for next round
		tlk->update_partials(tlk, idCentral,
							 idRight, idRight,
							 idLeft, idLeft);
		lnl2 = spr_score;
	} while (1);
	free_Optimizer(opt_3bl);
	tlk->tripod = false;
	tlk->use_upper = false;
	free_Parameters(parameters);
	
//	SingleTreeLikelihood_update_one_node(tlk, centralNode);
//	SingleTreeLikelihood_update_one_node(tlk, centralNode->left);
//	SingleTreeLikelihood_update_one_node(tlk, centralNode->right);
	
	return spr_score;
}

double spr_optimize_bl_openmp( struct TopologyOptimizer * opt ){

	double lnl = opt->model->logP(opt->model);
	
	if(opt->failures >= opt->max_failures) return lnl;
	
	unsigned nthreads = opt->threads;
	
	Model** pool = malloc(sizeof(Model*)*nthreads);
	Model** tlks = malloc(sizeof(Model*)*nthreads);
	pool[0] = opt->model;
	tlks[0] = opt->tlk;
	
	if(nthreads > 1){
		Hashtable* hash = new_Hashtable_string(10);
		hashtable_set_key_ownership( hash, false );
		hashtable_set_value_ownership( hash, false );
		for (size_t i = 1; i < nthreads; i++) {
			pool[i] = opt->model->clone(opt->model, hash);
			tlks[i] = Hashtable_get(hash, opt->tlk->name); // not ref++
			Hashtable_empty(hash);
		}
		free_Hashtable(hash);
	}
	
    int nNodes = Tree_node_count(((SingleTreeLikelihood*)(opt->tlk->obj))->tree);
	
    bool failed = false;
	int max_distance = opt->max_distance;
    
    double *scores = dvector(nNodes);
    int *prunes    = ivector(nNodes);
    int *grafts    = ivector(nNodes);
    double *d1 = dvector(nNodes);
	double *d2 = dvector(nNodes);
	double *d3 = dvector(nNodes);
	
    for ( int i = 0; i < nNodes; i++ ) {
        scores[i] = INFINITY;
    }
	
	if (opt->verbosity > 0) {
    	printf("\nStart SPR optimization LnL: %f\n\n", lnl);
	}
	
    opt->moves = 0;
	
    int count = 0;
    
    int rounds = 0;
	
	bool use_rescaling = SingleTreeLikelihood_rescaling(tlks[0]->obj);
	time_t start_time, end_time;
	double diff_time;
	
    while ( !failed && rounds < 1 ) {
		int potential = 0;
        count = 0;
        scores[count] = INFINITY;
        int tot = 0;
		
		time(&start_time);
		
		for (size_t i = 0; i < nthreads; i++) {
			Tree_init_depth(((SingleTreeLikelihood*)tlks[i]->obj)->tree);
		}
        
#pragma omp parallel for num_threads(nthreads)
        for ( int i = 0; i < nNodes; i++ ) {
            int tid = 0;
#if defined (_OPENMP)
            tid = omp_get_thread_num();
#endif
            SingleTreeLikelihood *tlk = tlks[tid]->obj;
			
			if (!use_rescaling && SingleTreeLikelihood_rescaling(tlk)) {
				SingleTreeLikelihood_use_rescaling(tlk, false);
				SingleTreeLikelihood_update_all_nodes(tlk);
			}
			
            Tree *tree = tlk->tree;
            Node *root = Tree_root(tree);
			Node* right = Node_right(Tree_root(tree));
            double local_score = INFINITY;
            int local_prune;
            int local_graft;
			double local_d1;
			double local_d2;
			double local_d3;
            
            // The node we want to prune
            Node *prune = Tree_node(tree, i);
            
            // we don't prune the root node and its right child
            if ( Node_isroot(prune) || ( Node_isroot(Node_parent(prune)) && prune == right) ) continue;
            
            for ( int j = 0; j < nNodes; j++ ) {
                
                // same node
                if( i == j ) continue;
                
                Node *graft = Tree_node(tree, j);
                
                
                //if( abs(Node_depth(prune)-Node_depth(graft)) > max_distance )continue;
                
                // we don't graft on the root node and its right child
                if ( Node_isroot(graft) || ( Node_isroot(Node_parent(graft)) && graft == right) ) continue;
                
                if( Node_isroot( Node_parent(prune) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(graft) == prune ) continue;
                    
                    // The sibling of prune is never evaluated (see for loop).
                    // If we get there it means that the sibling of prune has a parent other than root
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(graft) ) == root ) continue;
                    
                    // If we get there it means that the sibling of prune has a grand parent other than root
                    
                    // NNI case
                    // Skip if graft is the grand-nephew of prune
                    // The equivalent NNIs will prune the grand newphews and regraft them on the reachable son of root
                    if( Node_parent(Node_parent( Node_parent(graft)) ) == root ) continue;
                    
                }
                else if( Node_isroot( Node_parent(graft) ) ){
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft ) continue;
                    
                    // Adjacent edge (because of the fake rooting)
                    // testing root == root
                    if( Node_parent( Node_parent(prune) ) == root ) continue;
                    
                    // Adjacent edge
                    // Child j of the root and its nephiew i
                    if( Node_parent( Node_parent(prune)) != NULL ){
                        if( Node_parent(Node_parent( Node_parent(prune)) ) == root ) continue;
                    }
                    
                }
                else {
                    
                    // Adjacent edge
                    // siblings: no change in topology
                    if( Node_sibling(prune) == graft ) continue;
                    
                    // Adjacent edge
                    // parent-child: no change in topology
                    if( Node_parent(prune) == graft || Node_parent(graft) == prune ) continue;
                    
                    // NNI case I
                    // Skip when prune is the grand child of graft (seprated by one edge)
                    // NNI #1: prune with its uncle graft
                    // NNI #2: sibling of prune with its uncle graft
                    if( Node_sibling( Node_parent(prune) ) == graft || Node_parent( Node_parent(prune) ) == graft ) continue;
                    
                    // NNI case II
                    // Skip when prune is the grand father of graft (seprated by one edge)
                    // The equivalent NNIs are taken care of in the NNIs described in case I
                    if( Node_parent( Node_parent(graft) ) == prune ) continue;
                    
                    // NNI case III
                    // Avoid repeating the same NNI
                    // We always prune a newphew and graft it on its uncle, not the other way around
                    if( Node_sibling(prune) == Node_parent(graft) ) continue;
                    
                    
                    // NNI case IV
                    // For cousins with root grand parent: only 1 node (out of 4: its sibling and 2 cousins) is pruned and regrafted
                    if( Node_isroot( Node_parent( Node_parent(prune)) ) && Node_isroot( Node_parent( Node_parent(graft) )) ){
                        if( i > j || prune->parent->left == prune )continue;
                    }
                }
                int d = Node_graph_distance(prune, graft);
                
                if( d > max_distance ) continue;

#ifdef DEBUG_TOPOLOGY_SPR
				printf("%s [%d] %s [%d] ",graft->name, Node_depth(graft), prune->name, Node_depth(prune));
				Tree_print_newick(stdout, tree, true);
				printf("\n");
#endif
				Tree_save(tree);
				
                // SPR
                Node* centralNode = SPR_move(tree, prune, graft);
				
				double spr_score = optimize_tripod(pool[tid], centralNode);
				
#ifdef DEBUG_TOPOLOGY_SPR
                    printf("%s %s lnl %f (%f) %d %s ",graft->name, prune->name, spr_score, lnl, d, (spr_score > lnl ? "*" : ""));
				Tree_print_newick(stdout, tree, true);
				printf("\n");
#endif
				
                if ( (isinf(local_score) && spr_score > lnl) || (!isinf(local_score) && spr_score > local_score)  ) {
                    local_score = spr_score;
                    local_prune = Node_id(prune);
                    local_graft = Node_id(graft);
					local_d1 = Node_distance(centralNode);
					local_d2 = Node_distance(centralNode->left);
					local_d3 = Node_distance(centralNode->right);
                }
				tot++;
				
				Tree_revert(tree);
            }
            
            
            if (!isinf(local_score)) {
#pragma omp critical
                {
					if (isinf(scores[count])) {
						potential++;
					}
                    scores[count] = local_score;
                    prunes[count] = local_prune;
                    grafts[count] = local_graft;
					d1[count] = local_d1;
					d2[count] = local_d2;
					d3[count] = local_d3;
                    
                    count++;
                }
            }
        }
		
		
		if (opt->verbosity > 0) {
			time(&end_time);
			diff_time = difftime(end_time, start_time);
        	printf("%d potential SPRs out of %d [%f]\n", potential, tot, diff_time);
		}
		
        if( count != 0 ){
			double_int_pair_t** pairs = malloc(count*sizeof(double_int_pair_t*));
			for (int i = 0; i < count; i++) {
				pairs[i] = malloc(sizeof(double_int_pair_t));
				pairs[i]->index = i;
				pairs[i]->value = scores[i];
			}
			
            if(count > 1){
				qsort(pairs, count, sizeof(struct double_int_pair_t*), cmp_double_int_pair_desc);
            }
            
            SingleTreeLikelihood* tlk = opt->tlk->obj;
            Tree *tree = tlk->tree;
			Node* root = Tree_root(tree);
			
//			Tree_print_newick(stdout, tree, true);
//			printf("\n");
			
			int c = 0;
            // apply one by one until the score increases
            for ( int i = 0; i < count; i++ ) {
                Node *prune = Tree_node(tree, prunes[pairs[i]->index]);
                Node *graft = Tree_node(tree, grafts[pairs[i]->index]);
                
//                printf("%s %s %f\n", prune->name, graft->name, scores[i]);
//				Tree_print_newick(stdout, tree, true);
//				printf("\n");
				
                if( prune == graft ) continue;
                
                if( Node_isroot(prune) || Node_isroot(graft)  ) continue;
                
                // Adjacent edge
                // siblings: no change in topology
                if( Node_sibling(prune) == graft ) continue;
                
                
                Node *prune_parent = Node_parent(prune);
                Node *graft_parent = Node_parent(graft);
                
                // Adjacent edge
                // parent-child: no change in topology
                if( prune_parent == graft || graft_parent == prune ) continue;
				
				// Adjacent edge
				// parent-child: no change in topology
				if( Node_isroot( prune_parent ) && (graft_parent == prune || Node_parent(graft_parent) == root) ) continue;
				if( Node_isroot( graft_parent ) && (prune_parent == graft || Node_parent(prune_parent) == root) ) continue;
				
				Tree_save(tree);
				
				Node* centralNode = SPR_move(tree, prune, graft);

				Node_set_distance(centralNode, d1[pairs[i]->index]);
				Node_set_distance(centralNode->left, d2[pairs[i]->index]);
				Node_set_distance(centralNode->right, d3[pairs[i]->index]);
				
				double spr_score = scores[pairs[i]->index];
				
				// no need to reoptimize the first one since it will give the same estimates
				if(i != 0){
					spr_score = optimize_tripod(pool[0], centralNode);
				}
				
#ifdef DEBUG_TOPOLOGY_SPR
                printf("Parsimony %f (%f) d(prune,graft)=%d max %f [%s -> %s]\n", spr_score, scores[pairs[i]->index], ds[pairs[i]->index], lnl, Node_name(prune), Node_name(graft));
				Tree_print_newick(stdout, tree, true);
				printf("\n");
#endif
                // skip moves that increase the score
                if ( spr_score < lnl ) {
					Tree_revert(tree);
                }
                else {
#ifdef DEBUG_TOPOLOGY_SPR
                    printf("Parsimony %f (%f) d(prune,graft)=%d max %f\n", spr_score, scores[i], ds[i], lnl);
#endif
                    SingleTreeLikelihood_update_all_nodes(tlks[0]->obj);
                    for( int j = 1; j < opt->threads; j++ ){
						SingleTreeLikelihood* tlkj = (SingleTreeLikelihood*)tlks[j]->obj;
                        SPR_move(tlkj->tree, Tree_node(tlkj->tree, prunes[pairs[i]->index]), Tree_node(tlkj->tree, grafts[pairs[i]->index]));
                        SingleTreeLikelihood_update_all_nodes(tlkj);
                    }
                    
                    lnl = spr_score;
					c++;
					opt->moves++;
                }
            }

			if (opt->verbosity > 0) {
				time(&end_time);
				diff_time = difftime(end_time, start_time);
				printf("%d SPRs LnL: %f [%f]\n", c, lnl, diff_time);
			}
			
			for (int i = 0; i < count; i++) {
				free(pairs[i]);
			}
			free(pairs);
        }
        else {
            failed = true;
        }
        rounds++;
    }

    free(prunes);
    free(grafts);
    free(scores);
	free(d1);
	free(d2);
	free(d3);
	
	
	for ( int i = 1; i < opt->threads; i++ ) {
		pool[i]->free(pool[i]);
	}
	
	free(pool);
	free(tlks);
	
    opt->best_lnl = lnl;
	
	if (opt->moves == 0) {
		opt->failures++;
	}
	else{
		opt->failures = 0;
	}
	
    return lnl;
}


