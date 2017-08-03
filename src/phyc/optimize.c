/*
 *  optimimize.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/12/10.
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


#include "optimize.h"

#include <math.h>
#include <float.h>
#include <assert.h>
#include <signal.h>

#include "treelikelihood.h"
#include "parameters.h"
#include "matrix.h"
#include "tree.h"
#include "node.h"
#include "optimizer.h"
#include "mathconstant.h"
#include "random.h"

#include "exponential.h"
#include "lognormal.h"

#include "optapprox.h"
#include "topologyopt.h"

#include "treeio.h"



#pragma mark Private functions declaration

static double optimize_brent_heights_all_threshold( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param );
static double _brent_optimize_height_threshold( Parameters *params, double *grad, void *data );


static double _brent_optimize_frequency( Parameters *params, double *grad, void *data );
static double _brent_optimize_relative_rate( Parameters *params, double *grad, void *data );
static double _brent_optimize_sm_rates( Parameters *params, double *grad, void *data );
static double _brent_optimize_rate( Parameters *params, double *grad, void *data );
static double _brent_optimize_height( Parameters *params, double *grad, void *data );

static double _brent_optimize_scale_heights_rates( Parameters *params, double *grad, void *data );
static double _brent_optimize_scale_heights_rates2( Parameters *params, double *grad, void *data );
static double _brent_optimize_scale_heights( Parameters *params, double *grad, void *data );
static double _brent_optimize_scale_root_height_rate( Parameters *params, double *grad, void *data );
static double _brent_optimize_slide_height( Parameters *params, double *grad, void *data );
static double _brent_optimize_scale_heights_polynomial( Parameters *params, double *grad, void *data );

static double _brent_optimize_add_heights( Parameters *params, double *grad, void *data );// can add negative number too

static double _brent_optimize_branch_length_strict( Parameters *params, double *grad, void *data );


static double _cg_optimize_branch_length( Parameters *params, double *grad, void *data );
static double _cg_optimize_frequencies( Parameters *params, double *grad, void *data );
static double _cg_optimize_rel_rates( Parameters *params, double *grad, void *data );
static double _cg_optimize_rates( Parameters *params, double *grad, void *data );
static double _cg_optimize_all( Parameters *params, double *grad, void *data );
double _cg_optimize_sm_rates( Parameters *params, double *grad, void *data );
static double _cg_optimize_heightsRate( Parameters *params, double *grad, void *data );
static double _cg_optimize_heightsRate2( Parameters *params, double *grad, void *data );

static double optimize_brent_frequencies_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param );
static double optimize_brent_relative_rate_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param );
static double optimize_brent_sm_rates_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param );

static double optimize_scale_root_rate_height_strict( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk );
static double optimize_scale_rate_height_all( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk );
static double optimize_scale_heights_polynomial( SingleTreeLikelihood *tlk, Optimizer *optimizer, BrentData *data, Parameters *coefficient, double templk );
static double optimize_brent_heights_strict_local( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param );


static double _cg_optimize_scale_rate_heights_polynomial( Parameters *params, double *grad, void *data );



#pragma mark -
#pragma mark Objective Functions


double standard_loglikelihood_brent( void *data ){
    SingleTreeLikelihood *tlk = ((BrentData*)data)->tlk;
    return tlk->calculate(tlk);
}

double standard_loglikelihood_brent2( void *data ){
    SingleTreeLikelihood *tlk = ((BrentData*)data)->tlk;
    double prior = 0;
    for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
        Node *n = Tree_node(tlk->tree,i);
        if(n == Tree_root(tlk->tree) || n == Node_right(Tree_root(tlk->tree)) ) continue;
        prior += log(10) - 10*Node_distance(n);
    }
    return tlk->calculate(tlk)+prior;
}

// For calculation using lower and upper likelihoods
double standard_loglikelihood_upper_brent( void *data ){
    BrentData *brent = (BrentData*)data;
    SingleTreeLikelihood *tlk = brent->tlk;
    return tlk->calculate_upper(tlk, Tree_node(tlk->tree, brent->index_param));
}

double standard_loglikelihood_upper_brent2( void *data ){
    BrentData *brent = (BrentData*)data;
    SingleTreeLikelihood *tlk = brent->tlk;
    double prior = 0;
    for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
        Node *n = Tree_node(tlk->tree,i);
        if(n == Tree_root(tlk->tree) || n == Node_right(Tree_root(tlk->tree)) ) continue;
        prior += log(10) - 10*Node_distance(n);
    }
    return tlk->calculate_upper(tlk, Tree_node(tlk->tree, brent->index_param))+prior;
}

double standard_loglikelihood_multivariate( void *data ){
    SingleTreeLikelihood *tlk = ((MultivariateData*)data)->tlk;
    return tlk->calculate(tlk);
}


#pragma mark -
#pragma mark Misc

#ifdef TIMETEST

static void get_earliest(Node *node, double *earliest ){
    if(Node_isleaf(node)){
        *earliest = dmax(*earliest, Node_height(node));
    }
    else {
        get_earliest(Node_left(node), earliest);
        get_earliest(Node_right(node), earliest);
    }
}

// from height to timeParameter
void set_times(Node *node){
    double earliest,t;
    Node *kid = NULL;
    
    if( Node_isroot(node)){
        Node_set_t(node, Node_height(node));
    }
    if( !Node_isleaf(node)){
        kid = Node_left(node);
        if( !Node_isleaf(kid) ){
            earliest = 0;
            get_earliest(kid, &earliest);
            kid->timeParameter->value = (Node_height(kid) - earliest)/(Node_height(node) - earliest);
            //printf("%s scaler %e kid height %e parent height %e diff %e\n", Node_name(kid), kid->timeParameter->value, Node_height(kid), Node_height(node), (Node_height(node)-Node_height(kid)) );
        }
        set_times(Node_left(node));
        
        kid = Node_right(node);
        if( !Node_isleaf(kid) ){
            earliest = 0;
            get_earliest(kid, &earliest);
            kid->timeParameter->value = (Node_height(kid) - earliest)/(Node_height(node) - earliest);
        }
        set_times(Node_right(node));
    }
}

static void strict_heights_to_parameters( Node *node, Parameters *params ){
    
    if( !Node_isleaf(node) ){
        strict_heights_to_parameters(Node_left(node), params);
        strict_heights_to_parameters(Node_right(node), params);
    }
    
    if( !Node_isroot(node) && Node_right(Node_parent(node)) == node ) Parameters_add(params, node->timeParameter);
}
static void strict_parameters_to_heights( Node *node, Parameters *params, int *index ){
    
    
    if( !Node_isleaf(node) ){
        strict_parameters_to_heights(Node_left(node), params, index);
        strict_parameters_to_heights(Node_right(node), params, index);
    }
    
    
    if( !Node_isroot(node) &&  Node_right(Node_parent(node)) == node ){
        Node_set_t(node, Parameters_value(params, *index));
        (*index)++;
    }
}

void calculate_constrained(Node *node, BranchModel *bm){
    
    if( !Node_isleaf(node) ){
        calculate_constrained(Node_left(node), bm);
        calculate_constrained(Node_right(node), bm);
    }
    
    // left node is relative to right node
    if( !Node_isroot(node) && Node_left(Node_parent(node)) == node ){
        Node *sibling = Node_sibling(node);
        double bl = 0;
        while ( !Node_isleaf(sibling) ) {
            bl += Node_t(sibling);
            sibling = Node_right(sibling);
        }
        bl += Node_t(sibling);
        
        Node *node2 = node;
        while ( !Node_isleaf(node2) ) {
            node2 = Node_right(node2);
            bl -= Node_t(node2);
        }
        bl *= bm->get(bm, node);
        bl -= bm->get(bm, node) * (Node_height(node2) - Node_height(sibling));
        bl /= bm->get(bm, node);
        //printf("%s bl %e t %e diff %e\n",Node_name(node), bl, Node_t(node),Node_t(node)-bl);
        Node_set_t(node,bl);
    }
}
#endif

static double _correct_height( double h, double coef ){
    return (-1+sqrt(1+h*4*coef))/(coef*2);
}

static void _scale_heights_poly( Node *node, const double scaler ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( Node_isroot(node) || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
		_scale_heights_poly(node->left, scaler);
		_scale_heights_poly(node->right, scaler);
		
		double height = Node_height( node );
        double corrected_height = _correct_height(height, scaler);
		double max_height_son = dmax( Node_height( Node_left(node)), Node_height( Node_right(node) ) );
		double bl = (height  - max_height_son)/(height/corrected_height) + max_height_son;
		Node_set_height(node, bl);
        //printf("%f %f %f\n",height, (height/corrected_height), bl);
		
	}
}


// never called on the root
static void _get_lower_bound2( const Node *theNode, Node *node, double *lower, double threshold ){
    if( node == NULL )return;
    
    if( node!= theNode && (Node_distance(node) > threshold ) ){
        *lower = dmax( Node_height(node), *lower );
    }
    else {
        _get_lower_bound2(theNode, Node_left(node),  lower, threshold);
        _get_lower_bound2(theNode, Node_right(node), lower, threshold);
    }
    
}
// never called on the root
static void _get_lower_bound3( const Node *theNode, Node *node, double *lower, double threshold ){
    if( node == NULL )return;
    
    if( !Node_isleaf(node) ){
        if( Node_distance(Node_left(node)) <= threshold  || Node_distance(Node_right(node)) <= threshold ){
            *lower = dmax( Node_height(Node_left(node)), *lower );
            *lower = dmax( Node_height(Node_right(node)), *lower );
            
            if( Node_distance(Node_left(node)) <= threshold )  _get_lower_bound3(theNode, Node_left(node),  lower, threshold);
            if( Node_distance(Node_right(node)) <= threshold ) _get_lower_bound3(theNode, Node_right(node), lower, threshold);
        }
    }
    
}

// When we call this function node is the node that is going to also change its descendents if their distances are under the threshold.
// The function finds recursively the first node that has a child that is above the threshold and save its height or the hieght of the sibling
// which ever is the highest
static double _get_lower_bound( Node *node, double threshold, int *level ){
    
    double right = Node_height(Node_right(node));
    double left  = Node_height(Node_left(node));
    
    // Both chirldren are under the threshold
    if( (!Node_isleaf(Node_left(node))  && Node_distance(Node_left(node)) <= threshold )
       && ( !Node_isleaf(Node_right(node)) && Node_distance(Node_right(node)) <= threshold ) ){
        left  = _get_lower_bound(Node_left(node), threshold, level);
        right = _get_lower_bound(Node_right(node), threshold, level);
        (*level)++;
    }
    // Only left node is under the threshold and it is higher than its sibling. Visit the left node descendents
    else if( (!Node_isleaf(Node_left(node))  && Node_distance(Node_left(node)) <= threshold) && (Node_height(Node_left(node)) > Node_height(Node_right(node)) ) ){
        left  = _get_lower_bound(Node_left(node), threshold, level);
        (*level)++;
    }
    // Only righ node is under the threshold and it is higher than its sibling. Visit the right node descendents
    else if( (!Node_isleaf(Node_right(node))  && Node_distance(Node_right(node)) <= threshold) && (Node_height(Node_right(node)) > Node_height(Node_left(node)) ) ){
        right = _get_lower_bound(Node_right(node), threshold, level);
        (*level)++;
    }
    return dmax(left, right);
    
}

// When we call this function node is the node that is going to also change its descendents if their distances are under the threshold.
// Leaves are fixed so we do not modify them
static void _move_and_update_nodes( SingleTreeLikelihood *tlk, Node *node, double height, double threshold, int *level ){
    if( node == NULL ) return;
    
    if( !Node_isleaf(node) ){
        if( !Node_isleaf(Node_left(node))  && Node_distance(Node_left(node))  <= threshold ) _move_and_update_nodes(tlk, Node_left(node),  height, threshold, level);
        if( !Node_isleaf(Node_right(node)) && Node_distance(Node_right(node)) <= threshold ) _move_and_update_nodes(tlk, Node_right(node), height, threshold, level);
        
        /*if( (!Node_isleaf(Node_left(node))  && Node_distance(Node_left(node)) <= threshold )
           && ( !Node_isleaf(Node_right(node)) && Node_distance(Node_right(node)) <= threshold ) ){
            _move_and_update_nodes(tlk, Node_left(node),  height, threshold, level);
            _move_and_update_nodes(tlk, Node_right(node), height, threshold, level);
            (*level)++;
        }*/
        
        Node_set_height(node, height);
        SingleTreeLikelihood_update_three_nodes(tlk, node );
    }
    
//    if( !Node_isroot(node) ){
//        if( !Node_isleaf(Node_left(node)) && Node_distance(Node_left(node)) <= threshold ){
//            _move_and_update_nodes(tlk, Node_left(node), height, threshold);
//        }
//        if( !Node_isleaf(Node_right(node)) && Node_distance(Node_right(node)) <= threshold ){
//            _move_and_update_nodes(tlk, Node_right(node), height, threshold);
//        }
//    }
//    else{
//        if( !Node_isleaf(node) && Node_time_elapsed(Node_left(node)) < 1 ){
//            _move_and_update_nodes(tlk, Node_left(node), height, threshold);
//        }
//        if( !Node_isleaf(node) && Node_time_elapsed(Node_right(node)) < 1 ){
//            _move_and_update_nodes(tlk, Node_right(node), height, threshold);
//        }
//    }
    
//    Node_set_height(node, height);
//    SingleTreeLikelihood_update_three_nodes(tlk, node );
}

static void _get_max_down( Node *node, double *value ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( Node_isroot(node) || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
		//fprintf(stderr, "%s\n", node->name);
		if( Node_isleaf(Node_left(node)) ) {
			*value = dmin(*value, node->height->value -  Node_left(node)->height->value);
		}
		if( Node_isleaf(Node_right(node)) ) {
			*value = dmin(*value, node->height->value -  Node_right(node)->height->value);
		}
		_get_max_down(node->left,  value);
		_get_max_down(node->right, value);
	}
	else if( Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr) ){
		*value = dmin(*value, Node_parent(node)->height->value - node->height->value);
	}
}

// Used for expanding/shrinking (+/-) the tree
// Works for hetero and homochronous
// If there is no calibration the whole tree will be scaled
// otherwise the subtrees below each calibrated node will not be scaled
void get_adjustment_bounds( Tree *tree, double *down, double *up ){
	Constraint *cnstr = Tree_root(tree)->height->cnstr;
	
	if ( Constraint_upper_fixed( cnstr) ) {
		*up = Constraint_fupper(cnstr) - Node_height(Tree_root(tree));
	}
	else {
		*up = Node_height(Tree_root(tree))/2;
	}
	
	*down = Node_height(Tree_root(tree)) - Constraint_flower(cnstr);
	
	_get_max_down(Tree_root(tree), down);
	//*down -= 0.1;
	*down *= -1;
}

// Used for scaling (*) the tree
// Works for hetero and homochronous
// If the tree is calibrated then there should be NO fixed calibration
// It should be needed only for nodes above the calibrations AND if the root is not fixed
void set_scale_bounds( Tree *tree, Parameters *scaler ){
	const double low = 0.08;
	const double high = 1.5;
	
	double min = DBL_MIN;
	double max = DBL_MAX;
	bool isCalibrated = false;
	
	Node **nodes = Tree_nodes(tree);
	Parameter *p = NULL;
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		p = nodes[i]->height;
		if ( Constraint_fixed(p->cnstr)) continue;
		
		if( Constraint_lower_fixed(p->cnstr) ){
			min = dmax(min, Constraint_flower(p->cnstr)/p->value);
			isCalibrated = true;
		}
		if( Constraint_upper_fixed(p->cnstr) ){
			//max = dmin(max, p->value/Constraint_fupper(p->cnstr));
			max = dmin(max, Constraint_fupper(p->cnstr)/p->value);
			isCalibrated = true;
		}
	}
	if (isCalibrated) {
		Parameters_set_bounds(scaler, 0, min, max);
	}
	else {
		Parameters_set_bounds(scaler, 0, low, high);
	}
	
}

#pragma mark -
#pragma mark Optimization

bool optimization_interrupted = false;

static void _signal_callback_handler( int signum ) {
	//printf("Caught signal %d\n",signum);
	if ( optimization_interrupted == true ) {
		exit(SIGINT);
	}
	optimization_interrupted = true;
}

double optimize_singletreelikelihood( SingleTreeLikelihood *stlk ){
    
    
    if( stlk->opt.interruptible ){
        signal(SIGINT, _signal_callback_handler);
        optimization_interrupted = false;
	}
    
	Node **nodes = Tree_nodes(stlk->tree);
	
	
	BrentData *data_brent = new_BrentData( stlk );
	Parameters *oneparameter = new_Parameters(1);
    
    
	MultivariateData *data_multivariate = NULL;
    
	Optimizer *opt_freq = NULL;
	
	Optimizer *opt_rel_rate = NULL;
    
    Optimizer *opt_all = NULL;
	
	Optimizer *opt_sm = NULL;
	
	Optimizer *opt_bl = NULL;
	
	Optimizer *opt_rate = NULL;
	MultivariateData *multivariate_data_rate = NULL;
	
	Optimizer *opt_height = NULL;
	MultivariateData *multivariate_data_height = NULL;
	
	
	// scale only heights
	Parameters *height_scaler = NULL;
	Optimizer *opt_scale_height = NULL;
	
	// scale heights and rate
	Parameters *height_scaler_height_rate = NULL;
	Optimizer *opt_scale_height_rate = NULL;
	
	Parameters *height_shift = NULL;
	Optimizer *opt_shift_height = NULL;
	
	bool use_scaler_height      = false;
	bool use_scaler_height_rate = false;
    
	OptConfig opt = stlk->opt;
    
    Parameters *all_params = new_Parameters(10);
    
    
    if( opt.bl.optimize ){
        opt_all = new_Optimizer( OPT_CG_PR );
        //all_params = new_Parameters( Tree_node_count(stlk->tree)-1 + Parameters_count_optimizable(stlk->sm->m->freqs, NULL ) + Parameters_count_optimizable(stlk->sm->m->rates, NULL ));
        if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
        opt_set_data(opt_all, data_multivariate);
        opt_set_objective_function(opt_all, _cg_optimize_all);
        opt_set_tolfx(opt_all, opt.freqs.tolfx);
        
        Node *right = Node_right( Tree_root(stlk->tree) );
        Node *left  = Node_left( Tree_root(stlk->tree) );
        Node_set_distance(left, Node_distance(right)+Node_distance(left) );
        Node_set_distance(right, 0);
        Parameter_set_fixed(right->distance, true);
        SingleTreeLikelihood_update_one_node(stlk, right);
        SingleTreeLikelihood_update_one_node(stlk, left);
    }
    
	
	// Init frequencies of the substitution model
	if( opt.freqs.optimize ){
		int n = Parameters_count_optimizable(stlk->sm->m->freqs, NULL );
        
		
		if( n == 0 ){
			opt.freqs.optimize = false;
		}
		else{
            Parameters_add_parameters(all_params, stlk->sm->m->freqs);
            
			opt_freq = new_Optimizer( opt.freqs.method );
			opt_set_max_iteration(opt_freq, opt.freqs.max_iteration);
			
			if( opt.freqs.method == OPT_BRENT ){
				opt_set_data(opt_freq, data_brent );
				opt_set_objective_function(opt_freq, _brent_optimize_frequency );
				opt_set_tolx(opt_freq, opt.freqs.tolx);
			}
			else if( opt.freqs.method == OPT_CG_PR || opt.freqs.method == OPT_CG_FR ){
				if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
				opt_set_data(opt_freq, data_multivariate);
				opt_set_objective_function(opt_freq, _cg_optimize_frequencies);
				opt_set_tolfx(opt_freq, opt.freqs.tolfx);
			}
			else{
				fprintf(stderr, "Optmization not supported: frequencies %d %d\n", opt.freqs.method, OPT_BRENT  );
				exit(1);
			}
		}
	}
	
	// Init relative rates of the subsitution model
	if( opt.relative_rates.optimize ){
		int n = Parameters_count_optimizable(stlk->sm->m->rates, NULL );
		
		if( n == 0 ){
			opt.relative_rates.optimize = false;
		}
		else{
            Parameters_add_parameters(all_params, stlk->sm->m->rates);
            //Parameters_print(all_params);
            
			opt_rel_rate = new_Optimizer( opt.relative_rates.method );
			opt_set_max_iteration(opt_rel_rate, opt.relative_rates.max_iteration );
			
			if( n == 1 || opt.relative_rates.method == OPT_BRENT ){
				opt_set_data(opt_rel_rate, data_brent );
				opt_set_objective_function(opt_rel_rate, _brent_optimize_relative_rate );
                opt_set_tolx(opt_rel_rate, opt.relative_rates.tolx);
			}
            
			else if( opt.relative_rates.method == OPT_CG_PR || opt.relative_rates.method == OPT_CG_FR ){
				if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
				opt_set_data(opt_rel_rate, data_multivariate);
				opt_set_objective_function(opt_rel_rate, _cg_optimize_rel_rates);
				//opt_set_tolfx(opt_rel_rate, opt.relative_rates.tolfx);
                opt_set_tolfx(opt_rel_rate, 0.001);
			}
			else{
				fprintf(stderr, "Optmization not supported: relative rates\n" );
				exit(1);
			}
		}
        
	}
		
	
	// Init branch length optimization
	if( opt.bl.optimize ){
        
//        // fix the right node, set it to 0 and add its distance to the left node
//        // we should unfix it at the end
//        Node *right = Node_right( Tree_root(stlk->tree) );
//        Node *left  = Node_left( Tree_root(stlk->tree) );
//        Node_set_distance(left, Node_distance(right)+Node_distance(left) );
//        Node_set_distance(right, 0);
//        Parameter_set_fixed(right->distance, true);
//        SingleTreeLikelihood_update_one_node(stlk, right);
//        SingleTreeLikelihood_update_one_node(stlk, left);
        
        //MARK: check this
        // not sure what it is not set somewhere else
        Parameter_set_fixed( Tree_root(stlk->tree)->distance, true );
        
        int count = 0;
        for ( int i = 0; i < Tree_node_count(stlk->tree); i++ ) {
            //if( Parameter_value(nodes[i]->distance) < 1e-5) Parameter_set_fixed(nodes[i]->distance, true);
            if( !Parameter_fixed(nodes[i]->distance) ){
                //Parameters_add(all_params, nodes[i]->distance);
                count++;
            }
        }
		
		if( count == 0 ){
			opt.bl.optimize = false;
		}
		else{
            opt_bl = new_Optimizer( opt.bl.method );
            
            if(opt.bl.method == OPT_BRENT ){
                opt_set_data(opt_bl, data_brent);
                opt_set_objective_function(opt_bl, optimize_brent_branch_length);
                opt_set_max_iteration(opt_bl, opt.bl.max_iteration);
                opt_set_tolx( opt_bl, opt.bl.tolx);
            }
            else if( opt.bl.method == OPT_CG_PR || opt.bl.method == OPT_CG_FR ){
				if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
				opt_set_data(opt_bl, data_multivariate);
				opt_set_objective_function(opt_bl, _cg_optimize_branch_length);
				opt_set_tolfx(opt_bl, opt.bl.tolfx);
                
                printf("should check this\n");
                exit(1);
                //all_params = new_Parameters( Tree_node_count(stlk->tree)-1);
                
                for ( int i = 0; i < Tree_node_count(stlk->tree); i++ ) {
                    if( Parameter_fixed(nodes[i]->distance) ) continue;
                    Parameters_set(all_params, i, nodes[i]->distance);
                }
			}
            else {
                error("Optmization not supported: branch length\n" );
            }
        }
        
        
	}
    
    // Init rate heterogeneity
    opt_algorithm algo = OPT_BRENT;
    //opt_algorithm algo = OPT_CG_PR;
    if(Parameters_count(stlk->sm->rates) > 0){
        
        //Parameters_add(all_params, stlk->sm->rates);
        opt_sm = new_Optimizer( algo );
        
        if(algo == OPT_BRENT ){
            opt_set_data(opt_sm, data_brent );
            opt_set_objective_function(opt_sm, _brent_optimize_sm_rates);
        }
        else {
            if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
            opt_set_data(opt_sm, data_multivariate);
            opt_set_objective_function(opt_sm, _cg_optimize_sm_rates);
            opt_set_tolfx(opt_sm, opt.relative_rates.tolfx);
            
        }
        
    }
    //Parameters_print(all_params);
	
	// Init rate optimization
    int n_optim_rates = 0;
	if( opt.rates.optimize && stlk->bm != NULL ){
        n_optim_rates = Parameters_count_optimizable(stlk->bm->rates, NULL );
		
		if( n_optim_rates == 0 ){
			opt.rates.optimize = false;
		}
		else{
			opt_rate = new_Optimizer( opt.rates.method );
			if( opt.rates.method == OPT_BRENT ){
				opt_set_data(opt_rate, data_brent);
				opt_set_objective_function(opt_rate, _brent_optimize_rate);
                opt_set_max_iteration(opt_rate, opt.rates.max_iteration);
				opt_set_tolx(opt_rate, opt.rates.tolx);
			}
			else if( opt.rates.method == OPT_CG_PR || opt.rates.method == OPT_CG_FR ){
				error("Optmization not supported: rates\n" );
				multivariate_data_rate = new_MultivariateData( stlk, NULL); // should not be null
				opt_set_data(opt_rate, multivariate_data_rate);
				opt_set_objective_function(opt_rate, _cg_optimize_rates);
				opt_set_tolfx(opt_rate, opt.rates.tolfx);
			}
			else{
				error("Optmization not supported: rates\n" );
			}
			opt_set_max_iteration(opt_rate, opt.rates.max_iteration);
		}
		
	}
	else opt.rates.optimize = false;
	
	// Init height optimization
	if( opt.heights.optimize ){
        int n = 0;
		for ( int i = 0; i < Tree_node_count(stlk->tree); i++ ) {
			if ( !Parameter_fixed(nodes[i]->height) ){
				n++;
			}
		}
		
		if( n == 0 ){
			opt.heights.optimize = false;
		}
		else {
			opt_height = new_Optimizer( opt.heights.method );
			opt_set_max_iteration(opt_height, opt.heights.max_iteration);
			
			if( opt.heights.method == OPT_BRENT ){
				opt_set_data(opt_height, data_brent);
				opt_set_objective_function(opt_height, _brent_optimize_height_threshold);
				opt_set_tolx(opt_height, opt.heights.tolx);
			}
			else {
				error("Optimization not supported: heights\n" );
			}
            
			
			bool isRootFixed = Parameter_fixed(Tree_root(stlk->tree)->height );
			
            data_brent->backup = dvector( Tree_node_count(stlk->tree) + Parameters_count(stlk->bm->rates));
            
			if ( !isRootFixed ) {
				
				height_scaler = new_Parameters(1);
				height_shift = new_Parameters(1);
				
				// Scale heights
				Parameters_add( height_scaler, new_Parameter("scale.height", 0.99, new_Constraint(0.000001, 5) ) );
				
				opt_scale_height= new_Optimizer( OPT_BRENT );
				opt_set_data(opt_scale_height, data_brent);
				opt_set_objective_function(opt_scale_height, _brent_optimize_scale_heights );
				opt_set_max_iteration(opt_scale_height, 100);
				opt_set_tolx(opt_scale_height, 0.001);
				
				// Expand the whole tree
				Parameters_add( height_shift, new_Parameter("shift.height", 0.01, new_Constraint(0, 10) ) );
				
				opt_shift_height = new_Optimizer( OPT_BRENT );
				opt_set_data(opt_shift_height, data_brent);
				opt_set_objective_function(opt_shift_height, _brent_optimize_add_heights );
				opt_set_max_iteration(opt_shift_height, 100);
				opt_set_tolx(opt_shift_height, 0.001);
				
                
				
				use_scaler_height = true;
			}
            
            if ( opt.rates.optimize ) {
                height_scaler_height_rate = new_Parameters(1);
                // Scale heights and rates
                Parameters_add( height_scaler_height_rate, new_Parameter("scale.rate.height", 0.99, new_Constraint(0.000001, 5) ) );
                
                opt_scale_height_rate = new_Optimizer( OPT_BRENT );
                opt_set_data(opt_scale_height_rate, data_brent);
                opt_set_objective_function(opt_scale_height_rate, _brent_optimize_scale_heights_rates );
                opt_set_max_iteration(opt_scale_height_rate, 100);
                opt_set_tolx(opt_scale_height_rate, 0.001);
                
                use_scaler_height_rate = true;
            }
			
		}
		
	}
    
#ifdef TIMETEST
    Parameters *heightsRateParams = NULL;
    Optimizer *opt_heightsRates = NULL;
    
    if( opt.heights.optimize){
        int count = Tree_node_count(stlk->tree)-1;
        if(opt.rates.optimize){
            count += Parameters_count(stlk->bm->rates);
        }
        
        heightsRateParams = new_Parameters(count);
        for ( int i = 0; i < Tree_node_count(stlk->tree); i++) {
            Node *node = Tree_node(stlk->tree, i);
            if( Node_isroot(node)){
                double earliest = 0;
                get_earliest(node,&earliest);
                Parameter_set_bounds(node->timeParameter, earliest, Parameter_upper(node->height));
                Parameter_set_value(node->timeParameter, Node_height(node));
            }
            if( !Node_isleaf(node)){
                Parameters_add(heightsRateParams, node->timeParameter);
            }
        }
        set_times(Tree_root(stlk->tree));
        //strict_heights_to_parameters(Tree_root(stlk->tree),heightsRateParams);
        
        if(opt.rates.optimize){
            for ( int i = 0; i < Parameters_count(stlk->bm->rates); i++ ) {
                Parameters_add(heightsRateParams, Parameters_at(stlk->bm->rates, i));
            }
        }
        
        opt_heightsRates = new_Optimizer( OPT_CG_PR );
        if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
        opt_set_data(opt_heightsRates, data_multivariate);
        opt_set_objective_function(opt_heightsRates, _cg_optimize_heightsRate2);
        opt_set_tolfx(opt_heightsRates, 0.000001);
    }
#endif
    
    if( stlk->approx == TREELIKELIHOOD_APPROXIMATION_NONE ){
        data_brent->f = standard_loglikelihood_brent;
    }
    else if( stlk->approx == TREELIKELIHOOD_APPROXIMATION_HESSIAN_DIAGONAL ){
        data_brent->f = loglikelihood_approximation_diag;
    }
    else if( stlk->approx == TREELIKELIHOOD_APPROXIMATION_HESSIAN ){
        data_brent->f = loglikelihood_approximation_full;
    }
    else if( stlk->approx == TREELIKELIHOOD_APPROXIMATION_DERIVATIVE_DIAGONAL ){
        data_brent->f = loglikelihood_approximation_diag;
    }
	
    if( opt.verbosity > 0 && stlk->approx > 0 ){
        fprintf(stdout, "Optimization using approximation %d\n", stlk->approx);
    }
    
    
	if ( (opt.heights.optimize || opt.rates.optimize) && opt.bl.optimize ){
		fprintf(stderr, "Cannot optimize rate and branch length parameters togther\n" );
		exit(1);
	}
	
	bool isCalibrated = Tree_is_calibrated( stlk->tree );
    if( isCalibrated ){
        data_brent->threshold = -1;
    }
	//boolean isRootRangeCalibrated = Constraint_lower_fixed(Tree_root(stlk->tree)->height->cnstr) || Constraint_upper_fixed(Tree_root(stlk->tree)->height->cnstr);
    
	int max_rounds = opt.max_rounds;
    
	//double lnl = stlk->calculate(stlk);
    double lnl = data_brent->f(data_brent);
    
	
	if( !opt.freqs.optimize && !opt.relative_rates.optimize && !opt.rates.optimize && !opt.heights.optimize && !opt.pinv.optimize && !opt.gamma.optimize && !opt.bl.optimize ){
		max_rounds = -1;
	}
	
    TopologyOptimizer *topology = NULL;
    int topo_failure = 0;
    
	if ( opt.topology_optimize ) {
        topology = new_TopologyOptimizer( stlk, opt.topology_alogrithm );
        TopologyOptimizer_set_nthreads(topology, opt.topology_threads);
    }
    
    
    
	double fret = lnl;
	
	//if(opt.heights.optimize) opt.verbosity = 3;
    
	if ( opt.verbosity > 1 ) {
		fprintf(stdout, "Optimize branch lengths %s\n", (opt.bl.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize frequencies    %s\n", (opt.freqs.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize gamma          %s\n", (opt.gamma.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize pinv           %s\n", (opt.pinv.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize relative rates %s\n", (opt.relative_rates.optimize ? "yes" : "no"));
		
		fprintf(stdout, "Optimize heights        %s\n", (opt.heights.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize rates          %s\n", (opt.rates.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize topology       %s\n", (opt.topology_optimize ? "yes" : "no"));
		fprintf(stdout, "\n");
		
		if ( opt.rates.optimize ) {
			for (int i = 0; i < BranchModel_n_rate(stlk->bm); i++) {
				fprintf(stdout, "Rate %f fret = %f\n", Parameters_value(stlk->bm->rates, i), fret);
			}
		}
		if ( opt.heights.optimize ) {
			fprintf(stdout, "Root height %f\n", Node_height(Tree_root(stlk->tree)));
			Parameter_print(Tree_root(stlk->tree)->height);
		}
	}
	
	int failed_scaling  = 0;
	int failed_scaling2 = 0;
	int failed_shift    = 0;
	
	int iter_bl = (max_rounds > 1 ? 2 : 1);
    
    
    if ( stlk->opt.verbosity ==  1 ) {
        fprintf(stdout, "LnL: %f\n", lnl);
        if ( opt.rates.optimize ) {
            if ( Parameters_count(stlk->bm->rates) == 1 ) {
                fprintf(stdout, "Rate        rate: %f\n", Parameters_value(stlk->bm->rates, 0) );
            }
        }
        if( opt.heights.optimize ){
            fprintf(stdout, "Root height height: %f\n", Node_height(Tree_root(stlk->tree)));
        }
    }

#ifdef DEBUG_OPT_CLOCK
    FILE *testFile = NULL;
    if( stlk->bm != NULL ){
        testFile = fopen("opt.trees","w");
        assert(testFile);
        Tree_print_nexus_header_figtree(testFile, stlk->tree);
    }
#endif
    
    time_t start_time;
    time_t end_time;
    time(&start_time);
    
    bool optimize_freqs = opt.freqs.optimize;
    bool optimize_relative_rates = opt.relative_rates.optimize;
    
    //FIXME: remove
    data_brent->threshold = 1e-5;
    data_brent->threshold = -1;
    bool do_local = (stlk->bm != NULL && Parameters_count(stlk->bm->rates) == 1);
    
    if(false ){
        opt_optimize( opt_all, all_params, &lnl);
        lnl = -lnl;
    }
    else
        for ( int rounds = 0; rounds < max_rounds; rounds++ ) {
		//fret = 0;
		if ( opt.verbosity > 1  ) {
			fprintf(stdout, "========================================\n");
			fprintf(stdout, "Round %d lk = %f\n\n", rounds, lnl);
			if ( opt.rates.optimize ) {
				for (int i = 0; i < Parameters_count(stlk->bm->rates); i++) {
					fprintf(stdout, "Rate %f\n", Parameters_value(stlk->bm->rates, i) );
				}
			}
			if ( opt.heights.optimize ) {
				fprintf(stdout, "Root height %f\n", Node_height(Tree_root(stlk->tree)));
			}
			
		}
		
		if ( stlk->opt.verbosity > 0 ) {
			fprintf(stdout, "\n");
		}
        
//        double status = opt_optimize( opt_all, all_params, &fret);
//        fret = -fret;
//        if ( stlk->opt.verbosity > 0 ) {
//            fprintf(stdout, "All   LnL: %f {%f}\n", fret, fret-lnl );
//        }
        
		// Branch length
		if( opt.bl.optimize ){

			if( opt.bl.method == OPT_BRENT ){
                
                function_f f = data_brent->f;
                if( stlk->use_upper ){
                    f = data_brent->f;
                    data_brent->f = standard_loglikelihood_upper_brent;
                }
                
                //iter_bl = 1;
				fret = optimize_brent_branch_length_all(stlk, opt_bl, data_brent, oneparameter, iter_bl);
                
				data_brent->f = f;
                
                if( stlk->use_upper ){
                    SingleTreeLikelihood_update_all_nodes(stlk);
                    fret = f(data_brent);
                }
				
			}
            else {
				double status = opt_optimize( opt_all, all_params, &fret);
				if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");
                fret = -fret;
                
                if( status != OPT_SUCCESS ){
                    SingleTreeLikelihood_update_all_nodes(stlk);
                }
			}
            
            if ( fret - lnl < 1 || opt.topology_optimize ) {
                iter_bl = 1;
            }
            
            if ( stlk->opt.verbosity > 0 ) {
                fprintf(stdout, "Branch lengths LnL: %f tree length %f {%f}", fret, Tree_length(stlk->tree), fret-lnl );
                
                if( stlk->opt.verbosity > 1 ){
                    time(&end_time);
                    printf(" %f", difftime(end_time, start_time));
                }
                printf("\n");
            }
		}

//        if( opt.freqs.optimize  && opt.relative_rates.optimize ){
//            double status = opt_optimize( opt_all, all_params, &fret);
//            if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");
//            fret = -fret;
//            if( status != OPT_SUCCESS ){
//                stlk->sm->m->need_update = true;
//                SingleTreeLikelihood_update_all_nodes(stlk);
//            }
//            if ( stlk->opt.verbosity > 0 ) {
//				stlk->sm->m->update_frequencies(stlk->sm->m);
//				fprintf(stdout, "Frequencies    LnL: %f [", fret);
//				for ( int i = 0; i < stlk->sm->nstate; i++ ) {
//					fprintf( stdout, "%s%f", (i==0?"":","), stlk->sm->m->_freqs[i] );
//				}
//				fprintf(stdout, "]  {%f}\n", fret-lnl );
//                
//                fprintf(stdout, "Relative rates LnL: %f [", fret);
//                for ( int i = 0; i < Parameters_count(stlk->sm->m->rates); i++ ) {
//                    fprintf( stdout, "%s%f", (i==0?"":","), Parameters_value(stlk->sm->m->rates,i) );
//                }
//                fprintf(stdout, "]  {%f} [%d]\n", fret-lnl,opt.relative_rates.method );
//			}
//        }
//        
//        else
        if( optimize_freqs  && optimize_relative_rates ){
            double templnl;
            int count = 0;
            do {
                templnl = fret;
                // Frequencies
                if( opt.freqs.method == OPT_BRENT ){
                    fret = optimize_brent_frequencies_all(stlk, opt_freq, data_brent, oneparameter);
                }
                else {
                    double status = opt_optimize( opt_freq, stlk->sm->m->freqs, &fret);
                    if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");
                    fret = -fret;
                    if( status != OPT_SUCCESS ){
                        stlk->sm->m->need_update = true;
                        SingleTreeLikelihood_update_all_nodes(stlk);
                    }
                }
                //fprintf(stdout, "Frequencies    LnL: %f\n", fret);
                // Relative rates
                if( opt.relative_rates.method == OPT_BRENT ){
                    fret = optimize_brent_relative_rate_all(stlk, opt_rel_rate, data_brent, oneparameter);
                }
                else{
                    double status = opt_optimize( opt_rel_rate, stlk->sm->m->rates, &fret);
                    if( status == OPT_ERROR ) error("OPT.relative_rates No SUCCESS!!!!!!!!!!!!\n");
                    fret = -fret;
                    if( status != OPT_SUCCESS ){
                        stlk->sm->m->need_update = true;
                        SingleTreeLikelihood_update_all_nodes(stlk);
                    }
                }
                //fprintf(stdout, "Relative rates LnL: %f\n", fret);
                count++;
            }
            while ( fret - templnl > 0.01 );
            
            if ( stlk->opt.verbosity > 0 ) {
				stlk->sm->m->update_frequencies(stlk->sm->m);
				fprintf(stdout, "Frequencies    LnL: %f [", fret);
				for ( int i = 0; i < stlk->sm->nstate; i++ ) {
					fprintf( stdout, "%s%f", (i==0?"":","), stlk->sm->m->_freqs[i] );
				}
				fprintf(stdout, "]  {%f}\n", fret-lnl );
                
                fprintf(stdout, "Relative rates LnL: %f [", fret);
                for ( int i = 0; i < Parameters_count(stlk->sm->m->rates); i++ ) {
                    fprintf( stdout, "%s%f", (i==0?"":","), Parameters_value(stlk->sm->m->rates,i) );
                }
                fprintf(stdout, "]  {%f} [%s] #%d", fret-lnl, OPT_ALGORITHMS[opt.relative_rates.method], count );
                
                if( stlk->opt.verbosity > 1 ){
                    time(&end_time);
                    printf(" %f", difftime(end_time, start_time));
                }
                printf("\n");
			}
            if(rounds > 1) optimize_freqs = optimize_relative_rates = false;
        }
		
		else {
            // Frequencies
            if( optimize_freqs ){
                if( opt.freqs.method == OPT_BRENT ){
                    fret = optimize_brent_frequencies_all(stlk, opt_freq, data_brent, oneparameter);
                }
                else {
                    double status = opt_optimize( opt_freq, stlk->sm->m->freqs, &fret);
                    if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");
                    fret = -fret;
                    if( status != OPT_SUCCESS ){
                        stlk->sm->m->need_update = true;
                        SingleTreeLikelihood_update_all_nodes(stlk);
                    }
                }
                
                if ( stlk->opt.verbosity > 0 ) {
                    stlk->sm->m->update_frequencies(stlk->sm->m);
                    fprintf(stdout, "Frequencies    LnL: %f [", fret);
                    for ( int i = 0; i < stlk->sm->nstate; i++ ) {
                        fprintf( stdout, "%s%f", (i==0?"":","), stlk->sm->m->_freqs[i] );
                    }
                    fprintf(stdout, "]  {%f} [%s]", fret-lnl, OPT_ALGORITHMS[opt.freqs.method] );
                    
                    if( stlk->opt.verbosity > 1 ){
                        time(&end_time);
                        printf(" %f", difftime(end_time, start_time));
                    }
                    printf("\n");
                }
            }
            
            // Relative rates
            if( optimize_relative_rates ){
                if( opt.relative_rates.method == OPT_BRENT ){
                    fret = optimize_brent_relative_rate_all(stlk, opt_rel_rate, data_brent, oneparameter);
                }
                else{
                    double status = opt_optimize( opt_rel_rate, stlk->sm->m->rates, &fret);
                    if( status == OPT_ERROR ) error("OPT.relative_rates No SUCCESS!!!!!!!!!!!!\n");
                    fret = -fret;
                    if( status != OPT_SUCCESS ){
                        stlk->sm->m->need_update = true;
                        SingleTreeLikelihood_update_all_nodes(stlk);
                    }
                }
                
                if ( stlk->opt.verbosity > 0 ) {
                    fprintf(stdout, "Relative rates LnL: %f [", fret);
                    for ( int i = 0; i < Parameters_count(stlk->sm->m->rates); i++ ) {
                        fprintf( stdout, "%s%f", (i==0?"":","), Parameters_value(stlk->sm->m->rates,i) );
                    }
                    fprintf(stdout, "]  {%f} [%s]", fret-lnl, OPT_ALGORITHMS[opt.relative_rates.method] );
                    
                    if( stlk->opt.verbosity > 1 ){
                        time(&end_time);
                        printf(" %f", difftime(end_time, start_time));
                    }
                    printf("\n");
                }
            }
        }
        
//        if( opt.bl.optimize ){
//            
//			if( opt.bl.method == OPT_BRENT ){
//                
//                function_f f = data_brent->f;
//                if( stlk->use_upper ){
//                    f = data_brent->f;
//                    data_brent->f = standard_loglikelihood_upper_brent;
//                }
//                
//                //iter_bl = 1;
//				fret = optimize_brent_branch_length_all(stlk, opt_bl, data_brent, oneparameter, iter_bl);
//                
//				data_brent->f = f;
//                
//                if( stlk->use_upper ){
//                    SingleTreeLikelihood_update_all_nodes(stlk);
//                    fret = f(data_brent);
//                }
//				
//			}
//            else {
//				double status = opt_optimize( opt_all, all_params, &fret);
//				if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");
//                fret = -fret;
//                
//                if( status != OPT_SUCCESS ){
//                    SingleTreeLikelihood_update_all_nodes(stlk);
//                }
//			}
//            
//            
//            
//            if ( fret - lnl < 1 ) {
//                iter_bl = 1;
//            }
//            
//            if ( stlk->opt.verbosity > 0 ) {
//                fprintf(stdout, "Branch lengths LnL: %f tree length %f {%f}\n", fret, Tree_length(stlk->tree), fret-lnl );
//            }
//		}
        
		// Rate heterogeneity
        if( Parameters_count(stlk->sm->rates) > 0 ){
            if(algo == OPT_BRENT){
                fret = optimize_brent_sm_rates_all(stlk, opt_sm, data_brent, oneparameter);
            }
            else{
                double status= opt_optimize( opt_sm, stlk->sm->rates, &fret);
                
                if( status == OPT_ERROR ) error("OPT.SM.RATES No SUCCESS!!!!!!!!!!!!\n");
                fret = -fret;
            }
            
            bool flaggy = false;
            if ( stlk->opt.verbosity > 0 ) {
                for (int i = 0; i < Parameters_count(stlk->sm->rates); i++) {
                    if(strcmp(Parameters_name(stlk->sm->rates, i), "sitemodel.alpha") == 0){
                        fprintf(stdout, "Gamma          LnL: %f shape: %f  {%f}%s", fret, Parameters_value(stlk->sm->rates, i), fret-lnl, (i == Parameters_count(stlk->sm->rates)-1 ? "": "\n") );
                    }
                    else if(strcmp(Parameters_name(stlk->sm->rates, i), "sitemodel.pinv") == 0){
                        fprintf(stdout, "Prop invariant LnL: %f p: %f  {%f}%s", fret, Parameters_value(stlk->sm->rates, i), fret-lnl, (i == Parameters_count(stlk->sm->rates)-1 ? "": "\n") );
                    }
                    else{
                        flaggy = true;
                    }
                }
                if(flaggy){
                    fprintf(stdout,     "SiteModel      LnL: %f  {%f}", fret, fret-lnl );
                }
                
                if( stlk->opt.verbosity > 1 ){
                    time(&end_time);
                    printf(" %f", difftime(end_time, start_time));
                }
                printf("\n");
                //Parameters_print(stlk->sm->rates);
            }
        }
		
        //		if( opt.heights.optimize && use_brlen ){
        //			for (int i = 0; i < Parameters_count(ps_height); i++) {
        //				data_height->index_param = map_height[i];
        //				Parameters_set( oneparameter, 0,  Parameters_at(ps_height, i) );
        //				Node *node = nodes[ map_height[i] ];
        //
        //				if( Node_isroot(node) ){
        //					if ( !Constraint_fupper( Tree_root(stlk->tree)->height->cnstr ) ) {
        //						//Parameters_set_upper(oneparameter, 0, dmin(HEIGHT_MAX, Parameters_value(oneparameter, 0)*2) );
        //						Parameters_set_upper(oneparameter, 0, Parameters_value(oneparameter, 0)*2 );
        //					}
        //				}
        //
        //				double status = opt_optimize( opt_height, oneparameter, &fret);
        //				if( status == OPT_ERROR ) error("OPT.HEIGHTS No SUCCESS!!!!!!!!!!!!\n");
        //
        //				double value = Parameters_value(oneparameter, 0);
        //				Node_set_height( node,  value );
        //
        //				if( !Node_isroot(node) ){
        //					//value = Node_distance(node)/(Node_height(Node_parent(node)) - value);
        //					//stlk->bm->set(stlk->bm, map_height[i], value);
        //				}
        //				if ( !Node_isleaf(node) ) {
        //					Node *child = Node_left(node);
        //					value = Node_distance(child)/(Parameters_value(oneparameter, 0) - Node_height(child));
        //					stlk->bm->set(stlk->bm, child->postorder_idx, value);
        //					SingleTreeLikelihood_update_three_nodes(stlk, child->postorder_idx );
        //
        //					child = Node_right(node);
        //					value = Node_distance(child)/(Parameters_value(oneparameter, 0) - Node_height(child));
        //					stlk->bm->set(stlk->bm, child->postorder_idx, value);
        //					SingleTreeLikelihood_update_three_nodes(stlk, child->postorder_idx );
        //
        //				}
        //
        //				SingleTreeLikelihood_update_three_nodes(stlk, map_height[i] );
        //				Tree_constraint_heights(stlk->tree);
        //			}
        //
        //			if ( stlk->opt.verbosity > 0 ) {
        //				fprintf(stdout, "Root height    LnL: %f height: %f {%f}\n", -fret, Node_height(Tree_root(stlk->tree)), fabs(fret-lk) );
        //			}
        //		}
		
#ifdef TIMETEST
            if( opt.heights.optimize ){
                if( opt.rates.optimize ){
                    fret = optimize_brent_rate_all(stlk, opt_rate, data_brent, oneparameter);
                }
                double lnl2 = fret;
                double status = opt_optimize( opt_heightsRates, heightsRateParams, &fret);
                fret = -fret;
                if( status == OPT_ERROR ) error("Optimization heights and rates is not a success\n");
                
                if ( stlk->opt.verbosity > 0 ) {
                    fprintf(stdout, "Root height    LnL: %f height: %f {%f}\n", fret, Node_height(Tree_root(stlk->tree)), fret-lnl );
                }
                if( opt.rates.optimize ){
                    if(fabs(fret - lnl) > 0.1 ){
                        fret = optimize_brent_rate_all(stlk, opt_rate, data_brent, oneparameter);
                    }
                    if ( stlk->opt.verbosity > 0 ) {
                        if ( Parameters_count(stlk->bm->rates) == 1 ) {
                            fprintf(stdout, "Rate           LnL: %f rate: %f {%f}\n", fret, Parameters_value(stlk->bm->rates, 0), fret-lnl );
                        }
                        else {
                            fprintf(stdout, "Rates          LnL: %f {%f}\n", fret, fret-lnl );
                        }
                    }
                }
            }
#else
        
        // HEIGHT OPTIMIZATION
		if( opt.heights.optimize && !do_local ){
			function_f f = data_brent->f;
            if( stlk->use_upper ){
                f = data_brent->f;
                data_brent->f = standard_loglikelihood_upper_brent;
            }
            int n = 1;//(rounds == 0 && Parameters_count(stlk->bm->rates) == 1 ? 10 : 1 );
            for( int i = 0 ; i < n; i++){
                if( opt.heights.method == OPT_BRENT ){
                    SingleTreeLikelihood_update_all_nodes(stlk);
                    fret = optimize_brent_heights_all_threshold(stlk, opt_height, data_brent, oneparameter);
                    
                    // If threshold is used the likelihood can be lower than the previous iteration
                    if( fret < lnl ){
                        data_brent->threshold = -1;
                        lnl = fret;
                        n++;
                    }
                    
                    if ( stlk->opt.verbosity > 0 ) {
                        if ( stlk->opt.verbosity > 1 )fprintf(stderr, "---------------------------------------------------------------------\n");
                        fprintf(stdout, "Root height    LnL: %f height: %f {%f}\n", fret, Node_height(Tree_root(stlk->tree)), fret-lnl );
                    }
                    
                }
            }
            data_brent->f = f;
		}
        
        if ( stlk->bm != NULL ) {
            
            if( opt.heights.optimize && Parameters_count(stlk->bm->rates) == 1 && do_local){
                double templnl = fret;
                double *backup = dvector(Tree_node_count(stlk->tree));
                
                for( int i = 0 ; i < 5; i++){
                    Tree_heights_to_vector(stlk->tree, backup);
                    
                    fret = optimize_brent_heights_strict_local(stlk, opt_height, data_brent, oneparameter);
                    
                    if ( stlk->opt.verbosity > 1 ) printf("before %f new %f\n",templnl,fret);
                
                    if(fret < templnl ){
                        do_local = false;
                        Tree_vector_to_heights(backup, stlk->tree);
                        Tree_constraint_heights(stlk->tree);
                        SingleTreeLikelihood_update_all_nodes(stlk);
                        fret = templnl;
                        break;
                    }
                    else if( fret - templnl < 0.01  ){
                        do_local = false;
                        //templnl = fret;
                        break;
                    }
                    templnl = fret;
                }
                
                free(backup);
                if ( stlk->opt.verbosity > 0 ) {
                    if ( stlk->opt.verbosity > 1 )fprintf(stdout, "---------------------------------------------------------------------\n");
                    fprintf(stdout, " --- Root height    LnL: %f height: %f {%f}\n", fret, Node_height(Tree_root(stlk->tree)), fret-lnl );
                }
                
                if( !do_local ){
                    lnl = fret;
                    continue;
                }
            }
            
            // no constrained heights
            if( opt.heights.optimize && rounds == -1 ){
                double templk = fret;
                //            opt_set_objective_function(opt_scale_height_rate, _brent_optimize_scale_heights_polynomial);
                //
                //			fret = optimize_scale_heights_polynomial(stlk, opt_scale_height_rate, data_brent, height_scaler_height_rate, templk);
                //
                //            opt_set_objective_function(opt_scale_height_rate, _brent_optimize_scale_heights_rates );
                Parameters *params = new_Parameters(2);
                Parameter *coef = new_Parameter("coef", 0.01, new_Constraint(0.0001, 0.1));
                Parameters_add(params, coef);
                Parameters_add(params, Parameters_at(stlk->bm->rates, 0));
                
                Optimizer *opt_poly = new_Optimizer( OPT_CG_FR );
                opt_set_max_iteration(opt_poly, opt.relative_rates.max_iteration );
                if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);
                opt_set_data(opt_poly, data_multivariate);
                opt_set_objective_function(opt_poly, _cg_optimize_scale_rate_heights_polynomial);
                opt_set_tolfx(opt_poly, 0.001);
                
                double status = opt_optimize( opt_poly, params, &fret);
                if( status == OPT_ERROR ) error("OPT.poly No SUCCESS!!!!!!!!!!!!\n");
                fret = -fret;
                if( status != OPT_SUCCESS ){
                    SingleTreeLikelihood_update_all_nodes(stlk);
                }
                
                if ( opt.heights.optimize && opt.verbosity > 1 ){
                    fprintf(stdout, "---------------------------------------------------------------------\n");
                    if( fret > templk ){
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.POLY scaling success: old=%f new=%f scaler=%f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), fret-lnl );
                    }
                    else {
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.POLY scaling NO success: old=%f new=%f scaler=%f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), fret-lnl );
                        
                    }
                }
                free_Parameters_soft(params);
                free_Parameter(coef);
                free_Optimizer(opt_poly);
            }
            
            // STRICT CLOCK
            if( use_scaler_height_rate && Parameters_count(stlk->bm->rates) == 1 ){
                double templk = fret;
                
                // Scaling of heights and rate at the same time
                // seems to well on strict clocks only (tested with flu dataset)
                //if ( failed_scaling < 4 ) {
                fret = optimize_scale_rate_height_strict(stlk, opt_scale_height_rate, height_scaler_height_rate, templk);
                
                if ( opt.verbosity > 1 ){
                    fprintf(stdout, "---------------------------------------------------------------------\n");
                    if( fret > templk ){
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.RATES scaling success: old=%f new=%f scaler=%f height: %f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), Node_height(Tree_root(stlk->tree)), fret-lnl );
                    }
                    else {
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.RATES scaling NO success: old=%f new=%f scaler=%f height: %f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), Node_height(Tree_root(stlk->tree)), fret-lnl );
                        
                    }
                }
                if ( fret - templk < 0.01 ) {
                    failed_scaling++;
                }
                //}
                
                templk = fret;
                // Scaling of root height and its children if they have short branches and rate at the same time
                // Works very well when a child of the root is stuck to the root and the root cannot descrease its height
                
                if( !isCalibrated && false){
                    if( !Node_isleaf( Node_left(Tree_root(stlk->tree)) ) && !Node_isleaf( Node_right(Tree_root(stlk->tree)) ) ){
                        opt_set_objective_function(opt_scale_height_rate, _brent_optimize_scale_root_height_rate);
                        
                        fret = optimize_scale_root_rate_height_strict(stlk, opt_scale_height_rate, height_scaler_height_rate, templk);
                        
                        opt_set_objective_function(opt_scale_height_rate, _brent_optimize_scale_heights_rates );
                        
                        if ( opt.verbosity > 1 ){
                            fprintf(stdout, "---------------------------------------------------------------------\n");
                            if( fret > templk ){
                                fprintf(stdout, "OPT.SCALE.ROOT.HEIGHT.RATES scaling success: old=%f new=%f scaler=%f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), fret-lnl );
                            }
                            else {
                                fprintf(stdout, "OPT.SCALE.ROOT.HEIGHT.RATES scaling NO success: old=%f new=%f scaler=%f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), fret-lnl );
                                
                            }
                        }
                    }
                }
            }
            
            
            // Scaling of heights and rate at the same time
            // DISCRETE and LOCAL CLOCK
            if( use_scaler_height_rate && Parameters_count(stlk->bm->rates) > 1 ){//&& failed_scaling < -4 ) {
                double templk = fret;
                fret = optimize_scale_rate_height_all(stlk, opt_scale_height_rate, height_scaler_height_rate, templk);
                
                if ( opt.verbosity > 1 ){
                    fprintf(stdout, "---------------------------------------------------------------------\n");
                    if( fret > templk ){
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.RATES scaling success: old=%f new=%f scaler=%f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), fret-lnl );
                    }
                    else {
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.RATES scaling NO success: old=%f new=%f scaler=%f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), fret-lnl );
                        
                    }
                }
                if ( fret - templk < 0.01 ) {
                    failed_scaling++;
                }
            }
            
            
            
            
            // Scale heights only
            if( use_scaler_height && !isCalibrated ){//&& failed_scaling2 < -4 ){
                
                double templk = fret;
                double status;
                
                set_scale_bounds(stlk->tree, height_scaler);
                Parameters_set_value(height_scaler, 0, 1 );
                status = opt_optimize( opt_scale_height, height_scaler, &fret);
                if( status == OPT_ERROR ) error("OPT.SCALE.HEIGHTS No SUCCESS!!!!!!!!!!!!\n");
                
                fret = -fret;
                
                if( fret > templk ){
                    if ( opt.verbosity > 1 ){
                        fprintf(stdout, "---------------------------------------------------------------------\n");
                        fprintf(stdout, "OPT.SCALE.HEIGHTS               scaling success: old=%f new=%f scaler=%f [%f-%f] {%f}\n\n", templk, fret, Parameters_value( height_scaler, 0 ), Parameters_lower(height_scaler, 0), Parameters_upper(height_scaler, 0), fret-lnl );
                    }
                    SingleTreeLikelihood_scaler(stlk, Tree_root(stlk->tree), Parameters_value(height_scaler, 0));
                }
                else {
                    fret = templk;
                    if ( opt.verbosity > 1 ){
                        fprintf(stdout, "---------------------------------------------------------------------\n");
                        fprintf(stdout, "OPT.SCALE.HEIGHTS               scaling NO success: old=%f new=%f scaler=%f [%f-%f] {%f}\n\n", templk, fret, Parameters_value( height_scaler, 0 ), Parameters_lower(height_scaler, 0), Parameters_upper(height_scaler, 0), fret-lnl );
                    }
                }
                
                SingleTreeLikelihood_update_all_nodes(stlk);
                Tree_constraint_heights(stlk->tree);
                
                if ( fret - templk < 0.01 ) {
                    failed_scaling2++;
                }
            }
            
            // Add or substract time to tree
            // tree with calibrations, on the root too but not fixed
            // We don't scale if the root is fixed
            if( use_scaler_height ){//&& failed_shift < -4 ){
                double templk = fret;
                double status;
                double down;
                double up;
                get_adjustment_bounds(stlk->tree, &down, &up);
                
                Parameters_set_value(height_shift, 0, 0 );
                Parameters_set_bounds(height_shift, 0, down, up );
                
                status = opt_optimize( opt_shift_height, height_shift, &fret);
                if( status == OPT_ERROR ) error("OPT.SHIFT.HEIGHTS No SUCCESS!!!!!!!!!!!!\n");
                
                fret = -fret;
                
                if( fret > templk ){
                    if ( opt.verbosity > 1 ){
                        fprintf(stdout, "---------------------------------------------------------------------\n");
                        fprintf(stdout, "OPT.SHIFT.HEIGHTS               scaling success: old=%f new=%f shift=%f {%f}\n\n", templk, fret, Parameters_value( height_shift, 0 ), fret-lnl );
                    }
                    SingleTreeLikelihood_add_height(stlk, Tree_root(stlk->tree), Parameters_value(height_shift, 0));
                }
                else {
                    fret = templk;
                    if ( opt.verbosity > 1 ){
                        fprintf(stdout, "---------------------------------------------------------------------\n");
                        fprintf(stdout, "OPT.SHIFT.HEIGHTS               scaling NO success: old=%f new=%f shift=%f {%f}\n\n", templk, fret, Parameters_value( height_shift, 0 ), fret-lnl );
                    }
                }
                
                SingleTreeLikelihood_update_all_nodes(stlk);
                Tree_constraint_heights(stlk->tree);
                
                if ( fret - templk < 0.01 ) {
                    failed_shift++;
                }
            }
            
            // Rates in BranchModel
            if( opt.rates.optimize ){
                if( opt.rates.method == OPT_BRENT ){
                    fret = optimize_brent_rate_all(stlk, opt_rate, data_brent, oneparameter);
                }
                else{
                    double status = opt_optimize( opt_rate, stlk->bm->rates, &fret);
                    if( status == OPT_ERROR ) fprintf(stderr, "OPT.RATES No SUCCESS!!!!!!!!!!!!\n");
                    fret = -fret;
                }
                
                
                
                if ( stlk->opt.verbosity > 0 ) {
                    if ( stlk->opt.verbosity > 1 )fprintf(stdout, "---------------------------------------------------------------------\n");
                    if ( Parameters_count(stlk->bm->rates) > 0 ) {
                        fprintf(stdout, "Rate           LnL: %f rate: %f {%f}\n", fret, Parameters_value(stlk->bm->rates, 0), fret-lnl );
                    }
                    else {
                        fprintf(stdout, "Rates          LnL: %f {%f}\n", fret, fret-lnl );
                    }
                }
            }
            
            // Scaling of heights and rate at the same time
            // seems to well on strict clocks only (tested with flu dataset)
            // STRICT
            if( use_scaler_height_rate && Parameters_count(stlk->bm->rates) == 1 ){// && failed_scaling < 4 ) {
                
                double templk = fret;
                
                //opt_set_objective_function(opt_scale_height_rate, _brent_optimize_scale_heights_rates2 );
                fret = optimize_scale_rate_height_strict(stlk, opt_scale_height_rate, height_scaler_height_rate, templk);
                //opt_set_objective_function(opt_scale_height_rate, _brent_optimize_scale_heights_rates );
                if ( opt.verbosity > 1 ){
                    fprintf(stdout, "---------------------------------------------------------------------\n");
                    if( fret > templk ){
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.RATES scaling success: old=%f new=%f scaler=%f height: %f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), Node_height(Tree_root(stlk->tree)), fret-lnl );
                    }
                    else {
                        fprintf(stdout, "OPT.SCALE.HEIGHTS.RATES scaling NO success: old=%f new=%f scaler=%f height: %f {%f}\n\n", templk, fret, Parameters_value( height_scaler_height_rate, 0 ), Node_height(Tree_root(stlk->tree)), fret-lnl );
                        
                    }
                }
                if ( fret - templk < 0.01 ) {
                    failed_scaling++;
                }
            }
            
        }
		
		
        if( opt.topology_optimize && topo_failure < 2  ){
            //SingleTreeLikelihood_update_all_nodes(stlk);
            //fprintf(stdout, "Before topology    LnL: %f \n", stlk->calculate(stlk));
            
            fret = topology->optimize(topology);
            
            if ( stlk->opt.verbosity > 0 ) {
                fprintf(stdout, "Topology       LnL: %f NNIs: %d {%f} [%d]", fret, topology->moves, fret - lnl, topology->algorithm );
                
                if( stlk->opt.verbosity > 1 ){
                    time(&end_time);
                    printf(" %f", difftime(end_time, start_time));
                }
                printf("\n");
            }
            
            //if(topology->algorithm == TREE_SEARCH_NNI || topology->algorithm == TREE_SEARCH_NNNI) TopologyOptimizer_set_algorithm(topology, 1-topology->algorithm);
            
            if ( topology->moves == 0 ) {
                topo_failure++;
            }
            else {
                topo_failure = 0;
            }
        }
#endif
        
        if ( optimization_interrupted ) {
            fprintf(stderr, "Ctrl-c\n\n");
            fprintf(stderr, "0) To exit the program\n");
            fprintf(stderr, "1) Continue optimzation\n");
            fprintf(stderr, "2) Stop optimization\n");
            int answer = -1;
            
            while ( answer < 0 || answer > 2 ) {
                fprintf(stderr, "Choice: ");
                int r = scanf("%d",  &answer);
                if ( r == 0 ) {
                    answer = -1;
                }
            }
            
            switch (answer) {
                case 0:{
                    exit(SIGINT);
                    break;
                }
                case 1:{
                    optimization_interrupted = false;
                    break;
                }
                case 2:{
                    break;
                }
                default:{
                    break;
                }
            }
            
            if( optimization_interrupted ){
                lnl = fret;
                break;
            }
        }

        if(  fret - lnl < opt.precision ){//&& rounds > 2 ){
            if( opt.verbosity > 1 ){
                fprintf(stdout, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
                fprintf(stdout, "STOP LnL previous %f last %f round: %d\n\n", lnl, fret, rounds);
            }
            
//            if ( opt.heights.optimize){
//            SingleTreeLikelihood_scaler(stlk, Tree_root(stlk->tree), 0.6);
//            //Parameters_set_value(stlk->bm->rates, 0, Parameters_value(stlk->bm->rates, 0)/0.6);
//                Parameters_set_value(stlk->bm->rates, 0, 0.00146);
//            SingleTreeLikelihood_update_all_nodes(stlk);
//            lnl = stlk->calculate(stlk);
//            printf("lnl: %f\n", lnl);
//                opt.rates.optimize = false;
//            continue;
//            }
#ifdef TIMETEST
            lnl = fret;
            break;
#else
            if ( opt.heights.optimize && data_brent->threshold != -1 ) {
                
                double templk = fret;
                
                data_brent->threshold = -1;
                fret = optimize_brent_heights_all_threshold(stlk, opt_height, data_brent, oneparameter);
                
                if ( stlk->opt.verbosity > 0 ) {
                    if ( stlk->opt.verbosity > 1 )fprintf(stdout, "---------------------------------------------------------------------\n");
                    fprintf(stdout, "Root height    LnL: %f height: %f {%f} *\n", fret, Node_height(Tree_root(stlk->tree)), fret-templk );
                }
                
                if( fret-templk < opt.precision ){
                    lnl = templk;
                    break;
                }
            }
            else {
                lnl = fret;
                break;
            }
#endif
            
			
		}
		
		lnl = fret;
#ifdef DEBUG_OPT_CLOCK
        if( testFile != NULL && stlk->bm != NULL ){
            
            Node **nodes = Tree_nodes(stlk->tree);
            StringBuffer * buffer = new_StringBuffer(10);
            for ( int i = 0; i < Tree_node_count(stlk->tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    StringBuffer_empty(buffer);
                    StringBuffer_append_format(buffer, "%e", stlk->bm->get(stlk->bm,nodes[i]));
                    Node_set_annotation(nodes[i], "rate", buffer->c);
                    
                    if( stlk->bm->name == CLOCK_LOCAL ){
                        if ( stlk->bm->indicators[i] ) {
                            Node_set_annotation(nodes[i], "local", "1");
                        }
                    }
                    else if ( stlk->bm->name == CLOCK_DISCRETE ){
                        StringBuffer_empty(buffer);
                        StringBuffer_append_format(buffer, "%d", stlk->bm->map[ Node_id(nodes[i]) ]);
                        Node_set_annotation(nodes[i], "class", buffer->c);
                    }
                }
            }
            free_StringBuffer(buffer);
            
            fprintf(testFile, "tree TREE%d [&LnL=%f] = [&R] ", rounds, lnl);
            Tree_print_nexus_with_annotation(testFile, stlk->tree);
            fprintf(testFile, "\n");
            fflush(testFile);
        }
#endif
	}
#ifdef DEBUG_OPT_CLOCK
    if( testFile != NULL ){
        fprintf(testFile, "\nEnd;");
        fclose(testFile);
    }
#endif
    
    
    free_Parameters_soft(all_params);

#ifdef TIMETEST
    if ( opt_heightsRates != NULL) {
        free_Optimizer(opt_heightsRates);
        free_Parameters_soft(heightsRateParams);
    }
#endif
    
    free_BrentData(data_brent);
    if (data_multivariate != NULL ) free_MultivariateData( data_multivariate);
    
    free_Parameters_soft(oneparameter);
    
	if( opt.freqs.optimize ){
		free_Optimizer( opt_freq );
	}
	
	if( opt.relative_rates.optimize ){
		free_Optimizer( opt_rel_rate );
	}
    
    if( opt_all != NULL ){
        free_Optimizer(opt_all);
    }
    
    
    if( opt_sm != NULL ){
        free_Optimizer(opt_sm);
    }
	
	if ( opt.bl.optimize ){
		free_Optimizer( opt_bl );
    }
	
	if( opt.rates.optimize ){
		if (multivariate_data_rate != NULL ) free_MultivariateData( multivariate_data_rate);
		free_Optimizer( opt_rate );
	}
	
	if( opt.heights.optimize ){
		if (multivariate_data_height != NULL ) free_MultivariateData( multivariate_data_height );
		free_Optimizer( opt_height );
		
		if ( use_scaler_height ) {
			
			free_Parameters( height_scaler );
			free_Parameters( height_shift );
			
			free_Optimizer(opt_scale_height);
			free_Optimizer(opt_shift_height);
		}
		
		if ( use_scaler_height_rate ) {
			free_Parameters( height_scaler_height_rate );
			free_Optimizer(opt_scale_height_rate);
		}
	}
	if ( topology != NULL ) {
        free_TopologyOptimizer(topology);
    }
    
    
    //    if( stlk->bm!= NULL && Parameters_count(stlk->bm->rates) ==1 && stlk->hessian != NULL ){
    //        double ll = loglikelihood_approximation2(stlk);
    //        printf("Approx Optimize %f\n", ll);
    //        exit(0);
    //    }
	return lnl;
}

double optimize_singletreelikelihood2( SingleTreeLikelihood *stlk ){
	
	Node **nodes = Tree_nodes(stlk->tree);
	
    OptConfig opt = stlk->opt;
    
	MultivariateData *data_multivariate = new_MultivariateData( stlk, NULL);
    
    Optimizer *opt_all = new_Optimizer( OPT_CG_PR );
    
    Parameters *all_params = new_Parameters( Tree_node_count(stlk->tree)-1 + Parameters_count_optimizable(stlk->sm->m->freqs, NULL ) + Parameters_count_optimizable(stlk->sm->m->rates, NULL ));if( data_multivariate == NULL ) data_multivariate = new_MultivariateData( stlk, NULL);;
    
    opt_set_data(opt_all, data_multivariate);
    opt_set_objective_function(opt_all, _cg_optimize_all);
    opt_set_tolfx(opt_all, 0.00001);
    
    
    
	
	// Init frequencies of the substitution model
	if( opt.freqs.optimize ){
		int n = Parameters_count_optimizable(stlk->sm->m->freqs, NULL );
        if( n > 0 ) Parameters_add_parameters(all_params, stlk->sm->m->freqs);
        else opt.freqs.optimize = false;
	}
	
	// Init relative rates of the subsitution model
	if( opt.relative_rates.optimize ){
		int n = Parameters_count_optimizable(stlk->sm->m->rates, NULL );
        if( n > 0 ) Parameters_add_parameters(all_params, stlk->sm->m->rates);
        else opt.relative_rates.optimize = false;
	}
    
	if( opt.bl.optimize ){
        Node *right = Node_right( Tree_root(stlk->tree) );
        Node *left  = Node_left( Tree_root(stlk->tree) );
        Node_set_distance(left, Node_distance(right)+Node_distance(left) );
        Node_set_distance(right, 0);
        Parameter_set_fixed(right->distance, true);
        SingleTreeLikelihood_update_one_node(stlk, right);
        SingleTreeLikelihood_update_one_node(stlk, left);
        
        // not sure what it is not set somewhere else
        Parameter_set_fixed( Tree_root(stlk->tree)->distance, true );
        
        for ( int i = 0; i < Tree_node_count(stlk->tree); i++ ) {
            if( !Parameter_fixed(nodes[i]->distance) ){
                Parameters_add(all_params, nodes[i]->distance);
            }
        }
    }
    
    for (int i = 0; i < Parameters_count(stlk->sm->rates); i++) {
        Parameters_add_parameters(all_params, stlk->sm->rates);
        
	}
	
    //Parameters_print(all_params);
	
	
	int max_rounds = opt.max_rounds;
    
    double lnl = data_multivariate->f(data_multivariate);
	
    TopologyOptimizer *topology = NULL;
    int topo_failure = 0;
    
	if ( opt.topology_optimize ) {
        topology = new_TopologyOptimizer( stlk, opt.topology_alogrithm );
    }
    
    
    
	double fret = lnl;
	
	//if(opt.heights.optimize) opt.verbosity = 3;
    
	if ( opt.verbosity > 1 ) {
		fprintf(stdout, "Optimize branch lengths %s\n", (opt.bl.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize frequencies    %s\n", (opt.freqs.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize gamma          %s\n", (opt.gamma.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize pinv           %s\n", (opt.pinv.optimize ? "yes" : "no"));
		fprintf(stdout, "Optimize relative rates %s\n", (opt.relative_rates.optimize ? "yes" : "no"));
		fprintf(stdout, "\n");
	}
    
    if ( stlk->opt.verbosity ==  1 ) {
        fprintf(stdout, "LnL: %f\n", lnl);
    }
    
    time_t start_time;
    time(&start_time);
    
    
    for ( int rounds = 0; rounds < max_rounds; rounds++ ) {
        //fret = 0;
        if ( opt.verbosity > 1  ) {
            fprintf(stdout, "========================================\n");
            fprintf(stdout, "Round %d lk = %f\n\n", rounds, lnl);
            
        }
        
        if ( stlk->opt.verbosity > 0 ) {
            fprintf(stdout, "\n");
        }
        
        opt_optimize( opt_all, all_params, &fret);
        fret = -fret;
        if ( stlk->opt.verbosity > 0 ) {
            fprintf(stdout, "All   LnL: %f {%f}\n", fret, fret-lnl );
        }
        
        
        
        
        if( opt.topology_optimize && topo_failure < 2  ){
            //SingleTreeLikelihood_update_all_nodes(stlk);
            //fprintf(stdout, "Before topology    LnL: %f \n", stlk->calculate(stlk));
            
            fret = topology->optimize(topology);
            
            if ( stlk->opt.verbosity > 0 ) {
                fprintf(stdout, "Topology       LnL: %f NNIs: %d {%f} [%d]\n", fret, topology->moves, fret - lnl, topology->algorithm );
            }
            
            TopologyOptimizer_set_algorithm(topology, 1-topology->algorithm);
            
            if ( topology->moves == 0 ) {
                topo_failure++;
            }
            else {
                topo_failure = 0;
            }
        }
        
        
        if(  fret - lnl < opt.precision ){
            if( opt.verbosity > 1 ){
                fprintf(stdout, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
                fprintf(stdout, "STOP LnL previous %f last %f round: %d\n\n", lnl, fret, rounds);
            }
            
            lnl = fret;
            break;
        }
        
        lnl = fret;
    }
    
    free_MultivariateData( data_multivariate);
    
    if( opt_all != NULL ){
        free_Optimizer(opt_all);
    }
	
	
	if ( topology != NULL ) {
        free_TopologyOptimizer(topology);
    }
    
    
    //    if( stlk->bm!= NULL && Parameters_count(stlk->bm->rates) ==1 && stlk->hessian != NULL ){
    //        double ll = loglikelihood_approximation2(stlk);
    //        printf("Approx Optimize %f\n", ll);
    //        exit(0);
    //    }
	return lnl;
}



double optimize_brent_relative_rate_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = 0;
    double lnl2 = data->f(data);
    
    for (int j = 0; j < 3; j++){
        for (int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
            if( Parameters_fixed( stlk->sm->m->rates, i ) ){
                continue;
            }
            
            data->index_param = i;
            Parameters_set(param, 0, Parameters_at(stlk->sm->m->rates, i));
            
            status = opt_optimize( opt, param, &lnl);
            if( status == OPT_ERROR ) error("OPT.relative_rates No SUCCESS!!!!!!!!!!!!\n");
            
            stlk->sm->m->set_rate( stlk->sm->m, Parameters_value(param, 0), i );
            SingleTreeLikelihood_update_all_nodes(stlk);
            //printf("%f %f\n",Parameters_value(param, 0), -lnl);
        }
        
        if(-lnl - lnl2 < 0.001 ){
            break;
        }
        lnl2 = -lnl;
    }
	return -lnl;
}

double optimize_brent_frequencies_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = 0;
    for (int i = 0; i < Parameters_count(stlk->sm->m->freqs); i++) {
        if( Parameters_fixed( stlk->sm->m->freqs, i ) ){
            continue;
        }
        
        data->index_param = i;
        Parameters_set(param, 0, Parameters_at(stlk->sm->m->freqs, i));
        
        status = opt_optimize( opt, param, &lnl);
        if( status == OPT_ERROR ) error("OPT.frequencies No SUCCESS!!!!!!!!!!!!\n");
        
        stlk->sm->m->set_frequency( stlk->sm->m, Parameters_value(param, 0), i );
        SingleTreeLikelihood_update_all_nodes(stlk);
    }
	return -lnl;
}

double optimize_brent_sm_rates_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param ){
	opt_result status = OPT_NEED_CHECK;
    double lnl = 0;
    for (int i = 0; i < Parameters_count(stlk->sm->rates); i++) {
        
        data->index_param = i;
        Parameters_set(param, 0, Parameters_at(stlk->sm->rates, i));
        
        status = opt_optimize( opt, param, &lnl);
        if( status == OPT_ERROR ) error("OPT.SM.RATES No SUCCESS!!!!!!!!!!!!\n");
        
        stlk->sm->set_rate(stlk->sm, i, Parameters_value(param, 0));
        SingleTreeLikelihood_update_all_nodes(stlk);
    }
	return -lnl;
}


static void _optimize_brent_branch_length_all_aux( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, Node *node, double *lnl ){
    if( node == NULL ) return;
    
    _optimize_brent_branch_length_all_aux(stlk, opt, data, param, Node_left(node), lnl);
    _optimize_brent_branch_length_all_aux(stlk, opt, data, param, Node_right(node), lnl);
    
    if( !Node_isroot(node) ){
        if( Parameter_fixed( node->distance ) ){
            SingleTreeLikelihood_update_all_nodes(stlk);
        }
        else {
            data->index_param = Node_id(node);
            Parameters_set(param, 0, node->distance);
            //SingleTreeLikelihood_update_all_nodes(stlk);
            //SingleTreeLikelihood_update_one_node(stlk, node);
            //do{
                opt_result status = opt_optimize( opt, param, lnl);
                if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
            //}
            //while ( Parameters_is_at_boundry( param, 0, opt_tolx(opt) ) );
                
            Node_set_distance(node, Parameters_value(param, 0));
            //printf("%f %d %d %s\n",*lnl,data->index_param, opt_iterations(opt), Node_name(node));
            
            SingleTreeLikelihood_update_one_node(stlk, node);
            
        }
    }
}

double optimize_brent_branch_length_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, int iterations ){
//	opt_result status = OPT_NEED_CHECK;
//    Node **nodes = Tree_nodes(stlk->tree);
    
    double lnl = data->f(data); //stlk->calculate(stlk);
    double lnl_orginal = lnl;
    
    data->index_param = -1;
    
    Node *right = Node_right( Tree_root(stlk->tree) );
    Node *left  = Node_left( Tree_root(stlk->tree) );
    
    bool left_fixed  = Parameter_fixed(left->distance);
    bool right_fixed = Parameter_fixed(right->distance);
    
    // We don't need to optimize these branches together, at least for a reversible model
    // if at least 1 node is fixed we don't change these 2 branches
    if( !right_fixed && !left_fixed && stlk->sm->m->modeltype < NON_REVERSIBLE_DNA ){
        Node_set_distance(left, Node_distance(right)+Node_distance(left) );
        Node_set_distance(right, 0);
        Parameter_set_fixed(right->distance, true);
    }
    
    
    
    for (int j = 0; j < iterations; j++){
        SingleTreeLikelihood_update_all_nodes(stlk);
        
        _optimize_brent_branch_length_all_aux(stlk, opt, data, param, Tree_root(stlk->tree), &lnl);
        
        if( lnl-0.1 < lnl_orginal ){
            break;
        }
        lnl_orginal = lnl;
    }
    
    if( !right_fixed && !left_fixed && stlk->sm->m->modeltype < NON_REVERSIBLE_DNA ){
        Parameter_set_fixed(right->distance, false);
    }
    
    SingleTreeLikelihood_update_all_nodes(stlk);

    return -lnl;
}



double optimize_brent_rate_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param ){
    opt_result status = OPT_NEED_CHECK;
    double lnl = 0;
    for (int i = 0; i < BranchModel_n_rate(stlk->bm); i++) {
        if( Parameters_fixed(stlk->bm->rates, i) ){
            continue;
        }
        
        data->index_param = i;
        Parameters_set( param, 0,  Parameters_at(stlk->bm->rates, i ) );
        
        status = opt_optimize( opt, param, &lnl);
        if( status == OPT_ERROR ) fprintf(stderr, "OPT.RATES No SUCCESS!!!!!!!!!!!!\n");
        
        stlk->bm->set( stlk->bm, i, Parameters_value(param, 0) );
        SingleTreeLikelihood_update_all_nodes(stlk);
    }
    
    return -lnl;
}


static void _optimize_brent_height_all_post_order_aux( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, Node *node, double *lnl ){
    if( node == NULL ) return;
    
    _optimize_brent_height_all_post_order_aux(stlk, opt, data, param, Node_left(node), lnl);
    _optimize_brent_height_all_post_order_aux(stlk, opt, data, param, Node_right(node), lnl);
    
    if( Parameter_fixed( node->height ) ){
        SingleTreeLikelihood_update_all_nodes(stlk);
    }
    else {
        data->index_param = Node_id(node);
        Parameters_set(param, 0, node->height);
        
        if( Node_isroot(node ) ){
            bool root_upper_constrained = Constraint_fupper( Tree_root(stlk->tree)->height->cnstr );
            
            if ( !root_upper_constrained ) {
                Parameters_set_upper(param, 0, Parameters_value(param, 0)*2.0 );
            }
        }
        
        opt_result status = opt_optimize( opt, param, lnl);
        if( status == OPT_ERROR ) error("OPT.DISTANCE No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_height(node, Parameters_value(param, 0));
        Tree_constraint_heights(stlk->tree);
        
        SingleTreeLikelihood_update_three_nodes(stlk, node );
        
    }
}

// should not be used with approximation
double optimize_brent_height_all( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param ){

    double lnl = 0;
    
    _optimize_brent_height_all_post_order_aux(stlk, opt, data, param, Tree_root(stlk->tree), &lnl);
    
    SingleTreeLikelihood_update_all_nodes(stlk);
    
    return -lnl;
}


// should only be used without approximation
static void _optimize_brent_height_threshold_aux( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, Node *node, double *lnl ){
    if( node == NULL ) return;
    
    _optimize_brent_height_threshold_aux(stlk, opt, data, param, Node_left(node), lnl);
    _optimize_brent_height_threshold_aux(stlk, opt, data, param, Node_right(node), lnl);
    
    // Ignore fixed node (leaf are fixed)
    if( !Parameter_fixed( node->height ) ){
         // We ignore nodes that are <= threshold unless they are root or kids root
        if( ( Node_distance(node) <= data->threshold && !( Node_isroot( node ) || Node_isroot( Node_parent(node) ) ) ) ){
            SingleTreeLikelihood_update_all_nodes(stlk);
        }
        else {
            int level = 0;
            // This condition will be met if both its children are leaves but it sould not change anything
            if ( !Node_isroot( node ) && (Node_distance(Node_left(node)) <= data->threshold || Node_distance(Node_right(node)) <= data->threshold) ) {
                //double lower = -1;
                
                double lower = _get_lower_bound(node, data->threshold, &level);
                
                
                Parameter_set_lower(node->height, lower);
                //printf("++ %s %f l: %f r %f | %e %e\n", Node_name(node), Node_distance(node), Node_distance(Node_left(node)), Node_distance(Node_right(node)), lower, Node_height(Node_right(node)));
            }
            
            data->index_param = Node_id(node);
            Parameters_set( param, 0,  node->height );
            
            if( Node_isroot(node ) ){
                bool root_upper_constrained = Constraint_fupper( Tree_root(stlk->tree)->height->cnstr );
                
                if ( !root_upper_constrained ) {
                    Parameters_set_upper(param, 0, Parameters_value(param, 0)*2.0 );
                }
            }
            
            opt_result status = opt_optimize( opt, param, lnl);
            if( status == OPT_ERROR ) error("OPT.HEIGHTS No SUCCESS!!!!!!!!!!!!\n");
            int level2 = 0;
            _move_and_update_nodes(stlk, node, Parameters_value(param, 0), data->threshold, &level2);
            //printf("nde %s %d %d\n", node->name, level, level2);
            
        }
        Tree_constraint_heights(stlk->tree);
    }
    
}

double optimize_brent_heights_all_threshold( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param ){
    
    double lnl = 0;
    
    _optimize_brent_height_threshold_aux(stlk, opt, data, param, Tree_root(stlk->tree), &lnl);
    
    SingleTreeLikelihood_update_all_nodes(stlk);
    Tree_constraint_heights(stlk->tree);
    
    return -lnl;
}


double optimize_scale_rate_height_strict( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk ){
    
    double lnl = 0;
    double status;
    
    double *backup = dvector(Tree_node_count(tlk->tree));
    Tree_heights_to_vector(tlk->tree, backup);
    double rate = Parameters_value(tlk->bm->rates, 0);
    
    set_scale_bounds(tlk->tree, scaler);
    Parameters_set_upper(scaler, 0, 0.99);
    Parameters_set_value(scaler, 0, 0.9 );
    
    status = opt_optimize( optimizer, scaler, &lnl);
    if( status == OPT_ERROR ) error("OPT.SCALE.HEIGHTS.RATES No SUCCESS!!!!!!!!!!!!\n");
    
    lnl = -lnl;
    if( lnl > templk ){
        SingleTreeLikelihood_scaler(tlk, Tree_root(tlk->tree), Parameters_value(scaler, 0));
        Parameters_set_value(tlk->bm->rates, 0, Parameters_value(tlk->bm->rates,0)/Parameters_value(scaler, 0));
    }
    else {
        Tree_vector_to_heights(backup, tlk->tree);
        Parameters_set_value(tlk->bm->rates, 0, rate);
        lnl = templk;
    }
    SingleTreeLikelihood_update_all_nodes(tlk);
    Tree_constraint_heights(tlk->tree);
    
    free(backup);
    
    return lnl;
}

double optimize_scale_rate_height_strict2( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk ){
    
    double lnl = 0;
    double status;
    
    double *backup = dvector(Tree_node_count(tlk->tree));
    Tree_heights_to_vector(tlk->tree, backup);
    double rate = Parameters_value(tlk->bm->rates, 0);
    
    set_scale_bounds(tlk->tree, scaler);
    Parameters_set_upper(scaler, 0, 0.99);
    Parameters_set_value(scaler, 0, 0.9 );
    
    printf("lower scaler %f uppere scaler %f\n", Parameters_lower(scaler,0), Parameters_upper(scaler,0));
    
    status = opt_optimize( optimizer, scaler, &lnl);
    if( status == OPT_ERROR ) error("OPT.SCALE.HEIGHTS.RATES No SUCCESS!!!!!!!!!!!!\n");
    
    lnl = -lnl;
    if( lnl > templk ){
        SingleTreeLikelihood_scaler(tlk, Tree_root(tlk->tree), Parameters_value(scaler, 0));
        Parameters_set_value(tlk->bm->rates, 0, Parameters_value(tlk->bm->rates,0)*Parameters_value(scaler, 0));
    }
    else {
        Tree_vector_to_heights(backup, tlk->tree);
        Parameters_set_value(tlk->bm->rates, 0, rate);
        lnl = templk;
    }
    SingleTreeLikelihood_update_all_nodes(tlk);
    Tree_constraint_heights(tlk->tree);
    
    free(backup);
    
    return lnl;
}

double optimize_scale_root_rate_height_strict( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk ){

    double lnl = 0;
    double status;
    
    //set_scale_bounds(tlk->tree, scaler);
    Parameters_set_bounds(scaler, 0, 0.05, 0.99);
    Parameters_set_value(scaler, 0, 0.5 );
    
    double *backup = dvector(Tree_node_count(tlk->tree));
    Tree_heights_to_vector(tlk->tree, backup);
    double rate = Parameters_value(tlk->bm->rates, 0);
    
    Node *root = Tree_root(tlk->tree);
    Node *left  = Node_left(root);
    Node *right = Node_right(root);
    
    bool small_left  = (Node_time_elapsed(left) < 1);
    bool small_right = (Node_time_elapsed(right) < 1);
    
    status = opt_optimize( optimizer, scaler, &lnl);
    if( status == OPT_ERROR ) error("OPT.SCALE.HEIGHTS.RATES No SUCCESS!!!!!!!!!!!!\n");
    
    lnl = -lnl;
    
    //printf("=== %f %f %f\n", Parameters_value(scaler, 0), Node_time_elapsed(Node_left(Tree_root(tlk->tree))),Node_time_elapsed(Node_right(Tree_root(tlk->tree))));
    
    if( lnl > templk ){
        Node *root = Tree_root(tlk->tree);
        double dscaler = Parameters_value(scaler, 0);
        
        if ( dscaler > 1) {
            Node_set_height(root, Node_height( root )*dscaler);
            SingleTreeLikelihood_update_all_nodes(tlk);
        }
        else if( dscaler < 1 ){
            double height = Node_height( root );
            
            if(true && (small_left || small_right) ){
                double max_height_son = -1;
                
                if( small_left ){
                    max_height_son = dmax(Node_height( Node_left(left)) , Node_height( Node_right(left)));
                    //bl = (height  - max_height_son)*galer + max_height_son;
                    //Node_set_height( Node_left(root), bl);
                }
                if( small_right ){
                    max_height_son = dmax(max_height_son, dmax(Node_height( Node_left(right)) , Node_height( Node_right(right))));
                    //bl = (height  - max_height_son)*dscaler + max_height_son;
                    //Node_set_height( Node_right(root), bl);
                }
                
                double bl = (height  - max_height_son)*dscaler + max_height_son;
                if( small_left ) Node_set_height( Node_left(root), bl);
                if( small_right ) Node_set_height( Node_right(root), bl);
                Node_set_height( root, bl);
            }
            else {
                double max_height_son = dmax( Node_height( Node_left(root)), Node_height( Node_right(root) ) );
                double bl = (height  - max_height_son)*dscaler + max_height_son;
                Node_set_height(root, bl);
            }
            SingleTreeLikelihood_update_all_nodes(tlk);
        }
        else {
            exit(1);
        }
        
        Parameters_set_value(tlk->bm->rates, 0, Parameters_value(tlk->bm->rates,0)/dscaler);
    }
    else {
        Tree_vector_to_heights(backup, tlk->tree);
        Parameters_set_value(tlk->bm->rates, 0, rate);
        lnl = templk;
        
    }
    
    SingleTreeLikelihood_update_all_nodes(tlk);
    Tree_constraint_heights(tlk->tree);
    
    free(backup);
    
    return lnl;
}

double optimize_scale_rate_height_all( SingleTreeLikelihood *tlk, Optimizer *optimizer, Parameters *scaler, double templk ){

    double lnl = 0;
    double status;
    
    set_scale_bounds(tlk->tree, scaler);
    Parameters_set_upper(scaler, 0, 0.99);
    Parameters_set_value(scaler, 0, 0.9 );
    
    double *backup = dvector(Tree_node_count(tlk->tree));
    double *rates = dvector(Parameters_count(tlk->bm->rates));
    Tree_heights_to_vector(tlk->tree, backup);
    BranchModel_rates_to_vector(tlk->bm, rates);
    
    status = opt_optimize( optimizer, scaler, &lnl);
    if( status == OPT_ERROR ) error("OPT.SCALE.HEIGHTS.RATES No SUCCESS!!!!!!!!!!!!\n");
    
    lnl = -lnl;
    
    if( lnl > templk ){
        SingleTreeLikelihood_scaler(tlk, Tree_root(tlk->tree), Parameters_value(scaler, 0));
        for ( int i = 0; i < BranchModel_n_rate(tlk->bm); i++ ) {
            Parameters_set_value(tlk->bm->rates, i, Parameters_value(tlk->bm->rates, i)/Parameters_value(scaler, 0));
        }
    }
    else {
        Tree_vector_to_heights(backup, tlk->tree);
        BranchModel_vector_to_rates(tlk->bm, rates);
        lnl = templk;
        
    }
    SingleTreeLikelihood_update_all_nodes(tlk);
    
    Tree_constraint_heights(tlk->tree);
    
    free(rates);
    free(backup);
    
    return lnl;
}


// The age of old tips tend to grow more quickly than younger nodes
// Use a polynomial to "correct" it y=0.001x^2+x
double optimize_scale_heights_polynomial( SingleTreeLikelihood *tlk, Optimizer *optimizer, BrentData *data, Parameters *coefficient, double templk ){

    double lnl = 0;
    double status;
    
    double *backup = dvector(Tree_node_count(tlk->tree));
    Tree_heights_to_vector(tlk->tree, backup);
    
    Parameters_set_bounds(coefficient, 0, 0.0001, 0.1);
    Parameters_set_value(coefficient, 0, 0.01 );
    
    status = opt_optimize( optimizer, coefficient, &lnl);
    if( status == OPT_ERROR ) error("OPT.SCALE.HEIGHTS.RATES No SUCCESS!!!!!!!!!!!!\n");
    
    lnl = -lnl;
    
    if( lnl > templk ){
        _scale_heights_poly( Tree_root(tlk->tree), Parameters_value(coefficient, 0) );
        SingleTreeLikelihood_update_all_nodes(tlk);
        Tree_constraint_heights(tlk->tree);
        
        Optimizer *opt_height = new_Optimizer( OPT_BRENT );
        opt_set_max_iteration(opt_height, 100);
        
        opt_set_data(opt_height, data);
        opt_set_objective_function(opt_height, _brent_optimize_height);
        opt_set_tolx(opt_height, 0.001);
        Parameters *oneparameter = new_Parameters(1);
        lnl = optimize_brent_height_all(tlk, opt_height, data, oneparameter);
        free_Parameters_soft(oneparameter);
        free_Optimizer(opt_height);
        
        //        Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
        //        for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        //            if( !Parameter_fixed( nodes[i]->height ) ){
        //                //Node_set_height(nodes[i], _correct_height( Node_height(nodes[i]), Parameters_value(coefficient, 0)) );
        //                Node_set_height(nodes[i], _correct_height( Node_height(nodes[i]), Parameters_value(coefficient, 0)) );
        //            }
        //        }
    }
    else {
        Tree_vector_to_heights(backup, tlk->tree);
        lnl = templk;
        
    }
    SingleTreeLikelihood_update_all_nodes(tlk);
    Tree_constraint_heights(tlk->tree);

    free(backup);
    
    return lnl;
}



static void _count_internal_nodes( Node *node, int *count ){
    if( node == NULL ) return;
    
    _count_internal_nodes(node->left, count);
    _count_internal_nodes(node->right, count);
    if( !Node_isleaf(node) ){
        (*count)++;
    }
    
}

static void _opt_internal_nodes( Node *node, SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param, double *lnl2 ){
    if( node == NULL ) return;
    
    _opt_internal_nodes(node->left, stlk, opt, data, param, lnl2);
    _opt_internal_nodes(node->right, stlk, opt, data, param, lnl2);
    
    if( !Node_isleaf(node) ){
        data->index_param = Node_id(node);
        Parameters_set( param, 0,  node->height );
        
        double status = opt_optimize( opt, param, lnl2);
        if( status == OPT_ERROR ) error("OPT.HEIGHTS No SUCCESS!!!!!!!!!!!!\n");
        
        Node_set_height(node, Parameters_value(param, 0) );
        Tree_constraint_heights(stlk->tree);
        
        SingleTreeLikelihood_update_three_nodes(stlk, node );
    }
}

double optimize_brent_heights_strict_local( SingleTreeLikelihood *stlk, Optimizer *opt, BrentData *data, Parameters *param ){
   
    
    Node **nodes = Tree_get_nodes(stlk->tree, POSTORDER);
    int count = 0;
    
    for (int i = 0; i < Tree_node_count(stlk->tree); i++) {
        if( Parameter_fixed( nodes[i]->height ) ){
            continue;
        }
        count = 0 ;
        _count_internal_nodes(nodes[i], &count);
        
        // should check that taxa do not have the same sampling date
        if ( count < 6 ) {
            continue;
        }
        
        //int previous_node = stlk->node_id;
        if( Node_isroot(nodes[i]) ) TreeLikelihood_set_calculate(stlk, NULL);
        else TreeLikelihood_set_calculate(stlk, nodes[i]);
        double lnl2 = data->f(data);//stlk->calculate(stlk);
        //double temp = lnl2;
        
//        fprintf(stderr, "LnL at node %s %f ", Node_name(nodes[i]), lnl2);
//        
//        if( previous_node == Node_id(Node_left(nodes[i])) || previous_node == Node_id(Node_right(nodes[i])) ){
//            data->index_param = Node_id(nodes[i]);
//            Parameters_set( param, 0,  nodes[i]->height );
//            
//            double status = opt_optimize( opt, param, &lnl2);
//            if( status == OPT_ERROR ) error("OPT.HEIGHTS No SUCCESS!!!!!!!!!!!!\n");
//            
//            Node_set_height(nodes[i], Parameters_value(param, 0) );
//            Tree_constraint_heights(stlk->tree);
//            
//            SingleTreeLikelihood_update_three_nodes(stlk, nodes[i] );
//        }
//        else {
            _opt_internal_nodes(nodes[i], stlk, opt, data, param, &lnl2);
//        }
    
//        fprintf(stderr, " opt %f [%d] %s\n", -lnl2, ( previous_node == Node_id(Node_left(nodes[i])) || previous_node == Node_id(Node_right(nodes[i])) ), (lnl2 > temp ? "": "*"));
    }
    SingleTreeLikelihood_update_all_nodes(stlk);
    Tree_constraint_heights(stlk->tree);
    TreeLikelihood_set_calculate(stlk, NULL);
    return data->f(data);
}

/*****************************************************************************************************/
#pragma mark -
#pragma mark Brent

BrentData * new_BrentData( SingleTreeLikelihood *tlk){
	BrentData *data = (BrentData*)malloc( sizeof(BrentData) );
	assert(data);
	data->index_param = 0;
	data->tlk = tlk;
	data->backup = NULL;
    data->model = NULL;
    data->f = standard_loglikelihood_brent;
    
    data->threshold = 1e-5;
    //data->threshold = -1;
	return data;
}

void free_BrentData( BrentData *data ){
	if( data->backup != NULL ){
		free( data->backup);
	}
	free(data);
	data = NULL;
}


double _brent_optimize_frequency( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
	int index = mydata->index_param;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	stlk->sm->m->set_frequency( stlk->sm->m, Parameters_value(params, 0), index );
	SingleTreeLikelihood_update_all_nodes(stlk);
	
    return fabs(mydata->f(mydata));
}

double _brent_optimize_relative_rate( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	stlk->sm->m->set_rate( stlk->sm->m, Parameters_value(params, 0), mydata->index_param );
	SingleTreeLikelihood_update_all_nodes(stlk);
    double lnl = mydata->f(mydata);
    //printf("  %f %f\n",Parameters_value(params, 0),lnl);
    return fabs(lnl);
}

double _brent_optimize_sm_rates( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
    
    stlk->sm->set_rate(stlk->sm, mydata->index_param, Parameters_value(params, 0));
	SingleTreeLikelihood_update_all_nodes(stlk);
    return fabs(mydata->f(mydata));
}

double optimize_brent_branch_length( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
    Node *node = Tree_node(mydata->tlk->tree, mydata->index_param);
    Node_set_distance( node, Parameters_value(params, 0) );

    SingleTreeLikelihood_update_one_node( mydata->tlk, node);
    
    return fabs(mydata->f(mydata));
}


#pragma mark -
#pragma mark Clock

double _brent_optimize_rate( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
	int index = mydata->index_param;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	stlk->bm->set( stlk->bm, index, Parameters_value(params, 0) );
    
	// strict clock (map is not instanciated)
	if ( Parameters_count(stlk->bm->rates) == 1 ) {
		SingleTreeLikelihood_update_all_nodes(stlk);
	}
	else {
        Node **nodes = Tree_nodes(mydata->tlk->tree );
        
		unsigned *map = stlk->bm->map;
		
		for ( int i = 0; i < Tree_node_count(stlk->tree); i++ ) {
            if( Node_isroot(nodes[i]) ) continue;
            
			if ( map[i] == index ) {
				SingleTreeLikelihood_update_three_nodes(stlk, nodes[i]);
			}
		}
	}
    
    return fabs(mydata->f(mydata));
}



double _brent_optimize_height( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
    
    Node *node = Tree_node(mydata->tlk->tree, mydata->index_param);
    Node_set_height( node, Parameters_value(params, 0) );
	
	SingleTreeLikelihood_update_three_nodes(mydata->tlk, node );
    
    return fabs(mydata->f(mydata));
}


double _brent_optimize_height_threshold( Parameters *params, double *grad, void *data ){
	BrentData *mydata = (BrentData*)data;
    Node *node = Tree_node(mydata->tlk->tree, mydata->index_param);
    int level = 0;
    _move_and_update_nodes(mydata->tlk, node, Parameters_value(params, 0), mydata->threshold, &level);
    //printf("== %e [%e %e]\n\n", Parameters_value(params, 0), Parameters_lower(params, 0), Parameters_upper(params, 0));
    

    return fabs(mydata->f(mydata));
}


void scale_heights_rates( SingleTreeLikelihood*tlk, double scaler){
    
    SingleTreeLikelihood_scaler(tlk, Tree_root(tlk->tree), scaler);
    
    // assume that rates are not constrained
    for ( int j = 0; j < Parameters_count(tlk->bm->rates); j++) {
        Parameters_set_value(tlk->bm->rates, j, Parameters_value(tlk->bm->rates, j)/scaler);
    }
    SingleTreeLikelihood_update_all_nodes(tlk);
}

// works for 1 (strict) and more rates
double _brent_optimize_scale_heights_rates( Parameters *params, double *grad, void *data ){
	BrentData *opt_data = (BrentData*)data;
	
	Node **nodes = Tree_nodes( opt_data->tlk->tree );
	int i = 0;
	for ( ; i < Tree_node_count(opt_data->tlk->tree); i++) {
		opt_data->backup[i] = Node_height(nodes[i]);
	}
	//printf("s: %f height: %f ",Parameters_value(params, 0), Node_height(Tree_root(opt_data->tlk->tree)));
	SingleTreeLikelihood_scaler(opt_data->tlk, Tree_root(opt_data->tlk->tree), Parameters_value(params, 0));
    //printf("%f rate: %f ", Node_height(Tree_root(opt_data->tlk->tree)),Parameters_value(opt_data->tlk->bm->rates, 0));
	
	int j = 0;
	// assume that rates are not constrained
	for ( ; j < Parameters_count(opt_data->tlk->bm->rates); i++, j++) {
		opt_data->backup[i] = Parameters_value(opt_data->tlk->bm->rates, j);
		Parameters_set_value(opt_data->tlk->bm->rates, j, opt_data->backup[i]/Parameters_value(params, 0));
	}
    
    // strict clock (map is not instanciated)
	if ( Parameters_count(opt_data->tlk->bm->rates) == 1 ) {
		SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
	}
    // we only update the branches that depend on the rate
	else {
		unsigned *map = opt_data->tlk->bm->map;
		// root does not have a rate
		for ( int i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
            if( Node_isroot(nodes[i]) ) continue;
            
			if ( map[i] == opt_data->index_param ) {
				SingleTreeLikelihood_update_three_nodes(opt_data->tlk, nodes[i]);
			}
		}
	}
	
	
	double lk = fabs(opt_data->f(opt_data));
    
    //printf("%f lnl: %f\n", Parameters_value(opt_data->tlk->bm->rates, 0), lk);
    
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
		Node_set_height(nodes[i], opt_data->backup[i]);
	}
	
	for ( j = 0; j < Parameters_count(opt_data->tlk->bm->rates); i++, j++ ) {
		Parameters_set_value( opt_data->tlk->bm->rates, j, opt_data->backup[i] );
	}
	
	return lk;
}

// works for 1 (strict) and more rates
double _brent_optimize_scale_heights_rates2( Parameters *params, double *grad, void *data ){
	BrentData *opt_data = (BrentData*)data;
	
	Node **nodes = Tree_nodes( opt_data->tlk->tree );
	int i = 0;
	for ( ; i < Tree_node_count(opt_data->tlk->tree); i++) {
		opt_data->backup[i] = Node_height(nodes[i]);
	}
	
	SingleTreeLikelihood_scaler(opt_data->tlk, Tree_root(opt_data->tlk->tree), Parameters_value(params, 0));
	
	int j = 0;
	// assume that rates are not constrained
	for ( ; j < Parameters_count(opt_data->tlk->bm->rates); i++, j++) {
		opt_data->backup[i] = Parameters_value(opt_data->tlk->bm->rates, j);
		Parameters_set_value(opt_data->tlk->bm->rates, j, opt_data->backup[i]*Parameters_value(params, 0));
	}
    
    // strict clock (map is not instanciated)
	if ( Parameters_count(opt_data->tlk->bm->rates) == 1 ) {
		SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
	}
    // we only update the branches that depend on the rate
	else {
		unsigned *map = opt_data->tlk->bm->map;
		// root does not have a rate
		for ( int i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
            if( Node_isroot(nodes[i]) ) continue;
            
			if ( map[i] == opt_data->index_param ) {
				SingleTreeLikelihood_update_three_nodes(opt_data->tlk, nodes[i]);
			}
		}
	}
	
	
	double lk = fabs(opt_data->f(opt_data));
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
		Node_set_height(nodes[i], opt_data->backup[i]);
	}
	
	for ( j = 0; j < Parameters_count(opt_data->tlk->bm->rates); i++, j++ ) {
		Parameters_set_value( opt_data->tlk->bm->rates, j, opt_data->backup[i] );
	}
	
	return lk;
}

// works for 1 (strict) and more rates
double _brent_optimize_scale_root_height_rate( Parameters *params, double *grad, void *data ){
	BrentData *opt_data = (BrentData*)data;
	
    Node *root  = Tree_root(opt_data->tlk->tree);
    Node *left  = Node_left(root);
    Node *right = Node_right(root);
    
    double root_height  = Node_height(root);
    double left_height  = Node_height(left);
    double right_height = Node_height(right);
    
    double dscaler = Parameters_value(params, 0);
    
    if ( dscaler > 1) {
        Node_set_height(root, root_height*dscaler);
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
    }
    else if( dscaler < 1 ){
        double height = Node_height( root );
        //printf("+++ %f root %f %f %f %f %f\n", dscaler, height, Node_height(Node_left(root)),Node_height(Node_right(root)), Node_time_elapsed(Node_left(root)),Node_time_elapsed(Node_right(root)));
        
        if(true && (Node_time_elapsed(left) < 1 || Node_time_elapsed(right) < 1) ){
            double max_height_son = -1;
            
            if( Node_time_elapsed(left) < 1 ){
                max_height_son = dmax(Node_height( Node_left(left)) ,Node_height( Node_right(left)));
                //bl = (height  - max_height_son)*dscaler + max_height_son;
                //Node_set_height( left, bl);
            }
            if( Node_time_elapsed(right) < 1 ){
                max_height_son = dmax(max_height_son, dmax(Node_height( Node_left(right)) , Node_height( Node_right(right))));
                //bl = (height  - max_height_son)*dscaler + max_height_son;
                //Node_set_height( right, bl);
            }
            double bl = (height  - max_height_son)*dscaler + max_height_son;
            if( Node_time_elapsed(left) < 1 ) Node_set_height( left, bl);
            if( Node_time_elapsed(right) < 1 ) Node_set_height( right, bl);
            Node_set_height(root, bl);
        }
        else {
            double max_height_son = dmax( Node_height( Node_left(root)), Node_height( Node_right(root) ) );
            double bl = (height  - max_height_son)*dscaler + max_height_son;
            Node_set_height(root, bl);
        }
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
    }
    
    Node **nodes = Tree_nodes(opt_data->tlk->tree );
    
	// assume that rates are not constrained
	for ( int j = 0; j < Parameters_count(opt_data->tlk->bm->rates); j++) {
		opt_data->backup[j] = Parameters_value(opt_data->tlk->bm->rates, j);
		Parameters_set_value(opt_data->tlk->bm->rates, j, opt_data->backup[j]/Parameters_value(params, 0));
	}
	
    // strict clock (map is not instanciated)
	if ( Parameters_count(opt_data->tlk->bm->rates) == 1 ) {
		SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
	}
    // we only update the branches that depend on the rate
	else {
		unsigned *map = opt_data->tlk->bm->map;
		// root does not have a rate
		for ( int i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
            if( Node_isroot(nodes[i]) ) continue;
            
			if ( map[i] == opt_data->index_param ) {
				SingleTreeLikelihood_update_three_nodes(opt_data->tlk, nodes[i]);
			}
		}
	}
    
    double lk = fabs(opt_data->f(opt_data));
    //printf("%f\n", lk);
	
	Node_set_height(root, root_height);
	Node_set_height(left, left_height);
	Node_set_height(right, right_height);
	
	for ( int j = 0; j < Parameters_count(opt_data->tlk->bm->rates); j++ ) {
		Parameters_set_value( opt_data->tlk->bm->rates, j, opt_data->backup[j] );
	}
    
    // strict clock (map is not instanciated)
	if ( Parameters_count(opt_data->tlk->bm->rates) == 1 ) {
		SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
	}
    // we only update the branches that depend on the rate
	else {
		unsigned *map = opt_data->tlk->bm->map;
		// root does not have a rate
		for ( int i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
            if( Node_isroot(nodes[i]) ) continue;
            
			if ( map[i] == opt_data->index_param ) {
				SingleTreeLikelihood_update_three_nodes(opt_data->tlk, nodes[i]);
			}
		}
	}
    
	return lk;
}


// Update nodes above calibrations only
static void _update_above_calibrations_only( SingleTreeLikelihood *tlk, Node *node ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( Node_isroot(node) || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
        _update_above_calibrations_only(tlk, Node_left(node));
        _update_above_calibrations_only(tlk, Node_left(node));
        SingleTreeLikelihood_update_three_nodes(tlk, node);
	}
}

double _brent_optimize_scale_heights( Parameters *params, double *grad, void *data ){
	BrentData *opt_data = (BrentData*)data;
	
	Node **nodes = Tree_nodes( opt_data->tlk->tree );
	int i = 0;
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++) {
		opt_data->backup[i] = Node_height(nodes[i]);
	}
	
	SingleTreeLikelihood_scaler(opt_data->tlk, Tree_root(opt_data->tlk->tree), Parameters_value(params, 0));
	SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
    
    double lk = fabs(opt_data->f(opt_data));
	
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++) {
		Node_set_height(nodes[i], opt_data->backup[i]);
	}
    
    _update_above_calibrations_only(opt_data->tlk, Tree_root(opt_data->tlk->tree));
    Tree_constraint_heights(opt_data->tlk->tree);
    
	return lk;
}

// can add negative number too
double _brent_optimize_add_heights( Parameters *params, double *grad, void *data ){
	BrentData *opt_data = (BrentData*)data;
	
	Node **nodes = Tree_nodes( opt_data->tlk->tree );
	int i = 0;
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++) {
		opt_data->backup[i] = Node_height(nodes[i]);
	}
	
	SingleTreeLikelihood_add_height( opt_data->tlk, Tree_root(opt_data->tlk->tree), Parameters_value(params, 0));
	
    double lk = fabs(opt_data->f(opt_data));
	
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++) {
		Node_set_height(nodes[i], opt_data->backup[i]);
	}
    
    _update_above_calibrations_only(opt_data->tlk, Tree_root(opt_data->tlk->tree));
    Tree_constraint_heights(opt_data->tlk->tree);
	
	return lk;
}

double _brent_optimize_scale_heights_polynomial( Parameters *params, double *grad, void *data ){
	BrentData *opt_data = (BrentData*)data;
	
	Node **nodes = Tree_nodes( opt_data->tlk->tree );
	int i = 0;
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++) {
		opt_data->backup[i] = Node_height(nodes[i]);
	}
	
    //	for ( int i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
    //        if( !Parameter_fixed( nodes[i]->height ) ){
    //            Node_set_height(nodes[i], _correct_height( Node_height(nodes[i]), Parameters_value(params, 0)) );
    //        }
    //    }
    _scale_heights_poly( Tree_root(opt_data->tlk->tree), Parameters_value(params, 0) );
    Tree_constraint_heights(opt_data->tlk->tree);
    SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
    
    Optimizer *opt_height = new_Optimizer( OPT_BRENT );
    opt_set_max_iteration(opt_height, 100);
    
    opt_set_data(opt_height, opt_data);
    opt_set_objective_function(opt_height, _brent_optimize_height);
    opt_set_tolx(opt_height, 0.001);
    Parameters *oneparameter = new_Parameters(1);
    double fret = optimize_brent_height_all(opt_data->tlk, opt_height, opt_data, oneparameter);
    free_Parameters_soft(oneparameter);
    free_Optimizer(opt_height);
	
	//double lk = fabs(opt_data->tlk->calculate(opt_data->tlk));
    double lk = fabs(opt_data->f(opt_data));
    
    fprintf(stderr, "root: %f coef: %f LnL: %f Fret: %f\n", Node_height(Tree_root(opt_data->tlk->tree)), Parameters_value(params, 0), lk, fret);
	
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++) {
		Node_set_height(nodes[i], opt_data->backup[i]);
	}
    
    SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
    Tree_constraint_heights(opt_data->tlk->tree);
    
	return lk;
}

// change the height of the node, its rate and the rates of the 2 lineages below
double _brent_optimize_slide_height( Parameters *params, double *grad, void *data ){
	BrentData *opt_data = (BrentData*)data;
	
    Node **nodes = Tree_nodes(opt_data->tlk->tree );
    Node *node = nodes[opt_data->index_param];
    
	opt_data->backup[0] = Node_height(node);
	opt_data->backup[1] = opt_data->tlk->bm->get(opt_data->tlk->bm, node);
	opt_data->backup[2] = opt_data->tlk->bm->get(opt_data->tlk->bm, Node_right(node));
	opt_data->backup[3] = opt_data->tlk->bm->get(opt_data->tlk->bm, Node_left(node));
	
    double scaler = Node_time_elapsed(node)/( Node_height(Node_parent(node))-Parameters_value(params, 0));
    
    Node_set_height(node, Parameters_value(params, 0));
    
	opt_data->tlk->bm->set(opt_data->tlk->bm, opt_data->index_param, opt_data->backup[1]*scaler);
	opt_data->tlk->bm->set(opt_data->tlk->bm, Node_right(node)->postorder_idx, opt_data->backup[2]/scaler);
	opt_data->tlk->bm->set(opt_data->tlk->bm, Node_left(node)->postorder_idx, opt_data->backup[3]/scaler);
    
	SingleTreeLikelihood_update_three_nodes(opt_data->tlk, node);
    
    double lk = fabs(opt_data->f(opt_data));
    
    //fprintf(stderr, "lnl  %f %f %f %f %f\n", -lk, opt_data->backup[0],opt_data->backup[1],opt_data->backup[2],opt_data->backup[3]);
    //fprintf(stderr, "scaler %f %f %f %f %f\n\n", scaler, Node_height(node), opt_data->tlk->bm->get(opt_data->tlk->bm, opt_data->index_param),opt_data->tlk->bm->get(opt_data->tlk->bm, Node_right(node)->postorder_idx),opt_data->tlk->bm->get(opt_data->tlk->bm, Node_left(node)->postorder_idx) );
    
	Node_set_height(node, opt_data->backup[0]);
	opt_data->tlk->bm->set(opt_data->tlk->bm, opt_data->index_param, opt_data->backup[1]);
	opt_data->tlk->bm->set(opt_data->tlk->bm, Node_right(node)->postorder_idx, opt_data->backup[2]);
	opt_data->tlk->bm->set(opt_data->tlk->bm, Node_left(node)->postorder_idx, opt_data->backup[3]);
    
    SingleTreeLikelihood_update_three_nodes(opt_data->tlk, node);
    
	return lk;
}




/*****************************************************************************************************
 * CONJUGATE GRADIENT
 */
#pragma mark -
#pragma mark Conjugate gradient

double _cg_optimize_frequencies( Parameters *params, double *grad, void *data ){
	MultivariateData *mydata = (MultivariateData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	for (int i = 0; i < stlk->sm->nstate-1; i++) {
		stlk->sm->m->set_frequency( stlk->sm->m, Parameters_value(params, i), i );
	}
	SingleTreeLikelihood_update_all_nodes(stlk);
	//double lk = stlk->calculate(stlk);
    double lk = mydata->f(mydata);
	
	if ( grad != NULL ) {
		for (int i = 0; i < stlk->sm->nstate-1; i++) {
			double h = SQRT_EPS*(fabs(Parameters_value(params,i)) + 1.0);
			
			double oldx = Parameters_value(params,i);
			
			stlk->sm->m->set_frequency( stlk->sm->m, oldx + h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			//double fxplus = fabs(stlk->calculate(stlk));
			double fxplus = fabs(mydata->f(mydata));
			
			stlk->sm->m->set_frequency( stlk->sm->m, oldx - h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			//double fxminus = fabs(stlk->calculate(stlk));
			double fxminus = fabs(mydata->f(mydata));
			
			stlk->sm->m->set_frequency( stlk->sm->m, oldx, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			
			// Centered first derivative
			grad[i] = (fxplus-fxminus)/(2.0*h);
			
			//fprintf(stderr, "cg_optimize_frequencies lk=%f fxplus=%f fxminus=%f grad=%f h=%f\n",fabs(lk), fxplus, fxminus, grad[i], h);
		}
	}
	return fabs(lk);
}

double _cg_optimize_branch_length( Parameters *params, double *grad, void *data ){
	MultivariateData *mydata = (MultivariateData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
	Node **nodes = Tree_nodes(stlk->tree);
    
	for (int i = 0; i < Tree_node_count(stlk->tree); i++) {
        Node_set_distance(nodes[i], Parameters_value(params, i));
	}
	SingleTreeLikelihood_update_all_nodes(stlk);
	//double lk = stlk->calculate(stlk);
    double lk = mydata->f(mydata);
	
	if ( grad != NULL ) {
		for (int i = 0; i < Tree_node_count(stlk->tree)-1; i++) {
			double h = SQRT_EPS*(fabs(Parameters_value(params,i)) + 1.0);
			
			double oldx = Parameters_value(params,i);
			
            Node_set_distance(nodes[i], oldx + h);
			SingleTreeLikelihood_update_all_nodes(stlk);
			//double fxplus = fabs(stlk->calculate(stlk));
			double fxplus = fabs(mydata->f(mydata));
			
			Node_set_distance(nodes[i], oldx - h);
			SingleTreeLikelihood_update_all_nodes(stlk);
			//double fxminus = fabs(stlk->calculate(stlk));
			double fxminus = fabs(mydata->f(mydata));
			
			Node_set_distance(nodes[i], oldx);
			SingleTreeLikelihood_update_all_nodes(stlk);
			
			// Centered first derivative
			grad[i] = (fxplus-fxminus)/(2.0*h);
			
            //fprintf(stderr, "_cg_optimize_branch_length\n");
			//fprintf(stderr, "lk=%f fxplus=%f fxminus=%f grad=%f h=%f\n",fabs(lk), fxplus, fxminus, grad[i], h);
		}
	}
	return fabs(lk);
}

double _cg_optimize_rel_rates( Parameters *params, double *grad, void *data ){
	MultivariateData *mydata = (MultivariateData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	for (int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
        stlk->sm->m->set_rate( stlk->sm->m, Parameters_value(params, i), i );
	}
	SingleTreeLikelihood_update_all_nodes(stlk);
	//double lk = stlk->calculate(stlk);
    double lk = mydata->f(mydata);
	
	if ( grad != NULL ) {
		for (int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
			double h = SQRT_EPS*(fabs(Parameters_value(params,i)) + 1.0);
			
			double oldx = Parameters_value(params,i);
			
            stlk->sm->m->set_rate( stlk->sm->m, oldx + h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			double fxplus = fabs(mydata->f(mydata));
			
            stlk->sm->m->set_rate( stlk->sm->m, oldx - h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			double fxminus = fabs(mydata->f(mydata));
			
            stlk->sm->m->set_rate( stlk->sm->m, oldx, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			
			// Centered first derivative
			grad[i] = (fxplus-fxminus)/(2.0*h);
			
			//fprintf(stderr, "cg_optimize_rates value=%f fxplus=%f fxminus=%f grad=%f h=%e\n", oldx, fxplus, fxminus, grad[i], h);
        }
        //fprintf(stderr, "\n");
	}
    //fprintf(stderr, "LnL %f\n",lk);
	
    return fabs(lk);
}

double _cg_optimize_sm_rates( Parameters *params, double *grad, void *data ){
	MultivariateData *mydata = (MultivariateData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	for (int i = 0; i < Parameters_count(stlk->sm->rates); i++) {
        Parameters_set_value(stlk->sm->rates, i, Parameters_value(params, i));
	}
    stlk->sm->need_update = true;
	SingleTreeLikelihood_update_all_nodes(stlk);

    double lk = mydata->f(mydata);
	
	if ( grad != NULL ) {
		for (int i = 0; i < Parameters_count(stlk->sm->rates); i++) {
			double h = SQRT_EPS*(fabs(Parameters_value(params,i)) + 1.0);
			
			double oldx = Parameters_value(params,i);
			
			Parameters_set_value(stlk->sm->rates, i, oldx + h);
            stlk->sm->need_update = true;
			SingleTreeLikelihood_update_all_nodes(stlk);
			double fxplus = fabs(mydata->f(mydata));
			
			Parameters_set_value(stlk->sm->rates, i, oldx - h);
            stlk->sm->need_update = true;
			SingleTreeLikelihood_update_all_nodes(stlk);
			//double fxminus = fabs(stlk->calculate(stlk));
			double fxminus = fabs(mydata->f(mydata));
			
			Parameters_set_value(stlk->sm->rates, i, oldx);
            stlk->sm->need_update = true;
			SingleTreeLikelihood_update_all_nodes(stlk);
			
			// Centered first derivative
			grad[i] = (fxplus-fxminus)/(2.0*h);
			
			//fprintf(stderr, "cg_optimize_frequencies oldx %f lk=%f fxplus=%f fxminus=%f grad=%f h=%e\n",oldx, fabs(lk), fxplus, fxminus, grad[i], h);
		}
	}
	
    return fabs(lk);
}

double _cg_optimize_all( Parameters *params, double *grad, void *data ){
	MultivariateData *mydata = (MultivariateData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
    //Node **nodes = Tree_nodes(stlk->tree);
	
    int index = 0;
	for ( int i = 0; i < Parameters_count(stlk->sm->m->freqs); i++) {
        if( Parameters_fixed(stlk->sm->m->freqs, i) ) continue;
		stlk->sm->m->set_frequency( stlk->sm->m, Parameters_value(params, index++), i );
	}

    for ( int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
        if( Parameters_fixed(stlk->sm->m->rates, i) ) continue;
		stlk->sm->m->set_rate( stlk->sm->m, Parameters_value(params, index++), i );
	}
    
//    for ( int i = 0; i < Tree_node_count(stlk->tree); i++) {
//        if( Parameter_fixed(nodes[i]->distance) ) continue;
//        Node_set_distance(nodes[i], Parameters_value(params, index++));
//        
//	}
//    if(stlk->sm->shape != NULL && !Parameter_fixed(stlk->sm->shape) ){
//        SiteModel_set_alpha( stlk->sm, Parameters_value(params, index++) );
//    }
//    
//    if(stlk->sm->pinv != NULL && !Parameter_fixed(stlk->sm->pinv) ){
//        SiteModel_set_pinv( stlk->sm, Parameters_value(params, index++) );
//    }
    
	SingleTreeLikelihood_update_all_nodes(stlk);
	
    double lk = -mydata->f(mydata);
    
    //printf("LnL %f\n", -lk);
	
	if ( grad != NULL ) {
        index = 0;
        
		for ( int i = 0; i < Parameters_count(stlk->sm->m->freqs); i++ ) {
            if( Parameters_fixed(stlk->sm->m->freqs, i) ) continue;
            
			double h = SQRT_EPS*(fabs(Parameters_value(params, index)) + 1.0);
			
			double oldx = Parameters_value(params, index);
			
			stlk->sm->m->set_frequency( stlk->sm->m, oldx + h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			double fxplus = fabs(mydata->f(mydata));
			
			stlk->sm->m->set_frequency( stlk->sm->m, oldx - h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			double fxminus = fabs(mydata->f(mydata));
			
			stlk->sm->m->set_frequency( stlk->sm->m, oldx, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			
			// Centered first derivative
			grad[index] = (fxplus-fxminus)/(2.0*h);
			
			index++;
		}
        
        for ( int i = 0; i < Parameters_count(stlk->sm->m->rates); i++ ) {
            if( Parameters_fixed(stlk->sm->m->rates, i) ) continue;
            
			//double h = SQRT_EPS*(fabs(Parameters_value(params, index)) + 1.0);
            double h = 0.00001;
			
			double oldx = Parameters_value(params, index);
			
            stlk->sm->m->set_rate( stlk->sm->m, oldx + h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			double fxplus = fabs(mydata->f(mydata));
			
            stlk->sm->m->set_rate( stlk->sm->m, oldx - h, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			double fxminus = fabs(mydata->f(mydata));
			
            stlk->sm->m->set_rate( stlk->sm->m, oldx, i );
			SingleTreeLikelihood_update_all_nodes(stlk);
			
			// Centered first derivative
			grad[index] = (fxplus-fxminus)/(2.0*h);
			
			index++;
		}
        
//        calculate_all_dlnl_dt( stlk, grad+index );
//        for ( int i = 0; i < Tree_node_count(stlk->tree); i++) {
//            if( Parameter_fixed(nodes[i]->distance) ) continue;
//            
//            grad[index] = -grad[index];
//            if( grad[index] > 1000 ){
//                //grad[index] = 0;
//            }
//            //printf("%s %e %f\n", nodes[i]->name, Node_distance(nodes[i]), grad[index]);
//            
//            
////            double h = SQRT_EPS*(fabs(Parameters_value(params, index)) + 1.0);
////            double oldx = Parameters_value(params, index);
////            
////            if( oldx > 1e-5 ){
////                
////                Node_set_distance(nodes[i], oldx + h);
////                SingleTreeLikelihood_update_one_node(stlk, nodes[i]);
////                //SingleTreeLikelihood_update_all_nodes(stlk);
////                double fxplus = fabs(mydata->f(mydata));
////                
////                Node_set_distance(nodes[i], oldx - h);
////                SingleTreeLikelihood_update_one_node(stlk, nodes[i]);
////                //SingleTreeLikelihood_update_all_nodes(stlk);
////                double fxminus = fabs(mydata->f(mydata));
////                
////                Node_set_distance(nodes[i], oldx);
////                SingleTreeLikelihood_update_one_node(stlk, nodes[i]);
////                //SingleTreeLikelihood_update_all_nodes(stlk);
////                
////                // Centered first derivative
////                grad[index] = (fxplus-fxminus)/(2.0*h);
////            }
////            else {
////                // +
////                Node_set_distance( nodes[i], oldx + h );
////                SingleTreeLikelihood_update_one_node(stlk, nodes[i]);
////                double p = fabs(mydata->f(mydata));
////
////                // ++
////                Node_set_distance( nodes[i], oldx + 2*h );
////                SingleTreeLikelihood_update_one_node(stlk, nodes[i]);
////                double pp = fabs(mydata->f(mydata));
////
////
////                Node_set_distance( nodes[i], oldx );
////                SingleTreeLikelihood_update_one_node(stlk, nodes[i]);
////
////                grad[index] = (pp + -2*p + fabs(lk))/(h*h);
////            }
////            printf("%s %e %f %f %d\n", Node_name(nodes[i]), oldx, grad[index], temp[index],(oldx > 1e-5));
//            
//            index++;
//		}
//        
//        if(stlk->sm->shape != NULL && !Parameter_fixed(stlk->sm->shape)){
//            //double h = SQRT_EPS*(fabs(Parameters_value(params, index)) + 1.0);
//            double h = 0.000001;
//			double oldx = Parameters_value(params,index);
//			
//            SiteModel_set_alpha( stlk->sm, oldx + h );
//			SingleTreeLikelihood_update_all_nodes(stlk);
//			double fxplus = -mydata->f(mydata);
//			
//			SiteModel_set_alpha( stlk->sm, oldx - h );
//			SingleTreeLikelihood_update_all_nodes(stlk);
//			double fxminus = -mydata->f(mydata);
//			
//			SiteModel_set_alpha( stlk->sm, oldx );
//			
//			// Centered first derivative
//			grad[index] = (fxplus-fxminus)/(2.0*h);
//            printf("gamma %e %f (%e) %f %f\n", oldx, grad[index], h, fxplus, fxminus );
//            
//            index++;
//        }
//
//        if(stlk->sm->pinv != NULL ){
//            double h = SQRT_EPS*(fabs(Parameters_value(params, index)) + 1.0);
//			
//			double oldx = Parameters_value(params, index);
//			
//            SiteModel_set_pinv( stlk->sm, oldx + h );
//			SingleTreeLikelihood_update_all_nodes(stlk);
//			double fxplus = fabs(mydata->f(mydata));
//			
//			SiteModel_set_pinv( stlk->sm, oldx - h );
//			SingleTreeLikelihood_update_all_nodes(stlk);
//			double fxminus = fabs(mydata->f(mydata));
//			
//			SiteModel_set_pinv( stlk->sm, oldx );
//			
//			// Centered first derivative
//			grad[index] = (fxplus-fxminus)/(2.0*h);
//            
//            index++;
//        }
        
	}
    
//    lk = mydata->f(mydata);
//    
//    printf("LnL %f\n\n", lk);
    SingleTreeLikelihood_update_all_nodes(stlk);
	return lk;
}

double _cg_optimize_rates( Parameters *params, double *grad, void *data ){
	MultivariateData *mydata = (MultivariateData*)data;
	SingleTreeLikelihood *stlk = mydata->tlk;
	
	for (int i = 0; i < Parameters_count(stlk->bm->rates); i++) {
		stlk->bm->set(stlk->bm, i, Parameters_value(params, i) );
	}
	SingleTreeLikelihood_update_all_nodes(stlk);
	//double lk = stlk->calculate(stlk);
    double lk = mydata->f(mydata);
	
	if ( grad != NULL ) {
		for (int i = 0; i < Parameters_count(stlk->bm->rates); i++) {
			double h = SQRT_EPS*(fabs(Parameters_value(params,i)) + 1.0);
			
			double oldx = Parameters_value(params,i);
			
			stlk->bm->set(stlk->bm, i, oldx + h );
			SingleTreeLikelihood_update_all_nodes(stlk);
			//double fxplus = fabs(stlk->calculate(stlk));
			double fxplus = fabs(mydata->f(mydata));
			
			stlk->bm->set(stlk->bm, i, oldx - h );
			SingleTreeLikelihood_update_all_nodes(stlk);
			//double fxminus = fabs(stlk->calculate(stlk));
			double fxminus = fabs(mydata->f(mydata));
			
			stlk->bm->set(stlk->bm, i, oldx );
			SingleTreeLikelihood_update_all_nodes(stlk);
			
			// Centered first derivative
			grad[i] = (fxplus-fxminus)/(2.0*h);
			
			//fprintf(stderr, "cg_optimize_rates lk=%f fxplus=%f fxminus=%f grad=%f h=%f\n",fabs(lk), fxplus, fxminus, grad[i], h);
		}
	}
	return fabs(lk);
}

#ifdef TIMETEST

void calculate_derivatives(MultivariateData *mydata, Node *node, double *grad, int *index){
    
    if( !Node_isleaf(node)){
        calculate_derivatives(mydata,Node_left(node), grad, index);
        calculate_derivatives(mydata,Node_right(node), grad, index);
    }
    
    if( !Node_isroot(node) && Node_right(Node_parent(node)) == node ){
        SingleTreeLikelihood *stlk = mydata->tlk;
        
        double oldx = node->timeParameter->value;
        
        double h = SQRT_EPS*(fabs(oldx) + 1.0);
        printf("x %E h %E\n",oldx, h);
        
        node->timeParameter->value = oldx + h;
        calculate_constrained(Tree_root(stlk->tree),stlk->bm);
        SingleTreeLikelihood_update_all_nodes(stlk);
        double fxplus = fabs(mydata->f(mydata));
        
        node->timeParameter->value = oldx - h;
        calculate_constrained(Tree_root(stlk->tree),stlk->bm);
        SingleTreeLikelihood_update_all_nodes(stlk);
        double fxminus = fabs(mydata->f(mydata));
        
        node->timeParameter->value = oldx;
        calculate_constrained(Tree_root(stlk->tree),stlk->bm);
        SingleTreeLikelihood_update_all_nodes(stlk);
        
        // Centered first derivative
        grad[*index] = (fxplus-fxminus)/(2.0*h);
        ++*index;
    }
}

double _cg_optimize_heightsRate( Parameters *params, double *grad, void *data ){
    MultivariateData *mydata = (MultivariateData*)data;
    SingleTreeLikelihood *stlk = mydata->tlk;
    int nodeCount = Tree_node_count(stlk->tree);
    
    int index = 0;
    printf("strict_parameters_to_heights\n");
    strict_parameters_to_heights(Tree_root(stlk->tree), params, &index);
    
    if( (nodeCount-1)/2 == Parameters_count(params) -1){
        stlk->bm->set(stlk->bm, 0, Parameters_value(params, index));
    }
    calculate_constrained(Tree_root(stlk->tree),stlk->bm);
    SingleTreeLikelihood_update_all_nodes(stlk);
    
    double lk = -mydata->f(mydata);
    
    printf("LnL %f\n", -lk);
    
    if ( grad != NULL ) {
        index = 0;
        printf("derivative heights\n");
        calculate_derivatives(mydata,Tree_root(stlk->tree), grad, &index);
        printf("derivative heights done\n");
        if( (nodeCount-1)/2 == Parameters_count(params) -1){
            double oldx = Parameters_value(params, index);
            double h = SQRT_EPS*(fabs(oldx) + 1.0);
            
            stlk->bm->set(stlk->bm, 0, oldx + h);
            calculate_constrained(Tree_root(stlk->tree),stlk->bm);
            SingleTreeLikelihood_update_all_nodes(stlk);
            double fxplus = fabs(mydata->f(mydata));
            
            stlk->bm->set(stlk->bm, 0, oldx - h);
            calculate_constrained(Tree_root(stlk->tree),stlk->bm);
            SingleTreeLikelihood_update_all_nodes(stlk);
            double fxminus = fabs(mydata->f(mydata));
            
            stlk->bm->set(stlk->bm, 0, oldx);
            calculate_constrained(Tree_root(stlk->tree),stlk->bm);
            SingleTreeLikelihood_update_all_nodes(stlk);
            
            // Centered first derivative
            grad[index] = (fxplus-fxminus)/(2.0*h);
        }
    }
    printf("done\n");
    SingleTreeLikelihood_update_all_nodes(stlk);
    return lk;
}

void set_heights(Node *node){
    double earliest,height;
    Node *kid = NULL;
    
    if( Node_isroot(node)){
        Node_set_height(node, node->timeParameter->value);
    }
    if( !Node_isleaf(node)){
        kid = Node_left(node);
        if( !Node_isleaf(kid) ){
            earliest = 0;
            get_earliest(kid, &earliest);
            height = earliest + (Node_height(node) - earliest)*kid->timeParameter->value;
            Node_set_height(kid, height);
            //printf("%s scaler %e kid height %e parent height %e diff %e\n", Node_name(kid), kid->timeParameter->value, Node_height(kid), Node_height(node), (Node_height(node)-Node_height(kid)) );
        }
        set_heights(Node_left(node));
        
        kid = Node_right(node);
        if( !Node_isleaf(kid) ){
            earliest = 0;
            get_earliest(kid, &earliest);
            height = earliest + (Node_height(node) - earliest)*kid->timeParameter->value;
            Node_set_height(kid, height);
        }
        set_heights(Node_right(node));
    }
}




double _cg_optimize_heightsRate2( Parameters *params, double *grad, void *data ){
    
    int index = 0;
    MultivariateData *mydata = (MultivariateData*)data;
    SingleTreeLikelihood *stlk = mydata->tlk;

    for ( int i = 0 ; i < Tree_node_count(stlk->tree); i++) {
        Node *node = Tree_node(stlk->tree, i);
        if( !Node_isleaf(node) ){
            Node_set_t(node, Parameters_value(params, index++));
        }
    }
    for ( int i = 0; index < Parameters_count(params); index++, i++) {
        stlk->bm->set(stlk->bm, i, Parameters_value(params, index));
    }
    
    set_heights(Tree_root(stlk->tree));
    SingleTreeLikelihood_update_all_nodes(stlk);
    
    double lk = -mydata->f(mydata);
    //printf("LnL %f\n", -lk);
    
    if ( grad != NULL ) {
        index = 0;
        for ( int i = 0; i < Tree_node_count(stlk->tree); i++) {
            Node *node = Tree_node(stlk->tree, i);
            
            if( !Node_isleaf(node)){
                double oldx = node->timeParameter->value;
                
                double h = SQRT_EPS*(fabs(oldx) + 1.0);
                //printf("x %E h %E\n",oldx, h);
                
                node->timeParameter->value = oldx + h;
                set_heights(Tree_root(stlk->tree));
                SingleTreeLikelihood_update_all_nodes(stlk);
                double fxplus = -mydata->f(mydata);
                
                node->timeParameter->value = oldx - h;
                set_heights(Tree_root(stlk->tree));
                SingleTreeLikelihood_update_all_nodes(stlk);
                double fxminus = -mydata->f(mydata);
                
                node->timeParameter->value = oldx;
                set_heights(Tree_root(stlk->tree));
                SingleTreeLikelihood_update_all_nodes(stlk);
                
                // Centered first derivative
                grad[index++] = (fxplus-fxminus)/(2.0*h);
            }
        }
        
        for ( int i = 0; index < Parameters_count(params); index++, i++) {
            double oldx = Parameters_value(params, index);
            double h = SQRT_EPS*(fabs(oldx) + 1.0);
            
            stlk->bm->set(stlk->bm, i, oldx + h);
            SingleTreeLikelihood_update_all_nodes(stlk);
            double fxplus = -mydata->f(mydata);
            
            stlk->bm->set(stlk->bm, i, oldx - h);
            SingleTreeLikelihood_update_all_nodes(stlk);
            double fxminus = -mydata->f(mydata);
            
            stlk->bm->set(stlk->bm, i, oldx);
            SingleTreeLikelihood_update_all_nodes(stlk);
            
            // Centered first derivative
            grad[index] = (fxplus-fxminus)/(2.0*h);
            //printf("rate %e %e g %e f+ %f f- %f\n", oldx,h,grad[index],fxplus,fxminus);
        }
    }
    
    SingleTreeLikelihood_update_all_nodes(stlk);
    return lk;
}
#endif


double _cg_optimize_scale_rate_heights_polynomial( Parameters *params, double *grad, void *data ){
	MultivariateData *opt_data = (MultivariateData*)data;
	
	Node **nodes = Tree_nodes( opt_data->tlk->tree );
	int i = 0;
    double *backup = dvector(Tree_node_count(opt_data->tlk->tree));
	for ( i = 0; i < Tree_node_count(opt_data->tlk->tree); i++) {
		backup[i] = Node_height(nodes[i]);
	}
    double rate = opt_data->tlk->bm->get(opt_data->tlk->bm, 0);
	
    //	for ( int i = 0; i < Tree_node_count(opt_data->tlk->tree); i++ ) {
    //        if( !Parameter_fixed( nodes[i]->height ) ){
    //            Node_set_height(nodes[i], _correct_height( Node_height(nodes[i]), Parameters_value(params, 0)) );
    //        }
    //    }
    _scale_heights_poly( Tree_root(opt_data->tlk->tree), Parameters_value(params, 0) );
    Tree_constraint_heights(opt_data->tlk->tree);
    SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
    
    opt_data->tlk->bm->set(opt_data->tlk->bm, 0, Parameters_value(params, 1));
    
    //Parameters_print(params);
    fprintf(stderr, "coef %f rate %f\n", Parameters_value(params, 0),Parameters_value(params, 1));
    
    double lk = fabs(opt_data->f(opt_data));
    
    Tree_vector_to_heights(backup, opt_data->tlk->tree);
    
    // Compute gradient if needed
	if ( grad != NULL ) {
        //double h = SQRT_EPS*(fabs(Parameters_value(params, 0)) + 1.0);
        double h = 0.001 * Parameters_value(params, 0);
        fprintf(stderr, "graad h %e",h);
        
        double oldx = Parameters_value(params, 0);
        
        _scale_heights_poly( Tree_root(opt_data->tlk->tree), oldx + h );
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
        Tree_constraint_heights(opt_data->tlk->tree);
        double fxplus = fabs(opt_data->f(opt_data));
        
        Tree_vector_to_heights(backup, opt_data->tlk->tree);
        _scale_heights_poly( Tree_root(opt_data->tlk->tree), oldx - h );
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
        Tree_constraint_heights(opt_data->tlk->tree);
        double fxminus = fabs(opt_data->f(opt_data));
        
        Tree_vector_to_heights(backup, opt_data->tlk->tree);
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
        
        grad[0] = (fxplus-fxminus)/(2.0*h);
        
        
        //h = SQRT_EPS*(fabs(Parameters_value(params, 1)) + 1.0);
        h = 0.001 * Parameters_value(params, 1);
        fprintf(stderr, " h2 %e ",h);
        
        oldx = Parameters_value(params, 1);
        
        opt_data->tlk->bm->set(opt_data->tlk->bm, 0, oldx + h);
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
        fxplus = fabs(opt_data->f(opt_data));
        
        opt_data->tlk->bm->set(opt_data->tlk->bm, 0, oldx - h);
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
        fxminus = fabs(opt_data->f(opt_data));
        
        opt_data->tlk->bm->set(opt_data->tlk->bm, 0, oldx );
        SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
        
        grad[1] = (fxplus-fxminus)/(2.0*h);
        fprintf(stderr, " %f %f\n ",fxplus, fxminus);
        fprintf(stderr, "grads %f %f\n", grad[0], grad[1]);
    }
    
    //fprintf(stderr, "root: %f coef: %f LnL: %f Fret: %f\n", Node_height(Tree_root(opt_data->tlk->tree)), Parameters_value(params, 0), lk, fret);//exit(0);
    fprintf(stderr, "root: %f coef: %f rate: %f LnL: %f\n\n", Node_height(Tree_root(opt_data->tlk->tree)), Parameters_value(params, 0), Parameters_value(params, 1), lk);
	
    Tree_vector_to_heights(backup, opt_data->tlk->tree);
    opt_data->tlk->bm->set(opt_data->tlk->bm, 0, rate);
    
    SingleTreeLikelihood_update_all_nodes(opt_data->tlk);
    Tree_constraint_heights(opt_data->tlk->tree);
    
    //fprintf(stderr, "LnL: %f\n\n", fabs(opt_data->f(opt_data)));
    
	return lk;
}


MultivariateData * new_MultivariateData( SingleTreeLikelihood *tlk, int *map ){
	MultivariateData *data = (MultivariateData*)malloc( sizeof(MultivariateData) );
	assert(data);
	data->map_params = map;
	data->tlk = tlk;
    data->f = standard_loglikelihood_multivariate;
    data->model = NULL;
	return data;
}

// model should be freed
void free_MultivariateData( MultivariateData *data ){
	free(data);
	data = NULL;
}

