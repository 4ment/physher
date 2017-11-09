/*
 *  optimizega.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 8/10/13.
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

#include "optimize.h"

#include "treelikelihood.h"
#include "matrix.h"
#include "ga.h"
#include "random.h"

static void   _strict_mutate( GA *ga, Individual *individual );
static double _strict_fitness( GA *ga, Individual *individual );


// array has to sum to 1
static int _wheel( double *array, int len ){
    double accum = 0.0;
    double rnum = random_double();
    int i = 0;
    for ( ; i < len; i++ ) {
        accum += array[i];
        if( accum >= rnum ) break;
    }
    return i;
}

static void _scale_rate_heights( SingleTreeLikelihood *tlk ){
    Tree *tree = tlk->tree;
	const double low = 0.5;
	const double high = 1.5;
	
	double min = DBL_MIN;
	double max = DBL_MAX;
	bool isCalibrated = false;
	
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	Parameter *p = NULL;
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		p = nodes[i]->height;
		if ( !Parameter_estimate(p)) continue;
		
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
    
    double scaler;
    
	if (isCalibrated) {
        scaler = random_double4(min, max);
	}
	else {
        scaler = random_double4(low, high);
	}
	Tree_scale_heights(Tree_root(tree), scaler);
    Parameters_set_value(tlk->bm->rates, 0, Parameters_value(tlk->bm->rates, 0)/scaler);
    SingleTreeLikelihood_update_all_nodes(tlk);
    Tree_constraint_heights(tree);
}

static void _scale_root( SingleTreeLikelihood *tlk ){
    Tree *tree = tlk->tree;
    Node *root = Tree_root(tree);
    double max_down = dmax(Node_height(root->right), Node_height(root->left));
    double scaler = random_double4(dmin(Node_height(root)/max_down, 0.75), 1.5);
    double newheight = scaler*Node_height(root);
    Node_set_height(root, newheight);
    SingleTreeLikelihood_update_one_node(tlk, root);
    Tree_constrain_height(root);
    Tree_constrain_neighbor_heights(root);
}

static void _change_one_node( SingleTreeLikelihood *tlk ){
    Tree *tree = tlk->tree;
    int index = random_int(Tree_tip_count(tree)-3);
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    int count = 0;
    Node *node = NULL;
    for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
        if( !Node_isleaf(nodes[i]) ){
            if( count == index){
                node = nodes[i];
                break;
            }
            count++;
        }
    }
    //Parameter_print(node->height);
    double newheight = (random_double() * (Parameter_upper(node->height) - Parameter_lower(node->height))) + Parameter_lower(node->height);
    //printf("new %f\n", newheight);
    Node_set_height(node, newheight);
    SingleTreeLikelihood_update_one_node(tlk, node);
    Tree_constrain_height(node);
    Tree_constrain_neighbor_heights(node);
}

static void _change_rate( SingleTreeLikelihood *tlk ){
    double scaler = random_double4(0.75, 1.5);
    Parameters_set_value(tlk->bm->rates, 0, Parameters_value(tlk->bm->rates, 0)*scaler);
    SingleTreeLikelihood_update_all_nodes(tlk);
}

static void _init_pop( GA *ga ){
    SingleTreeLikelihood *tlk = (SingleTreeLikelihood*)ga->data;
    
    Node **nodes = Tree_get_nodes( tlk->tree, POSTORDER );
    
    double *chromosome = ga->individuals[0]->chromosome;
    int index = 0;
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        if( Parameter_estimate( nodes[i]->height) ){
            chromosome[index++] = Node_height( nodes[i] );
        }
    }
    chromosome[index] = Parameters_value(tlk->bm->rates, 0);
    
    for ( int j = 1; j < ga->pop_size; j++) {
        ga->copy_chromosome( ga->individuals[j]->chromosome, chromosome, ga->chromosome_size);
        _strict_mutate(ga, ga->individuals[j]);
        //_strict_fitness(ga, ga->individuals[j]);
    }
}

double optimize_ga( SingleTreeLikelihood *tlk ){
    double lnl=0;
    
    Tree *tree = tlk->tree;
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    int variable_heights = 0;
    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        if ( Parameter_estimate( nodes[i]->height) ) {
            variable_heights++;
        }
    }

    GA *ga = new_GA( GA_ROULETTE, GA_CHROMOSOME_DOUBLE, 1000, variable_heights+1 );
    
    ga_set_nelites(ga, 10);
    ga_set_ngeneration(ga, 10000);
    
    ga->mutate = _strict_mutate;
    ga->fitness = _strict_fitness;
    ga->mate    = NULL;
    
    ga->data = tlk;
    ga->use_max_no_improvement = true;
    ga->max_no_improvement = 100;
	
    ga->feedback_computeall = true;
    ga->verbosity = 1;
    ga->feedback_sampling = 1;
    
    _init_pop(ga);
    
    ga_evolve(ga);
    
    lnl = ga->maxfitness;
    
    free_GA(ga);
    return lnl;
}

void _strict_mutate( GA *ga, Individual *individual ){
    
    double *chromosome = individual->chromosome;
    
    SingleTreeLikelihood *tlk = (SingleTreeLikelihood*)ga->data;
    
    Node **nodes = Tree_get_nodes( tlk->tree, POSTORDER );
    int index = 0;
	for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        if( Parameter_estimate( nodes[i]->height) ){
            Node_set_height( nodes[i], chromosome[index] );
            index++;
        }
	}
    Parameters_set_value(tlk->bm->rates, 0, chromosome[index]);
    Tree_constraint_heights(tlk->tree);

    double schedules[4] = {0.25,0.25,0.25,0.25};
    int s = _wheel(schedules, 4);
    //fprintf(stderr,"Operator %d\n",s);
    switch ( s ) {
        case 0:
           // _scale_root(tlk);
            break;
        case 1:
            _change_one_node(tlk);
            break;
        case 2:
            _scale_rate_heights(tlk);
            break;
        case 3:
            _change_rate(tlk);
            break;
        default:
            break;
    }
    
    index = 0;
	for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        if( Parameter_estimate( nodes[i]->height) ){
            chromosome[index] = Node_height( nodes[i] );
            index++;
        }
	}
    chromosome[index] = Parameters_value(tlk->bm->rates, 0);

}

double _strict_fitness( GA *ga, Individual *individual ){
    
    double *chromosome = individual->chromosome;
    
    SingleTreeLikelihood *tlk = (SingleTreeLikelihood*)ga->data;
    
    Node **nodes = Tree_get_nodes( tlk->tree, POSTORDER );
    int index = 0;
	for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        if( Parameter_estimate( nodes[i]->height) ){
            Node_set_height( nodes[i], chromosome[index] );
            index++;
        }
	}
    Parameters_set_value(tlk->bm->rates, 0, chromosome[index]);
    
    SingleTreeLikelihood_update_all_nodes( tlk );
    individual->fitness = tlk->calculate(tlk);
    
    //fprintf( stderr, "LnL=%f Height=%f Rate=%e\n", individual->fitness, Node_height(Tree_root(tlk->tree)), tlk->bm->get(tlk->bm, 0));
    if( individual->fitness > ga->maxfitness ){
        ga->maxfitness = individual->fitness;
        ga->copy_chromosome(ga->best_chromosome, chromosome, ga->chromosome_size);
    }
	
	return individual->fitness;
}


