/*
 *  optimizesa.c
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
#include <math.h>

#include "treelikelihood.h"
#include "optimize.h"
#include "random.h"
#include "matrix.h"

#define BOLTZMANN_CONSTANT (1.3806503e-23)


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

static bool _boltzmann_acceptance( double newfx, double oldfx, double temperature ) {
	if( newfx < oldfx ){
		return true;
	}
	double p = exp( (oldfx-newfx )/(BOLTZMANN_CONSTANT*temperature));
	return p > random_double();
}

static bool _acceptance( double newfx, double oldfx, double temperature ) {
    //printf("%f %f\n",newfx,oldfx);
	if( newfx <= oldfx ){
		return true;
	}
    double random = random_double();
	return exp( (oldfx-newfx )/2) > random;
}

void scale_rate_heights( SingleTreeLikelihood *tlk ){
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

void scale_root( SingleTreeLikelihood *tlk ){
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

void change_one_node( SingleTreeLikelihood *tlk ){
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
    double newheight = (random_double() * (Parameter_upper(node->height) - Parameter_lower(node->height))) + Parameter_lower(node->height);
    Node_set_height(node, newheight);
    SingleTreeLikelihood_update_one_node(tlk, node);
    Tree_constrain_height(node);
    Tree_constrain_neighbor_heights(node);
}

void change_rate( SingleTreeLikelihood *tlk ){
     double scaler = random_double4(0.75, 1.5);
    Parameters_set_value(tlk->bm->rates, 0, Parameters_value(tlk->bm->rates, 0)*scaler);
    SingleTreeLikelihood_update_all_nodes(tlk);
}

double optimize_sa( SingleTreeLikelihood *tlk ){
    double lnl=0;
    double temperature = 1000;
    double cooling_rate = 0.002;
    int max_iterations = 100000000;
    double newlnl;
    double bestlnl;
    double oldlnl = fabs(tlk->calculate(tlk));
    bestlnl = oldlnl;
    double schedules[4] = {0.25,0.25,0.25,0.25};
    //double schedules[4] = {0,0.5,0.0,0.5};
    
    double *heights = dvector(Tree_node_count(tlk->tree));
    Tree_heights_to_vector(tlk->tree, heights);
    double rate = tlk->bm->get(tlk->bm, 0);
    
    for ( int i = 0; i < max_iterations; i++ ) {
        int s = _wheel(schedules, 4);
        
        
        switch ( s ) {
            case 0:
                scale_root(tlk);
                break;
            case 1:
                change_one_node(tlk);
                break;
            case 2:
                scale_rate_heights(tlk);
                break;
            case 3:
                change_rate(tlk);
                break;
            default:
                break;
        }
        
        newlnl = fabs(tlk->calculate(tlk));
        
        
        if( _acceptance(newlnl, oldlnl, temperature) ){
            printf("%d] Accept type: %d LnL new %f old %f p=%f best %f [%f] root %f rate %f\n",i, s, newlnl, oldlnl, exp( (oldlnl-newlnl )/temperature), -bestlnl, temperature, Node_height(Tree_root(tlk->tree)),Parameters_value(tlk->bm->rates, 0) );
            oldlnl = newlnl;
            Tree_heights_to_vector(tlk->tree, heights);
            rate = tlk->bm->get(tlk->bm, 0);
        }
        else {
            //if(i%100 == 0)printf("Reject type: %d LnL new %f old %f best %f [%f]\n", s, newlnl, oldlnl, -bestlnl, temperature);
            Tree_vector_to_heights(heights, tlk->tree);
            tlk->bm->set(tlk->bm, 0, rate);
            SingleTreeLikelihood_update_all_nodes(tlk);
            Tree_constraint_heights(tlk->tree);
        }
        if( newlnl < bestlnl){
            bestlnl = newlnl;
        }
        temperature *= 1.0-cooling_rate;
        temperature = dmax(temperature, 10);
        
    }
    return lnl;
}
