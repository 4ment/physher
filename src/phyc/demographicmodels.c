/*
 *  demographicmodels.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/2/10.
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

#include "demographicmodels.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "tree.h"
#include "mstring.h"
#include "matrix.h"
#include "combinatorics.h"


#pragma mark Coalescent

void free_Coalescent( Coalescent *coalescent ){

    switch ( coalescent->type ) {
        case CONSTANT_DEMOGRAPHY:{
            free_ConstantCoalescent(coalescent);
            break;
        }
        default:
            assert(0);
    }
}

Coalescent * clone_Coalescent( const Coalescent *coal, Tree *tree ){
    Coalescent *coalescent = NULL;
    switch ( coal->type ) {
        case CONSTANT_DEMOGRAPHY:{
            coalescent = clone_ConstantCoalescent(coal, tree);
            break;
        }
        default:
            assert(0);
    }
    return coalescent;
}

#pragma mark -
#pragma mark Constant coalescent

static double _constant_calculate( Coalescent* coal );
static void _update_intervals( Coalescent* coal );

Coalescent * new_ConstantCoalescent( Tree *tree, double theta ){
	Coalescent *coal = (Coalescent*)malloc(sizeof(Coalescent));
    assert(coal);
    coal->p = new_Parameters(1);
    Parameters_move(coal->p, new_Parameter("theta", theta, new_Constraint(0, 200)));
    coal->tree = tree;
	coal->type = CONSTANT_DEMOGRAPHY;
    coal->calculate = _constant_calculate;
    coal->lineages = ivector(Tree_node_count(tree));
    coal->times = dvector(Tree_node_count(tree));
    coal->iscoalescent = bvector(Tree_node_count(tree));
    coal->n = Tree_node_count(tree);
    _update_intervals(coal);
	return coal;
}

void free_ConstantCoalescent( Coalescent *coal ){
	free_Parameters(coal->p);
    free(coal->lineages);
    free(coal->times);
    free(coal->iscoalescent);
	free(coal);
}

Coalescent * clone_ConstantCoalescent( const Coalescent *coal, Tree *tree ){
    Coalescent *newcoal = new_ConstantCoalescent( tree, Parameters_value(coal->p, 0) );
    return newcoal;
}

double get_demographic_constant( const Coalescent *coal, double time ){
	return Parameters_value(coal->p, 0);
}

double get_intensity_constant( const Coalescent *coal, double time ){
	return time/Parameters_value(coal->p, 0);
}

void _post_traversal( Node *node, double lower, double upper, int *lineages, bool *iscoalescent ){
    if( node == NULL ){
        return;
    }
    _post_traversal(node->left, lower, upper, lineages, iscoalescent);
    _post_traversal(node->right, lower, upper, lineages, iscoalescent);
    
    if( !Node_isleaf(node) && Node_height(node) >= upper ){
        if( Node_height( Node_left(node) ) <= lower ){
            //fprintf(stderr, "%s %f %f\n", Node_name(Node_left(node)), Node_height( Node_left(node) ) , lower);
            (*lineages)++;
        }
        if( Node_height( Node_right(node) ) <=lower ){
            //fprintf(stderr, "%s %f %f\n", Node_name(Node_right(node)), Node_height( Node_right(node) ) , lower);
            (*lineages)++;
        }
        if( Node_height(node) == upper ){
            *iscoalescent = true;
        }
    }
}

void _update_intervals( Coalescent* coal ){
    Node **nodes = Tree_get_nodes(coal->tree, POSTORDER);
    memset(coal->times, 0, Tree_node_count(coal->tree)*sizeof(double));
    memset(coal->iscoalescent, 0, Tree_node_count(coal->tree)*sizeof(bool));
    memset(coal->lineages, 0, Tree_node_count(coal->tree)*sizeof(int));
    
    int n = 0;
    for ( int i = 0; i < Tree_node_count(coal->tree); i++ ) {
        if( Node_isleaf(nodes[i]) && Node_height(nodes[i]) != 0.0 ){
            coal->times[n] = Node_height(nodes[i]);
            n++;
        }
    }
    qsort(coal->times, n, sizeof(double), qsort_asc_dvector);
    
    //print_dvector(coal->times, n);

    int len = Tree_tip_count(coal->tree);
    n = 1;
    for ( int i = 1; i < Tree_node_count(coal->tree)-1; i++ ) {
        if( Node_isleaf(nodes[i]) && Node_height(nodes[i]) != 0.0  ){
            if( coal->times[n] == coal->times[n-1] ){
                len--;
                memmove( coal->times+n, coal->times+n+1, len*sizeof(double) );
            }
            else {
                n++;
            }
        }
    }
   
    //print_dvector(coal->times, n);

    for ( int i = 0; i < Tree_node_count(coal->tree); i++ ) {
        if( !Node_isleaf(nodes[i]) ){
            coal->times[n] = Node_height(nodes[i]);
            n++;
        }
    }
    
    qsort(coal->times, n, sizeof(double), qsort_asc_dvector);
    //print_dvector(coal->times, n);
    
    _post_traversal(Tree_root(coal->tree), 0.0, coal->times[0], &coal->lineages[0], &coal->iscoalescent[0]);
    for ( int i = 1; i < n; i++ ) {
        _post_traversal(Tree_root(coal->tree), coal->times[i-1], coal->times[i], &coal->lineages[i], &coal->iscoalescent[i]);
    }
    //print_ivector(coal->lineages, n);

    for ( int i = n-1; i > 0; i-- ) {
        coal->times[i] = coal->times[i] - coal->times[i-1];
    }
    coal->n = n;
    coal->need_update = false;
}

double _constant_calculate( Coalescent* coal ){
    if ( coal->need_update ) {
        _update_intervals(coal);
    }
	double lnl = 0;
    double theta = Parameters_value(coal->p, 0);
 
	for( int i = 0; i< coal->n; i++  ){
		double lambda = choose(coal->lineages[i], 2)/theta;
        lnl -= lambda * coal->times[i];
        
        if( coal->iscoalescent[i] ){
            lnl -= log(lambda);
        }
        
	}
    return lnl;
}


exponentialdemography * new_exponentialdemography( double n0, double r ){
	exponentialdemography *cd = (exponentialdemography*)malloc(sizeof(exponentialdemography));	
	assert(cd);
    cd->name = String_clone("exponential.demography");	
	cd->type = CONSTANT_DEMOGRAPHY;
	cd->n0 = n0;
	cd->r = r;
	cd->lower = 0;
	cd->upper = INFINITY;
	return cd;
}

void free_exponentialdemography( exponentialdemography *demo ){
	free(demo->name);
	free(demo);
}

double get_demographic_exp( exponentialdemography * demo, double time ){
	if( demo->r == 0. ) return demo->n0;
	else return demo->n0 * exp( -time * demo->r);
}

double get_intensity_exp( exponentialdemography * demo, double time ){
	if( demo->r == 0. ) return time/demo->n0;
	else return exp( time * demo->r)/demo->n0/demo->r;
}
