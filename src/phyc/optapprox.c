/*
 *  optapprox.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 23/01/2014.
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

#include "optapprox.h"
#include "matrix.h"
#include "optimize.h"

// Use a Taylor expansion of the likelihood function using a diagonal variance matrix
double loglikelihood_approximation( SingleTreeLikelihood *tlk  ){
    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
    double lnl = 0;
    double d = 0;
    Node *sibling = NULL;
    Node *node = NULL;
    Node *parent = NULL;
    
    for ( int i = 0; i < Tree_node_count(tlk->tree)-2; i++ ) {
        sibling = Node_sibling(nodes[i]);
        node = nodes[i];
        parent = Node_parent(node);
        
        if( Node_isroot(parent) ){
            double root_rate = ( Parameters_value(tlk->bm->rates, tlk->bm->map[i])*Node_time_elapsed(sibling)  + Parameters_value(tlk->bm->rates, tlk->bm->map[sibling->postorder_idx])*Node_time_elapsed(node)  ) / (Node_time_elapsed(node)+Node_time_elapsed(sibling) );
            d = Node_distance(node)+ Node_distance(Node_sibling(node));
            d -= (Node_time_elapsed(node)    * ( Parameters_value(tlk->bm->rates, tlk->bm->map[i])    + root_rate )*0.5);
            d -= (Node_time_elapsed(sibling) * ( Parameters_value(tlk->bm->rates, tlk->bm->map[sibling->postorder_idx]) + root_rate )*0.5);
        }
        else {
            d = Node_distance(node) - (Node_time_elapsed(node) * ( Parameters_value(tlk->bm->rates, tlk->bm->map[i])+Parameters_value(tlk->bm->rates, tlk->bm->map[parent->postorder_idx]) )*0.5);
        }
        lnl += d*d*tlk->hessian[i]/2.0;
        //printf("d %f var %f\n",d,tlk->hessian[i]);
    }
    return lnl;
}

// Use a Taylor expansion of the likelihood function using a diagonal variance matrix
double loglikelihood_approximation2( SingleTreeLikelihood *tlk  ){
    Node **nodes = Tree_nodes(tlk->tree);
    double lnl = 0;
    double d = 0;
    Node *sibling = NULL;
    Node *node = NULL;
    Node *parent = NULL;
    
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        node = nodes[i];
        
        if ( Node_isroot(node) || (Node_isroot(Node_parent(node)) && Node_right(Node_parent(node)) == node )) {
            continue;
        }
        
        sibling = Node_sibling(nodes[i]);
        
        parent = Node_parent(node);
        //if( Node_distance(node) < 1e-7 ) continue;
        //if( Node_issticky(node) ) continue;
        if( tlk->hessian[i] < 0 ) {
        
            if( Node_isroot(parent) ){
                d  = Node_distance(node) + Node_distance(Node_sibling(node));
                d -= Node_time_elapsed(node)    * tlk->bm->get(tlk->bm, node);
                d -= Node_time_elapsed(sibling) * tlk->bm->get(tlk->bm, sibling);
            }
            else {
                d = Node_distance(node) - (Node_time_elapsed(node) * tlk->bm->get(tlk->bm, node));
            }
            //lnl += d*d*tlk->hessian[i]/2.0;
            lnl += d*d*tlk->hessian[i];
            //printf("%s d %f f'' %f bl %e LnL %f\n",Node_name(node), d,tlk->hessian[count], Node_distance(node), lnl);
        }
    }
    lnl *= 0.5;
    //printf("LnL %f LnL_mle %f\n",lnl,tlk->lnl_bl);
    
    return lnl + tlk->lnl_bl;
}

// Use a Taylor expansion of the likelihood function using a  variance-covariance matrix
double loglikelihood_approximation3( SingleTreeLikelihood *tlk  ){
    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
    double lnl = 0;
    double d = 0;
    Node *sibling = NULL;
    Node *node = NULL;
    Node *parent = NULL;
    
    int dim = Tree_node_count(tlk->tree)-2;
    double *bs = dvector(dim);
    
    for ( int i = 0; i < dim; i++ ) {
        sibling = Node_sibling(nodes[i]);
        node = nodes[i];
        parent = Node_parent(node);
        
        if( Node_isroot(parent) ){
            d = Node_distance(node)+ Node_distance(Node_sibling(node));
            d -= Node_time_elapsed(node)    * tlk->bm->get(tlk->bm, node);
            d -= Node_time_elapsed(sibling) * tlk->bm->get(tlk->bm, sibling);
        }
        else {
            d = Node_distance(node) - (Node_time_elapsed(node) * tlk->bm->get(tlk->bm, node));
        }
        bs[i] = d;
    }
    double *temp = Matrix_mult(bs, tlk->hessian, 1, dim, dim, dim);
    
    for ( int i = 0; i < dim; i++ ) {
        lnl += temp[i] * bs[i];
    }
    lnl *= 0.5;
    
    free(bs);
    free(temp);
    
    return lnl + tlk->lnl_bl;
}

#pragma mark -

double loglikelihood_approximation_diag( void *data  ){
    SingleTreeLikelihood *tlk = ((BrentData*)data)->tlk;
    return loglikelihood_approximation2(tlk);
}

double loglikelihood_approximation_full( void *data  ){
    SingleTreeLikelihood *tlk = ((BrentData*)data)->tlk;
    return loglikelihood_approximation3(tlk);
}

