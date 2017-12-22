/*
 *  nj.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/08/13.
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

#include "nj.h"

#include <stdio.h>
#include <assert.h>

#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif

#include "tree.h"
#include "node.h"
#include "matrix.h"


static Node * cluster( Node **nodes, double **matrix, int nTips );

struct _Tree * new_NJ2( const char **taxa, size_t dim, double **matrix ){
    assert(dim > 3);
    Node **nodes = (Node**)malloc(sizeof(Node*)*dim);
    assert(nodes);
    
    int i = 0;
    for ( ; i < dim; i++ ) {
        nodes[i] = new_Node(NULL, taxa[i], i);
    }
    
    Node *root = cluster(nodes, matrix, dim);
    free(nodes);
    
    return new_Tree2(root, false);
}


// i row index
// j column index
void net_divergence( Node **nodes, double **matrix, int N ){
    int j = 0;
    for( int i = 0; i < N; i++ ){
        if( nodes[i] != NULL ){
            matrix[i][i] = 0;
            // horizontal
            for( j = i+1; j < N; j++ ){
                if( nodes[j] != NULL ){
                    matrix[i][i] += matrix[j][i];
                }
            }
            
            // vertical
            for( j = 0; j < i; j++ ){
                if( nodes[j] != NULL ){
                    matrix[i][i] += matrix[i][j];
                }
            }
        }
    }
}


// Update distance matrix
void redo_distance( Node **nodes, double **matrix, int ii, int jj, int N ){
    int i = 0;
    // horizontal
    for( i = jj+1; i < N; i++ ){
        if( nodes[i] != NULL && i != ii ){
            if( jj < i && i < ii ){
                matrix[ii][i] = ( matrix[ii][i] + matrix[i][jj] - matrix[ii][jj])/2;
            }
            else{
                matrix[i][ii] = ( matrix[i][ii] + matrix[i][jj] -  matrix[ii][jj])/2;
            }
        }
    }
    
    // vertical
    for( i = 0; i < jj; i++ ){
        if( nodes[i] != NULL ){
            matrix[ii][i] = (matrix[ii][i] + matrix[jj][i] - matrix[ii][jj])/2;
        }
    }
    
}

// i and j are the indices of the nodes to cluster
void joinThem( Node **nodes, double **matrix, int i, int j, int NN){
    Node *node = new_Node(NULL, NULL, NN);
    
    Node_addChild(node, nodes[i]);
    Node_addChild(node, nodes[j]);
    
    Node_set_parent(nodes[i], node);
    Node_set_parent(nodes[j], node);
    
    double left_bl = (matrix[i][j]*0.5) + ((matrix[i][i] - matrix[j][j])/(2*(NN-2)));
    
    if( left_bl < 0 ){
        Node_set_distance(nodes[i], matrix[i][j]);
        Node_set_distance(nodes[j], BL_MIN);
    } else{
        Node_set_distance(nodes[i], left_bl);
        Node_set_distance(nodes[j], matrix[i][j] - left_bl);
        if( matrix[i][j] - left_bl < 0 ){
            Node_set_distance(nodes[j], BL_MIN);
            Node_set_distance(nodes[i], matrix[i][j]);
        }
    }
    
    nodes[i] = node;
}

Node * last_join( Node **nodes, double **matrix, int N){
    double min = INFINITY;
    int mini = 0;
    int minj = 0;

    int j = 0;
    for(int i = 0; i < N; i++ ){
        for( j = 0; j < i; j++ ){
            // look in the upper part of the matrix: rate corrected
            if( matrix[j][i] < min && nodes[i] != NULL && nodes[j] != NULL ){
                mini = i; //lower part
                minj = j; //lower part
                min = matrix[j][i];
            }
        }
    }
    
    Node *node = new_Node(NULL, NULL, 0);
    
    Node_addChild(node, nodes[minj]);
    Node_addChild(node, nodes[mini]);
    
    Node_set_parent(nodes[mini], node);
    Node_set_parent(nodes[minj], node);
    
    Node_set_distance(nodes[mini], matrix[mini][minj] * 0.5);
    Node_set_distance(nodes[minj], matrix[mini][minj] - Node_distance(nodes[mini]));    
    
    nodes[mini]= node;
    nodes[minj] = NULL;
    
    return node;
}

// Calculate the corrected distance matrix Q for nodes that are available for clustering
// matrix[i][i] and matrix[j][j] are net divergences for node i and j
// Q is stored in the lower triangular of matrix
// Q(i,j) = d(i,j) - (sum_k d(i,k) + sum_k d(j,k))
void rate_corrected_matrix( Node **nodes, double **matrix, int N, int NN ){
    int j = 0;
    int NN2 = NN-2;
    for( int i = 0; i < N; i++ ){
        for( j = 0; j < i; j++ ){
            if( nodes[i] != NULL && nodes[j] != NULL ){
                matrix[j][i] = matrix[i][j] - ((matrix[i][i] + matrix[j][j])/NN2);
            }
        }
    }
}

Node * cluster( Node **nodes, double **matrix, int nTips ){
    double min = 0.0;
    int mini,minj;
    
    int i = 0;
    int j = 0;
    int NN = nTips;
    while( NN - 2 != 0 ){
        net_divergence(nodes, matrix, nTips);
        rate_corrected_matrix(nodes, matrix, nTips, NN);
        
        min = INFINITY;
        mini = 0;
        minj = 0;
        for( i = 0; i < nTips; i++ ){
            for( j = 0; j < i; j++ ){
                if( matrix[j][i] < min &&  nodes[i] != NULL && nodes[j] != NULL ){
                    mini = i;
                    minj = j;
                    min  = matrix[j][i];
                }
            }
        }
        joinThem( nodes, matrix, mini, minj, NN);
        redo_distance( nodes, matrix, mini, minj, nTips);
        
        nodes[minj] = NULL;
        
        NN--;       
    }
    return last_join(nodes, matrix, nTips);
}

void findMinIndexes( double **matrix, int ncluster, double *r, int *alias, int *imin, int *jmin){
    int i,j;
    double sij;
    double min = INFINITY;
    *imin = 0;
    *jmin = 0;
    double denom = 1.0/(ncluster-2);
    for( i = 0; i < ncluster; i++ ){
        for( j = i+1; j < ncluster; j++ ){
            sij = matrix[ alias[i] ][ alias[j] ] - (r[i] + r[j] ) * denom;
            
            if( sij < min ){
                *imin = i;
                *jmin = j;
                min  = sij;
            }
        }
    }
}

struct _Tree * new_NJ( const char **taxa, size_t dim, double **matrix ){
    Node **nodes = (Node**)malloc(sizeof(Node*)*dim);
    assert(nodes);
    int *alias = ivector(dim);
    double *r = dvector(dim);
    
    int imin=0;
    int jmin=0;
    
    int i = 0;
    for ( ; i < dim; i++ ) {
        nodes[i] = new_Node(NULL, taxa[i], i);
        alias[i] = i;
    }
    int ncluster = dim;
    
    while( ncluster > 2 ){
        // calculate net divergence
        for ( int i = 0; i < ncluster; i++ ) {
            r[i] = 0;
            for ( int j = 0; j < ncluster; j++ ) {
                r[i] += matrix[ alias[i] ][ alias[j] ];
            }
        }
        
        findMinIndexes(matrix, ncluster, r,alias, &imin, &jmin);
        
        Node *node = new_Node(NULL, NULL, ncluster);
        
        Node *inode = nodes[alias[imin]];
        Node *jnode = nodes[alias[jmin]];
        
        Node_addChild(node, inode);
        Node_addChild(node, jnode);
        
        Node_set_parent(inode, node);
        Node_set_parent(jnode, node);
        
        double il = (matrix[ alias[imin] ][ alias[jmin] ] + (r[imin] - r[jmin])/(ncluster-2))*0.5;
        double jl = matrix[ alias[imin] ][ alias[jmin] ]-il;
        
        Node_set_distance(inode, dmax(0.0, il));
        Node_set_distance(jnode, dmax(0.0, jl));
        
        nodes[alias[imin]] = node;
        nodes[alias[jmin]] = NULL;
        
        int k = 0;
        for ( ; k < imin; k++) {
            matrix[alias[k]][alias[imin]] = matrix[alias[imin]][alias[k]] = (matrix[ alias[k] ][ alias[imin] ] + matrix[ alias[k] ][ alias[jmin] ] - matrix[ alias[imin] ][ alias[jmin] ]) * 0.5;
        }
        for ( k++; k < jmin; k++) {
            matrix[alias[k]][alias[imin]] = matrix[alias[imin]][alias[k]] = (matrix[ alias[k] ][ alias[imin] ] + matrix[ alias[k] ][ alias[jmin] ] - matrix[ alias[imin] ][ alias[jmin] ]) * 0.5;
        }
        for ( k++; k < ncluster; k++) {
            matrix[alias[k]][alias[imin]] = matrix[alias[imin]][alias[k]] = (matrix[ alias[k] ][ alias[imin] ] + matrix[ alias[k] ][ alias[jmin] ] - matrix[ alias[imin] ][ alias[jmin] ]) * 0.5;
        }
        
        if( ncluster-jmin-1 != 0 ){
            memmove(&alias[jmin], &alias[jmin+1], sizeof(int)*(ncluster-jmin-1));
        }
        ncluster--;
    }
    
    Node *node = new_Node(NULL, NULL, ncluster);
    
    Node *inode = nodes[ alias[0] ];
    Node *jnode = nodes[ alias[1] ];
    
    Node_addChild(node, inode);
    Node_addChild(node, jnode);
    
    Node_set_parent(inode, node);
    Node_set_parent(jnode, node);
    
    double l = dmax(0.0, matrix[ alias[0] ][ alias[1] ]*0.5);
    
    Node_set_distance(inode, l);
    Node_set_distance(jnode, l);
    
    free(nodes);
    free(alias);
    free(r);
    
    return new_Tree2(node, false);
}

void findMinIndexes_float( float **matrix, int ncluster, float *r, int *alias, int *imin, int *jmin){
    int i,j;
    float sij;
    float min = INFINITY;
    *imin = 0;
    *jmin = 0;
    float denom = 1.0/(ncluster-2);
    for( i = 0; i < ncluster; i++ ){
        for( j = i+1; j < ncluster; j++ ){
            sij = matrix[ alias[i] ][ alias[j] ] - (r[i] + r[j] ) * denom;
            
            if( sij < min ){
                *imin = i;
                *jmin = j;
                min  = sij;
            }
        }
    }
}

struct _Tree * new_NJ_float( const char **taxa, size_t dim, float **matrix ){
    Node **nodes = (Node**)malloc(sizeof(Node*)*dim);
    assert(nodes);
    int *alias = ivector(dim);
    float *r = fvector(dim);
    
    int imin=0;
    int jmin=0;
    
    int i = 0;
    for ( ; i < dim; i++ ) {
        nodes[i] = new_Node(NULL, taxa[i], i);
        alias[i] = i;
    }
    int ncluster = dim;
    
    while( ncluster > 2 ){
        // calculate net divergence
        for ( int i = 0; i < ncluster; i++ ) {
            r[i] = 0;
            for ( int j = 0; j < ncluster; j++ ) {
                r[i] += matrix[ alias[i] ][ alias[j] ];
            }
        }
        
        findMinIndexes_float(matrix, ncluster, r,alias, &imin, &jmin);
        
        Node *node = new_Node(NULL, NULL, ncluster);
        
        Node *inode = nodes[alias[imin]];
        Node *jnode = nodes[alias[jmin]];
        
        Node_addChild(node, inode);
        Node_addChild(node, jnode);
        
        Node_set_parent(inode, node);
        Node_set_parent(jnode, node);
        
        float il = (matrix[ alias[imin] ][ alias[jmin] ] + (r[imin] - r[jmin])/(ncluster-2))*0.5;
        float jl = matrix[ alias[imin] ][ alias[jmin] ]-il;
        
        Node_set_distance(inode, fmaxf(0.0, il));
        Node_set_distance(jnode, fmaxf(0.0, jl));
        
        nodes[alias[imin]] = node;
        nodes[alias[jmin]] = NULL;
        
        int k = 0;
        for ( ; k < imin; k++) {
            matrix[alias[k]][alias[imin]] = matrix[alias[imin]][alias[k]] = (matrix[ alias[k] ][ alias[imin] ] + matrix[ alias[k] ][ alias[jmin] ] - matrix[ alias[imin] ][ alias[jmin] ]) * 0.5;
        }
        for ( k++; k < jmin; k++) {
            matrix[alias[k]][alias[imin]] = matrix[alias[imin]][alias[k]] = (matrix[ alias[k] ][ alias[imin] ] + matrix[ alias[k] ][ alias[jmin] ] - matrix[ alias[imin] ][ alias[jmin] ]) * 0.5;
        }
        for ( k++; k < ncluster; k++) {
            matrix[alias[k]][alias[imin]] = matrix[alias[imin]][alias[k]] = (matrix[ alias[k] ][ alias[imin] ] + matrix[ alias[k] ][ alias[jmin] ] - matrix[ alias[imin] ][ alias[jmin] ]) * 0.5;
        }
        
        if( ncluster-jmin-1 != 0 ){
            memmove(&alias[jmin], &alias[jmin+1], sizeof(int)*(ncluster-jmin-1));
        }
        ncluster--;
    }
    
    Node *node = new_Node(NULL, NULL, ncluster);
    
    Node *inode = nodes[ alias[0] ];
    Node *jnode = nodes[ alias[1] ];
    
    Node_addChild(node, inode);
    Node_addChild(node, jnode);
    
    Node_set_parent(inode, node);
    Node_set_parent(jnode, node);
    
    float l = fmaxf(0.0, matrix[ alias[0] ][ alias[1] ]*0.5);
    
    Node_set_distance(inode, l);
    Node_set_distance(jnode, l);
    
    free(nodes);
    free(alias);
    free(r);
    
    return new_Tree2(node, false);
}


