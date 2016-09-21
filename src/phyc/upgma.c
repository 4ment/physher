/*
 *  upgma.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 20/01/2015.
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

#include "upgma.h"

#include <string.h>
#include <assert.h>

#include "matrix.h"

static void _findMinIndexes( double **matrix, int ncluster, int *alias, int *imin, int *jmin);

Tree * new_UPGMA( const Sequences *sequences, double **_matrix ){
    int nTips = sequences->size;
    
    Node **nodes = (Node**)malloc(nTips*sizeof(Node*));
    assert(nodes);
    int *alias = ivector(nTips);
    double *h = dvector(nTips);
    int *counts = ivector(nTips);
    
    int imin=0;
    int jmin=0;
    
    int i = 0;
    for ( ; i < nTips; i++ ) {
        nodes[i] = new_Node(NULL, sequences->seqs[i]->name, i);
        alias[i] = i;
        counts[i] = 1;
    }
    int ncluster = nTips;
    
    while( ncluster > 2 ){
        
        _findMinIndexes(_matrix, ncluster, alias, &imin, &jmin);
        
        Node *node = new_Node(NULL, NULL, ncluster);
        
        Node *inode = nodes[alias[imin]];
        Node *jnode = nodes[alias[jmin]];
        
        Node_addChild(node, inode);
        Node_addChild(node, jnode);
        
        Node_set_parent(inode, node);
        Node_set_parent(jnode, node);
        
        double l = dmax(0.0, _matrix[ alias[imin] ][ alias[jmin] ]*0.5);
        
        Node_set_distance(inode, l-h[alias[imin]]);
        Node_set_distance(jnode, l-h[alias[jmin]]);
        
        nodes[alias[imin]] = node;
        counts[alias[imin]] += counts[alias[jmin]];
        h[alias[imin]] = l;
        nodes[alias[jmin]] = NULL;
        
        int k = 0;
        for ( ; k < imin; k++) {
            _matrix[alias[k]][alias[imin]] = _matrix[alias[imin]][alias[k]] = (counts[alias[imin]]*_matrix[ alias[k] ][ alias[imin] ] + counts[alias[jmin]]*_matrix[ alias[k] ][ alias[jmin] ]) / (counts[alias[imin]]+counts[alias[jmin]]);
        }
        for ( k++; k < jmin; k++) {
            _matrix[alias[k]][alias[imin]] = _matrix[alias[imin]][alias[k]] = (counts[alias[imin]]*_matrix[ alias[k] ][ alias[imin] ] + counts[alias[jmin]]*_matrix[ alias[k] ][ alias[jmin] ]) / (counts[alias[imin]]+counts[alias[jmin]]);
        }
        for ( k++; k < ncluster; k++) {
            _matrix[alias[k]][alias[imin]] = _matrix[alias[imin]][alias[k]] = (counts[alias[imin]]*_matrix[ alias[k] ][ alias[imin] ] + counts[alias[jmin]]*_matrix[ alias[k] ][ alias[jmin] ]) / (counts[alias[imin]]+counts[alias[jmin]]);
        }
        if(ncluster-jmin-1 != 0 ){
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
    
    double l = dmax(0.0, _matrix[ alias[0] ][ alias[1] ]*0.5);
    
    Node_set_distance(inode, l-h[ alias[0] ]);
    Node_set_distance(jnode, l-h[ alias[1] ]);
    
    free(counts);
    free(nodes);
    free(alias);
    free(h);
    
    return new_Tree2(node, false);
}

void _findMinIndexes( double **matrix, int ncluster, int *alias, int *imin, int *jmin){
    double min = INFINITY;
    *imin = 0;
    *jmin = 0;
    
    for( int i = 0; i < ncluster; i++ ){
        for( int j = i+1; j < ncluster; j++ ){
            double sij = matrix[ alias[i] ][ alias[j] ];
            
            if( sij < min ){
                *imin = i;
                *jmin = j;
                min  = sij;
            }
        }
    }
}

void _findMinIndexes_float( float **matrix, int ncluster, int *alias, int *imin, int *jmin){
    float min = INFINITY;
    *imin = 0;
    *jmin = 0;
    
    for( int i = 0; i < ncluster; i++ ){
        for( int j = i+1; j < ncluster; j++ ){
            float sij = matrix[ alias[i] ][ alias[j] ];
            
            if( sij < min ){
                *imin = i;
                *jmin = j;
                min  = sij;
            }
        }
    }
}

Tree * new_UPGMA_float( const Sequences *sequences, float **_matrix ){
    int nTips = sequences->size;
    
    Node **nodes = (Node**)malloc(nTips*sizeof(Node*));
    assert(nodes);
    int *alias = ivector(nTips);
    float *h = fvector(nTips);
    int *counts = ivector(nTips);
    
    int imin=0;
    int jmin=0;
    
    int i = 0;
    for ( ; i < nTips; i++ ) {
        nodes[i] = new_Node(NULL, sequences->seqs[i]->name, i);
        alias[i] = i;
        counts[i] = 1;
    }
    int ncluster = nTips;
    
    while( ncluster > 2 ){
        
        _findMinIndexes_float(_matrix, ncluster, alias, &imin, &jmin);
        
        Node *node = new_Node(NULL, NULL, ncluster);
        
        Node *inode = nodes[alias[imin]];
        Node *jnode = nodes[alias[jmin]];
        
        Node_addChild(node, inode);
        Node_addChild(node, jnode);
        
        Node_set_parent(inode, node);
        Node_set_parent(jnode, node);
        
        float l = fmaxf(0.0, _matrix[ alias[imin] ][ alias[jmin] ]*0.5);
        
        Node_set_distance(inode, l-h[alias[imin]]);
        Node_set_distance(jnode, l-h[alias[jmin]]);
        
        nodes[alias[imin]] = node;
        counts[alias[imin]] += counts[alias[jmin]];
        h[alias[imin]] = l;
        nodes[alias[jmin]] = NULL;
        
        int k = 0;
        for ( ; k < imin; k++) {
            _matrix[alias[k]][alias[imin]] = _matrix[alias[imin]][alias[k]] = (counts[alias[imin]]*_matrix[ alias[k] ][ alias[imin] ] + counts[alias[jmin]]*_matrix[ alias[k] ][ alias[jmin] ]) / (counts[alias[imin]]+counts[alias[jmin]]);
        }
        for ( k++; k < jmin; k++) {
            _matrix[alias[k]][alias[imin]] = _matrix[alias[imin]][alias[k]] = (counts[alias[imin]]*_matrix[ alias[k] ][ alias[imin] ] + counts[alias[jmin]]*_matrix[ alias[k] ][ alias[jmin] ]) / (counts[alias[imin]]+counts[alias[jmin]]);
        }
        for ( k++; k < ncluster; k++) {
            _matrix[alias[k]][alias[imin]] = _matrix[alias[imin]][alias[k]] = (counts[alias[imin]]*_matrix[ alias[k] ][ alias[imin] ] + counts[alias[jmin]]*_matrix[ alias[k] ][ alias[jmin] ]) / (counts[alias[imin]]+counts[alias[jmin]]);
        }
        if(ncluster-jmin-1 != 0 ){
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
    
    float l = fmaxf(0.0, _matrix[ alias[0] ][ alias[1] ]*0.5);
    
    Node_set_distance(inode, l-h[ alias[0] ]);
    Node_set_distance(jnode, l-h[ alias[1] ]);
    
    free(counts);
    free(nodes);
    free(alias);
    free(h);
    
    return new_Tree2(node, false);
}

