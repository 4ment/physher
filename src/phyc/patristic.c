/*
 *  patristic.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/2/12.
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

#include "tree.h"
#include "matrix.h"
#include "hashtable.h"
#include "patristic.h"

double ** calculate_patristic( Tree *tree ){
	Node **nodes = Tree_get_nodes(tree, POSTORDER );
	Node *a = NULL;
	Node *b = NULL;
	double d = 0.0;
	
	int nTips = Tree_tip_count(tree);
	
	double **matrix = dmatrix( nTips, nTips );;
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if (! Node_isleaf(nodes[i]) ) continue;
		for ( int j = 0; j < Tree_node_count(tree); j++ ){
			if ( !Node_isleaf(nodes[j])) continue;
			a = nodes[i];
			b = nodes[j];
			d = 0.0;
			while ( a != b ) {
				if ( a->depth < b->depth ) {
					d += Node_distance(b);
					b = Node_parent(b);
				}
				else if( a->depth > b->depth ){
					d += Node_distance(a);
					a = Node_parent(a);
				}
				else {
					d += Node_distance(a);
					d += Node_distance(b);
					a = Node_parent(a);
					b = Node_parent(b);
				}
			}
			matrix[Node_class_id(nodes[i])][Node_class_id(nodes[j])] = matrix[Node_class_id(nodes[j])][Node_class_id(nodes[i])] = d;
		}
	}
	return matrix;
}


double ** calculate_patristic2( Tree **trees, int nTree, Hashtable *hash ){
	
	int nTips = Tree_tip_count(trees[0]);
	int nNodes = Tree_node_count(trees[0]);

	
	Node *a = NULL;
	Node *b = NULL;
	double d = 0.0;
	int *pi = NULL;
	int *pj = NULL;
	int len = (nTips * nTips - nTips)/2;
	
	double **matrix = dmatrix( nTree, len );
	for ( int i = 0; i < nTree; i++) {
		for ( int j = 0; j < len; j++) {
			matrix[i][j] = -1;
		}
	}
	
	//int *bb = ivector(len);
	
	Node **nodes = NULL;
	int pos = 0;
	for ( int index = 0; index <nTree; index++ ) {
		nodes = Tree_get_nodes( trees[index], POSTORDER );
		
		for ( int ii = 0; ii < nNodes; ii++ ) {
			if (! Node_isleaf(nodes[ii]) ) continue;
			for ( int jj = ii+1; jj < nNodes; jj++ ){
				if ( !Node_isleaf(nodes[jj])) continue;
				
				pi = Hashtable_get(hash, Node_name(nodes[ii]) );
				pj = Hashtable_get(hash, Node_name(nodes[jj]) );
				
				if ( pi == NULL || pj == NULL ) {
					fprintf(stderr, "i %s %s\n", Node_name(nodes[ii]), Node_name(nodes[jj]));
					print_Hashtable(stderr,hash);
                    exit(1);
				}
				//fprintf(stderr, "%s (%f) %s (%f) ",  Node_name(nodes[ii]), Node_distance(nodes[ii]), Node_name(nodes[jj]), Node_distance(nodes[jj]) );
				int i = *pi;
				int j = *pj;
								
				a = nodes[ii];
				b = nodes[jj];
				d = 0.0;
				while ( a != b ) {
					if ( a->depth < b->depth ) {
						d += Node_distance(b);
						b = Node_parent(b);
					}
					else if( a->depth > b->depth ){
						d += Node_distance(a);
						a = Node_parent(a);
					}
					else {
						d += Node_distance(a);
						d += Node_distance(b);
						a = Node_parent(a);
						b = Node_parent(b);
					}
				}
				//matrix[index][i*nTips + j] = d;
				//matrix[index][j*nTips + i] = d;
				if ( i > j ) {
					swap_int(&i, &j);
				}
				pos = 0;
				for ( int k = 0; k < i; k++) {
					pos += nTips-k-1;
				}
				pos += j-i-1;
				matrix[index][pos] = d;
				//bb[pos]++;
				assert(pos < len);
				//fprintf(stderr, " %f (%d)\n", d,pos );exit(0);
			}
		}
		//print_ivector(bb, len);
		//memset(bb, 0, len*sizeof(int));
	}
	
	//free(bb);
	return matrix;
}
