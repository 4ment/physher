/*
 *  rf.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 14/6/12.
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
#include <string.h>
#include <math.h>

#include "tree.h"
#include "matrix.h"
#include "hashtable.h"

#include "splitsystem.h"

// return Robinson Foulds distance
double rf_distance( bool **s1, bool **s2, int splitCount, int length ){	
	//int size = Tree_node_count(tree1) - Tree_tip_count(tree1) - 1;
	
	// number of splits in t1 missing in t2
	int fn = 0;
	for ( int i = 0; i < splitCount; i++ ){
		if ( !hasSplit(s2, s1[i], splitCount, length) ) fn++;
	}
	
	// number of splits in t2 missing in t1
	int fp = 0;
	for ( int i = 0; i < splitCount; i++ ) {
		if ( !hasSplit(s1, s2[i], splitCount, length) ) fp++;
	}
	
	return 0.5*((double) (fp + fn));
}

// normalized distances
double rf_norm_distance( bool **s1, bool **s2, int splitCount, int length ){
	return rf_distance(s1, s2, splitCount, length)/splitCount;
}


// Calculate Branch score
// Khuner 1994
double Branch_score( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length ){
    
    Node ** nodes1 = Tree_nodes(t1);
    Node ** nodes2 = Tree_nodes(t2);
    
    bool *done = bvector(splitCount);
    double *d1 = dvector(splitCount);
    double *d2 = dvector(splitCount);
    
    //    int j = 0;
    //    int jj = 0;
    //    for (int i = 0; i < Tree_node_count(t1); i++) {
    //        if( !Node_isroot(nodes1[i]) ){
    //            //if( !Node_isleaf(nodes1[i]) ) {
    //                d1[j] = Node_distance(nodes1[i]);
    //                j++;
    //            //}
    //        }
    //        if( !Node_isroot(nodes2[i]) ){
    //            //if( !Node_isleaf(nodes2[i]) ) {
    //                d2[jj] = Node_distance(nodes2[i]);
    //                jj++;
    //            //}
    //        }
    //    }
    double l1 = Node_distance(Node_left(Tree_root(t1)));
    double r1 = Node_distance(Node_right(Tree_root(t1)));
    double l2 = Node_distance(Node_left(Tree_root(t2)));
    double r2 = Node_distance(Node_right(Tree_root(t2)));
    
    // Always ignore the right branch unless it is a leaf
    Node *n = Node_right(Tree_root(t1));
    // Right is a leaf so we skip left node
    if( Node_isleaf(n) ){
        Node_set_distance(n, Node_distance(n)+Node_distance(Node_left(Tree_root(t1))));
        Node_set_distance(Node_left(Tree_root(t1)), 0);
        n = Node_left(Tree_root(t1));
    }
    else {
        Node_set_distance(Node_left(Tree_root(t1)), Node_distance(Node_left(Tree_root(t1)))+Node_distance(Node_right(Tree_root(t1))));
        Node_set_distance(Node_right(Tree_root(t1)), 0);
    }
    
    int j = 0;
    for (int i = 0; i < Tree_node_count(t1); i++) {
        if( !Node_isroot(nodes1[i]) && nodes1[i] != n && !Node_isleaf(nodes1[i]) ){
            d1[j++] = Node_distance(nodes1[i]);
        }
        
    }
    
    n = Node_right(Tree_root(t2));
    if( Node_isleaf(n) ){
        Node_set_distance(n, Node_distance(n)+Node_distance(Node_left(Tree_root(t2))));
        Node_set_distance(Node_left(Tree_root(t2)), 0);
        n = Node_left(Tree_root(t2));
    }
    else {
        Node_set_distance(Node_left(Tree_root(t2)), Node_distance(Node_left(Tree_root(t2)))+Node_distance(Node_right(Tree_root(t2))));
        Node_set_distance(Node_right(Tree_root(t2)), 0);
    }
    
    j = 0;
    for (int i = 0; i < Tree_node_count(t2); i++) {
        if( !Node_isroot(nodes2[i]) && nodes2[i] != n && !Node_isleaf(nodes2[i]) ){
            d2[j++] = Node_distance(nodes2[i]);
        }
        
    }
    
    double sum  = 0.0;
    
    // s1
    for ( int i = 0; i < splitCount; i++ ){
        int index = hasSplit2(s2, s1[i], splitCount, length);
        if ( index >= 0 ){
            //sum += pow( Node_distance(nodes1[i]) - Node_distance(nodes2[i]), 2);
            sum += pow( d1[i] - d2[index], 2);
            done[index] = true;
        }
        else{
            //sum += Node_distance(nodes1[i])*Node_distance(nodes1[i]);
            sum += d1[i] * d1[i];
            done[index] = false;
        }
    }
    
    // s2
    for ( int i = 0; i < splitCount; i++ ){
        if( !done[i] ){
            //sum += Node_distance(nodes2[i])*Node_distance(nodes2[i]);
            sum += d2[i] * d2[i];
        }
    }
    
    // branches leading to tips
    for (int i = 0; i < Tree_node_count(t1); i++) {
        if( !Node_isleaf(nodes1[i]) ) continue;
        
        for (int j = 0; j < Tree_node_count(t2); j++) {
            if( Node_isleaf(nodes2[i]) && strcmp(Node_name(nodes1[i]), Node_name(nodes2[i])) == 0 ){
                sum += pow( Node_distance(nodes1[i]) - Node_distance(nodes2[i]), 2);
                break;
            }
        }
    }
    
    //	for ( int i = 0; i < splitCount; i++ ){
    //
    //		if ( !hasSplit(s2, s1[i], splitCount, length) ){
    //            sum += Node_distance(nodes1[i])*Node_distance(nodes1[i]);
    //        }
    //        else {
    //            sum += pow( Node_distance(nodes1[i]) - Node_distance(nodes2[i]), 2);
    //        }
    //
    //
    //		if ( !hasSplit(s1, s2[i], splitCount, length) ){
    //            sum += Node_distance(nodes2[i])*Node_distance(nodes2[i]);
    //        }
    //        else {
    //            sum += pow( Node_distance(nodes1[i]) - Node_distance(nodes2[i]), 2);
    //        }
    //	}
    
    Node_set_distance(Node_left(Tree_root(t1)), l1);
    Node_set_distance(Node_right(Tree_root(t1)), r1);
    Node_set_distance(Node_left(Tree_root(t2)), l2);
    Node_set_distance(Node_right(Tree_root(t2)), r2);
    
    free(d1);
    free(d2);
    free(done);
    return sqrt(sum);
}

double K_tree_score( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length ){
    
    Node ** nodes1 = Tree_nodes(t1);
    Node ** nodes2 = Tree_nodes(t2);
    
    bool *done = bvector(splitCount);
    double *d1 = dvector(splitCount);
    double *d2 = dvector(splitCount);
    
    double l1 = Node_distance(Node_left(Tree_root(t1)));
    double r1 = Node_distance(Node_right(Tree_root(t1)));
    double l2 = Node_distance(Node_left(Tree_root(t2)));
    double r2 = Node_distance(Node_right(Tree_root(t2)));
    
    // Always ignore the right branch unless it is a leaf
    Node *n = Node_right(Tree_root(t1));
    // Right is a leaf so we skip left node
    if( Node_isleaf(n) ){
        Node_set_distance(n, Node_distance(n)+Node_distance(Node_left(Tree_root(t1))));
        Node_set_distance(Node_left(Tree_root(t1)), 0);
        n = Node_left(Tree_root(t1));
    }
    else {
        Node_set_distance(Node_left(Tree_root(t1)), Node_distance(Node_left(Tree_root(t1)))+Node_distance(Node_right(Tree_root(t1))));
        Node_set_distance(Node_right(Tree_root(t1)), 0);
    }
    
    int j = 0;
    for (int i = 0; i < Tree_node_count(t1); i++) {
        if( !Node_isroot(nodes1[i]) && nodes1[i] != n && !Node_isleaf(nodes1[i]) ){
            d1[j++] = Node_distance(nodes1[i]);
        }
        
    }
    
    n = Node_right(Tree_root(t2));
    if( Node_isleaf(n) ){
        Node_set_distance(n, Node_distance(n)+Node_distance(Node_left(Tree_root(t2))));
        Node_set_distance(Node_left(Tree_root(t2)), 0);
        n = Node_left(Tree_root(t2));
    }
    else {
        Node_set_distance(Node_left(Tree_root(t2)), Node_distance(Node_left(Tree_root(t2)))+Node_distance(Node_right(Tree_root(t2))));
        Node_set_distance(Node_right(Tree_root(t2)), 0);
    }
    
    j = 0;
    for (int i = 0; i < Tree_node_count(t2); i++) {
        if( !Node_isroot(nodes2[i]) && nodes2[i] != n && !Node_isleaf(nodes2[i]) ){
            d2[j++] = Node_distance(nodes2[i]);
        }
    }
    
    // Calculate K
    double K = 0;
    double sumd2 = 0;
    
    for ( int i = 0; i < splitCount; i++ ){
        int index = hasSplit2(s2, s1[i], splitCount, length);
        if ( index >= 0 ){
            K += d1[i]*d2[index];
        }
    }
    
    // branches leading to tips
    for (int i = 0; i < Tree_node_count(t1); i++) {
        if( !Node_isleaf(nodes1[i]) ) continue;
        
        for (int j = 0; j < Tree_node_count(t2); j++) {
            if( Node_isleaf(nodes2[i]) && strcmp(Node_name(nodes1[i]), Node_name(nodes2[i])) == 0 ){
                K += Node_distance(nodes1[i]) * Node_distance(nodes2[i]);
                break;
            }
        }
    }
    
    for (int i = 0; i < Tree_node_count(t2); i++) {
        if( !Node_isroot(nodes2[i]) ) sumd2 += Node_distance(nodes2[i])*Node_distance(nodes2[i]);
    }
    
    K /= sumd2;
    
    // Calculate K tree score
    double sum  = 0.0;
    
    //K = 1;
    
    // s1
    for ( int i = 0; i < splitCount; i++ ){
        int index = hasSplit2(s2, s1[i], splitCount, length);
        // partition present in both trees
        if ( index >= 0 ){
            sum += pow( d1[i] - K*d2[index], 2);
            done[index] = true;
        }
        // partition from tree 1 but not in tree 2
        else{
            sum += pow(d1[i],2);
            done[index] = false;
        }
    }
    
    // partition from tree 2 but not in tree 1
    for ( int i = 0; i < splitCount; i++ ){
        if( !done[i] ){
            sum += pow(K*d2[i],2);
        }
    }
    
    // branches leading to tips
    for (int i = 0; i < Tree_node_count(t1); i++) {
        if( !Node_isleaf(nodes1[i]) ) continue;
        
        for (int j = 0; j < Tree_node_count(t2); j++) {
            if( Node_isleaf(nodes2[i]) && strcmp(Node_name(nodes1[i]), Node_name(nodes2[i])) == 0 ){
                sum += pow( Node_distance(nodes1[i]) - K*Node_distance(nodes2[i]), 2);
                break;
            }
        }
    }
    
    Node_set_distance(Node_left(Tree_root(t1)), l1);
    Node_set_distance(Node_right(Tree_root(t1)), r1);
    Node_set_distance(Node_left(Tree_root(t2)), l2);
    Node_set_distance(Node_right(Tree_root(t2)), r2);
    
    free(d1);
    free(d2);
    free(done);
    return sqrt(sum);
}

// Caclculate K tree score
// Soria-Carrasco 2007: The K tree score: quantification of differences in the relative branch length and topology of phylogenetic trees
// d is not a distance d(t1,t2) != d(t2,t1)
double K_tree_score4( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length ){
    
    Node ** nodes1 = Tree_nodes(t1);
    Node ** nodes2 = Tree_nodes(t2);
    
    bool *done = bvector(splitCount);
    
    // Calculate K
    double num    = 0.0;
    double denom  = 0.0;
    
    // s1
    for ( int i = 0; i < splitCount; i++ ){
        if ( hasSplit(s2, s1[i], splitCount, length) ){
            num += pow( Node_distance(nodes1[i]) - Node_distance(nodes2[i]), 2);
        }
        denom += Node_distance(nodes2[i])*Node_distance(nodes2[i]);
    }
    
    double K = num/denom;
    
    // Calculate K tree score
    double sum  = 0.0;
    
    // s1
    for ( int i = 0; i < splitCount; i++ ){
        int index = hasSplit2(s2, s1[i], splitCount, length);
        if ( index >= 0 ){
            sum += pow( Node_distance(nodes1[i]) - K*Node_distance(nodes2[index]), 2);
            done[index] = true;
        }
        else{
            sum += Node_distance(nodes1[i])*Node_distance(nodes1[i]);
            done[index] = false;
        }
    }
    
    // s2
    for ( int i = 0; i < splitCount; i++ ){
        if( !done[i] ){
            sum += pow(K*Node_distance(nodes2[i]), 2);
        }
    }
    
    free(done);
    return sqrt(sum);
}

double K_tree_score_crap( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length ){
    
    Node ** nodes1 = Tree_nodes(t1);
    Node ** nodes2 = Tree_nodes(t2);
    
    // Calculate t2 to scale K
    double num    = 0.0;
    double denom  = 0.0;
    for ( int i = 0; i < splitCount; i++ ){
        
        if ( hasSplit(s2, s1[i], splitCount, length) ){
            num += Node_distance(nodes1[i]) * Node_distance(nodes2[i]);
        }
        else if ( hasSplit(s1, s2[i], splitCount, length) ){
            num += Node_distance(nodes1[i]) * Node_distance(nodes2[i]);
        }
        
        denom += Node_distance(nodes2[i])*Node_distance(nodes2[i]);
    }
    
    double K = num/denom;
    
    // Calculate K score
    double sum  = 0.0;
    for ( int i = 0; i < splitCount; i++ ){
        
        if ( !hasSplit(s2, s1[i], splitCount, length) ){
            sum += Node_distance(nodes1[i])*Node_distance(nodes1[i]);
        }
        else {
            sum += pow( Node_distance(nodes1[i]) - K*Node_distance(nodes2[i]), 2);
        }
        
        
        if ( !hasSplit(s1, s2[i], splitCount, length) ){
            sum += Node_distance(nodes2[i])*Node_distance(nodes2[i]);
        }
        else {
            sum += pow( Node_distance(nodes1[i]) - Node_distance(nodes2[i]), 2);
        }
    }
    
    return sqrt(sum);
}



void test( Tree *tree, Tree *tree2 ){
    Hashtable *hash = new_Hashtable_string(Tree_tip_count(tree));
    int count = 0;
    Node **nodes = Tree_nodes(tree);
    
    int nNodes = Tree_node_count(tree);
	int splitCount = nNodes - Tree_tip_count(tree) - 1;
    
    for ( int i = 0; i < nNodes; i++ ) {
        if ( Node_isleaf(nodes[i]) ) {
            Hashtable_add(hash, String_clone(Node_name(nodes[i])), new_Int(count++));
        }
    }
    bool **split = getSplitSystem(hash, tree);
    bool **split2 = getSplitSystem(hash, tree2);
    
    double rf = rf_distance(split, split2, splitCount, nNodes);
    double rf_norm = rf_norm_distance(split, split2, splitCount, nNodes);
    
    printf("rf %f\n", rf);
    printf("rf normalized %f\n", rf_norm);
    
    free_Hashtable(hash);
}
