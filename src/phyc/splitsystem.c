/*
 *  splitsystem.c
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

#include "splitsystem.h"

#include <string.h>
#include <assert.h>

#include "matrix.h"

void markNode( Hashtable *hash, Node *node, bool *split) {
    if ( Node_isleaf(node) ) {
        int *index = Hashtable_get(hash, Node_name(node) );
        assert(index);
        split[*index] = true;
    }
    else {
        markNode(hash, node->left, split);
        markNode(hash, node->right, split);
    }
}

void getSplit(Hashtable *hash, Node *internalNode, bool *split, int length ) {
    
    // make sure split is reset
    for (int i = 0; i < length; i++) {
        split[i] = false;
    }
    
    // mark all leafs downstream of the node
    if( !Node_isleaf(internalNode) ){
        markNode(hash, internalNode->left, split);
        markNode(hash, internalNode->right, split);
    }
    else {
        int *index = Hashtable_get(hash, Node_name(internalNode) );
        split[*index] = true;
    }
    
    // standardize split (i.e. first index is alway true)
    if ( !split[0] ) {
        for (int i = 0; i < length; i++) {
            split[i] = !split[i];
        }
    }
}

// hash associate an index (value) to each taxon (key)
bool ** getSplitSystem( Hashtable *hash, Tree *tree ){
    int nNodes = Tree_node_count(tree);
    int splitCount = nNodes - Tree_tip_count(tree) - 1;
    bool **splits = bmatrix(splitCount, nNodes);
    
    Node **nodes = Tree_nodes(tree);
    int j = 0;
    for (int i = 0; i < nNodes; i++) {
        if( !Node_isroot(nodes[i]) ){
            if( !Node_isleaf(nodes[i])) {
                getSplit(hash, nodes[i], splits[j], Tree_tip_count(tree));
                j++;
            }
        }
        
    }
    return splits;
}

bool ** getSplitSystemUnrooted( Hashtable *hash, Tree *tree ){
    int nNodes = Tree_node_count(tree);
    int splitCount = nNodes - Tree_tip_count(tree) - 2;
    bool **splits = bmatrix(splitCount, nNodes);
    
    // Always ignore the right branch unless it is a leaf
    Node *n = Node_right(Tree_root(tree));
    // Right is a leaf so we skip left node
    if( Node_isleaf(Node_right(Tree_root(tree))) ){
        n = Node_left(Tree_root(tree));
    }
    
    Node **nodes = Tree_nodes(tree);
    int j = 0;
    for (int i = 0; i < nNodes; i++) {
        if( !Node_isroot(nodes[i]) && nodes[i] != n ){
            if( !Node_isleaf(nodes[i])) {
                getSplit(hash, nodes[i], splits[j], Tree_tip_count(tree));
                j++;
            }
        }
        
    }
    return splits;
}

// Internal and external
bool ** getSplitSystemAll( Hashtable *hash, Tree *tree ){
    int nNodes = Tree_node_count(tree);
    int splitCount = nNodes - 1;
    bool **splits = bmatrix(splitCount, nNodes);
    
    Node **nodes = Tree_nodes(tree);
    int j = 0;
    for (int i = 0; i < nNodes; i++) {
        if( !Node_isroot(nodes[i]) ){
            getSplit(hash, nodes[i], splits[j], Tree_tip_count(tree));
            j++;
        }
        
    }
    return splits;
}


bool isSame( bool *s1, bool *s2, int length ) {
    bool reverse = s1[0] != s2[0];
    
    if ( reverse  ) {
        for (int i = 0; i < length; i++) {
            // splits not identical
            if (s1[i] == s2[i]) return false;
        }
    }
    else {
        return memcmp(s1, s2, length * sizeof(bool)) == 0;
    }
    
    return true;
}

bool hasSplit( bool **splits, bool *split, int splitCount, int length ) {
    for ( int i = 0; i < splitCount; i++ ) {
        if ( isSame(split, splits[i], length)) return true;
    }
    
    return false;
}

int hasSplit2( bool **splits, bool *split, int splitCount, int length ) {
    for ( int i = 0; i < splitCount; i++ ) {
        if ( isSame(split, splits[i], length)) return i;
    }
    
    return -1;
}
