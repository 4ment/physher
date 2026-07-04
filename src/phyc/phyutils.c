// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "phyutils.h"


void PhyUtils_sort_sequences( Tree *tree, Sequences *seqs, treeorder to ){
    int *order = ivector(Tree_tip_count(tree));
    Node **nodes = Tree_get_nodes(tree, to);
    
    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if( Node_isleaf( nodes[i] ) ){
            int index = Sequences_get_index(seqs, nodes[i]->name);
			order[index] = i;
		}
	}
    
    Sequences_sort_from_ivector(seqs, order, Tree_tip_count(tree));
    free(order);
}

