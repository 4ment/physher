/*
 *  phyutils.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/06/13.
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

