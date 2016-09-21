/*
 *  treesearch.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 7/08/13.
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

#ifndef PhyC_treesearch_h
#define PhyC_treesearch_h

#include "tree.h"

typedef enum {
    TREE_SEARCH_NONE = -1,
    TREE_SEARCH_NNI,
	TREE_SEARCH_NNNI,
    TREE_SEARCH_SPR,
    
	TREE_SEARCH_PARSIMONY_SPR,
	TREE_SEARCH_PARSIMONY_NNI
}tree_search_algorithm;





bool SPR_move( Tree *tree, Node *prune, Node *graft );

void NNI_move( Tree *tree, Node *node1, Node *node2 );

#endif
