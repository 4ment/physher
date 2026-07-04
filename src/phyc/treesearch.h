// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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





Node* SPR_move( Tree *tree, Node *prune, Node *graft );

void NNI_move( Tree *tree, Node *node1, Node *node2 );

#endif
