// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_parsimony_h
#define PhyC_parsimony_h

#include "sitepattern.h"
#include "tree.h"

#define INT_WEIGHT 1

struct _Parsimony;
typedef struct _Parsimony Parsimony;

typedef double (*calculate_parsimony)(Parsimony*);

struct _Parsimony{
    SitePattern *sp;
    Tree *tree;
    
    int8_t **stateSets;
    uint8_t **states;
    calculate_parsimony calculate;
    void (*reconstruct)(Parsimony*);
    bool update;
    bool *update_nodes;
    double score;
    
#ifdef INT_WEIGHT
    int32_t *weights;
    int32_t **local_scores;
    int32_t *scores;
#else
    double **local_scores;
    double *scores;
#endif
};

Parsimony * new_Parsimony( SitePattern *sp, Tree *tree );

void free_Parsimony( Parsimony *parsimony );

Model * new_ParsimonyModel(char* name, Parsimony* parsimony, Model* tree);

Model * new_ParsimonyModel_from_json(json_node*node, Hashtable*hash);

void Parsimony_update_node( Parsimony *parsimony, Node *node );

void Parsimony_update_all_nodes( Parsimony *parsimony );

#endif
