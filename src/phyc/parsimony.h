/*
 *  parsimony.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 9/04/2014.
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

void Parsimony_update_node( Parsimony *parsimony, Node *node );

void Parsimony_update_all_nodes( Parsimony *parsimony );

#endif
