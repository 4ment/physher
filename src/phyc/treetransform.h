//
//  treetransform.h
//  physher
//
//  Created by mathieu on 11/7/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef treetransform_h
#define treetransform_h

#include <stdio.h>

#include "parameters.h"
#include "tree.h"

typedef struct TreeTransform {
    Tree* tree;
    Parameters* parameters;
    double* lowers;
    unsigned* map;
    Node** map_to_node;
    Parameter* (*parameter_of_node)(struct TreeTransform*, Node*);
    Node* (*node_of_parameter)(struct TreeTransform*, Parameter*);
    double (*inverse_transform)(struct TreeTransform*, Node*);
    void (*update)(struct TreeTransform*);
    void (*update_lowers)(struct TreeTransform*);
    void (*jvp)(struct TreeTransform*, const double*, double*);
    double (*dlog_jacobian)(struct TreeTransform*, Node*);
    void (*log_jacobian_gradient)(struct TreeTransform*, double*);
} TreeTransform;

TreeTransform* new_HeightTreeTransform(Tree* tree);

void free_TreeTransform(TreeTransform* tt);

Model* new_TreeTransformModel(const char* name, TreeTransform* coalescent, Model* tree);

#endif /* treetransform_h */
