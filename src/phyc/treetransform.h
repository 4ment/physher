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

typedef enum tree_transform_t{
	TREE_TRANSFORM_RATIO_NAIVE=0,
	TREE_TRANSFORM_RATIO,
    TREE_TRANSFORM_SHIFT,
    TREE_TRANSFORM_PROPORTION
}tree_transform_t;

typedef struct TreeTransform {
    Tree* tree;
    Parameters* parameters;
    double* lowers;
    size_t tipCount;
	tree_transform_t parameterization;
    bool* unknownLeaves;
    // ratio slot for each unknown-age leaf, indexed by node id; slot == rank in
    // increasing node-id order. Gives a traversal-independent leaf->slot mapping
    // (node ids need not follow postorder, e.g. when assigned from a taxa list).
    size_t* unknownLeafSlot;
    double (*inverse_transform)(struct TreeTransform*, Node*);
    void (*update)(struct TreeTransform*);
    void (*update_lowers)(struct TreeTransform*);
	double (*log_jacobian)(struct TreeTransform*);
    void (*jvp)(struct TreeTransform*, const double*, double*);
    void (*jvp2)(struct TreeTransform*, const double*);
    double (*dlog_jacobian)(struct TreeTransform*, Node*);
    void (*log_jacobian_gradient)(struct TreeTransform*, double*);
} TreeTransform;

void TreeTransform_backward(TreeTransform* obj, Parameters* parameters, const double* ingrad);

TreeTransform* new_HeightTreeTransform(Tree* tree, tree_transform_t parameterization);

void free_TreeTransform(TreeTransform* tt);

void TreeTransform_initialize_from_heights(TreeTransform* tt);

void TreeTransform_add_tree_model(TreeTransform* tt, Tree* tree);

void TreeTransform_jvp_with_heights(TreeTransform *tt, const double* heights, const double *height_gradient, double *gradient);

void TreeTransformModel_add_tree_model(Model* self, Model* tree);

Model* new_TreeTransformModel(const char* name, TreeTransform* tt, Model* tree);

Model* new_TreeTransformModel_from_json(json_node* node, Hashtable* hash);

#endif /* treetransform_h */
