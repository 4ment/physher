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
    Parameters* ratios;
    Parameter* rootHeight;
    double* lowers;
    size_t tipCount;
	tree_transform_t parameterization;
    double (*inverse_transform)(struct TreeTransform*, Node*);
    void (*update)(struct TreeTransform*);
    void (*update_lowers)(struct TreeTransform*);
	double (*log_jacobian)(struct TreeTransform*);
    void (*jvp)(struct TreeTransform*, const double*, double*);
    double (*dlog_jacobian)(struct TreeTransform*, Node*);
    void (*log_jacobian_gradient)(struct TreeTransform*, double*);
} TreeTransform;

TreeTransform* new_HeightTreeTransform(Tree* tree, tree_transform_t parameterization);

void free_TreeTransform(TreeTransform* tt);

Model* new_TreeTransformModel(const char* name, TreeTransform* coalescent, Model* tree);

void TreeTransform_jvp_with_heights(TreeTransform *tt, const double* heights, const double *height_gradient, double *gradient);

#endif /* treetransform_h */
