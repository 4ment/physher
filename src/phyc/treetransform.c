// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "treetransform.h"

#include "assert.h"
#include "strings.h"

#include "matrix.h"


// Assign each unknown-age leaf a stable ratio slot equal to its rank in
// increasing node-id order, stored by node id. Every forward/backward routine
// maps a leaf to its slot through this table, so the slot layout does not depend
// on traversal order (node ids do not necessarily follow postorder).
static void _build_unknown_leaf_slots(TreeTransform* tt){
    free(tt->unknownLeafSlot);
    tt->unknownLeafSlot = NULL;
    if(tt->unknownLeaves == NULL) return;
    tt->unknownLeafSlot = malloc(sizeof(size_t)*tt->tipCount);
    size_t slot = 0;
    for(size_t i = 0; i < tt->tipCount; i++){
        // leaves have node id == class_id in [0, tipCount); unused entries are
        // never read but initialized for cleanliness
        tt->unknownLeafSlot[i] = tt->unknownLeaves[i] ? slot++ : 0;
    }
}


#pragma region Shift transform
// ======================= Shift transform =======================

void tree_transform_shift_update_heights(Node *node, const double *parameters) {
	// set height quietly because ratios already notified the tree (and other upstream listeners)
	if (!Node_isleaf(node)) {
        tree_transform_shift_update_heights(node->left, parameters);
        tree_transform_shift_update_heights(node->right, parameters);
        double shift = parameters[Node_class_id(node)];
        double height = fmax(Node_height(node->left), Node_height(node->right));
        Node_set_height_quietly(node, height + shift);
    }
}

void _tree_transform_shift_update(TreeTransform *tt) {
    const double* values = Parameter_values(Parameters_at(tt->parameters, 0));
    tree_transform_shift_update_heights(Tree_root(tt->tree), values);
}

double _height_tree_inverse_shift_transform(TreeTransform *tt, Node *node) {
    return Node_height(node) - fmax(Node_height(node->left), Node_height(node->right));
}

void node_transform_vjp_shift(TreeTransform *tt, const double *height_gradient, double *gradient){
	memset(gradient, 0, sizeof(double)*(tt->tipCount-1));
    Node** nodes = Tree_get_nodes(tt->tree, PREORDER);
    size_t nodeCount = Tree_node_count(tt->tree);
    double* adjoints = clone_dvector(height_gradient, nodeCount);

    for(size_t i = 0; i < nodeCount; i++){
        Node* node = nodes[i];
        if(!Node_isleaf(node)){
            double grad = adjoints[Node_id(node)];
            gradient[Node_class_id(node)] = grad;
            double hl = Node_height(node->left);
            double hr = Node_height(node->right);

            if(hl > hr){
                adjoints[Node_id(node->left)] += grad;
            }
            else if(hl < hr){
                adjoints[Node_id(node->right)] += grad;
            }
            else{
                adjoints[Node_id(node->left)] += grad * 0.5;
                adjoints[Node_id(node->right)] += grad * 0.5;
            }
        }
    }
    free(adjoints);
}

double _node_transform_log_jacobian_zero(TreeTransform *tt){
	return 0.0;
}

double _node_transform_dlog_jacobian_zero(struct TreeTransform *tt, Node *node) {
    return 0.0;
}

void _node_transform_log_jacobian_gradient_zero(struct TreeTransform *tt, double *gradient) {
	
}

static void _shift_backward(TreeTransform* obj, Parameters* parameters, const double* ingrad){
	double* reparamGradient = dvector(obj->tipCount - 1);
	obj->vjp(obj, ingrad, reparamGradient);
    Parameter* shifts = Parameters_at(obj->parameters, 0);
    // accumulate the gradient of constrained parameters
    for(size_t i = 0; i < obj->tipCount - 1; i++){
        shifts->grad[i] += reparamGradient[i];
    }

    Parameter* shiftsx = Parameters_depends(parameters, shifts);

    // apply chain rule for unconstrained parameters
    if(shiftsx != shifts){
        shifts->transform->backward(shifts->transform, reparamGradient);
    }
	free(reparamGradient);
}
#pragma endregion Shift transform

#pragma region Ji gradient implementation
// ======================= Ji's gradient implementation =======================

// Efficient calculation of the ratio and root height gradient, adpated from BEAST.
// Xiang et al. Scalable Bayesian divergence time estimation with ratio transformation (2021).
// https://arxiv.org/abs/2110.13298
static double _epoch_gradient_addition(const Node* node, const Node* child, const double* lowers, const double* ratios, const double* ratio_gradient) {
    if (Node_isleaf(child)) {
        return 0.0;
    }
    size_t node_id = Node_id(node);
    size_t child_id = Node_id(child);
    size_t child_class_id = Node_class_id(child);
	double node_ratio = ratios[Node_class_id(node)];
	double child_ratio = ratios[child_class_id];
    
    // child_id and node_id are in the same epoch
    if (lowers[node_id] == lowers[child_id]) {
        return ratio_gradient[child_class_id] * child_ratio / node_ratio;
    }
    // NOT the same epoch
    else {
        double height = Node_height(node);
        return ratio_gradient[child_class_id] * child_ratio / (height - lowers[child_id]) * (height - lowers[node_id])/node_ratio;
    }
}

static double _root_height_gradient(Tree* tree, const double* ratios, const double* gradient_height) {
    size_t tipCount = Tree_tip_count(tree);
	size_t gradientLength = tipCount - 1;
    double* multiplierArray = dvector(gradientLength);
    multiplierArray[Node_class_id(Tree_root(tree))] = 1.0;
    Node** nodes = Tree_get_nodes(tree, PREORDER);
    size_t nodeCount = Tree_node_count(tree);
    
    for(size_t i = 1; i < nodeCount; i++){
		if(!Node_isleaf(nodes[i])){
			size_t node_class_id = Node_class_id(nodes[i]);
			multiplierArray[node_class_id] = ratios[node_class_id] * multiplierArray[Node_class_id(nodes[i]->parent)];
		}
    }
    double sum = 0.0;
    for (int i = 0; i < gradientLength; i++) {
        sum += gradient_height[i + tipCount] * multiplierArray[i];
    }
    free(multiplierArray);
    return sum;
}

void _update_ratios_gradient(Tree* tree, const double* lowers, const double* ratios, const double* gradient_height, double* gradient){
    Node** nodes = Tree_get_nodes(tree, POSTORDER);
    for (int i = 0; i < Tree_node_count(tree)-1; i++) {
        if(Node_isleaf(nodes[i])) continue;
        int node_id = Node_id(nodes[i]);
        int node_class_id = Node_class_id(nodes[i]);
        double height = Node_height(nodes[i]);
        gradient[node_class_id] += (height - lowers[node_id])/ratios[node_class_id] * gradient_height[node_id];
        gradient[node_class_id] += _epoch_gradient_addition(nodes[i], nodes[i]->left, lowers, ratios, gradient);
        gradient[node_class_id] += _epoch_gradient_addition(nodes[i], nodes[i]->right, lowers, ratios, gradient);
    }
}


void node_transform_vjp_efficient(TreeTransform *tt, const double *height_gradient, double *gradient){
	memset(gradient, 0.0, sizeof(double)*(tt->tipCount-1));
    const double* ratios = Parameter_values(Parameters_at(tt->parameters, 0));
	_update_ratios_gradient(tt->tree, tt->lowers, ratios, height_gradient, gradient);
	gradient[Node_class_id(Tree_root(tt->tree))] = _root_height_gradient(tt->tree, ratios, height_gradient);
}

void _node_transform_log_jacobian_gradient_efficient(struct TreeTransform *tt, double *gradient) {
	size_t nodeCount = Tree_node_count(tt->tree);
	double *log_time = dvector(nodeCount);
	unsigned root_class_id = Node_class_id(Tree_root(tt->tree));
    const double* ratioValues = Parameter_values(Parameters_at(tt->parameters, 0));
	
	Node** nodes = Tree_nodes(tt->tree);
	for (size_t i = tt->tipCount; i < nodeCount; i++) {
		log_time[Node_id(nodes[i])] = 1.0 / (Node_height(nodes[i]) - tt->lowers[Node_id(nodes[i])]);
	}
	log_time[Node_id(Tree_root(tt->tree))] = 0.0;
	double* jac_gradient = dvector(tt->tipCount-1);
	memset(jac_gradient, 0.0, sizeof(double)*(tt->tipCount-1));
	_update_ratios_gradient(tt->tree, tt->lowers, ratioValues, log_time, jac_gradient);
	
    double rootHeightGrad = _root_height_gradient(tt->tree, ratioValues, log_time);

	if(gradient == NULL){
        Parameter* ratios = Parameters_at(tt->parameters, 0);
        Parameter* rootHeight = Parameters_at(tt->parameters, 1);
        rootHeight->grad[0] += rootHeightGrad;
        if(rootHeight->transform){
            rootHeight->transform->backward(rootHeight->transform, &rootHeightGrad);
        }
        double* gradient = dvector(tt->tipCount-1);
        for (size_t i = tt->tipCount; i < nodeCount-1; i++) {
            size_t node_class_id = Node_class_id(nodes[i]);
            gradient[node_class_id] += jac_gradient [node_class_id] - 1.0/ratioValues[node_class_id];
            ratios->grad[node_class_id] += gradient[node_class_id];
        }
        if(ratios->transform){
            ratios->transform->backward(ratios->transform, gradient);
        }
        free(gradient);
    }
	else{
        gradient[root_class_id] += rootHeightGrad;
        for (size_t i = tt->tipCount; i < nodeCount-1; i++) {
            size_t node_class_id = Node_class_id(nodes[i]);
            gradient[node_class_id] += jac_gradient [node_class_id] - 1.0/ratioValues[node_class_id];
        }
    }
	
	free(jac_gradient);
	free(log_time);
}

#pragma endregion Ji gradient implementation

#pragma region Proportion gradient transform
// ======================= Proportion gradient =======================

void node_transform_vjp_backprop(TreeTransform *tt, const double *height_gradient, double *gradient){
    size_t nodeCount = Tree_node_count(tt->tree);
    Parameter* proportions = Parameters_at(tt->parameters, 0);
    size_t propCount = Parameter_size(proportions);
    size_t offset = propCount - (tt->tipCount - 2);
	memset(gradient, 0.0, sizeof(double)*(propCount + 1));
    double *adjoints = clone_dvector(height_gradient, nodeCount);
    Node** nodes = Tree_get_nodes(tt->tree, POSTORDER);
    Node* root = Tree_root(tt->tree);

    for (size_t i = 0; i < nodeCount-1; i++) {
        Node* node = nodes[i];
        size_t nodeIndex = Node_id(node);
        size_t parentIndex = Node_id(node->parent);
        if (!Node_isleaf(node)) {
            size_t nodeClassIndex = Node_class_id(node);
	        gradient[offset + nodeClassIndex] = adjoints[nodeIndex] * (Node_height(node->parent) - tt->lowers[nodeIndex]);
            adjoints[parentIndex] += adjoints[nodeIndex] * Parameter_value_at(proportions, offset + nodeClassIndex);
        }
        else if(tt->unknownLeaves[nodeIndex]){
            size_t slot = tt->unknownLeafSlot[nodeIndex];
            gradient[slot] = adjoints[nodeIndex] * (Node_height(node->parent) - tt->lowers[nodeIndex]);
            adjoints[parentIndex] += adjoints[nodeIndex] * Parameter_value_at(proportions, slot);
        }
    }
    gradient[offset + Node_class_id(root)] = adjoints[Node_id(root)];
    free(adjoints);
}

void _node_transform_log_jacobian_gradient_backprop(struct TreeTransform *tt, double *gradient) {
	size_t nodeCount = Tree_node_count(tt->tree);
    size_t rootId = Node_id(Tree_root(tt->tree));
	double *adjoints = dvector(nodeCount);
    memset(adjoints, 0.0, sizeof(double)*nodeCount);
	Node** nodes = Tree_get_nodes(tt->tree, POSTORDER);
	for (size_t i = 0; i < nodeCount-1; i++) {
        Node* node = nodes[i];
        if (!Node_isleaf(node) || tt->unknownLeaves[Node_id(node)]) {
            size_t parentIndex = Node_id(nodes[i]->parent);
            adjoints[parentIndex] += 1.0 / (Node_height(node->parent) - tt->lowers[Node_id(node)]);
        }
	}

    Parameter* ratios = Parameters_at(tt->parameters, 0);
    size_t propCount = Parameter_size(ratios);
    size_t offset = propCount - (tt->tipCount - 2);

    if(gradient == NULL){
        Parameter* rootHeight = Parameters_at(tt->parameters, 1);
        double* gradient = dvector(propCount + 1);
        for (size_t i = 0; i < nodeCount-1; i++) {
            Node* node = nodes[i];
            size_t nodeIndex = Node_id(node);
            size_t nodeClassIndex = Node_class_id(node);
            size_t parentIndex = Node_id(node->parent);

            if (!Node_isleaf(node)) {
                size_t idx = offset + nodeClassIndex;
                gradient[idx] += adjoints[nodeIndex] * (Node_height(node->parent) - tt->lowers[nodeIndex]);
                ratios->grad[idx] += gradient[idx];
                adjoints[parentIndex] += adjoints[nodeIndex] * Parameter_value_at(ratios, idx);
            }
            else if(tt->unknownLeaves[nodeIndex]){
                size_t slot = tt->unknownLeafSlot[nodeIndex];
                gradient[slot] += adjoints[nodeIndex] * (Node_height(node->parent) - tt->lowers[nodeIndex]);
                ratios->grad[slot] += gradient[slot];
                adjoints[parentIndex] += adjoints[nodeIndex] * Parameter_value_at(ratios, slot);
            }
        }
        rootHeight->grad[0] += adjoints[rootId];
        if(rootHeight->transform){
            rootHeight->transform->backward(rootHeight->transform, adjoints + rootId);
        }
        if(ratios->transform){
            ratios->transform->backward(ratios->transform, gradient);
        }
        free(gradient);
    }
	else{
        for (size_t i = 0; i < nodeCount-1; i++) {
            Node* node = nodes[i];
            if (!Node_isleaf(node)) {
                size_t nodeIndex = Node_id(node);
                size_t nodeClassIndex = Node_class_id(node);
                size_t parentIndex = Node_id(node->parent);
                gradient[nodeClassIndex] += adjoints[nodeIndex] * (Node_height(node->parent) - tt->lowers[nodeIndex]);
                adjoints[parentIndex] += adjoints[nodeIndex] * Parameter_value_at(ratios, nodeClassIndex);
            }
        }
        gradient[Node_class_id(Tree_root(tt->tree))] += adjoints[rootId];
    }

	free(adjoints);
}

#pragma endregion Proportion gradient transform

#pragma region Proportion transform
// parameters: parameters we want to differentiate with respect to. They can be constrained or unconstrained.
// ingrad: gradient of a function wrt heights
static void _proportions_backward(TreeTransform* obj, Parameters* parameters, const double* ingrad){
	Parameter* ratios = Parameters_at(obj->parameters, 0);
    Parameter* root = Parameters_at(obj->parameters, 1);
    size_t propCount = Parameter_size(ratios);
    double* reparamGradient = dvector(propCount + 1);
	obj->vjp(obj, ingrad, reparamGradient);
    // accumulate the gradient of constrained parameters
    size_t i = 0;
    while(i < propCount){
        ratios->grad[i] += reparamGradient[i];
        i++;
    }
    root->grad[0] += reparamGradient[i];

    Parameter* ratiosx = Parameters_depends(parameters, ratios);
    Parameter* rootx = Parameters_depends(parameters, root);

    // apply chain rule for unconstrained parameters
    if(ratiosx != ratios){
        ratios->transform->backward(ratios->transform, reparamGradient);
    }
    
    // rootx is NULL if the root is fixed
    if(rootx != NULL && rootx != root){
        root->transform->backward(root->transform, reparamGradient + propCount);
    }

	free(reparamGradient);
}

double _node_transform_log_jacobian(TreeTransform *tt){
	double logP = 0.0;
	Node** nodes = Tree_nodes(tt->tree);
    size_t nodeCount1 = Tree_node_count(tt->tree)-1;
	for(size_t i = tt->tipCount; i < nodeCount1; i++){
		logP += log(Node_height(Node_parent(nodes[i])) - tt->lowers[Node_id(nodes[i])]);
	}
    for(size_t i = 0; i < tt->tipCount; i++){
        if(tt->unknownLeaves[Node_id(nodes[i])]){
		    logP += log(Node_height(Node_parent(nodes[i])) - tt->lowers[Node_id(nodes[i])]);
        }
	}
	return logP;
}

void tree_transform_update_heights(Node *node, const double *ratios, double rootHeight, double *lowers, size_t offset) {
	// set height quietly because ratios already notified the tree (and other upstream listeners)
	if (!Node_isleaf(node)) {
        if (Node_isroot(node)) {
            Node_set_height_quietly(node, rootHeight);
        } else {
            double lower = lowers[Node_id(node)];
            Node_set_height_quietly(node, lower + (Node_height(Node_parent(node)) - lower) * ratios[offset + Node_class_id(node)]);
        }
        tree_transform_update_heights(node->left, ratios, rootHeight, lowers, offset);
        tree_transform_update_heights(node->right, ratios, rootHeight, lowers, offset);
    }
}

// should be called when topology changes
void tree_transform_collect_lowers(Node *node, TreeTransform *tt, double *lowers) {
    if (!Node_isleaf(node)) {
        tree_transform_collect_lowers(node->left, tt, lowers);
        tree_transform_collect_lowers(node->right, tt, lowers);
        lowers[Node_id(node)] = fmax(lowers[Node_id(Node_left(node))], lowers[Node_id(Node_right(node))]);
        if (Node_isroot(node)) {
            Parameter* rootHeight = Parameters_at(tt->parameters, 1);
            Parameter_set_lower(rootHeight, lowers[Node_id(node)]);
            Constraint_set_flower(rootHeight->cnstr, lowers[Node_id(node)]);
            if(rootHeight->transform != NULL){
                rootHeight->transform->lower = lowers[Node_id(node)];
                // printf("%f\n", lowers[Node_id(node)]);
                
            }
            // printf("root_height %p %s %f %f\n", rootHeight, Parameter_name(rootHeight), Parameter_value(rootHeight), lowers[Node_id(node)]);
        }
    } else if(!tt->unknownLeaves[Node_id(node)]) {
        lowers[Node_id(node)] = Node_height(node);
    }
    else {
        lowers[Node_id(node)] = 0.0;
    }
}

void _tree_transform_update(TreeTransform *tt) {
    Parameter* ratios = Parameters_at(tt->parameters, 0);
    const double* ratioValues = Parameter_values(ratios);
    double rootHeight = Parameter_value(Parameters_at(tt->parameters, 1));
    size_t offset = Parameter_size(ratios) - (tt->tipCount - 2);
    tree_transform_update_heights(Tree_root(tt->tree), ratioValues, rootHeight, tt->lowers, offset);
    if(offset > 0){
        Node** nodes = Tree_nodes(tt->tree);
        for(size_t i = 0; i < tt->tipCount; i++){
            if(tt->unknownLeaves[nodes[i]->id]){
                size_t slot = tt->unknownLeafSlot[nodes[i]->id];
                Node_set_height_quietly(nodes[i], Node_height(Node_parent(nodes[i])) * ratioValues[slot]);
            }
        }
    }
}

void _tree_transform_update_lowers(TreeTransform *tt) {
    tree_transform_collect_lowers(Tree_root(tt->tree), tt, tt->lowers);
}

double _height_tree_inverse_transform(TreeTransform *tt, Node *node) {
    if (Node_isroot(node)) return Node_height(node);
    return (Node_height(node) - tt->lowers[Node_id(node)]) / (Node_height(Node_parent(node)) - tt->lowers[Node_id(node)]);
}

#pragma endregion Proportion transform

#pragma region Proportion gradient inefficient algorithm

void _node_transform_dlog_jacobian_aux(TreeTransform *tt, const Node *noderef, Node *node, const double* proportions, double *dlogP, double *descendant) {
    if (!Node_isleaf(node)) {
        if (!Node_isroot(node) && node != noderef) {
            descendant[Node_id(node)] = descendant[Node_id(node->parent)] * proportions[Node_class_id(node)];
        } else if (!Node_isroot(node)) {
            descendant[Node_id(node)] =
                Node_height(Node_parent(node)) - tt->lowers[Node_id(node)];
        } else {
            descendant[Node_id(node)] = 1;
        }
        _node_transform_dlog_jacobian_aux(tt, noderef, node->left, proportions, dlogP, descendant);
        _node_transform_dlog_jacobian_aux(tt, noderef, node->right, proportions, dlogP, descendant);

        if (!Node_isroot(node) && node != noderef) {
            *dlogP += descendant[Node_id(node->parent)] / (Node_height(Node_parent(node)) - tt->lowers[Node_id(node)]);
        }
    }
}

// Returns derivative of the log det of the Jacobian
double _node_transform_dlog_jacobian(struct TreeTransform *tt, Node *node) {
    size_t nodeCount = Tree_node_count(tt->tree);
    double *descendant = dvector(nodeCount);
    const double* proportions = Parameter_values(Parameters_at(tt->parameters, 0));
    double adj = 0;
    _node_transform_dlog_jacobian_aux(tt, node, node, proportions, &adj, descendant);
    free(descendant);
    return adj;
}

// Update gradient with gradient of the log det of the Jacobian
void _node_transform_log_jacobian_gradient(struct TreeTransform *tt, double *gradient) {
    size_t nodeCount = Tree_node_count(tt->tree);
    double *descendant = dvector(nodeCount);
    Parameter* proportions = Parameters_at(tt->parameters, 0);
    const double* proportionValues = Parameter_values(proportions);
    Node **nodes = Tree_get_nodes(tt->tree, POSTORDER);
    if(gradient == NULL){
        Parameter* rootHeight = Parameters_at(tt->parameters, 1);
        double* gradient = dvector(tt->tipCount-1);
        for (size_t i = 0; i < nodeCount; i++) {
            if (Node_isleaf(nodes[i])) continue;
            double adj = 0;
            _node_transform_dlog_jacobian_aux(tt, nodes[i], nodes[i], proportionValues, &adj, descendant);
            gradient[nodes[i]->class_id] = adj;
            proportions->grad[nodes[i]->class_id] += adj;
        }
        size_t rootId = Node_class_id(Tree_root(tt->tree));
        rootHeight->grad[0] += gradient[rootId];
        if(rootHeight->transform){
            rootHeight->transform->backward(rootHeight->transform, gradient + rootId);
        }
        if(proportions->transform){
            proportions->transform->backward(proportions->transform, gradient);
        }
        free(gradient);
    }
    else{
    for (size_t i = 0; i < nodeCount; i++) {
        if (Node_isleaf(nodes[i])) continue;
        double adj = 0;
        _node_transform_dlog_jacobian_aux(tt, nodes[i], nodes[i], proportionValues, &adj, descendant);
        gradient[nodes[i]->class_id] += adj;
    }
    }
    free(descendant);
}

void product_of_ratios(Node *node, const double *grad, const double *ratios, double prod, double *out) {
    if (!Node_isleaf(node)) {
        prod *= ratios[node->class_id];
        *out += grad[node->id] * prod;
        product_of_ratios(node->left, grad, ratios, prod, out);
        product_of_ratios(node->right, grad, ratios, prod, out);
    }
}

void node_transform_vjp(TreeTransform *tt, const double *height_gradient, double *gradient) {
    size_t nodeCount = Tree_node_count(tt->tree);
    Node **nodes = Tree_get_nodes(tt->tree, POSTORDER);
    const double* ratioValues = Parameter_values(Parameters_at(tt->parameters, 0));
    for (size_t i = 0; i < nodeCount; i++) {
        Node *node = nodes[i];
        if (!Node_isleaf(node)) {
            double dhi_dri = Node_isroot(node) ? 1 : Node_height(node->parent) - tt->lowers[node->id];
            double accum = height_gradient[node->id];
            product_of_ratios(node->left, height_gradient, ratioValues, 1.0, &accum);
            product_of_ratios(node->right, height_gradient, ratioValues, 1.0, &accum);
            gradient[Node_class_id(node)] = accum * dhi_dri;
        }
    }
}

void TreeTransform_vjp_with_heights(TreeTransform *tt, const double* heights, const double *height_gradient, double *gradient) {
	size_t nodeCount = Tree_node_count(tt->tree);
	Node **nodes = Tree_get_nodes(tt->tree, POSTORDER);
    const double* ratioValues = Parameter_values(Parameters_at(tt->parameters, 0));
	for (size_t i = 0; i < nodeCount; i++) {
		Node *node = nodes[i];
		if (!Node_isleaf(node)) {
			double dhi_dri = Node_isroot(node) ? 1 : heights[node->parent->class_id] - tt->lowers[node->id];
            double accum = height_gradient[node->id];
			product_of_ratios(node->left, height_gradient, ratioValues, 1.0, &accum);
			product_of_ratios(node->right, height_gradient, ratioValues, 1.0, &accum);
			gradient[Node_class_id(node)] = accum * dhi_dri;
		}
	}
}

#pragma endregion Proportion gradient inefficient algorithm

// parameters: parameters we want to differentiate with respect to. They can be constrained or unconstrained.
// ingrad: gradient of a function wrt heights
void TreeTransform_backward(TreeTransform* obj, Parameters* parameters, const double* ingrad){
    if(obj->parameterization == TREE_TRANSFORM_SHIFT){
        _shift_backward(obj, parameters, ingrad);
    }
    else {
        _proportions_backward(obj, parameters, ingrad);
    }
}

void TreeTransform_initialize_from_heights(TreeTransform* tt){
    Tree* atree = tt->tree;
    assert(atree != NULL);
    if(tt->parameterization == TREE_TRANSFORM_SHIFT){
        Node **nodes = Tree_get_nodes(atree, POSTORDER);
        Parameter* shifts = Parameters_at(tt->parameters, 0);
        for (int i = 0; i < Tree_node_count(atree); i++) {
            Node* node = nodes[i];
            if(!Node_isleaf(node)){
                double s = tt->inverse_transform(tt, node);
                Parameter_set_value_at_quietly(shifts, s, node->class_id);
            }
        }
    }
    else{
        if(tt->update_lowers != NULL){
            tt->update_lowers(tt);
        }
        Parameter* ratios = Parameters_at(tt->parameters, 0);
        Parameter* rootHeight = Parameters_at(tt->parameters, 1);
        // unknown-leaf ratios occupy the first `offset` slots, internal proportions follow
        size_t offset = Parameter_size(ratios) - (tt->tipCount - 2);
        Node **nodes = Tree_get_nodes(atree, PREORDER);
        for (int i = 0; i < Tree_node_count(atree); i++) {
            Node* node = nodes[i];
            if(Node_isroot(node)){
                Parameter_set_value_quietly(rootHeight, Node_height(node));
                // printf("%p %s %f %f = node %s %p\n", rootHeight, Parameter_name(rootHeight),Parameter_value(rootHeight), Parameter_lower(rootHeight), node->height->name, node->height);
            }
            else if(!Node_isleaf(node)){
                double s = tt->inverse_transform(tt, node);
                Parameter_set_value_at_quietly(ratios, s, offset + node->class_id);
            }
        }
        if(offset > 0){
            Node** allNodes = Tree_nodes(atree);
            for (size_t i = 0; i < tt->tipCount; i++) {
                if(tt->unknownLeaves[i]){
                    double s = tt->inverse_transform(tt, allNodes[i]);
                    // unknown leaves typically start at their lower bound (height ~ 0),
                    // which maps to ratio ~ 0; use an interior default so the leaf age
                    // is not pinned at the boundary for optimization
                    if(s <= 0.0 || s >= 1.0) s = 0.5;
                    Parameter_set_value_at_quietly(ratios, s, tt->unknownLeafSlot[i]);
                }
            }
        }
	}
    
    Parameters_at(tt->parameters, 0)->listeners->fire(Parameters_at(tt->parameters, 0)->listeners, NULL, NULL, -1);
}

TreeTransform *new_HeightTreeTransform(Tree *tree, tree_transform_t parameterization) {
    TreeTransform *tt = malloc(sizeof(TreeTransform));
    tt->tree = tree;
    tt->unknownLeaves = bvector(Tree_tip_count(tree));
    if(Tree_unknown_leaves(tree) != NULL){
        memcpy(tt->unknownLeaves, Tree_unknown_leaves(tree), sizeof(bool)*Tree_tip_count(tree));
    }
    tt->tipCount = Tree_tip_count(tree);
    tt->unknownLeafSlot = NULL;
    _build_unknown_leaf_slots(tt);
	tt->parameterization = parameterization;
    tt->parameters = new_Parameters(2);
    Parameters_set_name2(tt->parameters, "reparam");
    tt->lowers = dvector(Tree_node_count(tree));
    tt->update = _tree_transform_update;
    tt->update_lowers = _tree_transform_update_lowers;
    tt->inverse_transform = _height_tree_inverse_transform;
    tt->log_jacobian = _node_transform_log_jacobian;
    tt->dlog_jacobian = _node_transform_dlog_jacobian;

	if(parameterization == TREE_TRANSFORM_RATIO_NAIVE){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient;
		tt->vjp = node_transform_vjp;
	}
	else if(parameterization == TREE_TRANSFORM_RATIO){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient_efficient;
		tt->vjp = node_transform_vjp_efficient;
	}
    else if(parameterization == TREE_TRANSFORM_PROPORTION){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient_backprop;
		tt->vjp = node_transform_vjp_backprop;
	}
	else if(parameterization == TREE_TRANSFORM_SHIFT){
        tt->update = _tree_transform_shift_update;
        tt->update_lowers = NULL;
        tt->inverse_transform = _height_tree_inverse_shift_transform;
        tt->log_jacobian = _node_transform_log_jacobian_zero;
        tt->dlog_jacobian = _node_transform_dlog_jacobian_zero;
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient_zero;
		tt->vjp = node_transform_vjp_shift;
	}
	else{
		fprintf(stderr, "Node height reparameterization not recognized\n");
		exit(2);
	}

    // The shift parameterization has a distinct layout: a single parameter with one
    // slot per internal node (indexed by class id, root included) holding
    // height - max(child heights). There is no separate root-height parameter, so it
    // must not share the ratio/proportion construction below (which splits the root
    // out and would size the parameter one slot too small: tipCount-2 vs tipCount-1).
    if (parameterization == TREE_TRANSFORM_SHIFT) {
        size_t shiftCount = tt->tipCount - 1;
        double *shift_values = dvector(shiftCount);
        Node **shiftNodes = Tree_get_nodes(tree, PREORDER);
        for (int i = 0; i < Tree_node_count(tree); i++) {
            Node *node = shiftNodes[i];
            if (!Node_isleaf(node)) {
                shift_values[Node_class_id(node)] = tt->inverse_transform(tt, node);
            }
        }
        Parameter *shifts = new_Parameter_with_postfix2("shifts", "", shift_values,
                                                        shiftCount,
                                                        new_Constraint(0, INFINITY));
        Parameter_set_model(shifts, MODEL_TREE_TRANSFORM);
        Parameters_move(tt->parameters, shifts);
        free(shift_values);
        // shift updates derive heights from child heights, not from lowers, and there
        // is no root-height parameter for tree_transform_collect_lowers to constrain.
        return tt;
    }

    Node **nodes = Tree_get_nodes(tree, PREORDER);

    // Unknown leaf ages become free parameters reparameterized like internal nodes.
    // Their ratios occupy the first `numUnknown` slots of the ratios parameter, in
    // leaf-id order (matching _tree_transform_update); internal-node proportions
    // follow at offset + class_id.
    size_t numUnknown = 0;
    for (size_t i = 0; i < tt->tipCount; i++) {
        if (tt->unknownLeaves[i]) numUnknown++;
    }
    size_t offset = numUnknown;
    size_t propCount = (tt->tipCount - 2) + numUnknown;

    StringBuffer *buffer = new_StringBuffer(10);
    double *ratio_values = dvector(propCount);
    Parameter* ratios = NULL;
    Parameter* rootHeight = NULL;
    for (int i = 0; i < Tree_node_count(tree); i++) {
        Node *node = nodes[i];
        if (Node_isleaf(node)) continue;

        if (Node_isroot(node)) {
            StringBuffer_set_string(buffer, node->name);
            StringBuffer_append_string(buffer, ".reparam");
            // The constraint is updated during lowers collection
            rootHeight = new_Parameter(buffer->c, Node_height(node), new_Constraint(0, INFINITY));
            Parameter_set_model(rootHeight, MODEL_TREE_TRANSFORM);

        } else {
            // position == node->class_id; shifted past the unknown-leaf slots
            ratio_values[offset + Node_class_id(node)] = Node_height(node) / Node_height(Node_parent(node));
        }
    }
    // unknown leaves have lower == 0, so their ratio is height / parent_height
    if (numUnknown > 0) {
        Node** allNodes = Tree_nodes(tree);
        for (size_t i = 0; i < tt->tipCount; i++) {
            if (tt->unknownLeaves[i]) {
                ratio_values[tt->unknownLeafSlot[i]] = Node_height(allNodes[i]) / Node_height(Node_parent(allNodes[i]));
            }
        }
    }
    ratios = new_Parameter_with_postfix2("ratios", "", ratio_values, propCount, new_Constraint(0., 1.));
    Parameter_set_model(ratios, MODEL_TREE_TRANSFORM);
    Constraint_set_flower(ratios->cnstr, 1.e-8);
    Constraint_set_fupper(ratios->cnstr, 1.0 - 1.e-8);

    Parameters_move(tt->parameters, ratios);
    Parameters_move(tt->parameters, rootHeight);
    free(ratio_values);
    free_StringBuffer(buffer);

    tree_transform_collect_lowers(Tree_root(tree), tt, tt->lowers);
    return tt;
}

void TreeTransform_add_tree(TreeTransform* tt, Tree* tree){
    tt->tree = tree;
    tt->tipCount = Tree_tip_count(tree);
    tt->lowers = dvector(Tree_node_count(tree));
    tt->unknownLeaves = bvector(Tree_tip_count(tree));
	memcpy(tt->unknownLeaves, Tree_unknown_leaves(tree), sizeof(bool)*Tree_tip_count(tree));
    _build_unknown_leaf_slots(tt);
    if(tt->update_lowers != NULL){
        tt->update_lowers(tt);
    }
}

// lower and upper constraints have to specified in the json file
TreeTransform *new_HeightTreeTransform2(Parameters* parameters, tree_transform_t parameterization) {
    TreeTransform *tt = malloc(sizeof(TreeTransform));
    tt->tree = NULL;// tree;
    tt->tipCount = 0;//Tree_tip_count(tree);
	tt->parameterization = parameterization;
    tt->parameters = parameters;
    tt->unknownLeaves = NULL;
    tt->unknownLeafSlot = NULL;
    tt->lowers = NULL;//dvector(Tree_node_count(tree));
    tt->update = _tree_transform_update;
    tt->update_lowers = _tree_transform_update_lowers;
    tt->inverse_transform = _height_tree_inverse_transform;
    tt->log_jacobian = _node_transform_log_jacobian;
    tt->dlog_jacobian = _node_transform_dlog_jacobian;
	if(parameterization == TREE_TRANSFORM_RATIO_NAIVE){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient;
		tt->vjp = node_transform_vjp;
	}
	else if(parameterization == TREE_TRANSFORM_RATIO){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient_efficient;
		tt->vjp = node_transform_vjp_efficient;
	}
    else if(parameterization == TREE_TRANSFORM_PROPORTION){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient_backprop;
		tt->vjp = node_transform_vjp_backprop;
	}
	else if(parameterization == TREE_TRANSFORM_SHIFT){
        tt->update = _tree_transform_shift_update;
        tt->update_lowers = NULL;
        tt->inverse_transform = _height_tree_inverse_shift_transform;
        tt->log_jacobian = _node_transform_log_jacobian_zero;
        tt->dlog_jacobian = _node_transform_dlog_jacobian_zero;
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient_zero;
		tt->vjp = node_transform_vjp_shift;
	}
	else{
		fprintf(stderr, "Node height reparameterization not recognized\n");
		exit(2);
	}

    if(parameterization == TREE_TRANSFORM_PROPORTION || parameterization == TREE_TRANSFORM_RATIO || parameterization == TREE_TRANSFORM_RATIO_NAIVE){
        Parameter* ratios = Parameters_at(tt->parameters, 0);
        Constraint_set_flower(ratios->cnstr, 1.e-8);
        Constraint_set_fupper(ratios->cnstr, 1.0 - 1.e-8);
    }
    else if(parameterization == TREE_TRANSFORM_SHIFT){
        Parameter* ratios = Parameters_at(tt->parameters, 0);
        Constraint_set_flower(ratios->cnstr, 1.e-8);
    }

    for(size_t i = 0; i < Parameters_count(parameters); i++){
        Parameter* p = Parameters_at(parameters, i);
        Parameter_set_model(p, MODEL_TREE_TRANSFORM);
    }

    return tt;
}

Model *clone_HeightTreeTransform(Model *self, Hashtable *hash) {
    TreeTransform *tt = self->obj;
    Model *mtree = self->data;
    Model *mtreeclone = NULL;
    if (Hashtable_exists(hash, mtree->name)) {
        mtreeclone = Hashtable_get(hash, mtree->name);
        mtreeclone->ref_count++;  // it is decremented at the end using free
    } else {
        mtreeclone = mtree->clone(mtree, hash);
        Hashtable_add(hash, mtreeclone->name, mtreeclone);
    }

    TreeTransform *ttnew = malloc(sizeof(TreeTransform));
    ttnew->tipCount = tt->tipCount;
    ttnew->tree = mtree->obj;
    ttnew->parameters = new_Parameters(2);
    Parameters_set_name2(ttnew->parameters, Parameters_name2(tt->parameters));
    ttnew->lowers = clone_dvector(tt->lowers, Tree_node_count(ttnew->tree));
    ttnew->unknownLeaves = NULL;
    ttnew->unknownLeafSlot = NULL;
    if(tt->unknownLeaves != NULL){
        ttnew->unknownLeaves = bvector(tt->tipCount);
        memcpy(ttnew->unknownLeaves, tt->unknownLeaves, sizeof(bool)*tt->tipCount);
        _build_unknown_leaf_slots(ttnew);
    }
    ttnew->update = tt->update;
    ttnew->update_lowers = tt->update_lowers;
    ttnew->inverse_transform = tt->inverse_transform;
    ttnew->log_jacobian = tt->log_jacobian;
    ttnew->dlog_jacobian = tt->dlog_jacobian;
    ttnew->log_jacobian_gradient = tt->log_jacobian_gradient;
    ttnew->vjp = tt->vjp;

    Parameter* ratios = clone_Parameter(Parameters_at(tt->parameters, 0));
    Parameter* rootHeight = clone_Parameter(Parameters_at(tt->parameters, 1));

     Parameters_move(ttnew->parameters, ratios);
     Parameters_move(ttnew->parameters, rootHeight);

    Model *clone = new_Model(MODEL_TREE_TRANSFORM, self->name, ttnew);
    for (size_t i = 0; i < Parameters_count(ttnew->parameters); i++) {
        Parameters_at(ttnew->parameters, i)->listeners->add(Parameters_at(ttnew->parameters, i)->listeners, clone);
    }
    Hashtable_add(hash, clone->name, clone);
    return clone;
}

void free_TreeTransform(TreeTransform *tt) {
    free_Parameters(tt->parameters);
    free(tt->lowers);
    free(tt->unknownLeaves);
    free(tt->unknownLeafSlot);
    free(tt);
}

#pragma mark -
#pragma mark Model

// If one proportion parameter changes, all node heights below the corresponding node also need to be updated
// The tree does not really care about which node heights would change but the tree likelihood needs to know
static void _fire_below(Model *self, Parameter* parameter, Node *node){
    if(!Node_isleaf(node)) {
        _fire_below(self, parameter, node->left);
        _fire_below(self, parameter, node->right);
        self->listeners->fire(self->listeners, self, parameter, node->id);
    }
}

// Returns the (slot)-th unknown-age leaf in id order, matching the leaf-ratio
// slot layout used throughout the transform (slots [0, numUnknown)).
static Node* _tree_transform_unknown_leaf(TreeTransform* tt, size_t slot){
    Node** nodes = Tree_nodes(tt->tree);
    size_t index = 0;
    for(size_t i = 0; i < tt->tipCount; i++){
        if(tt->unknownLeaves[i]){
            if(index == slot) return nodes[i];
            index++;
        }
    }
    return NULL;
}

void _tree_transform_model_handle_change(Model *self, Model *model, Parameter* parameter, int index) {
    TreeTransform *tt = (TreeTransform *)self->obj;
    // if the root height changed (parameter at index 1), all node heights need to be updated
    if(index < 0 || parameter == Parameters_at(tt->parameters, 1)){
        self->listeners->fire(self->listeners, self, parameter, -1);
    }
    else{
        // unknown-leaf ratios occupy the first `offset` slots of the ratios parameter
        size_t offset = Parameter_size(Parameters_at(tt->parameters, 0)) - (tt->tipCount - 2);
        if((size_t)index < offset){
            // an unknown leaf age changed: only its own branch length is affected
            Node* leaf = _tree_transform_unknown_leaf(tt, index);
            self->listeners->fire(self->listeners, self, parameter, Node_id(leaf));
        }
        else{
            Node* node = Tree_node(tt->tree, tt->tipCount + (index - offset));
            _fire_below(self, parameter, node);
        }
    }
}

// If one shift parameter changes, all node heights above the corresponding node also need to be updated
// To be consistent we should fire above like the proportions fire below, but the tree likelihood will update the partials above anyway
static void _tree_transform_model_handle_change_shift(Model *self, Model *model, Parameter* parameter, int index) {
    if(index < 0){
        self->listeners->fire(self->listeners, self, parameter, -1);
    }
    else{
        TreeTransform *tt = (TreeTransform *)self->obj;
        Node* node = Tree_node(tt->tree, tt->tipCount + index);
        Node* root = Tree_root(tt->tree);
        while(node != root){
            /* This guarantees that branches above are updated.
            In the treelikelihood the sibling node also need to be updated since they are the start of the branches
            If node + is updated, then branches starting from node - also needs to be updated because h_*=max(h_+,h_-) + s_*.
            All branches need to be updated except the 2 branches starting from C and D.

              -*-
             |   |
             +   -
            | | | |
            A B C D
            */
            // treelikelihood marks both children of node as dirty
            self->listeners->fire(self->listeners, self, parameter, Node_id(node));
            node = node->parent;
        }
    }
}


static void _tree_transform_model_store(Model *self) {
    if(!self->stored){
        TreeTransform *tt = (TreeTransform *)self->obj;
        Parameters_store(tt->parameters);
        self->stored = true;
    }
}

static void _tree_transform_model_restore(Model *self) {
    if(self->stored){
        TreeTransform *tt = (TreeTransform *)self->obj;
        Parameters_restore(tt->parameters);
        self->stored = false;
    }
}

static void _tree_transform_model_accept(Model *self) {
    if(self->stored){
        TreeTransform *tt = (TreeTransform *)self->obj;
        Parameters_accept(tt->parameters);
        self->stored = false;
    }
}

static void _tree_transform_model_free(Model *self) {
    if (self->ref_count == 1) {
        TreeTransform *tt = (TreeTransform *)self->obj;
        free_TreeTransform(tt);
        free_Model(self);
    } else {
        self->ref_count--;
    }
}

static Model *_tree_transform_model_clone(Model *self, Hashtable *hash) {
    if (Hashtable_exists(hash, self->name)) {
        return Hashtable_get(hash, self->name);
    }

    return NULL;
}

static double _tree_transform_model_logP(Model *self) {
	TreeTransform* tt = self->obj;
	self->lp = tt->log_jacobian(tt);
	return self->lp;
}

Model *new_TreeTransformModel(const char *name, TreeTransform *tt, Model *tree) {
    Model *model = new_Model(MODEL_TREE_TRANSFORM, name, tt);
    Parameters_add_listener(tt->parameters, model);
    Parameters_add_parameters_recursively(model->parameters, tt->parameters);

    model->free = _tree_transform_model_free;
    model->clone = _tree_transform_model_clone;
    model->store = _tree_transform_model_store;
    model->restore = _tree_transform_model_restore;
    model->accept = _tree_transform_model_accept;
    model->update = _tree_transform_model_handle_change;

	model->logP = _tree_transform_model_logP;

    model->data = tree;//never used could also create a leak with circular references

    if(tt->parameterization == TREE_TRANSFORM_SHIFT){
        model->update = _tree_transform_model_handle_change_shift;
    }
    return model;
}

void TreeTransformModel_add_tree_model(Model* self, Model* tree){
    self->data = tree; //never used could also create a leak with circular references
    TreeTransform_add_tree(self->obj, tree->obj);
}

Model* new_TreeTransformModel_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"proportions", JSON_OPTIONAL, JSON_ANY},
	    {"ratios", JSON_OPTIONAL, JSON_ANY},
	    {"root_height", JSON_OPTIONAL, JSON_ANY},
	    {"shifts", JSON_OPTIONAL, JSON_ANY},
	    {"transform", JSON_OPTIONAL, JSON_ANY},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

    char* id = get_json_node_value_string(node, "id");
    char* transform_desc = get_json_node_value_string(node, "transform");
    // json_node* treeRef = get_json_node(node, "tree");
    TreeTransform *tt = NULL;
    StringBuffer* buffer = new_StringBuffer(10);

    json_node* ratios_node = get_json_node(node, "ratios");
    json_node* proportions_node = get_json_node(node, "proportions");
    json_node* root_height_node = get_json_node(node, "root_height");
    json_node* shifts_height_node = get_json_node(node, "shifts");

    if((proportions_node != NULL || ratios_node != NULL) && root_height_node != NULL){
        Parameters* parameters = new_Parameters(2);
        StringBuffer_set_string(buffer, id);
        StringBuffer_append_string(buffer, ".ratios");
        Parameters_set_name2(parameters,  buffer->c);

        if(ratios_node != NULL && ratios_node->node_type == MJSON_OBJECT){
            Parameter* ratios = new_Parameter_from_json(ratios_node, hash);
            Parameters_move(parameters, ratios);

            Parameter* rootHeight = new_Parameter_from_json(root_height_node, hash);
            Parameters_move(parameters, rootHeight);

            if(transform_desc != NULL && strcasecmp(transform_desc, "ratio_naive") == 0){
                tt = new_HeightTreeTransform2(parameters, TREE_TRANSFORM_RATIO_NAIVE);
            }
            else{
                tt = new_HeightTreeTransform2(parameters, TREE_TRANSFORM_RATIO);
            }
        }
        else if(proportions_node != NULL && proportions_node->node_type == MJSON_OBJECT){
            Parameter* proportions = new_Parameter_from_json(proportions_node, hash);
            Parameters_move(parameters, proportions);

            Parameter* rootHeight = new_Parameter_from_json(root_height_node, hash);
            Parameters_move(parameters, rootHeight);

            tt = new_HeightTreeTransform2(parameters, TREE_TRANSFORM_PROPORTION);
        }
    }
    else if(shifts_height_node != NULL){
        Parameters* parameters = new_Parameters(1);
        StringBuffer_set_string(buffer, id);
        StringBuffer_append_string(buffer, ".shifts");
        Parameters_set_name2(parameters,  buffer->c);

        if(shifts_height_node->node_type == MJSON_OBJECT){
            Parameter* shifts = new_Parameter_from_json(shifts_height_node, hash);
            Parameters_move(parameters, shifts);
            tt = new_HeightTreeTransform2(parameters, TREE_TRANSFORM_SHIFT);
        }
        
    }
    
    Hashtable_add(hash, Parameters_name2(tt->parameters), tt->parameters);
    for (size_t i = 0; i < Parameters_count(tt->parameters); i++) {
        Hashtable_add(hash, Parameters_name(tt->parameters, i), Parameters_at(tt->parameters, i));
    }

    free_StringBuffer(buffer);

    tt->unknownLeaves = NULL;

    return new_TreeTransformModel(id, tt, NULL);
}
