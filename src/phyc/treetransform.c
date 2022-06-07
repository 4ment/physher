//
//  treetransform.c
//  physher
//
//  Created by mathieu on 11/7/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "treetransform.h"

#include "matrix.h"

// Efficient calculation of the ratio and root height gradient, adpated from BEAST.
// Xiang et al. Scalable Bayesian divergence time estimation with ratio transformation (2021).
// https://arxiv.org/abs/2110.13298
static double _epoch_gradient_addition(const Node* node, const Node* child, const double* lowers, Parameters* ratios, const double* ratio_gradient) {
    if (Node_isleaf(child)) {
        return 0.0;
    }
    size_t node_id = Node_id(node);
    size_t child_id = Node_id(child);
    size_t child_class_id = Node_class_id(child);
	double node_ratio = Parameters_value(ratios, Node_class_id(node));
	double child_ratio = Parameters_value(ratios, child_class_id);
    
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

static double _root_height_gradient(Tree* tree, Parameters* ratios, const double* gradient_height) {
	size_t gradientLength = Tree_tip_count(tree) - 1;
    double* multiplierArray = dvector(gradientLength);
    multiplierArray[Node_class_id(Tree_root(tree))] = 1.0;
    Node** nodes = Tree_get_nodes(tree, PREORDER);
    size_t nodeCount = Tree_node_count(tree);
    
    for(size_t i = 1; i < nodeCount; i++){
		if(!Node_isleaf(nodes[i])){
			size_t node_class_id = Node_class_id(nodes[i]);
			double node_ratio = Parameters_value(ratios, node_class_id);
			multiplierArray[node_class_id] = node_ratio * multiplierArray[Node_class_id(nodes[i]->parent)];
		}
    }
    double sum = 0.0;
    for (int i = 0; i < gradientLength; i++) {
        sum += gradient_height[i] * multiplierArray[i];
    }
    free(multiplierArray);
    return sum;
}

void _update_ratios_gradient(Tree* tree, const double* lowers, Parameters* ratios, const double* gradient_height, double* gradient){
    Node** nodes = Tree_get_nodes(tree, POSTORDER);
    for (int i = 0; i < Tree_node_count(tree)-1; i++) {
        if(Node_isleaf(nodes[i])) continue;
        int node_id = Node_id(nodes[i]);
        int node_class_id = Node_class_id(nodes[i]);
        double height = Node_height(nodes[i]);
        gradient[node_class_id] += (height - lowers[node_id])/Parameters_value(ratios, node_class_id) * gradient_height[node_class_id];
        gradient[node_class_id] += _epoch_gradient_addition(nodes[i], nodes[i]->left, lowers, ratios, gradient);
        gradient[node_class_id] += _epoch_gradient_addition(nodes[i], nodes[i]->right, lowers, ratios, gradient);
    }
}


void node_transform_jvp_efficient(TreeTransform *tt, const double *height_gradient, double *gradient){
	memset(gradient, 0.0, sizeof(double)*(tt->tipCount-1));
	_update_ratios_gradient(tt->tree, tt->lowers, tt->parameters, height_gradient, gradient);
	gradient[Node_class_id(Tree_root(tt->tree))] = _root_height_gradient(tt->tree, tt->parameters, height_gradient);
}

void _node_transform_log_jacobian_gradient_efficient(struct TreeTransform *tt, double *gradient) {
	size_t nodeCount = Tree_node_count(tt->tree);
	double *log_time = dvector(tt->tipCount-1);
	unsigned root_class_id = Node_class_id(Tree_root(tt->tree));
	
	Node** nodes = Tree_nodes(tt->tree);
	for (size_t i = tt->tipCount; i < nodeCount; i++) {
		log_time[Node_class_id(nodes[i])] = 1.0 / (Node_height(nodes[i]) - tt->lowers[Node_id(nodes[i])]);
	}
	log_time[root_class_id] = 0.0;
	double* jac_gradient = dvector(tt->tipCount-1);
	memset(jac_gradient, 0.0, sizeof(double)*(tt->tipCount-1));
	_update_ratios_gradient(tt->tree, tt->lowers, tt->parameters, log_time, jac_gradient);
	
	
	gradient[root_class_id] += _root_height_gradient(tt->tree, tt->parameters, log_time);
	
	for (size_t i = tt->tipCount; i < nodeCount-1; i++) {
        size_t node_class_id = Node_class_id(nodes[i]);
        gradient[node_class_id] += jac_gradient [node_class_id] -1.0 / Parameters_value(tt->parameters, node_class_id);
    }
	free(jac_gradient);
	free(log_time);
}

double _node_transform_log_jacobian(TreeTransform *tt){
	double logP = 0.0;
	Node** nodes = Tree_nodes(tt->tree);
    size_t nodeCount1 = Tree_node_count(tt->tree)-1;
	for(size_t i = tt->tipCount; i < nodeCount1; i++){
		logP += log(Node_height(Node_parent(nodes[i])) - tt->lowers[Node_id(nodes[i])]);
	}
	return logP;
}

void tree_transform_update_heights(Node *node, Parameters *parameters, double *lowers) {
	// set height quietly because ratios already notified the tree (and other upstream listeners)
	if (!Node_isleaf(node)) {
        double s = Parameters_value(parameters, Node_class_id(node));
        if (Node_isroot(node)) {
            Node_set_height_quietly(node, s);
        } else {
            double lower = lowers[Node_id(node)];
            Node_set_height_quietly(node, lower + (Node_height(Node_parent(node)) - lower) * s);
        }
        tree_transform_update_heights(node->left, parameters, lowers);
        tree_transform_update_heights(node->right, parameters, lowers);
    }
}

// should be called when topology changes
void tree_transform_collect_lowers(Node *node, TreeTransform *tt, double *lowers) {
    if (!Node_isleaf(node)) {
        tree_transform_collect_lowers(node->left, tt, lowers);
        tree_transform_collect_lowers(node->right, tt, lowers);
        lowers[Node_id(node)] = dmax(lowers[Node_id(Node_left(node))], lowers[Node_id(Node_right(node))]);
        if (Node_isroot(node)) {
            Parameter *p = Parameters_at(tt->parameters, Node_class_id(node));
            Constraint_set_lower(p->cnstr, lowers[Node_id(node)]);
        }
    } else {
        lowers[Node_id(node)] = Node_height(node);
    }
}

void _tree_transform_update(TreeTransform *tt) {
    tree_transform_update_heights(Tree_root(tt->tree), tt->parameters, tt->lowers);
}

void _tree_transform_update_lowers(TreeTransform *tt) {
    tree_transform_collect_lowers(Tree_root(tt->tree), tt, tt->lowers);
}

double _height_tree_inverse_transform(TreeTransform *tt, Node *node) {
    if (Node_isroot(node)) return Node_height(node);
    return (Node_height(node) - tt->lowers[Node_id(node)]) / (Node_height(Node_parent(node)) - tt->lowers[Node_id(node)]);
}

void _node_transform_dlog_jacobian_aux(TreeTransform *tt, const Node *noderef, Node *node, double *dlogP, double *descendant) {
    if (!Node_isleaf(node)) {
        if (!Node_isroot(node) && node != noderef) {
            descendant[Node_id(node)] = descendant[Node_id(node->parent)] * Parameters_value(tt->parameters, Node_class_id(node));
        } else if (!Node_isroot(node)) {
            descendant[Node_id(node)] =
                Node_height(Node_parent(node)) - tt->lowers[Node_id(node)];
        } else {
            descendant[Node_id(node)] = 1;
        }
        _node_transform_dlog_jacobian_aux(tt, noderef, node->left, dlogP, descendant);
        _node_transform_dlog_jacobian_aux(tt, noderef, node->right, dlogP, descendant);

        if (!Node_isroot(node) && node != noderef) {
            *dlogP += descendant[Node_id(node->parent)] / (Node_height(Node_parent(node)) - tt->lowers[Node_id(node)]);
        }
    }
}

// Returns derivative of the log det of the Jacobian
double _node_transform_dlog_jacobian(struct TreeTransform *tt, Node *node) {
    size_t nodeCount = Tree_node_count(tt->tree);
    double *descendant = dvector(nodeCount);
    double adj = 0;
    _node_transform_dlog_jacobian_aux(tt, node, node, &adj, descendant);
    free(descendant);
    return adj;
}

// Update gradient with gradient of the log det of the Jacobian
void _node_transform_log_jacobian_gradient(struct TreeTransform *tt, double *gradient) {
    size_t nodeCount = Tree_node_count(tt->tree);
    double *descendant = dvector(nodeCount);
    Node **nodes = Tree_get_nodes(tt->tree, POSTORDER);
    for (size_t i = 0; i < nodeCount; i++) {
        if (Node_isleaf(nodes[i])) continue;
        double adj = 0;
        _node_transform_dlog_jacobian_aux(tt, nodes[i], nodes[i], &adj, descendant);
        gradient[nodes[i]->class_id] += adj;
    }
    free(descendant);
}

void product_of_ratios(Node *node, const double *grad, Parameters *reparams, double prod, double *out) {
    if (!Node_isleaf(node)) {
        double updatedProd = Parameters_value(reparams, node->class_id) * prod;
        *out += grad[node->class_id] * updatedProd;
        product_of_ratios(node->left, grad, reparams, updatedProd, out);
        product_of_ratios(node->right, grad, reparams, updatedProd, out);
    }
}

void node_transform_jvp(TreeTransform *tt, const double *height_gradient, double *gradient) {
    size_t nodeCount = Tree_node_count(tt->tree);
    Node **nodes = Tree_get_nodes(tt->tree, POSTORDER);
    for (size_t i = 0; i < nodeCount; i++) {
        Node *node = nodes[i];
        if (!Node_isleaf(node)) {
            double dhi_dri = 1;
            if (!Node_isroot(node)) {
                dhi_dri = Node_height(node->parent) - tt->lowers[node->id];
            }
            unsigned reparamID = Node_class_id(node);
            gradient[reparamID] = height_gradient[node->class_id];
            product_of_ratios(node->left, height_gradient, tt->parameters, 1, gradient + reparamID);
            product_of_ratios(node->right, height_gradient, tt->parameters, 1, gradient + reparamID);
            gradient[reparamID] *= dhi_dri;
        }
    }
}

void TreeTransform_jvp_with_heights(TreeTransform *tt, const double* heights, const double *height_gradient, double *gradient) {
	size_t nodeCount = Tree_node_count(tt->tree);
	Node **nodes = Tree_get_nodes(tt->tree, POSTORDER);
	for (size_t i = 0; i < nodeCount; i++) {
		Node *node = nodes[i];
		if (!Node_isleaf(node)) {
			double dhi_dri = 1;
			if (!Node_isroot(node)) {
				dhi_dri = heights[node->parent->class_id] - tt->lowers[node->id];
			}
			unsigned reparamID = Node_class_id(node);
			gradient[reparamID] = height_gradient[node->class_id];
			product_of_ratios(node->left, height_gradient, tt->parameters, 1, gradient + reparamID);
			product_of_ratios(node->right, height_gradient, tt->parameters, 1, gradient + reparamID);
			gradient[reparamID] *= dhi_dri;
		}
	}
}

TreeTransform *new_HeightTreeTransform(Tree *tree, tree_transform_t parameterization) {
    TreeTransform *tt = malloc(sizeof(TreeTransform));
    tt->tree = tree;
    tt->tipCount = Tree_tip_count(tree);
	tt->parameterization = parameterization;
    tt->parameters = new_Parameters(tt->tipCount - 1);
    tt->ratios = new_Parameters(tt->tipCount - 2);
    tt->rootHeight = NULL;
    Parameters_set_name2(tt->parameters, "reparam");
    tt->lowers = dvector(Tree_node_count(tree));
    tt->update = _tree_transform_update;
    tt->update_lowers = _tree_transform_update_lowers;
    tt->inverse_transform = _height_tree_inverse_transform;
    tt->log_jacobian = _node_transform_log_jacobian;
    tt->dlog_jacobian = _node_transform_dlog_jacobian;
	if(parameterization == TREE_TRANSFORM_RATIO_NAIVE){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient;
		tt->jvp = node_transform_jvp;
	}
	else if(parameterization == TREE_TRANSFORM_RATIO){
		tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient_efficient;
		tt->jvp = node_transform_jvp_efficient;
	}
	else{
		fprintf(stderr, "Node height reparameterization not recognized\n");
		exit(2);
	}

    Node **nodes = Tree_get_nodes(tree, PREORDER);

    StringBuffer *buffer = new_StringBuffer(10);
    int *indices = ivector(tt->tipCount - 1);
    size_t count = 0;
    for (int i = 0; i < Tree_node_count(tree); i++) {
        Node *node = nodes[i];
        if (Node_isleaf(node)) continue;

        StringBuffer_set_string(buffer, node->name);
        StringBuffer_append_string(buffer, ".reparam");
        Parameter *p = NULL;
        if (Node_isroot(node)) {
            // The constraint is updated during lowers collection
            p = new_Parameter(buffer->c, Node_height(node), new_Constraint(0, INFINITY));
            tt->rootHeight = p;
            p->refCount++;

        } else {
            p = new_Parameter(buffer->c, Node_height(node) / Node_height(Node_parent(nodes[i])), new_Constraint(0., 1.));
            Parameters_add(tt->ratios, p);
        }
        p->model = MODEL_TREE_TRANSFORM;
        Parameters_move(tt->parameters, p);
        p->id = Node_class_id(node);
        indices[count++] = p->id;
    }
    // we want the reparam parameters to be sorted according to their ids
    // p->id == node->class_id == position in tt->parameters
    Parameters_sort_from_ivector(tt->parameters, indices);
    free(indices);
    free_StringBuffer(buffer);

    tree_transform_collect_lowers(Tree_root(tree), tt, tt->lowers);
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
    ttnew->parameters = new_Parameters(tt->tipCount - 1);
    ttnew->ratios = new_Parameters(tt->tipCount - 2);
    ttnew->rootHeight = NULL;
    Parameters_set_name2(ttnew->parameters, Parameters_name2(tt->parameters));
    Parameters_set_name2(ttnew->ratios, Parameters_name2(tt->ratios));
    ttnew->lowers = clone_dvector(tt->lowers, Tree_node_count(ttnew->tree));
    ttnew->update = tt->update;
    ttnew->update_lowers = tt->update_lowers;
    ttnew->inverse_transform = tt->inverse_transform;
    ttnew->log_jacobian = tt->log_jacobian;
    ttnew->dlog_jacobian = tt->dlog_jacobian;
    ttnew->log_jacobian_gradient = tt->log_jacobian_gradient;
    ttnew->jvp = tt->jvp;

    for (size_t i = 0; i < Parameters_count(tt->parameters); i++) {
        Parameters_move(ttnew->parameters, clone_Parameter(Parameters_at(tt->parameters, i)));
    }
    for (size_t i = 0; i < Parameters_count(tt->parameters)-1; i++) {
        Parameters_add(ttnew->ratios, Parameters_at(ttnew->parameters, i));
    }
    ttnew->rootHeight = Parameters_at(ttnew->parameters, Parameters_count(tt->parameters)-1);
    Parameter_set_name(ttnew->rootHeight, Parameter_name(tt->rootHeight));
    ttnew->rootHeight->refCount++;

    Model *clone = new_Model(MODEL_TREE_TRANSFORM, self->name, ttnew);
    for (size_t i = 0; i < Parameters_count(ttnew->parameters); i++) {
        Parameters_at(ttnew->parameters, i)->listeners->add(Parameters_at(ttnew->parameters, i)->listeners, clone);
    }
    Hashtable_add(hash, clone->name, clone);
    return clone;
}

void free_TreeTransform(TreeTransform *tt) {
    free_Parameters(tt->parameters);
    free_Parameters(tt->ratios);
    free_Parameter(tt->rootHeight);
    free(tt->lowers);
    free(tt);
}

#pragma mark -
#pragma mark Model

void _tree_transform_model_handle_change(Model *self, Model *model, int index) {
    TreeTransform *tt = (TreeTransform *)self->obj;
    self->listeners->fire(self->listeners, self, tt->tipCount + index);
}

static void _tree_transform_model_store(Model *self) {
    TreeTransform *tt = (TreeTransform *)self->obj;
    for (int i = 0; i < Parameters_count(tt->parameters); i++) {
        Parameter_store(Parameters_at(tt->parameters, i));
    }
}

static void _tree_transform_model_restore(Model *self) {
    TreeTransform *tt = (TreeTransform *)self->obj;
    bool height_changed = false;
    Parameter *r = NULL;
    for (int i = 0; i < Parameters_count(tt->parameters); i++) {
        r = Parameters_at(tt->parameters, i);
        if (Parameter_changed(r)) {
            Parameter_restore_quietly(r);
            height_changed = true;
        }
    }

    if (height_changed) {
        r->listeners->fire_restore(r->listeners, self, r->id);
    }
}

void _tree_transform_model_handle_restore(Model *self, Model *model, int index) {
    TreeTransform *tt = (TreeTransform *)self->obj;
    self->listeners->fire_restore(self->listeners, self, tt->tipCount + index);
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
	return tt->log_jacobian(tt);
}

static double _tree_transform_model_dlogP(Model *self, const Parameter* p){
	TreeTransform* tt = self->obj;
	Node* node = Tree_node(tt->tree, p->id + tt->tipCount);
	return tt->dlog_jacobian(tt, node);
}

Model *new_TreeTransformModel(const char *name, TreeTransform *tt, Model *tree) {
    Model *model = new_Model(MODEL_TREE_TRANSFORM, name, tt);
    for (int i = 0; i < Parameters_count(tt->parameters); i++) {
        Parameters_at(tt->parameters, i)->listeners->add(Parameters_at(tt->parameters, i)->listeners, model);
    }

    model->free = _tree_transform_model_free;
    model->clone = _tree_transform_model_clone;
    model->store = _tree_transform_model_store;
    model->restore = _tree_transform_model_restore;
    model->update = _tree_transform_model_handle_change;
    model->handle_restore = _tree_transform_model_handle_restore;
	
	model->logP = _tree_transform_model_logP;
	model->dlogP = _tree_transform_model_dlogP;

    model->data = tree;
    return model;
}
