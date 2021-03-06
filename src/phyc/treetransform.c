//
//  treetransform.c
//  physher
//
//  Created by mathieu on 11/7/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "treetransform.h"

#include "matrix.h"

void tree_transform_update_heights(Node *node, Parameters *parameters,
                                   unsigned *map, double *lowers) {
    if (!Node_isleaf(node)) {
        double s = Parameters_value(parameters, map[Node_id(node)]);
        if (Node_isroot(node)) {
            Node_set_height(node, s);
        } else {
            double lower = lowers[Node_id(node)];
            Node_set_height(node, lower + (Node_height(Node_parent(node)) - lower) * s);
        }
        tree_transform_update_heights(node->left, parameters, map, lowers);
        tree_transform_update_heights(node->right, parameters, map, lowers);
    }
}

// should be called when topology changes
void tree_transform_collect_lowers(Node *node, TreeTransform *tt, double *lowers) {
    if (!Node_isleaf(node)) {
        tree_transform_collect_lowers(node->left, tt, lowers);
        tree_transform_collect_lowers(node->right, tt, lowers);
        lowers[Node_id(node)] = dmax(lowers[Node_id(Node_left(node))], lowers[Node_id(Node_right(node))]);
        if (Node_isroot(node)) {
            Parameter *p = Parameters_at(tt->parameters, tt->map[Node_id(node)]);
            Constraint_set_lower(p->cnstr, lowers[Node_id(node)]);
        }
    } else {
        lowers[Node_id(node)] = Node_height(node);
    }
}

void _tree_transform_update(TreeTransform *tt) {
    tree_transform_update_heights(Tree_root(tt->tree), tt->parameters, tt->map, tt->lowers);
}

void _tree_transform_update_lowers(TreeTransform *tt) {
    tree_transform_collect_lowers(Tree_root(tt->tree), tt, tt->lowers);
}

Parameter *_parameter_of_node(struct TreeTransform *tt, Node *node) {
    return Parameters_at(tt->parameters, tt->map[Node_id(node)]);
}

Node *_node_of_parameter(struct TreeTransform *tt, Parameter *p) {
    return tt->map_to_node[p->id];
}

double _height_tree_inverse_transform(TreeTransform *tt, Node *node) {
    if (Node_isroot(node)) return Node_height(node);
    return (Node_height(node) - tt->lowers[Node_id(node)]) / (Node_height(Node_parent(node)) - tt->lowers[Node_id(node)]);
}

void _node_transform_dlog_jacobian_aux(TreeTransform *tt, const Node *noderef, Node *node, double *dlogP, double *descendant) {
    if (!Node_isleaf(node)) {
        if (!Node_isroot(node) && node != noderef) {
            descendant[Node_id(node)] = descendant[Node_id(node->parent)] * Parameters_value(tt->parameters, tt->map[Node_id(node)]);
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

void product_of_ratios(Node *node, const double *grad, unsigned *map, Parameters *reparams, double prod, double *out) {
    if (!Node_isleaf(node)) {
        double updatedProd = Parameters_value(reparams, map[node->id]) * prod;
        *out += grad[node->class_id] * updatedProd;
        product_of_ratios(node->left, grad, map, reparams, updatedProd, out);
        product_of_ratios(node->right, grad, map, reparams, updatedProd, out);
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
            unsigned reparamID = tt->map[Node_id(node)];
            gradient[reparamID] = height_gradient[node->class_id];
            product_of_ratios(node->left, height_gradient, tt->map, tt->parameters, 1, gradient + reparamID);
            product_of_ratios(node->right, height_gradient, tt->map, tt->parameters, 1, gradient + reparamID);
            gradient[reparamID] *= dhi_dri;
        }
    }
}

TreeTransform *new_HeightTreeTransform(Tree *tree) {
    TreeTransform *tt = malloc(sizeof(TreeTransform));
    tt->tree = tree;
    tt->parameters = new_Parameters(Tree_tip_count(tree) + 1);
    Parameters_set_name2(tt->parameters, "reparam");
    tt->map = uivector(Tree_node_count(tree));
    tt->map_to_node = malloc(sizeof(Node *) * (Tree_tip_count(tree) - 1));
    tt->lowers = dvector(Tree_node_count(tree));
    tt->update = _tree_transform_update;
    tt->update_lowers = _tree_transform_update_lowers;
    tt->parameter_of_node = _parameter_of_node;
    tt->node_of_parameter = _node_of_parameter;
    tt->inverse_transform = _height_tree_inverse_transform;
    tt->dlog_jacobian = _node_transform_dlog_jacobian;
    tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient;
    tt->jvp = node_transform_jvp;

    Node **nodes = Tree_get_nodes(tree, PREORDER);

    StringBuffer *buffer = new_StringBuffer(10);
    int *indices = ivector(Tree_tip_count(tt->tree) - 1);
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

        } else {
            p = new_Parameter(buffer->c, Node_height(node) / Node_height(Node_parent(nodes[i])), new_Constraint(0.00001, 0.9999));
        }
        p->model = MODEL_TREE;
        Parameters_move(tt->parameters, p);
        p->id = Node_class_id(node);
        tt->map_to_node[p->id] = node;
        tt->map[Node_id(node)] = p->id;
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
    ttnew->tree = mtree->obj;
    ttnew->parameters = new_Parameters(Tree_tip_count(ttnew->tree) + 1);
    Parameters_set_name2(ttnew->parameters, "reparam");
    ttnew->map = clone_uivector(tt->map, Tree_node_count(ttnew->tree));
    ttnew->map_to_node = malloc(sizeof(Node *) * (Tree_tip_count(ttnew->tree) - 1));
    ttnew->lowers = clone_dvector(tt->lowers, Tree_node_count(ttnew->tree));
    ttnew->update = tt->update;
    ttnew->update_lowers = tt->update_lowers;
    ttnew->parameter_of_node = tt->parameter_of_node;
    ttnew->node_of_parameter = tt->node_of_parameter;
    ttnew->inverse_transform = tt->inverse_transform;
    ttnew->dlog_jacobian = tt->dlog_jacobian;
    ttnew->log_jacobian_gradient = tt->log_jacobian_gradient;
    ttnew->jvp = tt->jvp;

    for (size_t i = 0; i < Parameters_count(tt->parameters); i++) {
        Parameters_move(ttnew->parameters, clone_Parameter(Parameters_at(tt->parameters, i)));
    }

    Node **nodes = Tree_nodes(ttnew->tree);
    for (int i = 0; i < Tree_node_count(ttnew->tree); i++) {
        Node *node = nodes[i];
        if (!Node_isleaf(node)) {
            tt->map_to_node[node->class_id] = node;
        }
    }

    Model *clone = new_Model(MODEL_TREE_TRANSFORM, self->name, ttnew);
    Hashtable_add(hash, clone->name, clone);
    return clone;
}

void free_TreeTransform(TreeTransform *tt) {
    free_Parameters(tt->parameters);
    free(tt->map);
    free(tt->map_to_node);
    free(tt->lowers);
    free(tt);
}

#pragma mark -
#pragma mark Model

void _tree_transform_model_handle_change(Model *self, Model *model, int index) {
    TreeTransform *tt = (TreeTransform *)self->obj;
    self->listeners->fire(self->listeners, self, tt->map_to_node[index]->id);
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
    self->listeners->fire_restore(self->listeners, self, tt->map_to_node[index]->id);
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

    model->data = tree;
    return model;
}
