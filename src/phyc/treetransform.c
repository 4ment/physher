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
        gradient[nodes[i]->id] += adj;
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
            gradient[node->id] = height_gradient[node->class_id];
            product_of_ratios(node->left, height_gradient, tt->map, tt->parameters, 1, gradient + node->id);
            product_of_ratios(node->right, height_gradient, tt->map, tt->parameters, 1, gradient + node->id);
            gradient[node->id] *= dhi_dri;
        }
    }
}

TreeTransform *new_HeightTreeTransform(Tree *tree) {
    TreeTransform *tt = malloc(sizeof(TreeTransform));
    tt->tree = tree;
    tt->parameters = new_Parameters(Tree_tip_count(tree) + 1);
    Parameters_set_name2(tt->parameters, "reparam");
    tt->map = uivector(Tree_node_count(tree));
    tt->lowers = dvector(Tree_node_count(tree));
    tt->update = _tree_transform_update;
    tt->update_lowers = _tree_transform_update_lowers;
    tt->parameter_of_node = _parameter_of_node;
    tt->inverse_transform = _height_tree_inverse_transform;
    tt->dlog_jacobian = _node_transform_dlog_jacobian;
    tt->log_jacobian_gradient = _node_transform_log_jacobian_gradient;
    tt->jvp = node_transform_jvp;

    Node **nodes = Tree_get_nodes(tree, PREORDER);
    int count = 0;
    for (int i = 0; i < Tree_node_count(tree); i++) {
        Node *node = nodes[i];
        Parameter *p = NULL;
        if (Node_isroot(node)) {
            // The constraint is updated during lowers collection
            p = new_Parameter(node->name, Node_height(node), new_Constraint(0, INFINITY));
            p->model = MODEL_TREE;
            Parameters_move(tt->parameters, p);
            p->id = -Node_id(node) - 1;
            tt->map[Node_id(node)] = count++;
        } else if (!Node_isleaf(node)) {
            p = new_Parameter(node->name, Node_height(node) / Node_height(Node_parent(nodes[i])), new_Constraint(0.00001, 0.9999));
            p->model = MODEL_TREE;
            Parameters_move(tt->parameters, p);
            p->id = -Node_id(node) - 1;
            tt->map[Node_id(node)] = count++;
        }
    }

    tree_transform_collect_lowers(Tree_root(tree), tt, tt->lowers);
    return tt;
}

void free_TreeTransform(TreeTransform *tt) {
    free_Parameters(tt->parameters);
    free(tt->map);
    free(tt->lowers);
    free(tt);
}
