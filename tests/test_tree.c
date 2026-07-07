// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>
#include <stdlib.h>

#include "minunit.h"

#include "phyc/matrix.h"
#include "phyc/parameters.h"
#include "phyc/tree.h"
#include "phyc/treetransform.h"

char* test_tree_serial() {
    double dates[4] = {5.0, 3.0, 0.0, 1.0};
    char* taxa[4] = {"A", "B", "C", "D"};
    Model* model =
        new_TimeTreeModel_from_newick("(A:2,(B:1.5,(C:2,D:1):2.5):2.5);", taxa, dates);
    TreeModel_set_transform(model, TREE_TRANSFORM_PROPORTION);
    Tree* tree = model->obj;
    Parameters* ratios = get_reparams(tree);

    double expected_node_heights[7] = {5.0, 3.0, 0.0, 1.0, 2.0, 4.5, 7.0};
    Node** nodes = Tree_nodes(tree);
    size_t tipCount = Tree_tip_count(tree);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        if (!Node_isleaf(nodes[i]))
            mu_assert(fabs(expected_node_heights[i] - Node_height(nodes[i])) < 1.e-7,
                      "node height not matching");
    }
    double expected_bounds[7] = {5.0, 3.0, 0.0, 1.0, 1.0, 3.0, 5.0};
    double logP = model->logP(model);
    double logP_expected = log(expected_node_heights[5] - expected_bounds[4]) +
                           log(expected_node_heights[6] - expected_bounds[5]);
    mu_assert(fabs(logP - logP_expected) < 1.e-7, "jacobian not matching");

    // logP must be stored in the model's lp field before being returned
    mu_assert(model->lp == logP, "TreeModel logP not stored in lp");
    mu_assert(model->full_logP(model) == logP, "TreeModel full_logP not matching");
    mu_assert(model->lp == logP, "TreeModel full_logP not stored in lp");

    // the TreeTransformModel (model->data) computes the same log-Jacobian and
    // must likewise store it in its own lp field
    Model* mtt = (Model*)model->data;
    double logP_tt = mtt->logP(mtt);
    mu_assert(fabs(logP_tt - logP_expected) < 1.e-7, "TreeTransform jacobian not matching");
    mu_assert(mtt->lp == logP_tt, "TreeTransformModel logP not stored in lp");
    mu_assert(mtt->full_logP(mtt) == logP_tt, "TreeTransformModel full_logP not matching");
    mu_assert(mtt->lp == logP_tt, "TreeTransformModel full_logP not stored in lp");

    double grads[3];
    double expected_grads[3] = {0.0, 1.1428571939468384, 0.3571428656578064};
    Tree_node_transform_jacobian_gradient(tree, grads);
    for (size_t i = 0; i < Parameters_count(ratios); i++) {
        mu_assert(fabs(grads[i] - expected_grads[i]) < 1.e-7,
                  "gradient jacobian not matching");
    }

    model->free(model);
    return NULL;
}

// The shift parameterization must reserve a slot for every unknown-age leaf and
// drive that leaf's height from it, both in the forward map (update) and the
// reverse map (vjp). Here leaf C has an unknown date, so the shift parameter must
// have one slot per internal node plus one for C. We finite-difference an
// arbitrary linear objective f(heights) = sum w_i * height_i against the gradient
// produced by the transform's vjp to confirm the unknown-leaf slot participates
// correctly (including when the leaf becomes the taller child and propagates up).
char* test_tree_shift_unknown_leaf() {
    double dates[4] = {5.0, 3.0, -1.0, 1.0};  // C (index 2) is unknown
    char* taxa[4] = {"A", "B", "C", "D"};
    Model* model =
        new_TimeTreeModel_from_newick("(A:2,(B:1.5,(C:2,D:1):2.5):2.5);", taxa, dates);
    TreeModel_set_transform(model, TREE_TRANSFORM_SHIFT);
    Tree* tree = model->obj;
    Model* mtt = (Model*)model->data;
    TreeTransform* tt = mtt->obj;
    Parameter* shifts = Parameters_at(tt->parameters, 0);

    size_t tipCount = Tree_tip_count(tree);
    size_t nShift = Parameter_size(shifts);
    // one shift per internal node (tipCount - 1) plus one for the unknown leaf C
    mu_assert(nShift == (tipCount - 1) + 1,
              "shift parameter must reserve a slot for the unknown-age leaf");

    // Make the unknown leaf the taller child of its parent so its shift also
    // propagates into ancestor heights (a stronger test than an isolated slot).
    // The unknown-leaf slots occupy [0, numUnknown); here slot 0 is leaf C.
    Parameter_set_value_at_quietly(shifts, 5.0, 0);
    tt->update(tt);

    size_t nodeCount = Tree_node_count(tree);
    Node** nodes = Tree_nodes(tree);
    double* w = dvector(nodeCount);  // weight per node id
    for (size_t i = 0; i < nodeCount; i++) w[i] = 0.3 + 0.7 * sin(1.0 + (double)i);

    double* analytic = dvector(nShift);
    tt->vjp(tt, w, analytic);  // d f / d shift, via the reverse map

    const double eps = 1.e-6;
    for (size_t j = 0; j < nShift; j++) {
        double v0 = Parameter_value_at(shifts, j);

        Parameter_set_value_at_quietly(shifts, v0 + eps, j);
        tt->update(tt);
        double fp = 0.0;
        for (size_t i = 0; i < nodeCount; i++)
            fp += w[Node_id(nodes[i])] * Node_height(nodes[i]);

        Parameter_set_value_at_quietly(shifts, v0 - eps, j);
        tt->update(tt);
        double fm = 0.0;
        for (size_t i = 0; i < nodeCount; i++)
            fm += w[Node_id(nodes[i])] * Node_height(nodes[i]);

        Parameter_set_value_at_quietly(shifts, v0, j);
        tt->update(tt);

        double fd = (fp - fm) / (2.0 * eps);
        mu_assert(fabs(fd - analytic[j]) < 1.e-4,
                  "shift transform vjp does not match finite differences");
    }

    free(w);
    free(analytic);
    model->free(model);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_tree_serial);
    mu_run_test(test_tree_shift_unknown_leaf);
    return NULL;
}

RUN_TESTS(all_tests);