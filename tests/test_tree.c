#include <stdlib.h>

#include "minunit.h"

#include "phyc/parameters.h"
#include "phyc/tree.h"
#include "phyc/treetransform.h"

char* test_tree_serial() {
    double dates[4] = {5.0, 3.0, 0.0, 1.0};
    char* taxa[4] = {"A", "B", "C", "D"};
    Model* model =
        new_TimeTreeModel_from_newick("(A:2,(B:1.5,(C:2,D:1):2.5):2.5);", taxa, dates);
    TreeModel_set_transform(model, TREE_TRANSFORM_RATIO);
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

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_tree_serial);
    return NULL;
}

RUN_TESTS(all_tests);