// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "minunit.h"

#include "phyc/matrix.h"
#include "phyc/nj.h"
#include "phyc/node.h"
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

// Collect the unrooted tree as a map from split to branch length. A split is the
// set of leaf ids on one side of a branch, canonicalized by keeping whichever
// side does not contain leaf 0 so that the two orientations of a branch agree.
// The root is a degree-2 artifact of the rooted storage: its two branches are a
// single branch of the unrooted tree and their splits are complementary, hence
// canonicalize to the same key, so accumulating lengths per key merges them.
static void collect_splits(Tree* tree, uint64_t* splits, double* lengths,
                           size_t* count) {
    Node** nodes = Tree_get_nodes(tree, POSTORDER);
    size_t nodeCount = Tree_node_count(tree);
    const double* bls = Tree_branch_lengths(tree);
    uint64_t* below = calloc(nodeCount, sizeof(uint64_t));
    uint64_t all = 0;
    *count = 0;

    for (size_t i = 0; i < nodeCount; i++) {
        Node* n = nodes[i];
        below[Node_id(n)] = Node_isleaf(n)
                                ? (uint64_t)1 << Node_id(n)
                                : below[Node_id(Node_left(n))] | below[Node_id(Node_right(n))];
        if (Node_isleaf(n)) all |= (uint64_t)1 << Node_id(n);
    }

    for (size_t i = 0; i < nodeCount; i++) {
        Node* n = nodes[i];
        if (Node_isroot(n)) continue;
        uint64_t key = below[Node_id(n)];
        if (key & 1) key = all & ~key;  // side without leaf 0
        size_t j = 0;
        while (j < *count && splits[j] != key) j++;
        if (j == *count) {
            splits[j] = key;
            lengths[j] = 0.0;
            (*count)++;
        }
        lengths[j] += bls[Node_id(n)];
    }
    free(below);
}

// Rerooting must leave the unrooted tree untouched: same splits, same lengths.
// It is the operation that changes which node is the right child of the root,
// and therefore which node is pinned to a zero-length branch and how every other
// node maps onto an entry of the distance vector, so a bookkeeping mistake shows
// up here as branch lengths attached to the wrong split.
char* test_tree_reroot() {
    size_t tipCount = 6;
    size_t branchCount = 2 * tipCount - 3;
    Parameter* bls = new_Parameter_full("bl", 0.1, branchCount, NULL);
    Tree* tree = new_Tree(
        "(A:0.11,B:0.22,((C:0.33,D:0.44):0.55,(E:0.66,F:0.77):0.88):0.99);", bls, true);
    Tree_init_branch_lengths(tree);

    size_t nodeCount = Tree_node_count(tree);
    uint64_t expectedSplits[16];
    double expectedLengths[16];
    size_t expectedCount = 0;
    collect_splits(tree, expectedSplits, expectedLengths, &expectedCount);
    mu_assert(expectedCount == branchCount, "unrooted tree must have 2n-3 branches");

    uint64_t splits[16];
    double lengths[16];
    size_t count = 0;

    for (size_t i = 0; i < nodeCount; i++) {
        Tree_reroot(tree, Tree_node(tree, i));

        Node* root = Tree_root(tree);
        mu_assert(Tree_zero_branch_node(tree) == Node_right(root),
                  "the pinned node must be the right child of the root");
        mu_assert(Tree_branch_lengths(tree)[Node_id(Node_right(root))] == 0.0,
                  "the pinned branch must have length zero");

        collect_splits(tree, splits, lengths, &count);
        mu_assert(count == expectedCount, "rerooting changed the number of branches");
        for (size_t j = 0; j < expectedCount; j++) {
            size_t k = 0;
            while (k < count && splits[k] != expectedSplits[j]) k++;
            mu_assert(k < count, "rerooting changed the topology");
            mu_assert(fabs(lengths[k] - expectedLengths[j]) < 1.e-10,
                      "rerooting changed a branch length");
        }
    }

    free_Tree(tree);  // releases bls, which the tree took a reference to
    return NULL;
}


// Neighbour joining must recover an additive distance matrix exactly, and the
// tree it hands back must be wired to the branch-length Parameter it was given:
// one slot per branch of the unrooted tree (2n-3), the right child of the root
// pinned to a zero-length branch, and the two branches incident to the root
// merged into the single unrooted branch they represent.
char* test_tree_nj() {
    // Distances read off ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6,E:0.7); which NJ
    // recovers exactly because the matrix is additive.
    const char* taxa[5] = {"A", "B", "C", "D", "E"};
    double d[5][5] = {
        {0.0, 0.3, 1.4, 1.5, 1.1},
        {0.3, 0.0, 1.5, 1.6, 1.2},
        {1.4, 1.5, 0.0, 0.9, 1.7},
        {1.5, 1.6, 0.9, 0.0, 1.8},
        {1.1, 1.2, 1.7, 1.8, 0.0},
    };
    double* matrix[5] = {d[0], d[1], d[2], d[3], d[4]};

    size_t tipCount = 5;
    size_t branchCount = 2 * tipCount - 3;
    Parameter* bls = new_Parameter_full("bl", 0.1, branchCount, NULL);
    Tree* tree = new_NJ(taxa, tipCount, matrix, bls);

    mu_assert(Tree_tip_count(tree) == tipCount, "NJ tree has the wrong tip count");
    mu_assert(Tree_node_count(tree) == 2 * tipCount - 1,
              "NJ tree is not a rooted binary tree");

    Node* root = Tree_root(tree);
    mu_assert(Tree_zero_branch_node(tree) == Node_right(root),
              "the pinned node must be the right child of the root");
    mu_assert(Tree_branch_lengths(tree)[Node_id(Node_right(root))] == 0.0,
              "the pinned branch must have length zero");

    // leaf ids are handed out in postorder, so look them up by name
    size_t leafID[5];
    Node** nodes = Tree_nodes(tree);
    for (size_t t = 0; t < tipCount; t++) {
        size_t i = 0;
        while (i < Tree_node_count(tree) &&
               (!Node_isleaf(nodes[i]) || strcmp(Node_name(nodes[i]), taxa[t]) != 0)) {
            i++;
        }
        mu_assert(i < Tree_node_count(tree), "NJ dropped a taxon");
        leafID[t] = Node_id(nodes[i]);
    }

    // the 2n-3 splits of the original tree, keyed the way collect_splits keys
    // them (a bitmask over leaf ids, complemented so leaf 0 is always outside)
    uint64_t expectedSplits[7];
    double expectedLengths[7] = {0.1, 0.2, 0.4, 0.5, 0.7, 0.3, 0.6};
    for (size_t t = 0; t < tipCount; t++) {
        expectedSplits[t] = (uint64_t)1 << leafID[t];
    }
    expectedSplits[5] = ((uint64_t)1 << leafID[0]) | ((uint64_t)1 << leafID[1]);
    expectedSplits[6] = ((uint64_t)1 << leafID[2]) | ((uint64_t)1 << leafID[3]);
    uint64_t all = 0;
    for (size_t t = 0; t < tipCount; t++) all |= (uint64_t)1 << leafID[t];
    for (size_t j = 0; j < 7; j++) {
        if (expectedSplits[j] & 1) expectedSplits[j] = all & ~expectedSplits[j];
    }

    uint64_t splits[8];
    double lengths[8];
    size_t count = 0;
    collect_splits(tree, splits, lengths, &count);
    mu_assert(count == branchCount, "NJ tree must have 2n-3 branches");
    for (size_t j = 0; j < branchCount; j++) {
        size_t k = 0;
        while (k < count && splits[k] != expectedSplits[j]) k++;
        mu_assert(k < count, "NJ did not recover the topology");
        mu_assert(fabs(lengths[k] - expectedLengths[j]) < 1.e-10,
                  "NJ did not recover a branch length");
    }

    // every branch length lives in the Parameter that was handed in, at the slot
    // Tree_branch_index reports for its node
    double total = 0.0;
    for (size_t i = 0; i < Tree_node_count(tree); i++) {
        Node* n = Tree_node(tree, i);
        if (Node_isroot(n) || n == Tree_zero_branch_node(tree)) continue;
        double bl = Parameter_value_at(bls, Tree_branch_index(tree, n));
        mu_assert(bl == Tree_branch_lengths(tree)[Node_id(n)],
                  "branch length does not come from the branch-length parameter");
        mu_assert(bl >= BL_MIN, "NJ branch length below the parameter lower bound");
        total += bl;
    }
    mu_assert(fabs(total - 2.8) < 1.e-10, "NJ tree length not matching");

    // a TreeModel must listen to the parameter it was built with: changing an
    // entry has to invalidate the cached branch lengths
    Model* mtree = new_TreeModel("tree", tree);
    Node* leafA = Tree_node(tree, leafID[0]);
    size_t slot = Tree_branch_index(tree, leafA);
    Parameter_set_value_at(bls, 0.25, slot);
    mu_assert(fabs(Tree_branch_lengths(tree)[leafID[0]] - 0.25) < 1.e-10,
              "TreeModel does not listen to the branch-length parameter");

    mtree->free(mtree);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_tree_serial);
    mu_run_test(test_tree_shift_unknown_leaf);
    mu_run_test(test_tree_reroot);
    mu_run_test(test_tree_nj);
    return NULL;
}

RUN_TESTS(all_tests);