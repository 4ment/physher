// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>
#include <stdlib.h>

#include "minunit.h"

#include "phyc/datatype.h"
#include "phyc/node.h"
#include "phyc/parameters.h"
#include "phyc/parsimony.h"
#include "phyc/sequence.h"
#include "phyc/sitepattern.h"
#include "phyc/tree.h"

// Build a 4-taxon nucleotide alignment. The unrooted topology is ((A,B),(C,D)),
// so the columns below have hand-computable Fitch scores:
//   AAAA -> 0   AACC -> 1 (A|B vs C|D split)   GGTT -> 1 (same split)
//   GGGT -> 1 (autapomorphy on D, one step on every tree)
// giving a total parsimony score of 3 over 4 sites. The two informative columns
// both favour the A|B vs C|D split, so ((A,B),(C,D)) is the unique best tree.
static SitePattern* build_sitepattern() {
    Sequences* aln = new_Sequences(4);
    aln->datatype = new_NucleotideDataType();
    aln->aligned = true;
    Sequences_add(aln, new_Sequence("A", "AAGG"));
    Sequences_add(aln, new_Sequence("B", "AAGG"));
    Sequences_add(aln, new_Sequence("C", "ACTG"));
    Sequences_add(aln, new_Sequence("D", "ACTT"));
    aln->length = 4;
    SitePattern* sp = new_SitePattern(aln);
    free_Sequences(aln);
    return sp;
}

// An unrooted 4-taxon tree owns 2n-3 = 5 branch lengths; passing a size-5
// distance vector makes tree.c pin the root's right child (zero-branch node).
static Tree* build_tree(const char* newick) {
    Parameter* bl =
        new_Parameter_full("bl", 0.1, 5, new_Constraint(0.0, INFINITY));
    return new_Tree(newick, bl, false);
}

static double total_branch_length(Tree* tree) {
    double sum = 0.0;
    Node** nodes = Tree_nodes(tree);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        if (!Node_isroot(nodes[i])) sum += Node_distance(nodes[i]);
    }
    return sum;
}

static double total_weight(SitePattern* sp) {
    double n = 0.0;
    for (int i = 0; i < sp->count; i++) n += sp->weights[i];
    return n;
}

char* test_parsimony_score() {
    SitePattern* sp = build_sitepattern();
    Tree* tree = build_tree("((A,B),(C,D));");
    Parsimony* pars = new_Parsimony(sp, tree);

    double score = pars->calculate(pars);
    mu_assert(fabs(score - 3.0) < 1.e-9, "parsimony score not matching");

    // A second call must return the cached score without recomputing.
    mu_assert(fabs(pars->calculate(pars) - 3.0) < 1.e-9,
              "cached parsimony score not matching");

    // A rearrangement that separates the (A,B) clade increases the score: both
    // informative columns now conflict with the ((A,C),(B,D)) split (score 5).
    free_Parsimony(pars);
    Tree* tree2 = build_tree("((A,C),(B,D));");
    Parsimony* pars2 = new_Parsimony(sp, tree2);
    mu_assert(pars2->calculate(pars2) > 3.0,
              "suboptimal topology should score worse");

    free_Parsimony(pars2);
    free_Tree(tree);
    free_Tree(tree2);
    free_SitePattern(sp);
    return NULL;
}

char* test_parsimony_branch_lengths() {
    SitePattern* sp = build_sitepattern();
    Tree* tree = build_tree("((A,B),(C,D));");
    Parsimony* pars = new_Parsimony(sp, tree);

    double score = pars->calculate(pars);
    double nsites = total_weight(sp);

    // Every substitution in a most-parsimonious reconstruction is attributed to
    // exactly one branch, so the branch lengths (changes / nsites) must sum to
    // score / nsites. Use min_length 0 so unchanged branches contribute nothing.
    Parsimony_init_branch_lengths(pars, 0.0);

    double sum = total_branch_length(tree);
    mu_assert(fabs(sum * nsites - score) < 1.e-9,
              "sum of parsimony branch lengths must equal the score / nsites");

    // Every branch length is non-negative and finite.
    Node** nodes = Tree_nodes(tree);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        if (Node_isroot(nodes[i])) continue;
        double bl = Node_distance(nodes[i]);
        mu_assert(bl >= 0.0 && isfinite(bl), "branch length must be non-negative");
    }

    free_Parsimony(pars);
    free_Tree(tree);
    free_SitePattern(sp);
    return NULL;
}

char* test_parsimony_constant_sites() {
    // A fully constant alignment has parsimony score 0: no branch carries a
    // substitution, so every branch is set to the requested floor.
    Sequences* aln = new_Sequences(4);
    aln->datatype = new_NucleotideDataType();
    aln->aligned = true;
    Sequences_add(aln, new_Sequence("A", "AAAA"));
    Sequences_add(aln, new_Sequence("B", "AAAA"));
    Sequences_add(aln, new_Sequence("C", "AAAA"));
    Sequences_add(aln, new_Sequence("D", "AAAA"));
    aln->length = 4;
    SitePattern* sp = new_SitePattern(aln);
    free_Sequences(aln);

    Tree* tree = build_tree("((A,B),(C,D));");
    Parsimony* pars = new_Parsimony(sp, tree);

    mu_assert(fabs(pars->calculate(pars)) < 1.e-9,
              "constant alignment must have score 0");

    const double floor = 1.e-6;
    Parsimony_init_branch_lengths(pars, floor);
    Node** nodes = Tree_nodes(tree);
    for (int i = 0; i < Tree_node_count(tree); i++) {
        if (Node_isroot(nodes[i]) ||
            Tree_branch_index(tree, nodes[i]) == NODE_NO_BRANCH)
            continue;
        mu_assert(fabs(Node_distance(nodes[i]) - floor) < 1.e-12,
                  "zero-change branch must be set to the floor");
    }

    free_Parsimony(pars);
    free_Tree(tree);
    free_SitePattern(sp);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_parsimony_score);
    mu_run_test(test_parsimony_branch_lengths);
    mu_run_test(test_parsimony_constant_sites);
    return NULL;
}

RUN_TESTS(all_tests);
