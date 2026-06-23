//
//  test_coalescent.c
//  physher
//  Created by Mathieu Fourment on 25/04/2020.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//
#include <assert.h>
#include <ctype.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "minunit.h"
#include "phyc/demographicmodels.h"
#include "phyc/gradient.h"
#include "phyc/parameters.h"
#include "phyc/tree.h"
#include "phyc/treetransform.h"

struct gsl_data {
    Model* model;
    size_t index;
};

double f_coalescent(double x, void* params) {
    struct gsl_data* data = params;
    Model* model = data->model;
    Coalescent* coal = model->obj;
    if (data->index >= Parameters_count(coal->p)) {
        Parameters* reparams = get_reparams(coal->tree);
        Parameters_set_value(reparams, data->index - Parameters_count(coal->p), x);
    } else {
        Parameters_set_value(coal->p, data->index, x);
    }
    return model->logP(model);
}

double gsl_coalescent_dlogPdx(struct gsl_data* data, double value, double eps) {
    gsl_function F;
    F.function = &f_coalescent;
    F.params = data;
    double result, abserr;
    gsl_deriv_central(&F, value, eps, &result, &abserr);
    return result;
}

char* test_skyride() {
    double dates[4] = {0.0, 0.0, 0.0, 0.0};
    char* taxa[4] = {"a", "b", "c", "d"};
    Model* mtree =
        new_TimeTreeModel_from_newick("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;

    Parameters* ps = new_Parameters(3);
    Parameters* reparams = get_reparams(tree);
    Parameter* ratios = Parameters_at(reparams, 0);
    Parameter* rootHeight = Parameters_at(reparams, 1);
    Parameters_add_parameters(ps, reparams);

    double thetas[3] = {3., 10., 4.};
    Parameter* popSizes =
        new_Parameter2("popSizes", thetas, 3, new_Constraint(0, INFINITY));
    Parameters_add(ps, popSizes);
    Coalescent* coal = new_SkyrideCoalescent(tree, popSizes);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -11.48749174278204;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching");

    model->gradient(model, ps);

    double true_gradient[6] = {3., 0.2, 0.5, -10.2, -7.4, -0.5583333};
    for (size_t i = 0; i < 3; i++) {
        mu_assert(fabs(popSizes->grad[i] - true_gradient[i] / thetas[i]) < 0.00001,
                  "d.logP/d.theta not matching");
    }
    for (size_t i = 0; i < 2; i++) {
        mu_assert(fabs(ratios->grad[i] - true_gradient[i + 3]) < 0.00001,
                  "d.logP/d.ratios not matching");
    }
    mu_assert(fabs(rootHeight->grad[0] - true_gradient[5]) < 0.00001,
              "d.logP/d.root not matching");

    model->free(model);
    mtree->free(mtree);
    free_Parameters(ps);
    return NULL;
}

char* test_skygrid() {
    double dates[4] = {0.0, 0.0, 0.0, 0.0};
    char* taxa[4] = {"a", "b", "c", "d"};
    Model* mtree =
        new_TimeTreeModel_from_newick("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates);
    Tree* tree = mtree->obj;
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);

    Parameters* ps = new_Parameters(5);
    Parameters* reparams = get_reparams(tree);
    Parameter* ratios = Parameters_at(reparams, 0);
    Parameter* rootHeight = Parameters_at(reparams, 1);
    Parameters_add_parameters(ps, reparams);

    double thetas[5] = {3., 10., 4., 2., 3.};
    Parameter* popSizes =
        new_Parameter2("popSizes", thetas, 5, new_Constraint(0, INFINITY));
    Parameters_add(ps, popSizes);
    double cutoff = 10;
    int grid = 5;
    Coalescent* coal = new_GridCoalescent(tree, popSizes, grid, cutoff);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -11.8751856;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching");

    model->gradient(model, ps);

    // the true derivatives are wrt log(theta) so the Jacobian adjustment needs to be
    // removed
    double true_gradient[8] = {3.5, 0.75, 0.1250, 1.25, -0.333333, -6.0, -10.0, -0.75};
    for (size_t i = 0; i < 5; i++) {
        mu_assert(fabs(popSizes->grad[i] - true_gradient[i] / thetas[i]) < 0.00001,
                  "d.logP/d.theta not matching");
    }
    for (size_t i = 0; i < 2; i++) {
        mu_assert(fabs(ratios->grad[i] - true_gradient[i + 5]) < 0.00001,
                  "d.logP/d.ratios not matching");
    }
    mu_assert(fabs(rootHeight->grad[0] - true_gradient[7]) < 0.00001,
              "d.logP/d.root not matching");

    // call second time to check it accumulates
    model->gradient(model, ps);
    for (size_t i = 0; i < 5; i++) {
        mu_assert(
            fabs(popSizes->grad[i] / 2.0 - true_gradient[i] / thetas[i]) < 0.00001,
            "d.logP/d.theta not matching");
    }
    for (size_t i = 0; i < 2; i++) {
        mu_assert(fabs(ratios->grad[i] / 2.0 - true_gradient[i + 5]) < 0.00001,
                  "d.logP/d.ratios not matching");
    }
    mu_assert(fabs(rootHeight->grad[0] / 2.0 - true_gradient[7]) < 0.00001,
              "d.logP/d.root not matching");

    model->free(model);
    mtree->free(mtree);
    free_Parameters(ps);
    return NULL;
}

char* test_constant_with_transform(tree_transform_t transform) {
    double dates[4] = {0.0, 0.0, 0.0, 0.0};
    char* taxa[4] = {"a", "b", "c", "d"};
    Model* mtree =
        new_TimeTreeModel_from_newick("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates);
    Tree* tree = mtree->obj;

    Parameters* ps = new_Parameters(1);

    double value = 3;
    Parameter* N = new_Parameter("", value, new_Constraint(0, INFINITY));
    Coalescent* coal = new_ConstantCoalescent(tree, N);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -13.2958368660;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching");

    double trueHeightGradient[3] = {-1.0, -0.6666666666666667, -0.3333333333333333};
    Parameters_add(ps, N);
    Node** nodes = Tree_nodes(tree);
    for (size_t i = 0; i < Tree_node_count(tree); i++) {
        if (!Node_isleaf(nodes[i])) {
            Parameters_add(ps, nodes[i]->height);
        }
    }
    model->gradient(model, ps);

    mu_assert(fabs(N->grad[0] - 2.333333333333333) < 1.e-10,
              "test_constant: d.logP/d.theta not matching");

    for (size_t i = 0; i < 3; i++) {
        mu_assert(
            fabs(Parameters_at(ps, i + 1)->grad[0] - trueHeightGradient[i]) < 1.e-10,
            "test_constant: d.logP/d.time not matching");
    }

    value = 7;
    Parameter_set_value(N, value);
    logP = model->logP(model);
    logP2 = -10.1234447329;
    mu_assert(fabs(logP - logP2) < 1.e-10, "logP not matching at 7");

    double eps = 0.0001;
    Parameter_set_value(N, value);
    Parameters_zero_grad(ps);

    model->gradient(model, ps);

    double trueHeightGradient2[3] = {-0.42857142857142855, -0.2857142857142857,
                                     -0.14285714285714285};
    double trueThetaGradient = 0.18367346938775514;

    mu_assert(fabs(N->grad[0] - trueThetaGradient) < 1.e-10,
              "test_constant: d.logP/d.theta not matching at 7");

    for (size_t i = 0; i < 3; i++) {
        mu_assert(
            fabs(Parameters_at(ps, i + 1)->grad[0] - trueHeightGradient2[i]) < 1.e-10,
            "test_constant: d.logP/d.time not matching at 7");
    }

    // struct gsl_data data = {model, 0};
    // double dlogPdp2 = gsl_coalescent_dlogPdx(&data, value, eps);

    // mu_assert(fabs(N->grad[0] - dlogPdp2) < 0.00001, "logP not matching after
    // change");

    model->gradient(model, ps);
    mu_assert(fabs(N->grad[0] / 2.0 - trueThetaGradient) < 1.e-7,
              "test_constant: d.logP/d.theta not matching when accumulating");
    for (size_t i = 0; i < 3; i++) {
        mu_assert(fabs(Parameters_at(ps, i + 1)->grad[0] / 2.0 -
                       trueHeightGradient2[i]) < 1.e-10,
                  "test_constant: d.logP/d.time not matching when accumulating");
    }

    // use reparameterization
    TreeModel_set_transform(mtree, transform);

    Parameters_removeAll(ps);
    Parameters_add(ps, N);
    Parameters* reparams = get_reparams(tree);
    Parameter* ratios = Parameters_at(reparams, 0);
    Parameter* rootHeight = Parameters_at(reparams, 1);
    Parameters_add_parameters(ps, reparams);
    Parameters_zero_grad(ps);

    model->gradient(model, ps);

    double trueHeightGradient3[3] = {-2.571428571428571, -5.142857142857142,
                                     -0.3571428571428571};
    mu_assert(fabs(N->grad[0] - trueThetaGradient) < 1.e-10,
              "dlogP/d.theta not matching with reparameterization");
    for (size_t i = 0; i < 2; i++) {
        mu_assert(
            fabs(ratios->grad[i] - trueHeightGradient3[i]) < 1.e-10,
            "test_constant: d.logP/d.ratios not matching with reparameterization");
    }
    mu_assert(fabs(rootHeight->grad[0] - trueHeightGradient3[2]) < 1.e-10,
              "test_constant: d.logP/d.root not matching with reparameterization");

    model->free(model);
    mtree->free(mtree);
    free_Parameter(N);
    free_Parameters(ps);
    return NULL;
}

char* test_constant_ratios() {
    test_constant_with_transform(TREE_TRANSFORM_RATIO);
    return NULL;
}

char* test_constant_proportions_naive() {
    test_constant_with_transform(TREE_TRANSFORM_RATIO_NAIVE);
    return NULL;
}

char* test_constant_proportions() {
    test_constant_with_transform(TREE_TRANSFORM_PROPORTION);
    return NULL;
}

char* test_constant_data() {
    size_t intervalCount = 7;
    double times[7] = {0, 0, 0, 0, 2, 4, 6};
    double value = 3;
    bool coalescent[7] = {false, false, false, false, true, true, true};
    Parameter* N = new_Parameter("", value, new_Constraint(0, INFINITY));
    Coalescent* coal =
        new_ConstantCoalescent_with_data(N, times, coalescent, intervalCount);

    Model* model = new_CoalescentModel2("", coal, NULL, NULL);

    double logP = model->logP(model);
    double logP2 = -13.2958368660;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching");

    // logP/full_logP must be stored in the model's lp field before returning
    mu_assert(model->lp == logP, "CoalescentModel logP not stored in lp");
    mu_assert(model->full_logP(model) == logP, "CoalescentModel full_logP not matching");
    mu_assert(model->lp == logP, "CoalescentModel full_logP not stored in lp");

    value = 7;
    Parameter_set_value(N, value);
    logP = model->logP(model);
    logP2 = -10.1234447329;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching after change");

    Parameters* ps = new_Parameters(1);
    Parameters_add(ps, N);

    model->gradient(model, ps);

    struct gsl_data data = {model, 0};
    double eps = 0.0001;
    double dlogPdp2 = gsl_coalescent_dlogPdx(&data, value, eps);
    mu_assert(fabs(N->grad[0] - dlogPdp2) < 0.0001, "dlogPdp not matching");

    model->free(model);
    free_Parameter(N);
    free_Parameters(ps);
    return NULL;
}

/*char* test_constant_clone() {
    Tree* tree = new_Tree("(((a:2,b:2):4,c:6):6,d:12);", false);
    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("letree", tree);

    double value = 3;
    Parameter* N = new_Parameter("", value, new_Constraint(0, INFINITY));
    Coalescent* coal = new_ConstantCoalescent(tree, N);

    Model* model = new_CoalescentModel("", coal, mtree);

    Hashtable* hash2 = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash2, false);
    hashtable_set_value_ownership(hash2, false);

    Model* clone = model->clone(model, hash2);

    model->free(model);
    mtree->free(mtree);
    free_Parameter(N);
    clone->free(clone);
    free_Hashtable(hash2);
    return NULL;
}*/

char* test_piecewise_linear() {
    double dates[4] = {0.0, 0.0, 0.0, 0.0};
    char* taxa[4] = {"a", "b", "c", "d"};
    Model* mtree =
        new_TimeTreeModel_from_newick("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates);
    Tree* tree = mtree->obj;

    Parameters* ps = new_Parameters(5);
    Node** nodes = Tree_nodes(tree);
    for (size_t i = 0; i < Tree_node_count(tree); i++) {
        if (!Node_isleaf(nodes[i])) {
            Parameters_add(ps, nodes[i]->height);
        }
    }
    double thetas[5] = {3., 10., 4., 2., 3.};
    Parameter* popSizes = new_Parameter2("", thetas, 5, new_Constraint(0, INFINITY));
    Parameters_add(ps, popSizes);
    double cutoff = 10;
    int grid = 5;
    Coalescent* coal = new_PiecewiseLinearGridCoalescent(tree, popSizes, grid, cutoff);
    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -11.08185677776700117647;
    mu_assert(fabs(logP - logP2) < 0.00001, "test_piecewise_linear: logP not matching");

    model->gradient(model, ps);
    double true_gradient_thetas[5] = {0.32063498962941356, 0.11153798261181064,
                                      0.17750252451894566, 0.33669080273686075,
                                      0.06921832582596682};
    double true_gradient_heights[3] = {-0.6744186046511627, -0.375,
                                       -0.3333333333333333};

    for (size_t i = 0; i < 5; i++) {
        mu_assert(fabs(popSizes->grad[i] - true_gradient_thetas[i]) < 0.00001,
                  "test_piecewise_linear: d.logP/d.theta not matching");
    }

    for (size_t i = 0; i < 3; i++) {
        mu_assert(
            fabs(Parameters_at(ps, i)->grad[0] - true_gradient_heights[i]) < 0.00001,
            "test_piecewise_linear: d.logP/d.time not matching");
    }

    model->free(model);
    mtree->free(mtree);
    free_Parameters(ps);
    return NULL;
}

/*char* test_piecewise_linear2() {
    double dates[4] = {0.0, 0.0, 0.0, 0.0};
    char* taxa[4] = {"a", "b", "c", "d"};
    Model* mtree =
        new_TimeTreeModel_from_newick("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates);
    Tree* tree = mtree->obj;

    Parameters* ps = new_Parameters(3);
    double thetas[5] = {3., 10., 4.};
    Parameter* popSizes = new_Parameter2("", thetas, 3, new_Constraint(0, INFINITY));
    Parameters_add(ps, popSizes);
    double cutoff = 16;
    int grid = 3;
    Coalescent* coal =
        new_PiecewiseLinearGridCoalescent(tree, ps, grid, cutoff);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -11.353573548873296;
    mu_assert(fabs(logP - logP2) < 0.00001,
              "test_piecewise_linear2: logP not matching");

    // size_t gradient_size = Coalescent_initialize_gradient(
    //     model, GRADIENT_FLAG_COALESCENT_THETA | GRADIENT_FLAG_TREE_HEIGHTS);
    double true_gradient_thetas[3] = {0.7349363530597866, 0.07563028216195936,
                                      -0.055451574843527904};
    double true_gradient_heights[3] = {-0.8157894736842107, -0.34848484848484845,
                                       -0.0357142857142857};

    // double* gradient = Coalescent_gradient(model);
    model->gradient(model, ps);
    for (int i = 0; i < 3; i++) {
        printf("%f %f\n", gradient[i], true_gradient_thetas[i]);
        mu_assert(fabs(gradient[i] - true_gradient_thetas[i]) < 0.00001,
                  "test_piecewise_linear2: d.logP/d.time not matching");
    }

    for (int i = 0; i < 3; i++) {
        mu_assert(fabs(gradient[i + 3] - true_gradient_heights[i]) < 0.00001,
                  "test_piecewise_linear2: d.logP/d.theta not matching");
    }

    model->free(model);
    mtree->free(mtree);
    free_Parameters(ps);
    return NULL;
}*/

// Generic finite-difference check of a coalescent model's gradient wrt the
// reparameterized node-height ratios, the root height, and the model's
// population parameters (`extra`). Exercises trees with and without unknown-age
// leaves (`nratios` = tipCount-2 + number of unknown leaves).
static char* _fd_coal_model_gradient(Model* model, Model* mtree, Parameters* extra,
                                     size_t nratios, const char* label) {
    Tree* tree = mtree->obj;
    Parameters* reparams = get_reparams(tree);
    Parameter* ratios = Parameters_at(reparams, 0);
    Parameter* root = Parameters_at(reparams, 1);
    mu_assert(Parameter_size(ratios) == nratios, "unexpected number of ratios");

    Parameters* ps = new_Parameters(2);
    Parameters_add_parameters(ps, reparams);
    Parameters_add_parameters(ps, extra);

    size_t offset = nratios - (Tree_tip_count(tree) - 2);
    for (size_t i = 0; i < offset; i++) {
        mu_assert(Parameter_value_at(ratios, i) > 1.e-3 &&
                      Parameter_value_at(ratios, i) < 1.0 - 1.e-3,
                  "unknown-leaf ratio initialized at the boundary");
    }

    double lp0 = model->logP(model);
    mu_assert(!isnan(lp0) && !isinf(lp0), "coalescent logP not finite");

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    double g_ratio[16];
    for (size_t i = 0; i < nratios; i++) g_ratio[i] = ratios->grad[i];
    double g_root = root->grad[0];
    double g_extra[32];
    size_t nextra = 0;
    for (size_t c = 0; c < Parameters_count(extra); c++) {
        Parameter* p = Parameters_at(extra, c);
        for (size_t j = 0; j < Parameter_size(p); j++) g_extra[nextra++] = p->grad[j];
    }

    double h = 1.e-6;
    double worst = 0.0;

    for (size_t i = 0; i < nratios; i++) {
        double v0 = Parameter_value_at(ratios, i);
        Parameter_set_value_at(ratios, v0 + h, i);
        double lp = model->logP(model);
        Parameter_set_value_at(ratios, v0 - h, i);
        double lm = model->logP(model);
        Parameter_set_value_at(ratios, v0, i);
        double fd = (lp - lm) / (2.0 * h);
        double err = fabs(fd - g_ratio[i]) / (1.0 + fabs(g_ratio[i]));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-4, "coalescent ratio gradient does not match finite difference");
    }

    {
        double v0 = Parameter_value(root);
        double hh = h * (1.0 + fabs(v0));
        Parameter_set_value(root, v0 + hh);
        double lp = model->logP(model);
        Parameter_set_value(root, v0 - hh);
        double lm = model->logP(model);
        Parameter_set_value(root, v0);
        double fd = (lp - lm) / (2.0 * hh);
        double err = fabs(fd - g_root) / (1.0 + fabs(g_root));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-4, "coalescent root-height gradient does not match finite difference");
    }

    size_t idx = 0;
    for (size_t c = 0; c < Parameters_count(extra); c++) {
        Parameter* p = Parameters_at(extra, c);
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double v0 = Parameter_value_at(p, j);
            double hh = h * (1.0 + fabs(v0));
            Parameter_set_value_at(p, v0 + hh, j);
            double lp = model->logP(model);
            Parameter_set_value_at(p, v0 - hh, j);
            double lm = model->logP(model);
            Parameter_set_value_at(p, v0, j);
            double fd = (lp - lm) / (2.0 * hh);
            double err = fabs(fd - g_extra[idx]) / (1.0 + fabs(g_extra[idx]));
            if (err > worst) worst = err;
            idx++;
            mu_assert(err < 1.e-4, "coalescent population gradient does not match finite difference");
        }
    }

    printf("  [%s nratios=%zu] worst relative FD error = %.3e\n", label, nratios, worst);
    free_Parameters(ps);
    return NULL;
}

// Exponential coalescent gradient (n0, growth rate, ratios, root height), checked
// against finite differences with and without an unknown-age leaf.
// caterpillar tree (distinct internal heights 1.3, 3, 5.5, 8.5) to avoid the
// tied-height degeneracy that makes finite differences unreliable
#define COAL_FD_NEWICK "((((a:1.3,b:1.3):1.7,c:3):2.5,d:5.5):3,e:8.5);"

char* test_exponential_proportions_gradient_fd() {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    double dates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    Model* mtree = new_TimeTreeModel_from_newick(COAL_FD_NEWICK, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;
    Parameters* exp = new_Parameters(2);
    Parameters_add(exp, new_Parameter("n0", 4.0, new_Constraint(0, INFINITY)));
    Parameters_add(exp, new_Parameter("growth", 0.1, new_Constraint(-INFINITY, INFINITY)));
    Coalescent* coal = new_ExponentialCoalescent(tree, exp);
    Model* model = new_CoalescentModel("coalescent", coal, mtree);
    char* r = _fd_coal_model_gradient(model, mtree, exp, 3, "exponential");
    model->free(model);
    mtree->free(mtree);
    free_Parameters(exp);
    return r;
}

char* test_exponential_proportions_leaf_gradient_fd() {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    double dates[5] = {0.0, 0.0, -1.0, 0.0, 0.0};
    Model* mtree = new_TimeTreeModel_from_newick(COAL_FD_NEWICK, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;
    Parameters* exp = new_Parameters(2);
    Parameters_add(exp, new_Parameter("n0", 4.0, new_Constraint(0, INFINITY)));
    Parameters_add(exp, new_Parameter("growth", 0.1, new_Constraint(-INFINITY, INFINITY)));
    Coalescent* coal = new_ExponentialCoalescent(tree, exp);
    Model* model = new_CoalescentModel("coalescent", coal, mtree);
    char* r = _fd_coal_model_gradient(model, mtree, exp, 4, "exponential-leaf");
    model->free(model);
    mtree->free(mtree);
    free_Parameters(exp);
    return r;
}

// (COAL_FD_NEWICK is defined above; its heights also avoid the skygrid grid
// points at cutoff 9, grid 4 -> 2.25, 4.5, 6.75, 9.)

// Skygrid grid points sit at cutoff*i/(grid-1). For grid=4 the lines are at
// cutoff/3, 2*cutoff/3 and cutoff; the cutoffs below (10 and 7) keep them off the
// caterpillar coalescent heights (1.3, 3, 5.5, 8.5) so finite differences are
// valid. cutoff=10 is OLDER than the root (8.5) -> a trailing empty grid segment;
// cutoff=7 is YOUNGER than the root -> the oldest popSize covers the tail.
static char* _fd_skygrid(double cutoff, double* dates, size_t nratios, const char* label) {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    Model* mtree = new_TimeTreeModel_from_newick(COAL_FD_NEWICK, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;
    double thetas[4] = {3., 10., 4., 2.};
    Parameter* popSizes = new_Parameter2("popSizes", thetas, 4, new_Constraint(0, INFINITY));
    Parameters* extra = new_Parameters(1);
    Parameters_add(extra, popSizes);
    Coalescent* coal = new_GridCoalescent(tree, popSizes, 4, cutoff);
    Model* model = new_CoalescentModel("coalescent", coal, mtree);
    char* r = _fd_coal_model_gradient(model, mtree, extra, nratios, label);
    model->free(model);
    mtree->free(mtree);
    free_Parameters(extra);
    return r;
}

char* test_skygrid_cutoff_old_gradient_fd() {
    double dates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    return _fd_skygrid(10.0, dates, 3, "skygrid-cutoff>root");
}

char* test_skygrid_cutoff_old_leaf_gradient_fd() {
    double dates[5] = {0.0, 0.0, -1.0, 0.0, 0.0};
    return _fd_skygrid(10.0, dates, 4, "skygrid-cutoff>root-leaf");
}

char* test_skygrid_cutoff_young_gradient_fd() {
    double dates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    return _fd_skygrid(7.0, dates, 3, "skygrid-cutoff<root");
}

char* test_skygrid_cutoff_young_leaf_gradient_fd() {
    double dates[5] = {0.0, 0.0, -1.0, 0.0, 0.0};
    return _fd_skygrid(7.0, dates, 4, "skygrid-cutoff<root-leaf");
}

// Piecewise-linear grid coalescent: shares height_gradient_from_interval_gradient
// (the grid-point OOB fix), so exercise it the same way as skygrid -- grid points
// off the coalescent heights, cutoff both older and younger than the root.
static char* _fd_piecewise(double cutoff, double* dates, size_t nratios,
                           const char* label) {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    Model* mtree = new_TimeTreeModel_from_newick(COAL_FD_NEWICK, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;
    double thetas[4] = {3., 10., 4., 2.};
    Parameter* popSizes = new_Parameter2("popSizes", thetas, 4, new_Constraint(0, INFINITY));
    Parameters* extra = new_Parameters(1);
    Parameters_add(extra, popSizes);
    Coalescent* coal = new_PiecewiseLinearGridCoalescent(tree, popSizes, 4, cutoff);
    Model* model = new_CoalescentModel("coalescent", coal, mtree);
    char* r = _fd_coal_model_gradient(model, mtree, extra, nratios, label);
    model->free(model);
    mtree->free(mtree);
    free_Parameters(extra);
    return r;
}

char* test_piecewise_cutoff_old_gradient_fd() {
    double dates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    return _fd_piecewise(10.0, dates, 3, "piecewise-cutoff>root");
}

char* test_piecewise_cutoff_old_leaf_gradient_fd() {
    double dates[5] = {0.0, 0.0, -1.0, 0.0, 0.0};
    return _fd_piecewise(10.0, dates, 4, "piecewise-cutoff>root-leaf");
}

char* test_piecewise_cutoff_young_gradient_fd() {
    double dates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    return _fd_piecewise(7.0, dates, 3, "piecewise-cutoff<root");
}

char* test_piecewise_cutoff_young_leaf_gradient_fd() {
    double dates[5] = {0.0, 0.0, -1.0, 0.0, 0.0};
    return _fd_piecewise(7.0, dates, 4, "piecewise-cutoff<root-leaf");
}

char* test_skyride_proportions_gradient_fd() {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    double dates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    Model* mtree = new_TimeTreeModel_from_newick(COAL_FD_NEWICK, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;
    double thetas[4] = {3., 10., 4., 2.};
    Parameter* popSizes = new_Parameter2("popSizes", thetas, 4, new_Constraint(0, INFINITY));
    Parameters* extra = new_Parameters(1);
    Parameters_add(extra, popSizes);
    Coalescent* coal = new_SkyrideCoalescent(tree, popSizes);
    Model* model = new_CoalescentModel("coalescent", coal, mtree);
    char* r = _fd_coal_model_gradient(model, mtree, extra, 3, "skyride");
    model->free(model);
    mtree->free(mtree);
    free_Parameters(extra);
    return r;
}

char* test_skyride_proportions_leaf_gradient_fd() {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    double dates[5] = {0.0, 0.0, -1.0, 0.0, 0.0};
    Model* mtree = new_TimeTreeModel_from_newick(COAL_FD_NEWICK, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;
    double thetas[4] = {3., 10., 4., 2.};
    Parameter* popSizes = new_Parameter2("popSizes", thetas, 4, new_Constraint(0, INFINITY));
    Parameters* extra = new_Parameters(1);
    Parameters_add(extra, popSizes);
    Coalescent* coal = new_SkyrideCoalescent(tree, popSizes);
    Model* model = new_CoalescentModel("coalescent", coal, mtree);
    char* r = _fd_coal_model_gradient(model, mtree, extra, 4, "skyride-leaf");
    model->free(model);
    mtree->free(mtree);
    free_Parameters(extra);
    return r;
}

// Verify the constant-coalescent gradient (wrt reparameterized node-height
// ratios, root height, and theta) against central finite differences, using the
// PROPORTION transform. `nratios` is tipCount-2 plus the number of unknown-age
// leaves, so the same routine exercises trees with and without unknown leaves.
static char* _fd_coalescent_gradient(const char* newick, char** taxa, double* dates,
                                     size_t nratios) {
    Model* mtree = new_TimeTreeModel_from_newick(newick, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;

    Parameter* N = new_Parameter("theta", 4.0, new_Constraint(0, INFINITY));
    Coalescent* coal = new_ConstantCoalescent(tree, N);
    Model* model = new_CoalescentModel("coalescent", coal, mtree);

    Parameters* reparams = get_reparams(tree);
    Parameter* ratios = Parameters_at(reparams, 0);
    Parameter* root = Parameters_at(reparams, 1);
    mu_assert(Parameter_size(ratios) == nratios, "unexpected number of ratios");

    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, N);
    Parameters_add_parameters(ps, reparams);

    // unknown-leaf ratios occupy the first `offset` slots and should start
    // at an interior point (not pinned at the lower bound)
    size_t offset = nratios - (Tree_tip_count(tree) - 2);
    for (size_t i = 0; i < offset; i++) {
        mu_assert(Parameter_value_at(ratios, i) > 1.e-3 &&
                      Parameter_value_at(ratios, i) < 1.0 - 1.e-3,
                  "unknown-leaf ratio initialized at the boundary");
    }

    double lp0 = model->logP(model);
    mu_assert(!isnan(lp0) && !isinf(lp0), "coalescent logP not finite");

    Parameters_zero_grad(ps);
    model->gradient(model, ps);
    double g_N = N->grad[0];
    double g_root = root->grad[0];
    double g_ratio[16];
    for (size_t i = 0; i < nratios; i++) g_ratio[i] = ratios->grad[i];

    double h = 1.e-6;
    double worst = 0.0;

    for (size_t i = 0; i < nratios; i++) {
        double v0 = Parameter_value_at(ratios, i);
        Parameter_set_value_at(ratios, v0 + h, i);
        double lp = model->logP(model);
        Parameter_set_value_at(ratios, v0 - h, i);
        double lm = model->logP(model);
        Parameter_set_value_at(ratios, v0, i);
        double fd = (lp - lm) / (2.0 * h);
        double err = fabs(fd - g_ratio[i]) / (1.0 + fabs(g_ratio[i]));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-4,
                  "coalescent ratio gradient does not match finite difference");
    }

    {
        double v0 = Parameter_value(root);
        double hh = h * (1.0 + fabs(v0));
        Parameter_set_value(root, v0 + hh);
        double lp = model->logP(model);
        Parameter_set_value(root, v0 - hh);
        double lm = model->logP(model);
        Parameter_set_value(root, v0);
        double fd = (lp - lm) / (2.0 * hh);
        double err = fabs(fd - g_root) / (1.0 + fabs(g_root));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-4,
                  "coalescent root-height gradient does not match finite difference");
    }

    {
        double v0 = Parameter_value(N);
        double hh = h * (1.0 + fabs(v0));
        Parameter_set_value(N, v0 + hh);
        double lp = model->logP(model);
        Parameter_set_value(N, v0 - hh);
        double lm = model->logP(model);
        Parameter_set_value(N, v0);
        double fd = (lp - lm) / (2.0 * hh);
        double err = fabs(fd - g_N) / (1.0 + fabs(g_N));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-4,
                  "coalescent theta gradient does not match finite difference");
    }

    printf("  [coalescent nratios=%zu] worst relative FD error = %.3e\n", nratios, worst);

    model->free(model);
    mtree->free(mtree);
    free_Parameter(N);
    free_Parameters(ps);
    return NULL;
}

// Verify the log-det-Jacobian of the height reparameterization and its gradient.
// The tree model's logP is the log-det-Jacobian; mtree->gradient accumulates its
// gradient wrt the ratios and root height. Checked against finite differences.
static char* _fd_logdet_jacobian(const char* newick, char** taxa, double* dates,
                                 size_t nratios) {
    Model* mtree = new_TimeTreeModel_from_newick(newick, taxa, dates);
    TreeModel_set_transform(mtree, TREE_TRANSFORM_PROPORTION);
    Tree* tree = mtree->obj;

    Parameters* reparams = get_reparams(tree);
    Parameter* ratios = Parameters_at(reparams, 0);
    Parameter* root = Parameters_at(reparams, 1);
    mu_assert(Parameter_size(ratios) == nratios, "unexpected number of ratios");

    size_t offset = nratios - (Tree_tip_count(tree) - 2);
    for (size_t i = 0; i < offset; i++) {
        mu_assert(Parameter_value_at(ratios, i) > 1.e-3 &&
                      Parameter_value_at(ratios, i) < 1.0 - 1.e-3,
                  "unknown-leaf ratio initialized at the boundary");
    }

    double lj0 = mtree->logP(mtree);
    mu_assert(!isnan(lj0) && !isinf(lj0), "log-det-Jacobian not finite");

    Parameters_zero_grad(reparams);
    mtree->gradient(mtree, reparams);
    double g_root = root->grad[0];
    double g_ratio[16];
    for (size_t i = 0; i < nratios; i++) g_ratio[i] = ratios->grad[i];

    double h = 1.e-6;
    double worst = 0.0;

    for (size_t i = 0; i < nratios; i++) {
        double v0 = Parameter_value_at(ratios, i);
        Parameter_set_value_at(ratios, v0 + h, i);
        double lp = mtree->logP(mtree);
        Parameter_set_value_at(ratios, v0 - h, i);
        double lm = mtree->logP(mtree);
        Parameter_set_value_at(ratios, v0, i);
        double fd = (lp - lm) / (2.0 * h);
        double err = fabs(fd - g_ratio[i]) / (1.0 + fabs(g_ratio[i]));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-4,
                  "log-det-Jacobian ratio gradient does not match finite difference");
    }

    {
        double v0 = Parameter_value(root);
        double hh = h * (1.0 + fabs(v0));
        Parameter_set_value(root, v0 + hh);
        double lp = mtree->logP(mtree);
        Parameter_set_value(root, v0 - hh);
        double lm = mtree->logP(mtree);
        Parameter_set_value(root, v0);
        double fd = (lp - lm) / (2.0 * hh);
        double err = fabs(fd - g_root) / (1.0 + fabs(g_root));
        if (err > worst) worst = err;
        mu_assert(err < 1.e-4,
                  "log-det-Jacobian root-height gradient does not match finite difference");
    }

    printf("  [logdet-jacobian nratios=%zu] worst relative FD error = %.3e\n", nratios, worst);

    mtree->free(mtree);
    return NULL;
}

char* test_logdet_jacobian_gradient_fd() {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    double dates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    return _fd_logdet_jacobian("((a:4,b:4):4,((c:2,d:2):2,e:4):4);", taxa, dates, 3);
}

char* test_logdet_jacobian_leaf_gradient_fd() {
    char* taxa[5] = {"a", "b", "c", "d", "e"};
    // taxon c is unknown; its parent is an internal node, so c's Jacobian term
    // log(parent_height) contributes to that internal node's ratio gradient
    double dates[5] = {0.0, 0.0, -1.0, 0.0, 0.0};
    return _fd_logdet_jacobian("((a:4,b:4):4,((c:2,d:2):2,e:4):4);", taxa, dates, 4);
}

char* test_constant_proportions_gradient_fd() {
    char* taxa[4] = {"a", "b", "c", "d"};
    double dates[4] = {0.0, 1.0, 2.0, 4.0};
    return _fd_coalescent_gradient("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates, 2);
}

char* test_constant_proportions_leaf_gradient_fd() {
    char* taxa[4] = {"a", "b", "c", "d"};
    // taxon d has an unknown age: it becomes a free (reparameterized) parameter
    double dates[4] = {0.0, 1.0, 2.0, -1.0};
    return _fd_coalescent_gradient("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates, 3);
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_constant_proportions_naive);
    mu_run_test(test_constant_proportions);
    mu_run_test(test_constant_ratios);
    mu_run_test(test_constant_data);
    // // mu_run_test(test_constant_clone);
    mu_run_test(test_skyride);
    mu_run_test(test_skygrid);
    mu_run_test(test_piecewise_linear);
    // mu_run_test(test_piecewise_linear2);
    mu_run_test(test_constant_proportions_gradient_fd);
    mu_run_test(test_constant_proportions_leaf_gradient_fd);
    mu_run_test(test_logdet_jacobian_gradient_fd);
    mu_run_test(test_logdet_jacobian_leaf_gradient_fd);
    mu_run_test(test_skyride_proportions_gradient_fd);
    mu_run_test(test_skyride_proportions_leaf_gradient_fd);
    mu_run_test(test_exponential_proportions_gradient_fd);
    mu_run_test(test_exponential_proportions_leaf_gradient_fd);
    mu_run_test(test_skygrid_cutoff_old_gradient_fd);
    mu_run_test(test_skygrid_cutoff_old_leaf_gradient_fd);
    mu_run_test(test_skygrid_cutoff_young_gradient_fd);
    mu_run_test(test_skygrid_cutoff_young_leaf_gradient_fd);
    mu_run_test(test_piecewise_cutoff_old_gradient_fd);
    mu_run_test(test_piecewise_cutoff_old_leaf_gradient_fd);
    mu_run_test(test_piecewise_cutoff_young_gradient_fd);
    mu_run_test(test_piecewise_cutoff_young_leaf_gradient_fd);
    return NULL;
}

RUN_TESTS(all_tests);
