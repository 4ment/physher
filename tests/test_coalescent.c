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
    TreeModel_set_transform(mtree, TREE_TRANSFORM_RATIO);
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
    TreeModel_set_transform(mtree, TREE_TRANSFORM_RATIO);

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

char* test_constant() {
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
    TreeModel_set_transform(mtree, TREE_TRANSFORM_RATIO);

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

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_constant);
    mu_run_test(test_constant_data);
    // // mu_run_test(test_constant_clone);
    mu_run_test(test_skyride);
    mu_run_test(test_skygrid);
    mu_run_test(test_piecewise_linear);
    // mu_run_test(test_piecewise_linear2);
    return NULL;
}

RUN_TESTS(all_tests);
