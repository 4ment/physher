//
//  test_coalescent.c
//  physher
//  Created by Mathieu Fourment on 25/04/2020.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "minunit.h"

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "phyc/demographicmodels.h"
#include "phyc/parameters.h"
#include "phyc/tree.h"

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
        new_TreeModel_from_newick("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates);
    Tree* tree = mtree->obj;
    Tree_update_heights(tree);

    Parameters* ps = new_Parameters(3);
    double thetas[3] = {3., 10., 4.};
    for (int i = 0; i < 3; i++) {
        Parameter* p = new_Parameter("", log(thetas[i]), new_Constraint(0, INFINITY));
        Parameters_move(ps, p);
        p->id = i;
    }
    Coalescent* coal = new_SkyrideCoalescent(tree, ps, COALESCENT_THETA_LOG);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -11.48749174278204;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching");

    Parameters* reparams = get_reparams(tree);
    for (int i = 0; i < Parameters_count(reparams); i++) {
        Parameters_add(ps, Parameters_at(reparams, i));
    }

    Parameters_set_value(ps, 0, log(thetas[0]));
    model->prepare_gradient(model, ps);
    double true_gradient[6] = {3., 0.2, 0.5, -10.2, -7.4, -0.5583333};
    for (int i = 0; i < 3; i++) {
        double dlogPdp = model->dlogP(model, Parameters_at(ps, i));
        mu_assert(fabs(dlogPdp - true_gradient[i]) < 0.00001,
                  "d.logP/d.theta not matching");
    }
    // for(int i = 3; i < 6; i++){
    // 	double dlogPdp = model->dlogP(model, Parameters_at(ps, i));
    // 	printf("%f\n", dlogPdp);
    // 	mu_assert(fabs(dlogPdp - true_gradient[i]) < 0.00001, "d.logP/d.t not
    // matching");
    // }
    struct gsl_data data = {model, 0};
    // 	for(int i = 0; i < 3; i++){
    // 		data.index = i;
    // 		double dlogPdp2 = gsl_coalescent_dlogPdx(&data, log(thetas[i]), 0.0001);
    // 	}
    for (int i = 0; i < 3; i++) {
        data.index = i + 3;
        double dlogPdp2 =
            gsl_coalescent_dlogPdx(&data, Parameters_value(reparams, i), 0.0001);
    }
    model->free(model);
    mtree->free(mtree);
    free_Parameters(ps);
    return NULL;
}

char* test_skygrid() {
    double dates[4] = {0.0, 0.0, 0.0, 0.0};
    char* taxa[4] = {"a", "b", "c", "d"};
    Model* mtree =
        new_TreeModel_from_newick("(((a:2,b:2):4,c:6):6,d:12);", taxa, dates);
    Tree* tree = mtree->obj;

    Parameters* ps = new_Parameters(5);
    double thetas[5] = {3., 10., 4., 2., 3.};
    for (int i = 0; i < 5; i++) {
        Parameter* p = new_Parameter("", log(thetas[i]), new_Constraint(0, INFINITY));
        Parameters_move(ps, p);
        p->id = i;
    }
    double cutoff = 10;
    int grid = 5;
    Coalescent* coal = new_GridCoalescent(tree, ps, grid, cutoff, COALESCENT_THETA_LOG);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -11.8751856;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching");

    Parameters* reparams = get_reparams(tree);
    for (int i = 0; i < Parameters_count(reparams); i++) {
        Parameters_add(ps, Parameters_at(reparams, i));
    }

    Parameters_set_value(ps, 0, log(thetas[0]));
    model->prepare_gradient(model, ps);
    double true_gradient[8] = {3.5, 0.75, 0.1250, 1.25, -0.333333, -6.0, -10.0, -0.75};
    for (int i = 0; i < 5; i++) {
        double dlogPdp = model->dlogP(model, Parameters_at(ps, i));
        mu_assert(fabs(dlogPdp - true_gradient[i]) < 0.00001,
                  "d.logP/d.theta not matching");
    }
    // for(int i = 5; i < 8; i++){
    // 	double dlogPdp = model->dlogP(model, Parameters_at(ps, i));
    // 	mu_assert(fabs(dlogPdp - true_gradient[i]) < 0.00001, "d.logP/d.t not
    // matching");
    // }
    struct gsl_data data = {model, 0};
    for (int i = 0; i < 3; i++) {
        data.index = i + 5;
        double dlogPdp2 =
            gsl_coalescent_dlogPdx(&data, Parameters_value(reparams, i), 0.0001);
    }
    model->free(model);
    mtree->free(mtree);
    free_Parameters(ps);
    return NULL;
}

char* test_constant() {
    Tree* tree = new_Tree("(((a:2,b:2):4,c:6):6,d:12);", false);
    init_heights_from_distances(tree);

    Model* mtree = new_TreeModel("letree", tree);

    double value = 3;
    Parameter* N = new_Parameter("", value, new_Constraint(0, INFINITY));
    Coalescent* coal = new_ConstantCoalescent(tree, N);

    Model* model = new_CoalescentModel("", coal, mtree);

    double logP = model->logP(model);
    double logP2 = -13.2958368660;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching");

    value = 7;
    Parameter_set_value(N, value);
    logP = model->logP(model);
    logP2 = -10.1234447329;
    mu_assert(fabs(logP - logP2) < 0.00001, "logP not matching after change");

    double eps = 0.0001;
    Parameters* ps = new_Parameters(1);
    Parameters_add(ps, N);
    Parameter_set_value(N, value);
    model->prepare_gradient(model, ps);
    double dlogPdp = model->dlogP(model, N);
    struct gsl_data data = {model, 0};
    double dlogPdp2 = gsl_coalescent_dlogPdx(&data, value, eps);

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
    Parameters_set_value(ps, 0, value);
    model->prepare_gradient(model, ps);
    double eps = 0.0001;
    double dlogPdp = model->dlogP(model, N);
    struct gsl_data data = {model, 0};
    double dlogPdp2 = gsl_coalescent_dlogPdx(&data, value, eps);
    mu_assert(fabs(dlogPdp - dlogPdp2) < 0.0001, "dlogPdp not matching");

    model->free(model);
    free_Parameter(N);
    free_Parameters(ps);
    return NULL;
}

char* test_constant_clone() {
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
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_constant);
    mu_run_test(test_constant_data);
    mu_run_test(test_constant_clone);
    mu_run_test(test_skyride);
    mu_run_test(test_skygrid);
    return NULL;
}

RUN_TESTS(all_tests);
