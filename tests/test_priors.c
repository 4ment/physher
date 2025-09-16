#include <math.h>

#include "minunit.h"
#include "phyc/ctmcscale.h"
#include "phyc/distoneonx.h"
#include "phyc/gradient.h"
#include "phyc/matrix.h"
#include "phyc/parameters.h"

char* test_oneonx() {
    double value = 2.0;
    Parameter* x = new_Parameter("", value, new_Constraint(0, INFINITY));
    Parameters* parameters = new_Parameters(1);
    Parameters_move(parameters, x);

    DistributionModel* dm = new_OneOnXDistributionModel(parameters);
    Model* prior = new_DistributionModel2("id", dm);

    double trueLogP = -log(value);
    double trueGradient = -1.0 / value;

    double logP = prior->logP(prior);
    mu_assert(fabs(logP - trueLogP) < .000001, "logP not matching");

    prior->gradient(prior, parameters);
    mu_assert(fabs(x->grad[0] - trueGradient) < .000001, "gradient not matching");

    Parameters_zero_grad(parameters);

    value = 0.1;
    trueLogP = -log(value);
    trueGradient = -1.0 / value;

    Parameter_set_value(x, value);
    logP = prior->logP(prior);
    mu_assert(fabs(logP - trueLogP) < .000001, "logP not matching with 0.1");

    // derivative wrt to x
    prior->gradient(prior, parameters);
    mu_assert(fabs(x->grad[0] - trueGradient) < .000001,
              "gradient not matching with 0.1");

    // use reparameterization
    value = 0.01;
    trueLogP = -log(value);
    trueGradient = -1.0 / value;

    Parameter* parameter =
        new_Parameter("a", log(value), new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_Transform_with_parameter(NULL, 0, INFINITY, parameter);
    x->transform = transform;
    // Parameters_add(transform->parameter->listeners->parameters, x);
    transform->parameter->listeners->add_parameter(transform->parameter->listeners, x);
    Parameter_fire(parameter, -1);  // prior needs to be updated

    logP = prior->logP(prior);
    mu_assert(fabs(logP - trueLogP) < .000001, "logP not matching with transform");

    // derivative wrt to constrained parameter
    Parameters_zero_grad(parameters);
    prior->gradient(prior, parameters);
    mu_assert(fabs(x->grad[0] - trueGradient) < .000001,
              "gradient not matching with transform");

    // derivative wrt unconstrained parameter
    Parameters_removeAll(parameters);
    Parameters_add(parameters, parameter);
    prior->gradient(prior, parameters);
    mu_assert(fabs(parameter->grad[0] - trueGradient * value) < .000001,
              "gradient not matching wrt to unconstrained parameter");

    free_Parameters(parameters);
    prior->free(prior);

    return NULL;
}

char* test_ctmcscale() {
    Parameters* parameters = new_Parameters(1);
    double value = 0.001;
    Parameter* x = new_Parameter("", value, new_Constraint(0, INFINITY));
    Parameters_add(parameters, x);

    double dates[5] = {0.0, 1.0, 2.0, 3.0, 12.0};
    char* taxa[5] = {"A_0", "B_1", "C_2", "D_3", "E_12"};
    Model* mtree = new_TimeTreeModel_from_newick(
        "((((A_0:1.5,B_1:0.5):2.5,C_2:2):2,D_3:3):10,E_12:4);", taxa, dates);
    Tree* tree = mtree->obj;
    Parameters* nodeHeights = new_Parameters(4);
    for (size_t i = 0; i < Tree_node_count(tree); i++) {
        Node* node = Tree_node(tree, i);
        if (!Node_isleaf(node)) {
            Parameters_add(nodeHeights, node->height);
        }
    }
    Parameters_add_parameters(parameters, nodeHeights);

    DistributionModel* dm = new_CTMCScale_with_parameters(parameters, tree);
    Model* prior = new_CTMCScaleModel("ctmcscale", dm, mtree);

    double trueLogP = 4.475351922659342;
    double trueGradient = -525.5;
    double trueHeightGradient[4] = {0.018607843667268753, 0.018607843667268753,
                                    0.018607843667268753, 0.037215687334537506};

    double logP = prior->logP(prior);
    mu_assert(fabs(logP - trueLogP) < 1.e-8, "logP not matching");

    prior->gradient(prior, parameters);
    mu_assert(fabs(x->grad[0] - trueGradient) < 1.e-8, "x gradient not matching");
    for (size_t i = 0; i < 4; i++) {
        Parameter* height = Parameters_at(nodeHeights, i);
        mu_assert(fabs(height->grad[0] - trueHeightGradient[i]) < .000001,
                  "height gradient not matching");
    }

    // test CTMCModel_gradient used by C++ wrapper
    double* gradient = dvector(5);
    CTMCModel_gradient(prior, GRADIENT_FLAG_TREE_HEIGHTS | GRADIENT_FLAG_CLOCK_RATE,
                       gradient);
    for (size_t i = 0; i < Parameters_count(nodeHeights); i++) {
        mu_assert(fabs(gradient[i] - trueHeightGradient[i]) < 1.e-8,
                  "gradient not matching wrt to heights with CTMCModel_gradient");
    }
    mu_assert(fabs(gradient[4] - trueGradient) < 1.e-8,
              "gradient not matching wrt to x with CTMCModel_gradient");
    free(gradient);

    // use reparameterization
    value = 0.01;
    trueLogP = 3.094559376151536;
    trueGradient = -75.5;

    Parameter* parameter =
        new_Parameter("a", log(value), new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_Transform_with_parameter(NULL, 0, INFINITY, parameter);
    x->transform = transform;
    // Parameters_add(transform->parameter->listeners->parameters, x);
    transform->parameter->listeners->add_parameter(transform->parameter->listeners, x);
    Parameter_fire(parameter, -1);  // prior needs to be updated

    logP = prior->logP(prior);
    mu_assert(fabs(logP - trueLogP) < 1.e-8, "logP not matching with transform");

    // derivative wrt to constrained parameter
    Parameters_zero_grad(parameters);
    prior->gradient(prior, parameters);
    mu_assert(fabs(x->grad[0] - trueGradient) < 1.e-8,
              "gradient not matching with transform");

    // derivative wrt unconstrained parameter
    Parameters_zero_grad(parameters);
    Parameters_removeAll(parameters);
    Parameters_add(parameters, parameter);
    prior->gradient(prior, parameters);
    mu_assert(fabs(parameter->grad[0] - trueGradient * value) < 1.e-8,
              "gradient not matching wrt to unconstrained parameter");

    for (size_t i = 0; i < Parameters_count(nodeHeights); i++) {
        mu_assert(Parameters_at(nodeHeights, 0)->grad[0] == 0.0,
                  "gradient should be 0 for node heights");
    }

    free_Parameters(parameters);
    free_Parameters(nodeHeights);
    free_Parameter(x);
    prior->free(prior);

    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_oneonx);
    mu_run_test(test_ctmcscale);
    return NULL;
}

RUN_TESTS(all_tests);