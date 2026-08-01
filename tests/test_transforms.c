// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <stdlib.h>

#include "minunit.h"
#include "phyc/parameters.h"
#include "phyc/transforms.h"

char* test_exp() {
    // clang-format off
    /*
    Checked with pytorch

    z = torch.distributions.ExpTransform().inv(torch.tensor([0.1, 0.2, 0.3], dtype=torch.double))
    def f(x):
        return torch.distributions.ExpTransform()(x)
    torch.autograd.functional.jacobian(f, z).t().tolist()

    z.requires_grad = True
    x = torch.distributions.ExpTransform()(z)
    x.backward(torch.tensor([1.,2.,3.]))
    z.grad.tolist()

    x = torch.distributions.ExpTransform()(z)
    jac = torch.distributions.ExpTransform().log_abs_det_jacobian(z, x)
    jac.backward(torch.ones(3)))
    z.grad.tolist()
    */
    // clang-format on
    double initValues[4] = {0.2, 0.1, 0.4};
    double values[4] = {0.1, 0.2, 0.3};
    double unConstrainedValues[3];
    unConstrainedValues[0] = log(values[0]);
    unConstrainedValues[1] = log(values[1]);
    unConstrainedValues[2] = log(values[2]);
    Parameter* parameter = new_Parameter2("a", unConstrainedValues, 3,
                                          new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_Transform_with_parameter(NULL, 0, INFINITY, parameter);
    Parameter* positive =
        new_Parameter2("positive", initValues, 3, new_Constraint(0.0, INFINITY));
    positive->transform = transform;
    // Parameters_add(transform->parameter->listeners->parameters, positive);
    transform->parameter->listeners->add_parameter(transform->parameter->listeners, positive);
    const double* values2 = Parameter_values(positive);
    for (size_t i = 0; i < Parameter_size(positive); i++) {
        mu_assert(fabs(values2[i] - exp(unConstrainedValues[i])) < 1.e-7,
                  "exp: constrained values not matching");
    }

    const double ingrad[3] = {1., 2., 3.};
    const double trueGrad[3] = {0.10000000000000002, 0.4, 0.8999999999999999};
    transform->backward(transform, ingrad);
    for (size_t i = 0; i < 3; i++) {
        printf("%f %f\n", parameter->grad[i], trueGrad[i]);
        // mu_assert(fabs(parameter->grad[i] - exp(unConstrainedValues[i])*ingrad[i])
        // < 1.e-7,
        mu_assert(fabs(parameter->grad[i] - trueGrad[i]) < 1.e-7,
                  "exp: gradient not matching");
    }

    double trueLogDetJacobians[3] = {-2.3025850929940455, -1.6094379124341003,
                                     -1.2039728043259361};
    double trueLogDetJacobian =
        trueLogDetJacobians[0] + trueLogDetJacobians[1] + trueLogDetJacobians[2];
    double logDetJacobian = transform->log_det_jacobian(transform);
    mu_assert(fabs(trueLogDetJacobian - logDetJacobian) < 1.e-7,
              "exp: log det Jacobian not matching");

    double trueLogDetJacobianGradient[3] = {1.0, 1.0, 1.0};
    Parameter_zero_grad(parameter);
    transform->gradient_log_det_jacobian(transform);

    for (size_t i = 0; i < 3; i++) {
        // printf("%f %f\n", parameter->grad[i], trueLogDetJacobianGradient[i]);
        mu_assert(fabs(parameter->grad[i] - trueLogDetJacobianGradient[i]) < 1.e-7,
                  "exp: log det Jacobian gradient not matching");
    }
    return NULL;
}
// The optimizer works on the unconstrained parameter, so the soft (f) bounds set
// on the constrained parameter must be pushed through the transform onto the
// unconstrained leaf. y = log(x - lower) is increasing, so lower maps to lower
// and upper maps to upper.
char* test_exp_fbounds() {
    double unConstrainedValues[3] = {log(0.1), log(0.2), log(0.3)};
    double initValues[3] = {0.1, 0.2, 0.3};

    // lower = 0: y = log(x)
    Parameter* leaf = new_Parameter2("a", unConstrainedValues, 3,
                                     new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_Transform_with_parameter(NULL, 0, INFINITY, leaf);
    Parameter* positive =
        new_Parameter2("positive", initValues, 3, new_Constraint(0.0, INFINITY));
    positive->transform = transform;
    transform->parameter->listeners->add_parameter(transform->parameter->listeners,
                                                   positive);

    Parameter_set_fupper(positive, 10.0);
    mu_assert(fabs(Parameter_fupper(positive) - 10.0) < 1.e-7,
              "exp bounds: fupper not set on constrained parameter");
    mu_assert(fabs(Parameter_fupper(leaf) - log(10.0)) < 1.e-7,
              "exp bounds: fupper not transformed onto unconstrained parameter");

    Parameter_set_flower(positive, 0.01);
    mu_assert(fabs(Parameter_flower(positive) - 0.01) < 1.e-7,
              "exp bounds: flower not set on constrained parameter");
    mu_assert(fabs(Parameter_flower(leaf) - log(0.01)) < 1.e-7,
              "exp bounds: flower not transformed onto unconstrained parameter");

    // the hard bounds are untouched by the f setters
    mu_assert(Parameter_lower(positive) == 0.0 && isinf(Parameter_upper(positive)),
              "exp bounds: hard bounds of constrained parameter modified");
    mu_assert(isinf(Parameter_lower(leaf)) && isinf(Parameter_upper(leaf)),
              "exp bounds: hard bounds of unconstrained parameter modified");

    // flower == lower is the edge of the support: log(0) = -infinity
    Parameter_set_flower(positive, 0.0);
    mu_assert(isinf(Parameter_flower(leaf)) && Parameter_flower(leaf) < 0,
              "exp bounds: flower at the bound is not -infinity");

    // shifted support lower = 2: y = log(x - 2)
    double shiftedValues[2] = {log(1.0), log(3.0)};
    double shiftedInit[2] = {3.0, 5.0};
    Parameter* shiftedLeaf =
        new_Parameter2("b", shiftedValues, 2, new_Constraint(-INFINITY, INFINITY));
    Transform* shiftedTransform =
        new_Transform_with_parameter(NULL, 2.0, INFINITY, shiftedLeaf);
    Parameter* shifted =
        new_Parameter2("shifted", shiftedInit, 2, new_Constraint(2.0, INFINITY));
    shifted->transform = shiftedTransform;
    shiftedTransform->parameter->listeners->add_parameter(
        shiftedTransform->parameter->listeners, shifted);

    Parameter_set_flower(shifted, 2.5);
    Parameter_set_fupper(shifted, 10.0);
    mu_assert(fabs(Parameter_flower(shiftedLeaf) - log(0.5)) < 1.e-7,
              "exp bounds: shifted flower not transformed onto unconstrained parameter");
    mu_assert(fabs(Parameter_fupper(shiftedLeaf) - log(8.0)) < 1.e-7,
              "exp bounds: shifted fupper not transformed onto unconstrained parameter");

    return NULL;
}

char* test_sigmoid() {
    // clang-format off
    /*
    Checked with pytorch

    z = torch.tensor([0.1, -1.0], requires_grad = True)
    x = torch.distributions.SigmoidTransform()(z)
    x.backward(torch.tensor([1., 1.]))
    z.grad.tolist()

    x = torch.distributions.SigmoidTransform()(z)
    jac = torch.distributions.SigmoidTransform().log_abs_det_jacobian(z, x)
    jac.backward(torch.ones(2))
    z.grad.tolist()
    */
    // clang-format on
    double unConstrainedValues[2] = {0.1, -1.0};
    double constrainedValues[2] = {0.5249791741371155, 0.2689414322376251};
    Parameter* parameter = new_Parameter2("a", unConstrainedValues, 2,
                                          new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_Transform_with_parameter(NULL, 0, 1.0, parameter);
    Parameter* constrained =
        new_Parameter2("constrained", unConstrainedValues, 2, new_Constraint(0.0, 1.0));
    constrained->transform = transform;
    // Parameters_add(transform->parameter->listeners->parameters, constrained);
    transform->parameter->listeners->add_parameter(transform->parameter->listeners, constrained);
    const double* values2 = Parameter_values(constrained);
    for (size_t i = 0; i < Parameter_size(constrained); i++) {
        mu_assert(fabs(values2[i] - constrainedValues[i]) < 1.e-7,
                  "sigmoid: constrained values not matching");
    }

    const double ingrad[3] = {1., 1.};
    double trueGradient[2] = {0.24937604367733002, 0.1966119408607483};
    transform->backward(transform, ingrad);
    for (size_t i = 0; i < 2; i++) {
        mu_assert(fabs(parameter->grad[i] - trueGradient[i]) < 1.e-7,
                  "sigmoid: gradient not matching");
    }

    double trueLogDetJacobians[2] = {-1.3887933492660522, -1.6265232563018799};
    double trueLogDetJacobian = trueLogDetJacobians[0] + trueLogDetJacobians[1];
    double logDetJacobian = transform->log_det_jacobian(transform);
    mu_assert(fabs(logDetJacobian - trueLogDetJacobian) < 1.e-7,
              "sigmoid: log det Jacobian not matching");

    double trueLogDetJacobianGradient[2] = {-0.04995834827423096, 0.46211716532707214};
    Parameter_zero_grad(parameter);
    transform->gradient_log_det_jacobian(transform);

    for (size_t i = 0; i < 2; i++) {
        mu_assert(fabs(parameter->grad[i] - trueLogDetJacobianGradient[i]) < 1.e-7,
                  "sigmoid: log det Jacobian gradient not matching");
    }
    return NULL;
}

char* test_simplex() {
    // clang-format off
    /*
    Checked with pytorch

    z = torch.distributions.StickBreakingTransform().inv(torch.tensor([0.1, 0.2, 0.3, 0.4], dtype=torch.double))
    def f(x):
        return torch.distributions.StickBreakingTransform()(x)
    torch.autograd.functional.jacobian(f, z).t().tolist()

    z.requires_grad = True
    x = torch.distributions.StickBreakingTransform()(z)
    x.backward(torch.tensor([1.,2,3,4]))
    z.grad.tolist()

    x = torch.distributions.StickBreakingTransform()(z)
    jac = torch.distributions.StickBreakingTransform().log_abs_det_jacobian(z, x)
    jac.backward()
    z.grad.tolist()
    */
    // clang-format on
    double initValues[4] = {0.2, 0.1, 0.4, 0.3};
    double values[4] = {0.1, 0.2, 0.3, 0.4};
    double unConstrainedValues[3] = {-1.0986122886681093, -0.5596157879354225,
                                     -0.2876820724517808};
    Parameter* parameter = new_Parameter2("a", unConstrainedValues, 3,
                                          new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_SimplexTransform_with_parameter(NULL, parameter);
    Parameter* simplex =
        new_Parameter2("simplex", initValues, 4, new_Constraint(0.0, 1.0));
    simplex->transform = transform;
    // Parameters_add(transform->parameter->listeners->parameters, simplex);
    transform->parameter->listeners->add_parameter(transform->parameter->listeners, simplex);
    const double* values2 = Parameter_values(simplex);
    for (size_t i = 0; i < Parameter_size(simplex); i++) {
        mu_assert(fabs(values2[i] - values[i]) < 1.e-7,
                  "simplex: constrained values not matching");
    }

    double jacobian[12];
    double trueJacobian[3][4] = {
        {0.09000000000000002, -0.020000000000000004, -0.03000000000000001,
         -0.04000000000000001},
        {0.0, 0.15555555555555556, -0.06666666666666668, -0.0888888888888889},
        {0.0, 0.0, 0.17142857142857146, -0.17142857142857146}};
    transform->jacobian(transform, jacobian);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 4; j++) {
            mu_assert(fabs(jacobian[i * 4 + j] - trueJacobian[i][j]) < 1.e-7,
                      "simplex: jacobian not matching");
        }
    }

    const double ingrad[4] = {1., 2., 3., 4.};
    const double trueGradient[3] = {-0.2, -0.24444444444444446, -0.17142857142857149};
    transform->backward(transform, ingrad);
    for (size_t i = 0; i < 3; i++) {
        mu_assert(fabs(parameter->grad[i] - trueGradient[i]) < 1.e-7,
                  "simplex: gradient not matching");
    }

    double trueLogDetJacobian = -6.0322865416282365;
    double logDetJacobian = transform->log_det_jacobian(transform);
    mu_assert(fabs(trueLogDetJacobian - logDetJacobian) < 1.e-7,
              "simplex: log det Jacobian not matching");

    double trueLogDetJacobianGradient[3] = {0.5999999999999999, 0.33333333333333337,
                                            0.1428571428571428};
    Parameter_zero_grad(parameter);
    transform->gradient_log_det_jacobian(transform);

    for (size_t i = 0; i < 3; i++) {
        mu_assert(fabs(parameter->grad[i] - trueLogDetJacobianGradient[i]) < 1.e-7,
                  "log det Jacobian gradient not matching");
    }
    return NULL;
}

// The "proportions" simplex transform S -> X (pure stick-breaking from the
// break-fractions S in (0,1)^{K-1}), verified two ways:
//   (A) S fed by a raw box-constrained (0,1) leaf,
//   (B) S fed by a logit-backed unconstrained leaf (U -> S -> X composition),
// which must refresh X through the chain and backprop dL/dX all the way to U.
// Reference numbers computed with numpy (finite differences).
char* test_simplex_proportions() {
    double initValues[4] = {0.2, 0.1, 0.4, 0.3};
    double values[4] = {0.1, 0.2, 0.3, 0.4};
    // break-fractions of [0.1, 0.2, 0.3, 0.4]: s_k = x_k / (1 - sum_{j<k} x_j)
    double propValues[3] = {0.1, 0.22222222222222224, 0.4285714285714286};

    // (A) raw (0,1) leaf S
    Parameter* parameter =
        new_Parameter2("s", propValues, 3, new_Constraint(0.0, 1.0));
    Transform* transform =
        new_SimplexTransform_with_parameter("proportions", parameter);
    Parameter* simplex =
        new_Parameter2("simplex", initValues, 4, new_Constraint(0.0, 1.0));
    simplex->transform = transform;
    transform->parameter->listeners->add_parameter(transform->parameter->listeners,
                                                   simplex);
    const double* values2 = Parameter_values(simplex);
    for (size_t i = 0; i < Parameter_size(simplex); i++) {
        mu_assert(fabs(values2[i] - values[i]) < 1.e-7,
                  "simplex proportions: constrained values not matching");
    }

    double jacobian[12];
    double trueJacobian[3][4] = {
        {1.0, -0.22222222222222224, -0.33333333333333337, -0.4444444444444445},
        {0.0, 0.9, -0.3857142857142857, -0.5142857142857143},
        {0.0, 0.0, 0.7, -0.7}};
    transform->jacobian(transform, jacobian);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 4; j++) {
            mu_assert(fabs(jacobian[i * 4 + j] - trueJacobian[i][j]) < 1.e-7,
                      "simplex proportions: jacobian not matching");
        }
    }

    const double ingrad[4] = {1., 2., 3., 4.};
    const double trueGradient[3] = {-2.2222222222222223, -1.4142857142857141, -0.7};
    transform->backward(transform, ingrad);
    for (size_t i = 0; i < 3; i++) {
        mu_assert(fabs(parameter->grad[i] - trueGradient[i]) < 1.e-7,
                  "simplex proportions: gradient not matching");
    }

    double trueLogDetJacobian = -0.46203545959655873;  // log(0.9) + log(0.7)
    double logDetJacobian = transform->log_det_jacobian(transform);
    mu_assert(fabs(trueLogDetJacobian - logDetJacobian) < 1.e-7,
              "simplex proportions: log det Jacobian not matching");

    double trueLogDetJacobianGradient[3] = {-2.2222222222222223, -1.2857142857142858,
                                            0.0};
    Parameter_zero_grad(parameter);
    transform->gradient_log_det_jacobian(transform);
    for (size_t i = 0; i < 3; i++) {
        mu_assert(fabs(parameter->grad[i] - trueLogDetJacobianGradient[i]) < 1.e-7,
                  "simplex proportions: log det Jacobian gradient not matching");
    }

    // (B) S fed by a logit-backed unconstrained leaf: U --logit--> S --S->X--> X.
    // Seed U so that S matches propValues above.
    double u[3] = {logit(propValues[0]), logit(propValues[1]), logit(propValues[2])};
    Parameter* leaf =
        new_Parameter2("u", u, 3, new_Constraint(-INFINITY, INFINITY));
    Transform* logitT = new_Transform_with_parameter(NULL, 0.0, 1.0, leaf);
    Parameter* s = new_Parameter2("s2", propValues, 3, new_Constraint(0.0, 1.0));
    s->transform = logitT;
    leaf->listeners->add_parameter(leaf->listeners, s);
    Transform* sxT = new_SimplexTransform_with_parameter("proportions", s);
    Parameter* x = new_Parameter2("x2", initValues, 4, new_Constraint(0.0, 1.0));
    x->transform = sxT;
    s->listeners->add_parameter(s->listeners, x);

    // X must refresh through the whole U -> S -> X chain.
    const double* xv = Parameter_values(x);
    for (size_t i = 0; i < 4; i++) {
        mu_assert(fabs(xv[i] - values[i]) < 1.e-7,
                  "simplex proportions (chain): constrained values not matching");
    }

    // Backprop dL/dX -> dL/dS (into s->grad) -> dL/dU (into leaf->grad). The
    // additive logit offset does not affect ds/du, so leaf->grad must equal the
    // fused stan simplex backward for the same X (cf. test_simplex).
    double trueLeafGradient[3] = {-0.2, -0.24444444444444446, -0.17142857142857149};
    Parameter_zero_grad(leaf);
    Parameter_zero_grad(s);
    sxT->backward(sxT, ingrad);       // dL/dS into s->grad
    logitT->backward(logitT, s->grad);  // dL/dS -> dL/dU into leaf->grad
    for (size_t i = 0; i < 3; i++) {
        mu_assert(fabs(leaf->grad[i] - trueLeafGradient[i]) < 1.e-7,
                  "simplex proportions (chain): leaf gradient not matching");
    }
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_exp);
    mu_run_test(test_exp_fbounds);
    mu_run_test(test_sigmoid);
    mu_run_test(test_simplex);
    mu_run_test(test_simplex_proportions);
    return NULL;
}

RUN_TESTS(all_tests);
