#include <stdlib.h>

#include "minunit.h"
#include "phyc/parameters.h"

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
    Parameters_add(transform->parameter->listeners->parameters, positive);
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
    Parameters_add(transform->parameter->listeners->parameters, constrained);
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
    Parameters_add(transform->parameter->listeners->parameters, simplex);
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

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_exp);
    mu_run_test(test_sigmoid);
    mu_run_test(test_simplex);
    return NULL;
}

RUN_TESTS(all_tests);
