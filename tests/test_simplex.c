#include <stdlib.h>

#include "minunit.h"
#include "phyc/parameters.h"
#include "phyc/simplex.h"

char* test_simplex() {
    // clang-format off
    /*
    Checked with pytorch

    z = torch.distributions.StickBreakingTransform().inv(torch.tensor([0.1, 0.2, 0.3, 0.4], dtype=torch.double))
    def f(x):
        return torch.distributions.StickBreakingTransform()(x)
    torch.autograd.functional.jacobian(f, z).t().tolist()
    */
    // clang-format on
    double values[4] = {0.1, 0.2, 0.3, 0.4};
    double constrainedValues[3] = {-1.0986122886681093, -0.5596157879354225,
                                   -0.2876820724517808};
    Simplex* simplex = new_Simplex_with_values("a", values, 4);
    for (size_t i = 0; i < Parameters_count(simplex->parameters); i++) {
        mu_assert(fabs(Parameters_value(simplex->parameters, i) -
                       constrainedValues[i]) < 1.e-7,
                  "constrained values not matching");
    }

    double gradient[4] = {0., 0., 0., 0.};
    double trueGradient[3][4] = {
        {0.09000000000000002, -0.020000000000000004, -0.03000000000000001,
         -0.04000000000000001},
        {0.0, 0.15555555555555556, -0.06666666666666668, -0.0888888888888889},
        {0.0, 0.0, 0.17142857142857146, -0.17142857142857146}};
    for (size_t i = 0; i < 3; i++) {
        simplex->gradient(simplex, i, gradient);
        for (size_t j = 0; j < 4; j++) {
            mu_assert(fabs(gradient[j] - trueGradient[i][j]) < 1.e-7,
                      "gradient not matching");
        }
    }
    free_Simplex(simplex);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_simplex);
    return NULL;
}

RUN_TESTS(all_tests);