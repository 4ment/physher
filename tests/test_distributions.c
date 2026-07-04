//
//  test_distributions.c
//  physher
//  Created by Mathieu Fourment on 13/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include "minunit.h"
#include "phyc/distdirichlet.h"
#include "phyc/distexp.h"
#include "phyc/distgamma.h"
#include "phyc/distkumaraswamy.h"
#include "phyc/distlognormal.h"
#include "phyc/distmultinormal.h"
#include "phyc/distnormal.h"
#include "phyc/matrix.h"
#include "phyc/distweibull.h"
#include "phyc/parameters.h"

#pragma region Exponential Distribution

char* test_exponential_distribution_aux(Parameter* x, Parameter* lambda,
                                        double* dlogPdx, double* dlogPdlambda) {
    const double* lambdaValues = Parameter_values(lambda);
    const double* xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(1);
    Parameters_move(parameters, lambda);

    DistributionModel* dm = new_ExponentialDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_EXPONENTIAL_RATE);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if (Parameter_size(lambda) > 1) {
        for (size_t i = 0; i < Parameter_size(x); i++) {
            logP2 += log(gsl_ran_exponential_pdf(xValues[i], 1 / lambdaValues[i]));
        }
    } else {
        for (size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_exponential_pdf(xValues[i], 1 / lambdaValues[0]));
        }
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters* ps = new_Parameters(2);
    Parameters_add(ps, lambda);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for (size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for (size_t i = 0; i < Parameter_size(lambda); i++) {
        mu_assert(fabs(lambda->grad[i] - dlogPdlambda[i]) < 0.0001,
                  "dlogPdlambda not matching");
    }

    model->free(model);
    return NULL;
}

char* test_exponential_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double lambdaValue = 2.0;
    Parameter* lambda =
        new_Parameter("lambda", lambdaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-2.0, -2.0, -2.0};
    double dlogPdlambda = -1.0099999999999998;

    return test_exponential_distribution_aux(x, lambda, dlogPdx, &dlogPdlambda);
}

char* test_exponential_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double lambdaValues[] = {0.1, 1.0, 2.0};
    Parameter* lambda =
        new_Parameter2("lambda", lambdaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-0.1, -1.0, -2.0};
    double dlogPdlambda[3] = {8.0, 0.5, 0.49};

    return test_exponential_distribution_aux(x, lambda, dlogPdx, dlogPdlambda);
}

#pragma endregion

#pragma region Normal Distribution

char* test_normal_distribution_aux(Parameter* x, Parameter* mu, Parameter* sigma,
                                   double* dlogPdx, double* dlogPdmu,
                                   double* dlogPdsigma) {
    const double* muValues = Parameter_values(mu);
    const double* sigmaValues = Parameter_values(sigma);
    const double* xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, mu);
    Parameters_move(parameters, sigma);

    DistributionModel* dm = new_NormalDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if (Parameter_size(mu) > 1) {
        for (size_t i = 0; i < Parameter_size(x); i++) {
            logP2 +=
                log(gsl_ran_gaussian_pdf(xValues[i] - muValues[i], sigmaValues[i]));
        }
    } else {
        for (size_t i = 0; i < 3; i++) {
            logP2 +=
                log(gsl_ran_gaussian_pdf(xValues[i] - muValues[0], sigmaValues[0]));
        }
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, mu);
    Parameters_add(ps, sigma);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for (size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for (size_t i = 0; i < Parameter_size(mu); i++) {
        mu_assert(fabs(mu->grad[i] - dlogPdmu[i]) < 0.0001, "dlogPdmu not matching");
        mu_assert(fabs(sigma->grad[i] - dlogPdsigma[i]) < 0.0001,
                  "dlogPdsigma not matching");
    }

    model->free(model);
    return NULL;
}

char* test_normal_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValue = 2.0;
    double sigmaValue = 0.5;
    Parameter* mu = new_Parameter("mu", muValue, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter("sigma", sigmaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-0.0, 6.0, 7.96};
    double dlogPdmu = -13.96;
    double dlogPdsigma = 43.680800000000005;

    return test_normal_distribution_aux(x, mu, sigma, dlogPdx, &dlogPdmu, &dlogPdsigma);
}

char* test_normal_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValues[] = {1.0, 2.0, 3.0};
    double sigmaValues[] = {0.5, 0.1, 0.2};
    Parameter* mu =
        new_Parameter2("mu", muValues, 3, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma =
        new_Parameter2("sigma", sigmaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-4.0, 149.99999999999997, 74.75};
    double dlogPdmu[3] = {4.0, -149.99999999999997, -74.75};
    double dlogPdsigma[3] = {6.0, 2239.999999999999, 1112.5124999999996};

    return test_normal_distribution_aux(x, mu, sigma, dlogPdx, dlogPdmu, dlogPdsigma);
}

// Central finite-difference of model logP wrt element i of parameter p.
static double fd_logP_grad(Model* model, Parameter* p, size_t i) {
    const double h = 1e-5;
    double v = Parameter_value_at(p, i);
    Parameter_set_value_at(p, v + h, i);
    double lpPlus = model->logP(model);
    Parameter_set_value_at(p, v - h, i);
    double lpMinus = model->logP(model);
    Parameter_set_value_at(p, v, i);
    return (lpPlus - lpMinus) / (2.0 * h);
}

// Exercise a Normal distribution whose dm->x holds several Parameters (multiple
// random variables sharing the same distribution), for both modes:
//   mode 1 (prior):      scalar mu/sigma shared across every x element
//   mode 2 (variational): mu/sigma sized to the total number of x elements
// Analytic gradients (accumulated into the leaves) are checked against finite
// differences of logP.
char* test_normal_distribution_multi_aux(Parameters* xs, Parameter* mu,
                                         Parameter* sigma) {
    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, mu);
    Parameters_add(parameters, sigma);

    DistributionModel* dm = new_NormalDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    // logP must equal the sum of per-element Gaussian densities.
    double logP = model->logP(model);
    double logP2 = 0;
    size_t index = 0;
    size_t scalar = (Parameter_size(mu) == 1);
    const double* muValues = Parameter_values(mu);
    const double* sigmaValues = Parameter_values(sigma);
    for (size_t k = 0; k < Parameters_count(xs); k++) {
        Parameter* x = Parameters_at(xs, k);
        const double* xValues = Parameter_values(x);
        for (size_t j = 0; j < Parameter_size(x); j++) {
            size_t idx = scalar ? 0 : index;
            logP2 += log(gsl_ran_gaussian_pdf(xValues[j] - muValues[idx],
                                              sigmaValues[idx]));
            index++;
        }
    }
    mu_assert(fabs(logP - logP2) < 1e-10, "multi-x logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);

    Parameters* ps = new_Parameters(2 + Parameters_count(xs));
    Parameters_add(ps, mu);
    Parameters_add(ps, sigma);
    for (size_t k = 0; k < Parameters_count(xs); k++) {
        Parameters_add(ps, Parameters_at(xs, k));
    }

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    // Compare every accumulated leaf gradient with its finite difference.
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "multi-x gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// mode 1: two x parameters of differing sizes share one scalar mu and sigma.
char* test_normal_distribution_multi_prior() {
    double x0v[] = {1.5, 2.5};
    double x1v[] = {2.0};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(-INFINITY, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(-INFINITY, INFINITY));

    Parameter* mu = new_Parameter("mu", 2.0, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter("sigma", 1.0, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_normal_distribution_multi_aux(xs, mu, sigma);
}

// mode 2: two x parameters whose elements each have their own mu and sigma.
char* test_normal_distribution_multi() {
    double x0v[] = {1.5, 2.5};
    double x1v[] = {2.0};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(-INFINITY, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(-INFINITY, INFINITY));

    double muValues[] = {1.0, 2.0, 3.0};
    double sigmaValues[] = {1.0, 1.5, 0.8};
    Parameter* mu =
        new_Parameter2("mu", muValues, 3, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma =
        new_Parameter2("sigma", sigmaValues, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_normal_distribution_multi_aux(xs, mu, sigma);
}

// Build a positive parameter x = exp(y) backed by an unconstrained leaf y, wired
// exactly like new_Parameter_from_json: the constrained parameter carries the
// transform and listens to its leaf. Returns the constrained parameter; the leaf
// is reachable via p->transform->parameter.
static Parameter* new_exp_transformed_parameter(const char* name, double y) {
    Parameter* leaf = new_Parameter2(name, &y, 1, new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_Transform_with_parameter(NULL, 0, INFINITY, leaf);
    double x;
    transform->get(transform, &x);
    Parameter* p = new_Parameter2(name, &x, 1, new_Constraint(0, INFINITY));
    p->transform = transform;
    transform->parameter->listeners->add_parameter(transform->parameter->listeners, p);
    return p;
}

// mode 1 with mu and sigma each carrying an exp-transform: the distribution holds
// the constrained mu/sigma, but the optimizer differentiates the unconstrained
// leaves. This drives the transform backward() path with a shared scalar
// hyperparameter over multiple x — the case the mode-1 fix targets. Gradients of
// the unconstrained leaves are checked against finite differences.
char* test_normal_distribution_multi_prior_transformed() {
    double x0v[] = {1.5, 2.5};
    double x1v[] = {2.0};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(-INFINITY, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(-INFINITY, INFINITY));

    Parameter* mu = new_exp_transformed_parameter("mu", log(1.5));      // mu = 1.5
    Parameter* sigma = new_exp_transformed_parameter("sigma", log(1.0));  // sigma = 1

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, mu);
    Parameters_add(parameters, sigma);

    DistributionModel* dm = new_NormalDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP (transformed mu/sigma): %f\n", logP);

    // The optimizer sees the unconstrained leaves, so the gradient flows through
    // the exp transform's backward() into leaf->grad.
    Parameter* muUnc = mu->transform->parameter;
    Parameter* sigmaUnc = sigma->transform->parameter;
    Parameters* ps = new_Parameters(4);
    Parameters_add(ps, muUnc);
    Parameters_add(ps, sigmaUnc);
    Parameters_add(ps, x0);
    Parameters_add(ps, x1);

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "transformed multi-x gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}
#pragma endregion

#pragma region Half Normal Distribution

char* test_half_normal_distribution_aux(Parameter* x, Parameter* sigma, double* dlogPdx,
                                        double* dlogPdsigma) {
    const double* sigmaValues = Parameter_values(sigma);
    const double* xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(1);
    Parameters_move(parameters, sigma);

    DistributionModel* dm = new_HalfNormalDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if (Parameter_size(sigma) > 1) {
        for (size_t i = 0; i < Parameter_size(x); i++) {
            logP2 += log(gsl_ran_gaussian_pdf(xValues[i], sigmaValues[i]) * 2.0);
        }
    } else {
        for (size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_gaussian_pdf(xValues[i], sigmaValues[0]) * 2.0);
        }
    }
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    mu_assert(logP2 == logP, "logP not matching");
    Parameters* ps = new_Parameters(2);
    Parameters_add(ps, sigma);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for (size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for (size_t i = 0; i < Parameter_size(sigma); i++) {
        mu_assert(fabs(sigma->grad[i] - dlogPdsigma[i]) < 0.0001,
                  "dlogPdsigma not matching");
    }

    model->free(model);
    return NULL;
}

char* test_half_normal_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double sigmaValue = 0.5;
    Parameter* sigma = new_Parameter("sigma", sigmaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-8.0, -2.0, -0.04};
    double dlogPdsigma = 28.000799999999998;

    return test_half_normal_distribution_aux(x, sigma, dlogPdx, &dlogPdsigma);
}

char* test_half_normal_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double sigmaValues[] = {0.5, 0.1, 0.2};
    Parameter* sigma =
        new_Parameter2("sigma", sigmaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-8.0, -49.99999999999999, -0.24999999999999997};
    double dlogPdsigma[3] = {30.0, 239.99999999999991, -4.9875};

    return test_half_normal_distribution_aux(x, sigma, dlogPdx, dlogPdsigma);
}
#pragma endregion

#pragma region Lognormal Distribution

char* test_lognormal_distribution_aux(Parameter* x, Parameter* mu, Parameter* sigma,
                                      double* dlogPdx, double* dlogPdmu,
                                      double* dlogPdsigma) {
    const double* muValues = Parameter_values(mu);
    const double* sigmaValues = Parameter_values(sigma);
    const double* xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, mu);
    Parameters_move(parameters, sigma);

    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_LOGNORMAL_MU_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if (Parameter_size(mu) > 1) {
        for (size_t i = 0; i < Parameter_size(x); i++) {
            logP2 +=
                log(gsl_ran_lognormal_pdf(xValues[i], muValues[i], sigmaValues[i]));
        }
    } else {
        for (size_t i = 0; i < 3; i++) {
            logP2 +=
                log(gsl_ran_lognormal_pdf(xValues[i], muValues[0], sigmaValues[0]));
        }
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, mu);
    Parameters_add(ps, sigma);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for (size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for (size_t i = 0; i < Parameter_size(mu); i++) {
        mu_assert(fabs(mu->grad[i] - dlogPdmu[i]) < 0.0001, "dlogPdmu not matching");
        mu_assert(fabs(sigma->grad[i] - dlogPdsigma[i]) < 0.0001,
                  "dlogPdsigma not matching");
    }

    model->free(model);
    return NULL;
}

char* test_lognormal_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValue = 2.0;
    double sigmaValue = 0.5;
    Parameter* mu = new_Parameter("mu", muValue, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter("sigma", sigmaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {2.113705638880109, 19.545177444479563, 2542.068074395236};
    double dlogPdmu = -42.42068074395236;
    double dlogPdsigma = 414.7134337096188;

    return test_lognormal_distribution_aux(x, mu, sigma, dlogPdx, &dlogPdmu,
                                           &dlogPdsigma);
}

char* test_lognormal_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double muValues[] = {1.0, 2.0, 3.0};
    double sigmaValues[] = {0.5, 0.1, 0.2};
    Parameter* mu =
        new_Parameter2("mu", muValues, 3, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma =
        new_Parameter2("sigma", sigmaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {0.11370563888010943, 536.629436111989, 18912.925464970223};
    double dlogPdmu[3] = {-1.2274112777602189, -269.3147180559945, -190.12925464970223};
    double dlogPdsigma[3] = {-1.2467307776135135, 7243.041736157981, 7224.826694730264};

    return test_lognormal_distribution_aux(x, mu, sigma, dlogPdx, dlogPdmu,
                                           dlogPdsigma);
}

char* test_lognormal_distribution_prior_mean_stdev() {
    Parameter* m = new_Parameter("mean", 2.0, new_Constraint(0, INFINITY));
    Parameter* stdev = new_Parameter("stdev", 0.5, new_Constraint(0, INFINITY));
    double mValue = Parameter_value(m);
    double stdevValue = Parameter_value(stdev);

    Parameter* x1 = new_Parameter("x1", 0.5, new_Constraint(0, INFINITY));
    double x3Values[] = {2.0, 0.5, 0.01};
    Parameter* x3 = new_Parameter2("x3", x3Values, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x3);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, m);
    Parameters_move(parameters, stdev);

    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_LOGNORMAL_MEAN_STDEV);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    for (size_t i = 0; i < 3; i++) {
        double sigma = sqrt(log(1.0 + (stdevValue * stdevValue) / (mValue * mValue)));
        double mu = log(mValue) - 0.5 * sigma * sigma;
        logP2 += log(gsl_ran_lognormal_pdf(x3Values[i], mu, sigma));
    }
    mu_assert(logP2 == logP, "logP not matching");
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, m);
    Parameters_add(ps, stdev);
    Parameters_add(ps, x3);
    model->gradient(model, ps);
    double dlogPdx[3] = {-0.7499999999999998, 42.73370751301173, 8589.547081367036};
    double dlogPdm = -292.9289997760196;
    double dlogPdstdev = 954.1913499637261;

    mu_assert(fabs(x3->grad[0] - dlogPdx[0]) < 0.0001, "dlogPdx1 not matching");
    mu_assert(fabs(x3->grad[1] - dlogPdx[1]) < 0.0001, "dlogPdx2 not matching");
    mu_assert(fabs(x3->grad[2] - dlogPdx[2]) < 0.0001, "dlogPdx3 not matching");
    mu_assert(fabs(m->grad[0] - dlogPdm) < 0.0001, "dlogPdm not matching");
    mu_assert(fabs(stdev->grad[0] - dlogPdstdev) < 0.0001, "dlogPdstdev not matching");

    printf("Gradients:\n");
    for (size_t i = 0; i < 3; i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    model->free(model);
    return NULL;
}

// Exercise a LogNormal distribution whose dm->x holds several Parameters, for
// both modes and both parameterizations. Analytic gradients (accumulated into the
// leaves) are checked against finite differences of logP.
char* test_lognormal_distribution_multi_aux(Parameters* xs, Parameter* mu,
                                            Parameter* sigma,
                                            distribution_parameterization param) {
    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, mu);
    Parameters_add(parameters, sigma);

    DistributionModel* dm =
        new_LogNormalDistributionModel_with_parameters(parameters, xs, param);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP: %f\n", logP);

    Parameters* ps = new_Parameters(2 + Parameters_count(xs));
    Parameters_add(ps, mu);
    Parameters_add(ps, sigma);
    for (size_t k = 0; k < Parameters_count(xs); k++) {
        Parameters_add(ps, Parameters_at(xs, k));
    }

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    // Compare every accumulated leaf gradient with its finite difference.
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "multi-x lognormal gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// mode 1: two x parameters of differing sizes share one scalar mu and sigma.
char* test_lognormal_distribution_multi_prior() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    Parameter* mu = new_Parameter("mu", 0.5, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma = new_Parameter("sigma", 0.6, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_lognormal_distribution_multi_aux(xs, mu, sigma,
                                                 DISTRIBUTION_LOGNORMAL_MU_SIGMA);
}

// mode 2: two x parameters whose elements each have their own mu and sigma.
char* test_lognormal_distribution_multi() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    double muValues[] = {0.1, 0.5, 1.0};
    double sigmaValues[] = {0.5, 0.8, 0.3};
    Parameter* mu =
        new_Parameter2("mu", muValues, 3, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma =
        new_Parameter2("sigma", sigmaValues, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_lognormal_distribution_multi_aux(xs, mu, sigma,
                                                 DISTRIBUTION_LOGNORMAL_MU_SIGMA);
}

// mode 1 with mu and sigma each carrying an exp-transform: the optimizer
// differentiates the unconstrained leaves, driving the transform backward() path
// with a shared scalar hyperparameter over multiple x.
char* test_lognormal_distribution_multi_prior_transformed() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    Parameter* mu = new_exp_transformed_parameter("mu", log(1.5));        // mu = 1.5
    Parameter* sigma = new_exp_transformed_parameter("sigma", log(0.6));  // sigma = 0.6

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, mu);
    Parameters_add(parameters, sigma);

    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_LOGNORMAL_MU_SIGMA);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP (transformed mu/sigma): %f\n", logP);

    // The optimizer sees the unconstrained leaves, so the gradient flows through
    // the exp transform's backward() into leaf->grad.
    Parameter* muUnc = mu->transform->parameter;
    Parameter* sigmaUnc = sigma->transform->parameter;
    Parameters* ps = new_Parameters(4);
    Parameters_add(ps, muUnc);
    Parameters_add(ps, sigmaUnc);
    Parameters_add(ps, x0);
    Parameters_add(ps, x1);

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "transformed multi-x lognormal gradient does not match FD");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// mean/stdev parameterization, mode 1: two x share one scalar mean and stdev.
char* test_lognormal_distribution_multi_prior_mean_stdev() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    Parameter* mean = new_Parameter("mean", 2.0, new_Constraint(0, INFINITY));
    Parameter* stdev = new_Parameter("stdev", 0.5, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_lognormal_distribution_multi_aux(xs, mean, stdev,
                                                 DISTRIBUTION_LOGNORMAL_MEAN_STDEV);
}

// mean/stdev parameterization, mode 2: two x whose elements each have their own
// mean and stdev (exercises the vector branches).
char* test_lognormal_distribution_multi_mean_stdev() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    double meanValues[] = {1.5, 2.0, 1.0};
    double stdevValues[] = {0.5, 0.8, 0.3};
    Parameter* mean = new_Parameter2("mean", meanValues, 3, new_Constraint(0, INFINITY));
    Parameter* stdev =
        new_Parameter2("stdev", stdevValues, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_lognormal_distribution_multi_aux(xs, mean, stdev,
                                                 DISTRIBUTION_LOGNORMAL_MEAN_STDEV);
}

#pragma endregion

#pragma region Gamma Distribution

char* test_gamma_distribution_aux(Parameter* x, Parameter* alpha, Parameter* beta,
                                  double* dlogPdx, double* dlogPdalpha,
                                  double* dlogPdbeta,
                                  distribution_parameterization parameterization) {
    const double* alphaValues = Parameter_values(alpha);
    const double* betaValues = Parameter_values(beta);
    const double* xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, alpha);
    Parameters_move(parameters, beta);

    DistributionModel* dm =
        new_GammaDistributionModel_with_parameters(parameters, xs, parameterization);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if (Parameter_size(alpha) > 1) {
        for (size_t i = 0; i < Parameter_size(x); i++) {
            double betaValue = parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE
                                   ? 1.0 / betaValues[i]
                                   : betaValues[i];
            logP2 += log(gsl_ran_gamma_pdf(xValues[i], alphaValues[i], betaValue));
        }
    } else {
        double betaValue = parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE
                               ? 1.0 / betaValues[0]
                               : betaValues[0];
        for (size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_gamma_pdf(xValues[i], alphaValues[0], betaValue));
        }
    }
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    mu_assert(logP2 == logP, "logP not matching");
    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, alpha);
    Parameters_add(ps, beta);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for (size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for (size_t i = 0; i < Parameter_size(alpha); i++) {
        mu_assert(fabs(alpha->grad[i] - dlogPdalpha[i]) < 0.0001,
                  "dlogPdalpha not matching");
        mu_assert(fabs(beta->grad[i] - dlogPdbeta[i]) < 0.0001,
                  "dlogPdbeta not matching");
    }

    model->free(model);
    return NULL;
}

char* test_gamma_distribution_shape_rate_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValue = 2.0;
    double betaValue = 0.5;
    Parameter* alpha = new_Parameter("alpha", alphaValue, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter("beta", betaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {0.0, 1.5, 99.5};
    double dlogPdalpha = -7.952964732963327;
    double dlogPdbeta = 9.49;

    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, &dlogPdalpha,
                                       &dlogPdbeta, DISTRIBUTION_GAMMA_SHAPE_RATE);
}

char* test_gamma_distribution_shape_rate() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValues[] = {1.0, 2.0, 3.0};
    double betaValues[] = {0.5, 0.1, 0.2};
    Parameter* alpha =
        new_Parameter2("alpha", alphaValues, 3, new_Constraint(0, INFINITY));
    Parameter* beta =
        new_Parameter2("beta", betaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-0.5, 1.9, 199.8};
    double dlogPdalpha[3] = {0.5772156649015329, -3.4185166086524577,
                             -7.137392433520659};
    double dlogPdbeta[3] = {0.0, 19.5, 14.99};
    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, dlogPdalpha, dlogPdbeta,
                                       DISTRIBUTION_GAMMA_SHAPE_RATE);
}

char* test_gamma_distribution_shape_scale() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValues[] = {1.0, 2.0, 3.0};
    double betaValues[] = {0.5, 0.1, 0.2};
    Parameter* alpha =
        new_Parameter2("alpha", alphaValues, 3, new_Constraint(0, INFINITY));
    Parameter* beta =
        new_Parameter2("beta", betaValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-2.0, -8.0, 195.0};
    double dlogPdalpha[3] = {1.9635100260214235, 1.1866535773356337,
                             -3.9185166086524577};
    double dlogPdbeta[3] = {6.0, 30.0, -14.75};
    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, dlogPdalpha, dlogPdbeta,
                                       DISTRIBUTION_GAMMA_SHAPE_SCALE);
}

char* test_gamma_distribution_shape_scale_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double alphaValue = 2.0;
    double betaValue = 0.5;
    Parameter* alpha = new_Parameter("alpha", alphaValue, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter("beta", betaValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-1.5, 0.0, 98.0};
    double dlogPdalpha = -3.794081649603656;
    double dlogPdbeta = -1.9600000000000009;

    return test_gamma_distribution_aux(x, alpha, beta, dlogPdx, &dlogPdalpha,
                                       &dlogPdbeta, DISTRIBUTION_GAMMA_SHAPE_SCALE);
}

// Exercise a Gamma distribution whose dm->x holds several Parameters, for both
// modes. Analytic gradients (accumulated into the leaves) are checked against
// finite differences of logP.
char* test_gamma_distribution_multi_aux(Parameters* xs, Parameter* alpha,
                                        Parameter* beta,
                                        distribution_parameterization parameterization) {
    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, alpha);
    Parameters_add(parameters, beta);

    DistributionModel* dm = new_GammaDistributionModel_with_parameters(
        parameters, xs, parameterization);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP: %f\n", logP);

    Parameters* ps = new_Parameters(2 + Parameters_count(xs));
    Parameters_add(ps, alpha);
    Parameters_add(ps, beta);
    for (size_t k = 0; k < Parameters_count(xs); k++) {
        Parameters_add(ps, Parameters_at(xs, k));
    }

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "multi-x gamma gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// mode 1: two x parameters of differing sizes share one scalar alpha and beta.
char* test_gamma_distribution_multi_prior() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    Parameter* alpha = new_Parameter("alpha", 2.0, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter("beta", 1.5, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_gamma_distribution_multi_aux(xs, alpha, beta,
                                             DISTRIBUTION_GAMMA_SHAPE_RATE);
}

// mode 2: two x parameters whose elements each have their own alpha and beta.
char* test_gamma_distribution_multi() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    double alphaValues[] = {2.0, 1.5, 3.0};
    double betaValues[] = {1.0, 2.0, 1.5};
    Parameter* alpha =
        new_Parameter2("alpha", alphaValues, 3, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter2("beta", betaValues, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_gamma_distribution_multi_aux(xs, alpha, beta,
                                             DISTRIBUTION_GAMMA_SHAPE_SCALE);
}

// mode 1 with alpha and beta each carrying an exp-transform: the optimizer
// differentiates the unconstrained leaves, driving the transform backward() path
// with a shared scalar hyperparameter over multiple x.
char* test_gamma_distribution_multi_prior_transformed() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    Parameter* alpha = new_exp_transformed_parameter("alpha", log(2.0));
    Parameter* beta = new_exp_transformed_parameter("beta", log(1.5));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, alpha);
    Parameters_add(parameters, beta);

    DistributionModel* dm = new_GammaDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_GAMMA_SHAPE_RATE);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP (transformed alpha/beta): %f\n", logP);

    Parameter* alphaUnc = alpha->transform->parameter;
    Parameter* betaUnc = beta->transform->parameter;
    Parameters* ps = new_Parameters(4);
    Parameters_add(ps, alphaUnc);
    Parameters_add(ps, betaUnc);
    Parameters_add(ps, x0);
    Parameters_add(ps, x1);

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "transformed multi-x gamma gradient does not match FD");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// Validate the vendored derivative of the regularized lower incomplete gamma
// P(a, x) wrt a against a central finite difference of GSL's gsl_sf_gamma_inc_P,
// spanning both evaluation regimes (x < a+1 series and x >= a+1 continued fraction).
char* test_gamma_p_grad_a() {
    double as[] = {0.5, 1.0, 2.0, 5.0, 10.0};
    double xs[] = {0.1, 0.5, 1.0, 3.0, 8.0, 15.0};
    double h = 1e-6;
    for (size_t ia = 0; ia < sizeof(as) / sizeof(as[0]); ia++) {
        for (size_t ix = 0; ix < sizeof(xs) / sizeof(xs[0]); ix++) {
            double a = as[ia], x = xs[ix];
            double analytic = gamma_p_grad_a(a, x);
            double fd = (gsl_sf_gamma_inc_P(a + h, x) - gsl_sf_gamma_inc_P(a - h, x)) /
                        (2.0 * h);
            printf("dP/da(a=%.1f, x=%.1f): %g (fd %g)\n", a, x, analytic, fd);
            mu_assert(fabs(analytic - fd) < 1e-5 * (1.0 + fabs(fd)),
                      "gamma_p_grad_a does not match finite difference");
        }
    }
    return NULL;
}

// Validate the implicit reparameterization gradient. For each element the
// reparameterized sample satisfies P(a, z) = u with z = x/theta and u fixed, so
// dx/da = theta * d/da Pinv(u; a, 1) and dx/dbeta = -z/beta^2 (rate). These are
// independent oracles (CDF inversion / closed form) for the rgradient output.
char* test_gamma_rgradient() {
    double xVals[] = {1.5, 3.0, 0.5};
    Parameter* x = new_Parameter2("x", xVals, 3, new_Constraint(0, INFINITY));
    double aVals[] = {2.0, 1.5, 3.0};
    double bVals[] = {1.0, 2.0, 1.5};  // rate
    Parameter* alpha = new_Parameter2("alpha", aVals, 3, new_Constraint(0, INFINITY));
    Parameter* beta = new_Parameter2("beta", bVals, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);
    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, alpha);
    Parameters_move(parameters, beta);

    DistributionModel* dm = new_GammaDistributionModel_with_parameters(
        parameters, xs, DISTRIBUTION_GAMMA_SHAPE_RATE);
    Model* model = new_DistributionModel2("dist", dm);

    double upstream[] = {0.7, -1.3, 0.4};  // downstream dL/dx_i
    Parameters* ps = new_Parameters(2);
    Parameters_add(ps, alpha);
    Parameters_add(ps, beta);
    Parameters_zero_grad(ps);
    for (size_t i = 0; i < 3; i++) x->grad[i] = upstream[i];

    dm->rgradient(dm);

    double h = 1e-4;
    for (size_t i = 0; i < 3; i++) {
        double a = aVals[i], bta = bVals[i], theta = 1.0 / bta;
        double z = xVals[i] / theta;
        double u = gsl_cdf_gamma_P(z, a, 1.0);
        double dxda = theta *
                      (gsl_cdf_gamma_Pinv(u, a + h, 1.0) -
                       gsl_cdf_gamma_Pinv(u, a - h, 1.0)) /
                      (2.0 * h);
        double expA = upstream[i] * dxda;
        double expB = upstream[i] * (-z / (bta * bta));
        printf("alpha[%zu] %g (inv %g)  beta[%zu] %g (exact %g)\n", i,
               alpha->grad[i], expA, i, beta->grad[i], expB);
        mu_assert(fabs(alpha->grad[i] - expA) < 1e-3 * (1.0 + fabs(expA)),
                  "gamma rgradient dalpha does not match CDF-inversion FD");
        mu_assert(fabs(beta->grad[i] - expB) < 1e-8 * (1.0 + fabs(expB)),
                  "gamma rgradient dbeta does not match closed form");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

#pragma endregion

#pragma region Weibull Distribution

char* test_weibull_distribution_aux(Parameter* x, Parameter* scale, Parameter* shape,
                                    double* dlogPdx, double* dlogPdscale,
                                    double* dlogPdshape) {
    const double* scaleValues = Parameter_values(scale);
    const double* shapeValues = Parameter_values(shape);
    const double* xValues = Parameter_values(x);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, scale);
    Parameters_move(parameters, shape);

    DistributionModel* dm =
        new_WeibullDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    if (Parameter_size(scale) > 1) {
        for (size_t i = 0; i < Parameter_size(x); i++) {
            logP2 += log(gsl_ran_weibull_pdf(xValues[i], scaleValues[i], shapeValues[i]));
        }
    } else {
        for (size_t i = 0; i < 3; i++) {
            logP2 += log(gsl_ran_weibull_pdf(xValues[i], scaleValues[0], shapeValues[0]));
        }
    }
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    mu_assert(logP2 == logP, "logP not matching");
    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, scale);
    Parameters_add(ps, shape);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for (size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }

    for (size_t i = 0; i < Parameter_size(scale); i++) {
        mu_assert(fabs(scale->grad[i] - dlogPdscale[i]) < 0.0001,
                  "dlogPdscale not matching");
        mu_assert(fabs(shape->grad[i] - dlogPdshape[i]) < 0.0001,
                  "dlogPdshape not matching");
    }

    model->free(model);
    return NULL;
}

char* test_weibull_distribution_prior() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double scaleValue = 0.5;
    double shapeValue = 2.0;
    Parameter* scale = new_Parameter("scale", scaleValue, new_Constraint(0, INFINITY));
    Parameter* shape = new_Parameter("shape", shapeValue, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-15.5, -2.0, 99.92};
    double dlogPdscale = 56.001599999999996;
    double dlogPdshape = -23.204873613024333;

    return test_weibull_distribution_aux(x, scale, shape, dlogPdx, &dlogPdscale,
                                         &dlogPdshape);
}

char* test_weibull_distribution() {
    double xValues[] = {2.0, 0.5, 0.01};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, INFINITY));

    double scaleValues[] = {1.0, 2.0, 3.0};
    double shapeValues[] = {0.5, 0.1, 0.2};
    Parameter* scale =
        new_Parameter2("scale", scaleValues, 3, new_Constraint(0, INFINITY));
    Parameter* shape =
        new_Parameter2("shape", shapeValues, 3, new_Constraint(0, INFINITY));

    double dlogPdx[3] = {-0.6035533905932737, -1.9741101126592249,
                         -86.39154343676122};
    double dlogPdscale[3] = {0.20710678118654757, -0.006472471835193794,
                             -0.04536152187746261};
    double dlogPdshape[3] = {1.712889037091398, 9.820544975847271,
                             1.119016197373924};

    return test_weibull_distribution_aux(x, scale, shape, dlogPdx, dlogPdscale,
                                         dlogPdshape);
}

// Exercise a Weibull distribution whose dm->x holds several Parameters, for both
// modes. Analytic gradients (accumulated into the leaves) are checked against
// finite differences of logP.
char* test_weibull_distribution_multi_aux(Parameters* xs, Parameter* scale,
                                          Parameter* shape) {
    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, scale);
    Parameters_add(parameters, shape);

    DistributionModel* dm =
        new_WeibullDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP: %f\n", logP);

    Parameters* ps = new_Parameters(2 + Parameters_count(xs));
    Parameters_add(ps, scale);
    Parameters_add(ps, shape);
    for (size_t k = 0; k < Parameters_count(xs); k++) {
        Parameters_add(ps, Parameters_at(xs, k));
    }

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "multi-x weibull gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// mode 1: two x parameters of differing sizes share one scalar scale and shape.
char* test_weibull_distribution_multi_prior() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    Parameter* scale = new_Parameter("scale", 1.5, new_Constraint(0, INFINITY));
    Parameter* shape = new_Parameter("shape", 2.0, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_weibull_distribution_multi_aux(xs, scale, shape);
}

// mode 2: two x parameters whose elements each have their own scale and shape.
char* test_weibull_distribution_multi() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    double scaleValues[] = {1.0, 2.0, 1.5};
    double shapeValues[] = {2.0, 1.5, 3.0};
    Parameter* scale =
        new_Parameter2("scale", scaleValues, 3, new_Constraint(0, INFINITY));
    Parameter* shape =
        new_Parameter2("shape", shapeValues, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    return test_weibull_distribution_multi_aux(xs, scale, shape);
}

// mode 1 with scale and shape each carrying an exp-transform: the optimizer
// differentiates the unconstrained leaves, driving the transform backward() path
// with a shared scalar hyperparameter over multiple x.
char* test_weibull_distribution_multi_prior_transformed() {
    double x0v[] = {1.5, 2.0};
    double x1v[] = {0.8};
    Parameter* x0 = new_Parameter2("x0", x0v, 2, new_Constraint(0, INFINITY));
    Parameter* x1 = new_Parameter2("x1", x1v, 1, new_Constraint(0, INFINITY));

    Parameter* scale = new_exp_transformed_parameter("scale", log(1.5));
    Parameter* shape = new_exp_transformed_parameter("shape", log(2.0));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);

    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, scale);
    Parameters_add(parameters, shape);

    DistributionModel* dm =
        new_WeibullDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP (transformed scale/shape): %f\n", logP);

    Parameter* scaleUnc = scale->transform->parameter;
    Parameter* shapeUnc = shape->transform->parameter;
    Parameters* ps = new_Parameters(4);
    Parameters_add(ps, scaleUnc);
    Parameters_add(ps, shapeUnc);
    Parameters_add(ps, x0);
    Parameters_add(ps, x1);

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "transformed multi-x weibull gradient does not match FD");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

#pragma region Kumaraswamy Distribution

char* test_kumaraswamy_distribution_aux(Parameter* x, Parameter* a, Parameter* b,
                                        double expectedLogP, double* dlogPdx,
                                        double* dlogPda, double* dlogPdb) {
    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);

    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, a);
    Parameters_move(parameters, b);

    DistributionModel* dm =
        new_KumaraswamyDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    printf("LogP: %f, expected: %f\n", logP, expectedLogP);
    mu_assert(fabs(logP - expectedLogP) < 0.0001, "logP not matching");

    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, a);
    Parameters_add(ps, b);
    Parameters_add(ps, x);
    model->gradient(model, ps);

    printf("Gradients:\n");
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            printf(" %f", p->grad[j]);
        }
        printf("\n");
    }

    for (size_t i = 0; i < Parameter_size(x); i++) {
        mu_assert(fabs(x->grad[i] - dlogPdx[i]) < 0.0001, "dlogPdx not matching");
    }
    for (size_t i = 0; i < Parameter_size(a); i++) {
        mu_assert(fabs(a->grad[i] - dlogPda[i]) < 0.0001, "dlogPda not matching");
        mu_assert(fabs(b->grad[i] - dlogPdb[i]) < 0.0001, "dlogPdb not matching");
    }

    model->free(model);
    return NULL;
}

char* test_kumaraswamy_distribution_prior() {
    double xValues[] = {0.2, 0.5, 0.9};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, 1));

    double aValue = 2.0;
    double bValue = 3.0;
    Parameter* a = new_Parameter("a", aValue, new_Constraint(0, INFINITY));
    Parameter* b = new_Parameter("b", bValue, new_Constraint(0, INFINITY));

    double expectedLogP = -1.0111377485550814;
    double dlogPdx[3] = {4.166666666667, -0.666666666667, -17.836257309942};
    double dlogPda = 0.586609365998;
    double dlogPdb = -0.989235273794;

    return test_kumaraswamy_distribution_aux(x, a, b, expectedLogP, dlogPdx,
                                             &dlogPda, &dlogPdb);
}

char* test_kumaraswamy_distribution() {
    double xValues[] = {0.2, 0.5, 0.9};
    Parameter* x = new_Parameter2("x", xValues, 3, new_Constraint(0, 1));

    double aValues[] = {2.0, 0.5, 4.0};
    double bValues[] = {3.0, 1.5, 0.8};
    Parameter* a = new_Parameter2("a", aValues, 3, new_Constraint(0, INFINITY));
    Parameter* b = new_Parameter2("b", bValues, 3, new_Constraint(0, INFINITY));

    double expectedLogP = 0.6061456320728635;
    double dlogPdx[3] = {4.166666666667, -2.207106781187, 5.029175147814};
    double dlogPda[3] = {-0.975318086398, 2.143555481454, 0.104437661531};
    double dlogPdb[3] = {0.292511338813, -0.561280510633, 0.182595638456};

    return test_kumaraswamy_distribution_aux(x, a, b, expectedLogP, dlogPdx,
                                             dlogPda, dlogPdb);
}

#pragma endregion

#pragma region Dirichlet Distribution

// Check the Dirichlet gradient wrt both the concentration alpha and the simplex x
// against finite differences of logP. alpha and x are plain (untransformed)
// parameters so the analytic gradients land directly in their leaves; GSL's
// dirichlet_lnpdf uses log(x_i) without renormalizing, so its FD matches the
// ambient gradient component-wise.
char* test_dirichlet_distribution() {
    double xVals[] = {0.2, 0.3, 0.5};
    Parameter* x = new_Parameter2("x", xVals, 3, new_Constraint(0, INFINITY));
    double aVals[] = {2.0, 3.0, 1.5};
    Parameter* alpha = new_Parameter2("alpha", aVals, 3, new_Constraint(0, INFINITY));

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);
    Parameters* parameters = new_Parameters(1);
    Parameters_move(parameters, alpha);

    DistributionModel* dm =
        new_DirichletDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = gsl_ran_dirichlet_lnpdf(3, aVals, xVals);
    printf("LogP: %f, LogP2: %f\n", logP, logP2);
    mu_assert(fabs(logP - logP2) < 1e-10, "dirichlet logP not matching");

    Parameters* ps = new_Parameters(2);
    Parameters_add(ps, alpha);
    Parameters_add(ps, x);
    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "dirichlet gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

#pragma endregion

#pragma region Multivariate Normal Distribution

// Independent oracle for the MVN log density with a lower-triangular Cholesky
// factor L (packed row-major), computed by hand:
//   logP = -k/2 log(2pi) - sum_i log(L_ii) - 1/2 || L^{-1}(x-mu) ||^2
static double mvn_logP_oracle(const double* mu, const double* Lpacked,
                              const double* x, size_t dim) {
    double* r = dvector(dim);
    double* w = dvector(dim);
    for (size_t i = 0; i < dim; i++) r[i] = x[i] - mu[i];
    // Forward substitution L w = r, reading L in packed lower-triangular order.
    double logdet = 0;
    for (size_t i = 0; i < dim; i++) {
        double s = r[i];
        size_t base = i * (i + 1) / 2;
        for (size_t j = 0; j < i; j++) s -= Lpacked[base + j] * w[j];
        double Lii = Lpacked[base + i];
        w[i] = s / Lii;
        logdet += log(Lii);
    }
    double quad = 0;
    for (size_t i = 0; i < dim; i++) quad += w[i] * w[i];
    free(r);
    free(w);
    return -0.5 * dim * (log(2.0) + log(M_PI)) - logdet - 0.5 * quad;
}

// logP against the hand-rolled oracle and every gradient (x, mu, packed L)
// against a finite difference of logP.
char* test_multivariate_normal_distribution() {
    size_t dim = 3;
    double muValues[] = {0.5, -0.2, 1.0};
    // L (lower-triangular Cholesky), packed row-major: rows [1.0][0.3,0.8][0.1,-0.2,0.5]
    double LValues[] = {1.0, 0.3, 0.8, 0.1, -0.2, 0.5};
    double xValues[] = {0.4, 0.1, 0.9};

    Parameter* mu =
        new_Parameter2("mu", muValues, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* L =
        new_Parameter2("L", LValues, 6, new_Constraint(-INFINITY, INFINITY));
    Parameter* x =
        new_Parameter2("x", xValues, dim, new_Constraint(-INFINITY, INFINITY));

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);
    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, mu);
    Parameters_move(parameters, L);

    DistributionModel* dm =
        new_MultivariateNormalDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = mvn_logP_oracle(muValues, LValues, xValues, dim);
    printf("MVN LogP: %f, oracle: %f\n", logP, logP2);
    mu_assert(fabs(logP - logP2) < 1e-10, "MVN logP not matching oracle");

    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, mu);
    Parameters_add(ps, L);
    Parameters_add(ps, x);
    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "MVN gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// A full covariance matrix (dim*dim) is Cholesky-decomposed internally; check logP
// against the oracle fed the corresponding L, and x/mu gradients against FD.
char* test_multivariate_normal_distribution_full_cov() {
    size_t dim = 2;
    double muValues[] = {0.5, -0.2};
    // Sigma = [[1.0, 0.3],[0.3, 0.5]] -> L = [[1.0],[0.3, sqrt(0.5-0.09)]]
    double sigmaValues[] = {1.0, 0.3, 0.3, 0.5};
    double xValues[] = {0.4, 0.1};
    double Lpacked[] = {1.0, 0.3, sqrt(0.5 - 0.09)};

    Parameter* mu =
        new_Parameter2("mu", muValues, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma =
        new_Parameter2("sigma", sigmaValues, 4, new_Constraint(-INFINITY, INFINITY));
    Parameter* x =
        new_Parameter2("x", xValues, dim, new_Constraint(-INFINITY, INFINITY));

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);
    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, mu);
    Parameters_move(parameters, sigma);

    DistributionModel* dm =
        new_MultivariateNormalDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = mvn_logP_oracle(muValues, Lpacked, xValues, dim);
    printf("MVN (full cov) LogP: %f, oracle: %f\n", logP, logP2);
    mu_assert(fabs(logP - logP2) < 1e-10, "MVN full-cov logP not matching oracle");

    // Only x and mu gradients are defined for the full-covariance parameterization.
    Parameters* ps = new_Parameters(2);
    Parameters_add(ps, mu);
    Parameters_add(ps, x);
    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "MVN full-cov gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// Reparameterization gradient: x = mu + L z. Given downstream dL/dx (in x->grad)
// and the standard-normal draw z (stashed in dm->tempx by rsample), we expect
//   dL/dmu_a  = dL/dx_a
//   dL/dL_ab  = dL/dx_a * z_b   (packed lower-triangular, a >= b)
char* test_multivariate_normal_rgradient() {
    size_t dim = 3;
    double muValues[] = {0.5, -0.2, 1.0};
    double LValues[] = {1.0, 0.3, 0.8, 0.1, -0.2, 0.5};
    double xValues[] = {0.0, 0.0, 0.0};

    Parameter* mu =
        new_Parameter2("mu", muValues, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* L =
        new_Parameter2("L", LValues, 6, new_Constraint(-INFINITY, INFINITY));
    Parameter* x =
        new_Parameter2("x", xValues, dim, new_Constraint(-INFINITY, INFINITY));

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, x);
    Parameters* parameters = new_Parameters(2);
    Parameters_move(parameters, mu);
    Parameters_move(parameters, L);

    DistributionModel* dm =
        new_MultivariateNormalDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double z[] = {0.7, -1.1, 0.3};
    double upstream[] = {0.9, -0.4, 1.2};  // downstream dL/dx
    for (size_t i = 0; i < dim; i++) {
        dm->tempx[i] = z[i];
        x->grad[i] = upstream[i];
    }
    Parameters_zero_grad(parameters);

    dm->rgradient(dm);

    for (size_t a = 0; a < dim; a++) {
        mu_assert(fabs(mu->grad[a] - upstream[a]) < 1e-12,
                  "MVN rgradient dmu does not match");
    }
    size_t idx = 0;
    for (size_t a = 0; a < dim; a++) {
        for (size_t b = 0; b <= a; b++) {
            double expected = upstream[a] * z[b];
            printf("dL/dL[%zu,%zu]: %f (expected %f)\n", a, b, L->grad[idx], expected);
            mu_assert(fabs(L->grad[idx] - expected) < 1e-12,
                      "MVN rgradient dL does not match");
            idx++;
        }
    }

    model->free(model);
    return NULL;
}

// Several observations sharing one MVN(mu, L): logP is the sum of per-observation
// densities, and mu/L gradients accumulate over observations while each x_k has its
// own gradient. All checked against finite differences of logP.
char* test_multivariate_normal_multi() {
    size_t dim = 3;
    double muValues[] = {0.5, -0.2, 1.0};
    double LValues[] = {1.0, 0.3, 0.8, 0.1, -0.2, 0.5};
    double x0v[] = {0.4, 0.1, 0.9};
    double x1v[] = {-0.3, 0.7, 0.2};

    Parameter* mu =
        new_Parameter2("mu", muValues, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* L =
        new_Parameter2("L", LValues, 6, new_Constraint(-INFINITY, INFINITY));
    Parameter* x0 =
        new_Parameter2("x0", x0v, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* x1 =
        new_Parameter2("x1", x1v, dim, new_Constraint(-INFINITY, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);
    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, mu);
    Parameters_add(parameters, L);

    DistributionModel* dm =
        new_MultivariateNormalDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = 0;
    for (size_t k = 0; k < Parameters_count(xs); k++) {
        logP2 += mvn_logP_oracle(muValues, LValues,
                                 Parameter_values(Parameters_at(xs, k)), dim);
    }
    printf("MVN multi LogP: %f, oracle: %f\n", logP, logP2);
    mu_assert(fabs(logP - logP2) < 1e-10, "MVN multi logP not matching oracle");

    Parameters* ps = new_Parameters(2 + Parameters_count(xs));
    Parameters_add(ps, mu);
    Parameters_add(ps, L);
    for (size_t k = 0; k < Parameters_count(xs); k++) {
        Parameters_add(ps, Parameters_at(xs, k));
    }
    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        printf("dlogP/d%s:", Parameter_name(p));
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double fd = fd_logP_grad(model, p, j);
            printf(" %f (fd %f)", p->grad[j], fd);
            mu_assert(fabs(p->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                      "MVN multi gradient does not match finite difference");
        }
        printf("\n");
    }

    free_Parameters(ps);
    model->free(model);
    return NULL;
}

// Reparameterization gradient with several observations: mu and L gradients sum
// over observations, each using its own stashed z_k.
char* test_multivariate_normal_multi_rgradient() {
    size_t dim = 3;
    double muValues[] = {0.5, -0.2, 1.0};
    double LValues[] = {1.0, 0.3, 0.8, 0.1, -0.2, 0.5};
    double zero[] = {0.0, 0.0, 0.0};

    Parameter* mu =
        new_Parameter2("mu", muValues, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* L =
        new_Parameter2("L", LValues, 6, new_Constraint(-INFINITY, INFINITY));
    Parameter* x0 =
        new_Parameter2("x0", zero, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* x1 =
        new_Parameter2("x1", zero, dim, new_Constraint(-INFINITY, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);
    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, mu);
    Parameters_add(parameters, L);

    DistributionModel* dm =
        new_MultivariateNormalDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double z[2][3] = {{0.7, -1.1, 0.3}, {0.2, 0.5, -0.9}};
    double up[2][3] = {{0.9, -0.4, 1.2}, {-0.5, 0.8, 0.1}};
    for (size_t k = 0; k < 2; k++) {
        Parameter* xk = Parameters_at(xs, k);
        for (size_t i = 0; i < dim; i++) {
            dm->tempx[k * dim + i] = z[k][i];
            xk->grad[i] = up[k][i];
        }
    }
    Parameters_zero_grad(parameters);

    dm->rgradient(dm);

    for (size_t a = 0; a < dim; a++) {
        double expected = up[0][a] + up[1][a];
        mu_assert(fabs(mu->grad[a] - expected) < 1e-12,
                  "MVN multi rgradient dmu does not match");
    }
    size_t idx = 0;
    for (size_t a = 0; a < dim; a++) {
        for (size_t b = 0; b <= a; b++) {
            double expected = up[0][a] * z[0][b] + up[1][a] * z[1][b];
            printf("dL/dL[%zu,%zu]: %f (expected %f)\n", a, b, L->grad[idx], expected);
            mu_assert(fabs(L->grad[idx] - expected) < 1e-12,
                      "MVN multi rgradient dL does not match");
            idx++;
        }
    }

    model->free(model);
    return NULL;
}

// Central finite-difference of the entropy wrt element i of parameter p.
static double fd_entropy_grad(DistributionModel* dm, Parameter* p, size_t i) {
    const double h = 1e-6;
    double v = Parameter_value_at(p, i);
    Parameter_set_value_at(p, v + h, i);
    double hp = dm->entropy(dm);
    Parameter_set_value_at(p, v - h, i);
    double hm = dm->entropy(dm);
    Parameter_set_value_at(p, v, i);
    return (hp - hm) / (2.0 * h);
}

// Entropy against the closed form k/2(1+log2pi) + sum_i log L_ii (times the number
// of observations), and its gradient wrt mu (zero) and L against finite differences.
char* test_multivariate_normal_entropy() {
    size_t dim = 3;
    double muValues[] = {0.5, -0.2, 1.0};
    double LValues[] = {1.0, 0.3, 0.8, 0.1, -0.2, 0.5};  // diag: 1.0, 0.8, 0.5
    double x0v[] = {0.4, 0.1, 0.9};
    double x1v[] = {-0.3, 0.7, 0.2};

    Parameter* mu =
        new_Parameter2("mu", muValues, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* L =
        new_Parameter2("L", LValues, 6, new_Constraint(-INFINITY, INFINITY));
    Parameter* x0 =
        new_Parameter2("x0", x0v, dim, new_Constraint(-INFINITY, INFINITY));
    Parameter* x1 =
        new_Parameter2("x1", x1v, dim, new_Constraint(-INFINITY, INFINITY));

    Parameters* xs = new_Parameters(2);
    Parameters_move(xs, x0);
    Parameters_move(xs, x1);
    Parameters* parameters = new_Parameters(2);
    Parameters_add(parameters, mu);
    Parameters_add(parameters, L);

    DistributionModel* dm =
        new_MultivariateNormalDistributionModel_with_parameters(parameters, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double H = dm->entropy(dm);
    double logdet = log(1.0) + log(0.8) + log(0.5);
    double Href = 2 * (0.5 * dim * (1.0 + log(2.0) + log(M_PI)) + logdet);
    printf("MVN entropy: %f, oracle: %f\n", H, Href);
    mu_assert(fabs(H - Href) < 1e-10, "MVN entropy not matching closed form");

    Parameters_zero_grad(parameters);
    dm->gradient_entropy(dm, parameters);

    // mu gradient is zero; L gradient matches finite differences of the entropy.
    for (size_t i = 0; i < dim; i++) {
        mu_assert(fabs(mu->grad[i]) < 1e-12, "MVN entropy dmu should be zero");
    }
    for (size_t j = 0; j < Parameter_size(L); j++) {
        double fd = fd_entropy_grad(dm, L, j);
        printf("dH/dL[%zu]: %f (fd %f)\n", j, L->grad[j], fd);
        mu_assert(fabs(L->grad[j] - fd) < 1e-4 * (1.0 + fabs(fd)),
                  "MVN entropy dL does not match finite difference");
    }

    model->free(model);
    return NULL;
}

#pragma endregion

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_exponential_distribution);
    mu_run_test(test_exponential_distribution_prior);

    mu_run_test(test_normal_distribution);
    mu_run_test(test_normal_distribution_prior);
    mu_run_test(test_normal_distribution_multi);
    mu_run_test(test_normal_distribution_multi_prior);
    mu_run_test(test_normal_distribution_multi_prior_transformed);

    mu_run_test(test_half_normal_distribution);
    mu_run_test(test_half_normal_distribution_prior);

    mu_run_test(test_lognormal_distribution);
    mu_run_test(test_lognormal_distribution_prior);
    mu_run_test(test_lognormal_distribution_prior_mean_stdev);
    mu_run_test(test_lognormal_distribution_multi);
    mu_run_test(test_lognormal_distribution_multi_prior);
    mu_run_test(test_lognormal_distribution_multi_prior_transformed);
    mu_run_test(test_lognormal_distribution_multi_prior_mean_stdev);
    mu_run_test(test_lognormal_distribution_multi_mean_stdev);

    mu_run_test(test_gamma_distribution_shape_rate);
    mu_run_test(test_gamma_distribution_shape_rate_prior);
    mu_run_test(test_gamma_distribution_shape_scale);
    mu_run_test(test_gamma_distribution_shape_scale_prior);
    mu_run_test(test_gamma_distribution_multi);
    mu_run_test(test_gamma_distribution_multi_prior);
    mu_run_test(test_gamma_distribution_multi_prior_transformed);
    mu_run_test(test_gamma_p_grad_a);
    mu_run_test(test_gamma_rgradient);

    mu_run_test(test_weibull_distribution);
    mu_run_test(test_weibull_distribution_prior);
    mu_run_test(test_weibull_distribution_multi);
    mu_run_test(test_weibull_distribution_multi_prior);
    mu_run_test(test_weibull_distribution_multi_prior_transformed);

    mu_run_test(test_kumaraswamy_distribution);
    mu_run_test(test_kumaraswamy_distribution_prior);

    mu_run_test(test_dirichlet_distribution);

    mu_run_test(test_multivariate_normal_distribution);
    mu_run_test(test_multivariate_normal_distribution_full_cov);
    mu_run_test(test_multivariate_normal_rgradient);
    mu_run_test(test_multivariate_normal_multi);
    mu_run_test(test_multivariate_normal_multi_rgradient);
    mu_run_test(test_multivariate_normal_entropy);

    return NULL;
}

RUN_TESTS(all_tests);
