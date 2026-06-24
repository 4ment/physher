//
//  test_compound.c
//  physher
//
//  Gradient tests for the mixture (weighted) compound model.
//

#include <math.h>

#include "minunit.h"
#include "phyc/compoundmodel.h"
#include "phyc/distmodel.h"
#include "phyc/distnormal.h"
#include "phyc/matrix.h"
#include "phyc/parameters.h"

// Build a 2-component normal mixture over the shared scalar parameter x:
//   logP = logsumexp_i( log w_i + log N(x; mu_i, sigma_i) )
// The returned model owns its components; the parameters passed in are shared
// (refcounted) so the caller can keep manipulating them.
static Model* build_normal_mixture(Parameter* x, Parameter* mu0, Parameter* sigma0,
                                   Parameter* mu1, Parameter* sigma1,
                                   Parameter* weights) {
    Parameters* xs0 = new_Parameters(1);
    Parameters_add(xs0, x);
    Parameters* par0 = new_Parameters(2);
    Parameters_add(par0, mu0);
    Parameters_add(par0, sigma0);
    DistributionModel* dm0 = new_NormalDistributionModel_with_parameters(
        par0, xs0, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* m0 = new_DistributionModel2("c0", dm0);

    Parameters* xs1 = new_Parameters(1);
    Parameters_add(xs1, x);
    Parameters* par1 = new_Parameters(2);
    Parameters_add(par1, mu1);
    Parameters_add(par1, sigma1);
    DistributionModel* dm1 = new_NormalDistributionModel_with_parameters(
        par1, xs1, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* m1 = new_DistributionModel2("c1", dm1);

    CompoundModel* cm = new_CompoundModel();
    cm->add(cm, m0);
    cm->add(cm, m1);
    cm->weights = weights;
    weights->refCount++;
    Model* model = new_CompoundModel2("mixture", cm);

    // cm->add and new_CompoundModel2 hold their own refs now.
    m0->free(m0);
    m1->free(m1);
    return model;
}

// Central finite-difference gradient of model->logP wrt every entry of ps.
static void finite_difference_gradient(Model* model, const Parameters* ps,
                                       double* out) {
    const double eps = 1e-6;
    size_t index = 0;
    for (size_t j = 0; j < Parameters_count(ps); j++) {
        Parameter* p = Parameters_at(ps, j);
        for (size_t k = 0; k < Parameter_size(p); k++) {
            double v = Parameter_value_at(p, k);
            Parameter_set_value_at(p, v + eps, k);
            double plus = model->logP(model);
            Parameter_set_value_at(p, v - eps, k);
            double minus = model->logP(model);
            Parameter_set_value_at(p, v, k);
            out[index++] = (plus - minus) / (2.0 * eps);
        }
    }
}

char* test_compound_mixture_gradient() {
    Parameter* x = new_Parameter("x", 0.7, new_Constraint(-INFINITY, INFINITY));
    Parameter* mu0 = new_Parameter("mu0", -1.0, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma0 = new_Parameter("sigma0", 0.8, new_Constraint(0, INFINITY));
    Parameter* mu1 = new_Parameter("mu1", 1.5, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma1 = new_Parameter("sigma1", 1.2, new_Constraint(0, INFINITY));

    double weightValues[] = {0.3, 0.7};
    Parameter* weights =
        new_Parameter2("weights", weightValues, 2, new_Constraint(0, INFINITY));

    Model* model =
        build_normal_mixture(x, mu0, sigma0, mu1, sigma1, weights);

    // logP must equal the analytic logsumexp of the two weighted components.
    double logN0 = -0.5 * log(2.0 * M_PI) - log(0.8) -
                   0.5 * pow((0.7 - (-1.0)) / 0.8, 2.0);
    double logN1 = -0.5 * log(2.0 * M_PI) - log(1.2) -
                   0.5 * pow((0.7 - 1.5) / 1.2, 2.0);
    double a = log(0.3) + logN0;
    double b = log(0.7) + logN1;
    double expected = (a > b ? a : b) + log1p(exp(-fabs(a - b)));
    double logP = model->logP(model);
    printf("logP: %f expected: %f\n", logP, expected);
    mu_assert(fabs(logP - expected) < 1e-9, "mixture logP not matching logsumexp");

    // logP/full_logP must be stored in the model's lp field before returning
    mu_assert(model->lp == logP, "CompoundModel logP not stored in lp");
    mu_assert(model->full_logP(model) == logP, "CompoundModel full_logP not matching");
    mu_assert(model->lp == logP, "CompoundModel full_logP not stored in lp");

    Parameters* ps = new_Parameters(6);
    Parameters_add(ps, x);
    Parameters_add(ps, mu0);
    Parameters_add(ps, mu1);
    Parameters_add(ps, sigma0);
    Parameters_add(ps, sigma1);
    Parameters_add(ps, weights);

    size_t size = Parameters_size(ps);
    double* fd = dvector(size);
    finite_difference_gradient(model, ps, fd);

    // Analytic gradient (caller zeroes first, per the convention).
    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    printf("Analytic vs finite-difference gradient:\n");
    size_t index = 0;
    for (size_t j = 0; j < Parameters_count(ps); j++) {
        Parameter* p = Parameters_at(ps, j);
        for (size_t k = 0; k < Parameter_size(p); k++) {
            printf("d/d%s[%zu]: %f vs %f\n", Parameter_name(p), k, p->grad[k],
                   fd[index]);
            mu_assert(fabs(p->grad[k] - fd[index]) < 1e-5,
                      "mixture gradient not matching finite differences");
            index++;
        }
    }

    free(fd);
    free_Parameters(ps);
    model->free(model);
    free_Parameter(x);
    free_Parameter(mu0);
    free_Parameter(sigma0);
    free_Parameter(mu1);
    free_Parameter(sigma1);
    free_Parameter(weights);
    return NULL;
}

// The mixture gradient must accumulate into ps (it is itself a submodel of the
// posterior compound, whose other submodels have already written into these
// buffers). It must not zero ps.
char* test_compound_mixture_gradient_accumulates() {
    Parameter* x = new_Parameter("x", 0.2, new_Constraint(-INFINITY, INFINITY));
    Parameter* mu0 = new_Parameter("mu0", 0.0, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma0 = new_Parameter("sigma0", 1.0, new_Constraint(0, INFINITY));
    Parameter* mu1 = new_Parameter("mu1", 2.0, new_Constraint(-INFINITY, INFINITY));
    Parameter* sigma1 = new_Parameter("sigma1", 0.5, new_Constraint(0, INFINITY));

    double weightValues[] = {0.4, 0.6};
    Parameter* weights =
        new_Parameter2("weights", weightValues, 2, new_Constraint(0, INFINITY));

    Model* model = build_normal_mixture(x, mu0, sigma0, mu1, sigma1, weights);

    Parameters* ps = new_Parameters(3);
    Parameters_add(ps, x);
    Parameters_add(ps, mu0);
    Parameters_add(ps, weights);

    size_t size = Parameters_size(ps);

    // Reference gradient from a clean (zeroed) state.
    Parameters_zero_grad(ps);
    model->gradient(model, ps);
    double* clean = dvector(size);
    size_t index = 0;
    for (size_t j = 0; j < Parameters_count(ps); j++) {
        Parameter* p = Parameters_at(ps, j);
        for (size_t k = 0; k < Parameter_size(p); k++) {
            clean[index++] = p->grad[k];
        }
    }

    // Seed a non-zero baseline (as a sibling model would have left behind) and
    // re-run: the result must be baseline + clean gradient, proving ps was not
    // zeroed.
    double* baseline = dvector(size);
    index = 0;
    for (size_t j = 0; j < Parameters_count(ps); j++) {
        Parameter* p = Parameters_at(ps, j);
        for (size_t k = 0; k < Parameter_size(p); k++) {
            baseline[index] = 0.5 + 0.1 * index;
            p->grad[k] = baseline[index];
            index++;
        }
    }
    model->gradient(model, ps);

    index = 0;
    for (size_t j = 0; j < Parameters_count(ps); j++) {
        Parameter* p = Parameters_at(ps, j);
        for (size_t k = 0; k < Parameter_size(p); k++) {
            printf("d/d%s[%zu]: %f (baseline %f + clean %f = %f)\n",
                   Parameter_name(p), k, p->grad[k], baseline[index], clean[index],
                   baseline[index] + clean[index]);
            mu_assert(fabs(p->grad[k] - (baseline[index] + clean[index])) < 1e-9,
                      "mixture gradient did not accumulate into ps");
            index++;
        }
    }

    free(clean);
    free(baseline);
    free_Parameters(ps);
    model->free(model);
    free_Parameter(x);
    free_Parameter(mu0);
    free_Parameter(sigma0);
    free_Parameter(mu1);
    free_Parameter(sigma1);
    free_Parameter(weights);
    return NULL;
}

// Same mixture, but the weights are a simplex parameter and the optimizer sees
// the unconstrained (stick-breaking) parameter. The weights gradient must be
// routed through transform->backward into the unconstrained parameter's grad
// (the weightsx != cm->weights branch).
char* test_compound_mixture_gradient_simplex() {
    const size_t K = 3;
    Parameter* x = new_Parameter("x", 0.3, new_Constraint(-INFINITY, INFINITY));
    double muVals[3] = {-1.0, 0.5, 2.0};
    double sigVals[3] = {0.7, 1.0, 0.5};
    Parameter* mu[3];
    Parameter* sigma[3];

    CompoundModel* cm = new_CompoundModel();
    for (size_t i = 0; i < K; i++) {
        char name[16];
        snprintf(name, 16, "mu%zu", i);
        mu[i] = new_Parameter(name, muVals[i], new_Constraint(-INFINITY, INFINITY));
        snprintf(name, 16, "sigma%zu", i);
        sigma[i] = new_Parameter(name, sigVals[i], new_Constraint(0, INFINITY));

        Parameters* xs = new_Parameters(1);
        Parameters_add(xs, x);
        Parameters* par = new_Parameters(2);
        Parameters_add(par, mu[i]);
        Parameters_add(par, sigma[i]);
        DistributionModel* dm = new_NormalDistributionModel_with_parameters(
            par, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
        Model* m = new_DistributionModel2(name, dm);
        cm->add(cm, m);
        m->free(m);
    }

    // Unconstrained stick-breaking parameter (dim K-1) feeding a simplex (dim K).
    double zVals[2] = {0.1, -0.4};
    Parameter* z = new_Parameter2("weights.unres", zVals, K - 1,
                                  new_Constraint(-INFINITY, INFINITY));
    Transform* transform = new_SimplexTransform_with_parameter(NULL, z);
    double simplexInit[3] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    Parameter* simplex =
        new_Parameter2("weights", simplexInit, K, new_Constraint(0.0, 1.0));
    simplex->transform = transform;
    transform->parameter->listeners->add_parameter(transform->parameter->listeners,
                                                    simplex);

    cm->weights = simplex;
    simplex->refCount++;
    Model* model = new_CompoundModel2("mixture", cm);

    // The optimizer holds the unconstrained parameter, not the simplex.
    mu_assert(Parameters_depends(new_Parameters(0), simplex) == NULL,
              "sanity: simplex not in empty set");

    Parameters* ps = new_Parameters(K + 2);
    Parameters_add(ps, x);
    for (size_t i = 0; i < K; i++) Parameters_add(ps, mu[i]);
    Parameters_add(ps, z);

    mu_assert(Parameters_depends(ps, simplex) == z,
              "weights must resolve to the unconstrained parameter");

    size_t size = Parameters_size(ps);
    double* fd = dvector(size);
    finite_difference_gradient(model, ps, fd);

    Parameters_zero_grad(ps);
    model->gradient(model, ps);

    printf("Analytic vs finite-difference gradient (simplex weights):\n");
    size_t index = 0;
    for (size_t j = 0; j < Parameters_count(ps); j++) {
        Parameter* p = Parameters_at(ps, j);
        for (size_t k = 0; k < Parameter_size(p); k++) {
            printf("d/d%s[%zu]: %f vs %f\n", Parameter_name(p), k, p->grad[k],
                   fd[index]);
            mu_assert(fabs(p->grad[k] - fd[index]) < 1e-5,
                      "simplex mixture gradient not matching finite differences");
            index++;
        }
    }

    free(fd);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_compound_mixture_gradient);
    mu_run_test(test_compound_mixture_gradient_accumulates);
    mu_run_test(test_compound_mixture_gradient_simplex);
    return NULL;
}

RUN_TESTS(all_tests);
