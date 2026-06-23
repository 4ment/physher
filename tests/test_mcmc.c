#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "minunit.h"
#include "phyc/distmodel.h"
#include "phyc/distnormal.h"
#include "phyc/matrix.h"
#include "phyc/mcmc.h"

// Defined in mcmc.c / operator.c but not exposed in the public headers.
extern void run(MCMC* mcmc);
extern bool operator_slider(Operator* op, double* logHR);

// Build a standard-normal prior N(0,1) on a single random variable x.
static Model* make_normal_prior(Parameters* xs) {
    Parameters* params = new_Parameters(2);
    Parameters_move(params, new_Parameter("mu", 0.0, new_Constraint(-INFINITY, INFINITY)));
    Parameters_move(params, new_Parameter("sigma", 1.0, new_Constraint(0, INFINITY)));
    DistributionModel* dm = new_NormalDistributionModel_with_parameters(
        params, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    Model* model = new_DistributionModel2("prior", dm);
    free_Parameters(params);  // dm holds its own references
    return model;
}

// A proposal that corrupts the state and *then* reports failure, exercising the
// "propose returned false" path in run().
static bool propose_fail_after_mutate(Operator* op, double* logHR) {
    Parameter_set_value(Parameters_at(op->x, 0), 999.0);
    *logHR = 0;
    return false;
}

static MCMC make_mcmc(Model* model, Operator** ops, size_t op_count,
                      size_t length, gsl_rng* rng) {
    MCMC mcmc;
    memset(&mcmc, 0, sizeof(MCMC));
    mcmc.model = model;
    mcmc.operators = ops;
    mcmc.operator_count = op_count;
    mcmc.chain_length = length;
    mcmc.chain_temperature = -1;  // plain posterior path in _calculate_logP
    mcmc.log_count = 0;
    mcmc.logs = NULL;
    mcmc.verbose = 0;
    mcmc.interruptible = false;
    mcmc.generalized = false;
    mcmc.bf = false;
    mcmc.tuning_frequency = 1;
    mcmc.rng = rng;
    return mcmc;
}

// Fix #4: a failed proposal must roll back any partial mutation and be counted
// as a rejection, never accepted.
char* test_mcmc_failed_proposal_restores_state() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 42);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, new_Parameter("x", 0.5, new_Constraint(-INFINITY, INFINITY)));
    Parameter* x = Parameters_at(xs, 0);
    Model* model = make_normal_prior(xs);

    Operator* op = calloc(1, sizeof(Operator));
    op->weight = 1.0;
    op->x = xs;
    op->propose = propose_fail_after_mutate;
    op->optimize = NULL;
    op->rng = rng;

    Operator* ops[1] = {op};
    MCMC mcmc = make_mcmc(model, ops, 1, 10, rng);
    run(&mcmc);

    mu_assert(op->accepted_count == 0,
              "failed proposals must never be accepted");
    mu_assert(op->rejected_count == 10,
              "every failed proposal must be counted as one rejection");
    mu_assert(fabs(Parameter_value(x) - 0.5) < 1e-12,
              "state must be restored to its pre-proposal value");
    mu_assert(model->stored == false,
              "store flag must be cleared after the restore");

    free(op);
    model->free(model);
    free_Parameters(xs);
    gsl_rng_free(rng);
    return NULL;
}

// End-to-end: a symmetric slider targeting N(0,1) must produce a chain whose
// empirical mean and variance match the prior. Exercises the accept, reject,
// store, restore and Metropolis-ratio paths of run().
char* test_mcmc_samples_normal_prior() {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 7);

    Parameters* xs = new_Parameters(1);
    Parameters_move(xs, new_Parameter("x", 0.0, new_Constraint(-INFINITY, INFINITY)));
    Parameter* x = Parameters_at(xs, 0);
    Model* model = make_normal_prior(xs);

    Operator* op = calloc(1, sizeof(Operator));
    op->weight = 1.0;
    op->x = xs;
    op->propose = operator_slider;
    op->optimize = NULL;  // fixed proposal width for a reproducible target
    op->parameters = dvector(1);
    op->parameters[0] = 3.0;  // slider window
    op->all = false;
    op->rng = rng;

    Operator* ops[1] = {op};
    // One MH step per run() call keeps the model state between calls, letting
    // us read a sample after each step without wiring up a logger.
    MCMC mcmc = make_mcmc(model, ops, 1, 1, rng);

    const size_t burnin = 10000;
    const size_t nsamples = 40000;
    double sum = 0.0, sum2 = 0.0;
    size_t accepted = 0;
    for (size_t i = 0; i < burnin + nsamples; i++) {
        run(&mcmc);
        accepted += op->accepted_count;  // counts reset to 0 each run() call
        if (i >= burnin) {
            double v = Parameter_value(x);
            sum += v;
            sum2 += v * v;
        }
    }

    double mean = sum / nsamples;
    double var = sum2 / nsamples - mean * mean;
    double accept_rate = (double)accepted / (burnin + nsamples);

    mu_assert(accept_rate > 0.05 && accept_rate < 0.95,
              "chain should be mixing (non-degenerate acceptance rate)");
    mu_assert(fabs(mean) < 0.1, "empirical mean should match prior mean 0");
    mu_assert(fabs(var - 1.0) < 0.15, "empirical variance should match prior variance 1");

    free(op->parameters);
    free(op);
    model->free(model);
    free_Parameters(xs);
    gsl_rng_free(rng);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_mcmc_failed_proposal_restores_state);
    mu_run_test(test_mcmc_samples_normal_prior);
    return NULL;
}

RUN_TESTS(all_tests);
