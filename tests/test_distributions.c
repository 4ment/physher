//
//  test_clone.c
//  physher
//  Created by Mathieu Fourment on 13/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
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

#include "phyc/distbeta.h"
#include "phyc/distcauchy.h"
#include "phyc/distdirichlet.h"
#include "phyc/distexp.h"
#include "phyc/distgamma.h"
#include "phyc/distkumaraswamy.h"
#include "phyc/distlognormal.h"
#include "phyc/distnormal.h"
#include "phyc/matrix.h"
#include "phyc/parameters.h"
#include "phyc/random.h"
#include "phyc/simplex.h"

Parameters** create_1_parameters_dim(double* p, size_t dim) {
    Parameters** ps = malloc(sizeof(Parameters*));
    ps[0] = new_Parameters_with_name("dim0", dim);
    for (int i = 0; i < dim; i++)
        Parameters_move(ps[0],
                        new_Parameter("first", p[i], new_Constraint(0, INFINITY)));
    return ps;
}

Parameters** create_2_parameters_dim(const double* first, const double* second,
                                     size_t dim) {
    Parameters** ps = malloc(sizeof(Parameters*) * 2);
    ps[0] = new_Parameters_with_name("dim0", dim);
    ps[1] = new_Parameters_with_name("dim1", dim);
    for (size_t i = 0; i < dim; i++) {
        Parameters_move(ps[0],
                        new_Parameter("first", first[i], new_Constraint(0, INFINITY)));
        Parameters_move(
            ps[1], new_Parameter("second", second[i], new_Constraint(0, INFINITY)));
    }
    return ps;
}

Parameters* create_x_dim(const double* x, size_t dim) {
    Parameters* xs = new_Parameters_with_name("x", dim);
    for (size_t i = 0; i < dim; i++)
        Parameters_move(xs, new_Parameter("x", x[i], new_Constraint(0, INFINITY)));
    return xs;
}

double gsl_1p_logP(double (*lpdf)(const double, const double), const double* xs,
                   size_t xdim, const double* p, size_t pdim) {
    double logP = 0;
    if (pdim == 1) {
        for (size_t i = 0; i < xdim; i++) {
            logP += log(lpdf(xs[i], p[0]));
        }
    } else {
        for (size_t i = 0; i < xdim; i++) {
            logP += log(lpdf(xs[i], p[i]));
        }
    }
    return logP;
}

double gsl_2p_logP(double (*lpdf)(const double, const double, const double),
                   const double* xs, size_t xdim, const double* p, size_t pdim) {
    double logP = 0;
    if (pdim == 1) {
        for (size_t i = 0; i < xdim; i++) {
            logP += log(lpdf(xs[i], p[0], p[1]));
        }
    } else {
        for (size_t i = 0; i < xdim; i++) {
            logP += log(lpdf(xs[i], p[i], p[i + pdim]));
        }
    }

    return logP;
}

double gsl_1p_dlogPdx(double (*lpdf)(const double, const double), const double* xs,
                      size_t xdim, size_t idx, const double* p, size_t pdim,
                      double eps) {
    double p0 = p[0];
    if (pdim > 1) {
        p0 = p[idx];
    }
    return (log(lpdf(xs[idx] + eps, p0)) - log(lpdf(xs[idx] - eps, p0))) / (2 * eps);
}

double gsl_1p_d2logPdx(double (*lpdf)(const double, const double), const double* xs,
                       size_t xdim, size_t idx, const double* p, size_t pdim,
                       double eps) {
    double p0 = p[0];
    if (pdim > 1) {
        p0 = p[idx];
    }

    double logP = log(lpdf(xs[idx], p0));
    double pp = log(lpdf(xs[idx] + eps, p0));
    double mm = log(lpdf(xs[idx] - eps, p0));
    return (pp + mm - 2 * logP) / (eps * eps);
}

double gsl_2p_dlogPdx(double (*lpdf)(const double, const double, const double),
                      const double* xs, size_t xdim, size_t idx, const double* p,
                      size_t pdim, double eps) {
    double p0 = p[0];
    double p1 = p[1];
    if (pdim > 1) {
        p0 = p[idx];
        p1 = p[idx + pdim];
    }
    return (log(lpdf(xs[idx] + eps, p0, p1)) - log(lpdf(xs[idx] - eps, p0, p1))) /
           (2 * eps);
}

double gsl_1p_dlogPdp(double (*lpdf)(const double, const double), const double* xs,
                      size_t xdim, const double* p, size_t pdim, size_t idx,
                      double eps) {
    double dlogP = 0;
    if (pdim > 1)
        return (log(lpdf(xs[idx], p[idx] + eps)) - log(lpdf(xs[idx], p[idx] - eps))) /
               (2 * eps);
    for (int i = 0; i < xdim; i++) {
        dlogP +=
            (log(lpdf(xs[i], p[0] + eps)) - log(lpdf(xs[idx], p[0] - eps))) / (2 * eps);
    }
    return dlogP;
}

double gsl_2p_dlogPdp(double (*lpdf)(const double, const double, const double),
                      const double* xs, size_t xdim, const double* p, size_t pdim,
                      size_t dim, size_t idx, double eps) {
    double dlogP = 0;
    double p0 = p[0];
    double p1 = p[1];
    if (dim == 0) {
        if (pdim > 1)
            return (log(lpdf(xs[idx], p[idx] + eps, p[idx + pdim])) -
                    log(lpdf(xs[idx], p[idx] - eps, p[idx + pdim]))) /
                   (2 * eps);

        for (int i = 0; i < xdim; i++) {
            dlogP += (log(lpdf(xs[i], p0 + eps, p1)) - log(lpdf(xs[i], p0 - eps, p1))) /
                     (2 * eps);
        }
    } else {
        if (pdim > 1)
            return (log(lpdf(xs[idx], p[idx], p[idx + pdim] + eps)) -
                    log(lpdf(xs[idx], p[idx], p[idx + pdim] - eps))) /
                   (2 * eps);

        for (int i = 0; i < xdim; i++) {
            dlogP += (log(lpdf(xs[i], p0, p1 + eps)) - log(lpdf(xs[i], p0, p1 - eps))) /
                     (2 * eps);
        }
    }
    return dlogP;
}

double gsl_gamma_logP(const double* xs, size_t xdim, const double* p, size_t pdim,
                      bool rate) {
    if (rate) {
        double* p2 = clone_dvector(p, 2 * pdim);
        for (int i = 0; i < pdim; i++) p2[i + pdim] = 1.0 / p[i + pdim];
        double logP = gsl_2p_logP(&gsl_ran_gamma_pdf, xs, xdim, p2, pdim);
        free(p2);
        return logP;
    }
    return gsl_2p_logP(&gsl_ran_gamma_pdf, xs, xdim, p, pdim);  // gsl shape scale
}

double gsl_gamma_dlogPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                         size_t pdim, double eps, bool rate) {
    if (rate) {
        double* p2 = clone_dvector(p, 2 * pdim);
        for (int i = 0; i < pdim; i++) p2[i + pdim] = 1.0 / p[i + pdim];
        double dlogP = gsl_2p_dlogPdx(&gsl_ran_gamma_pdf, xs, xdim, idx, p2, pdim, eps);
        free(p2);
        return dlogP;
    } else {
        return gsl_2p_dlogPdx(&gsl_ran_gamma_pdf, xs, xdim, idx, p, pdim, eps);
    }
}

double gsl_gamma_dlogPdp(const double* xs, size_t xdim, const double* p, size_t pdim,
                         size_t dim, size_t idx, double eps) {
    return gsl_2p_dlogPdp(gsl_ran_gamma_pdf, xs, xdim, p, pdim, dim, idx, eps);
}

double gsl_lognormal_logP(const double* xs, size_t xdim, const double* p, size_t pdim) {
    return gsl_2p_logP(&gsl_ran_lognormal_pdf, xs, xdim, p, pdim);
}

double gsl_lognormal_dlogPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                             size_t pdim, double eps) {
    return gsl_2p_dlogPdx(&gsl_ran_lognormal_pdf, xs, xdim, idx, p, pdim, eps);
}

double gsl_lognormal_dlogPdp(const double* xs, size_t xdim, const double* p,
                             size_t pdim, size_t dim, size_t idx, double eps) {
    return gsl_2p_dlogPdp(gsl_ran_lognormal_pdf, xs, xdim, p, pdim, dim, idx, eps);
}

double gsl_beta_logP(const double* xs, size_t xdim, const double* p, size_t pdim) {
    return gsl_2p_logP(&gsl_ran_beta_pdf, xs, xdim, p, pdim);
}

double gsl_beta_dlogPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                        size_t pdim, double eps) {
    return gsl_2p_dlogPdx(&gsl_ran_beta_pdf, xs, xdim, idx, p, pdim, eps);
}

double gsl_beta_dlogPdp(const double* xs, size_t xdim, const double* p, size_t pdim,
                        size_t dim, size_t idx, double eps) {
    return gsl_2p_dlogPdp(gsl_ran_beta_pdf, xs, xdim, p, pdim, dim, idx, eps);
}

double gsl_kumaraswamy_logP(const double* xs, size_t xdim, const double* p,
                            size_t pdim) {
    double logP = 0;
    if (pdim == 1) {
        for (size_t i = 0; i < xdim; i++) {
            logP += log(p[1] * p[0] * pow(xs[i], p[0] - 1.0) *
                        pow(1.0 - pow(xs[i], p[0]), p[1] - 1.0));
        }
    } else {
        for (size_t i = 0; i < xdim; i++) {
            logP += log(p[i + pdim] * p[i] * pow(xs[i], p[i] - 1.0) *
                        pow(1.0 - pow(xs[i], p[i]), p[i + pdim] - 1.0));
        }
    }

    return logP;
}

double f_kumaraswamy(double x, void* params) {
    double* p = params;
    return log(p[0] * p[1] * pow(x, p[0] - 1.0) * pow(1.0 - pow(x, p[0]), p[1] - 1.0));
}

double gsl_kumaraswamy_dlogPdx(const double* xs, size_t xdim, size_t idx,
                               const double* p, size_t pdim, double eps) {
    gsl_function F;
    double params[2] = {p[0], p[1]};
    F.function = &f_kumaraswamy;
    F.params = params;
    double result, abserr;
    gsl_deriv_central(&F, xs[idx], eps, &result, &abserr);
    return result;
}

double gsl_kumaraswamy_dlogPdp(const double* xs, size_t xdim, const double* p,
                               size_t pdim, size_t idx, double eps) {
    return 0;
}

double* center_xs(const double* xs, size_t xdim, const double* p, size_t pdim) {
    double* newxs = clone_dvector(xs, xdim);
    if (pdim == 1)
        for (int i = 0; i < xdim; i++) newxs[i] -= p[0];
    else
        for (int i = 0; i < xdim; i++) newxs[i] -= p[i];
    return newxs;
}

// MARK: Normal
double gsl_normal_logP(const double* xs, size_t xdim, const double* p, size_t pdim,
                       bool sigma) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double* newp = clone_dvector(p + pdim, pdim);
    double result;
    if (!sigma) {
        for (int i = 0; i < pdim; i++) newp[i] = sqrt(1.0 / p[pdim + i]);
        result = gsl_1p_logP(&gsl_ran_gaussian_pdf, newxs, xdim, newp, pdim);
    } else {
        result = gsl_1p_logP(&gsl_ran_gaussian_pdf, newxs, xdim, newp, pdim);
    }
    free(newp);
    free(newxs);
    return result;
}

double gsl_normal_dlogPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                          size_t pdim, double eps, bool sigma) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double* newp = clone_dvector(p + pdim, pdim);
    double dlogPdx;
    if (!sigma) {
        for (int i = 0; i < pdim; i++) newp[i] = sqrt(1.0 / p[pdim + i]);
        dlogPdx =
            gsl_1p_dlogPdx(&gsl_ran_gaussian_pdf, newxs, xdim, idx, newp, pdim, eps);
    } else {
        dlogPdx =
            gsl_1p_dlogPdx(&gsl_ran_gaussian_pdf, newxs, xdim, idx, newp, pdim, eps);
    }
    free(newxs);
    free(newp);
    return dlogPdx;
}

double gsl_normal_d2logPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                           size_t pdim, double eps, bool sigma) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double* newp = clone_dvector(p + pdim, pdim);
    double dlogPdx;
    if (!sigma) {
        for (int i = 0; i < pdim; i++) newp[i] = sqrt(1.0 / p[pdim + i]);
        dlogPdx =
            gsl_1p_d2logPdx(&gsl_ran_gaussian_pdf, newxs, xdim, idx, newp, pdim, eps);
    } else {
        dlogPdx =
            gsl_1p_d2logPdx(&gsl_ran_gaussian_pdf, newxs, xdim, idx, newp, pdim, eps);
    }
    free(newxs);
    free(newp);
    return dlogPdx;
}

double gsl_normal_dlogPdp(const double* xs, size_t xdim, const double* p, size_t pdim,
                          size_t idx, double eps, bool sigma) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double* newp = clone_dvector(p + pdim, pdim);
    double result;
    if (!sigma) {
        for (int i = 0; i < pdim; i++) newp[i] = sqrt(1.0 / p[pdim + i]);
        result = gsl_1p_dlogPdp(&gsl_ran_gaussian_pdf, xs, xdim, newp, pdim, idx, eps);
    } else {
        result = gsl_1p_dlogPdp(&gsl_ran_gaussian_pdf, xs, xdim, newp, pdim, idx, eps);
    }
    return result;
}

// MARK: HalfNormal
double gsl_halfnormal_logP(const double* xs, size_t xdim, const double* p, size_t pdim,
                           bool sigma) {
    return gsl_normal_logP(xs, xdim, p, pdim, sigma) / 2.0;
}

double gsl_halfnormal_dlogPdx(const double* xs, size_t xdim, size_t idx,
                              const double* p, size_t pdim, double eps, bool sigma) {
    return gsl_normal_dlogPdx(xs, xdim, idx, p, pdim, eps, sigma);
}

double gsl_halfnormal_d2logPdx(const double* xs, size_t xdim, size_t idx,
                               const double* p, size_t pdim, double eps, bool sigma) {
    return gsl_normal_d2logPdx(xs, xdim, idx, p, pdim, eps, sigma);
}

double gsl_halfnormal_dlogPdp(const double* xs, size_t xdim, const double* p,
                              size_t pdim, size_t idx, double eps, bool sigma) {
    return gsl_normal_dlogPdp(xs, xdim, p, pdim, idx, eps, sigma);
}

// MARK: Cauchy

double gsl_cauchy_logP(const double* xs, size_t xdim, const double* p, size_t pdim) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double result = gsl_1p_logP(&gsl_ran_cauchy_pdf, newxs, xdim, p + pdim, pdim);
    free(newxs);
    return result;
}

double gsl_cauchy_dlogPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                          size_t pdim, double eps) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double dlogPdx =
        gsl_1p_dlogPdx(&gsl_ran_cauchy_pdf, newxs, xdim, idx, p + pdim, pdim, eps);
    free(newxs);
    return dlogPdx;
}

double gsl_cauchy_d2logPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                           size_t pdim, double eps) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double dlogPdx =
        gsl_1p_d2logPdx(&gsl_ran_cauchy_pdf, newxs, xdim, idx, p + pdim, pdim, eps);
    free(newxs);
    return dlogPdx;
}

double gsl_cauchy_dlogPdp(const double* xs, size_t xdim, const double* p, size_t pdim,
                          size_t idx, double eps) {
    double* newxs = center_xs(xs, xdim, p, pdim);
    double dlogPd =
        gsl_1p_dlogPdp(&gsl_ran_cauchy_pdf, newxs, xdim, p + pdim, pdim, idx, eps);
    free(newxs);
    return dlogPd;
}

double gsl_exponential_logP(const double* xs, size_t xdim, const double* p, size_t pdim,
                            bool ismean) {
    double* newp = clone_dvector(p, pdim);
    if (!ismean) {
        for (int i = 0; i < pdim; i++) newp[i] = 1.0 / p[i];
    }
    double res = gsl_1p_logP(&gsl_ran_exponential_pdf, xs, xdim, newp, pdim);
    free(newp);
    return res;
}

double gsl_exponential_dlogPdx(const double* xs, size_t xdim, size_t idx,
                               const double* p, size_t pdim, double eps, bool ismean) {
    double* newp = clone_dvector(p, pdim);
    if (!ismean) {
        for (int i = 0; i < pdim; i++) newp[i] = 1.0 / p[i];
    }
    double res =
        gsl_1p_dlogPdx(&gsl_ran_exponential_pdf, xs, xdim, idx, newp, pdim, eps);
    free(newp);
    return res;
}

double gsl_exponential_d2logPdx(const double* xs, size_t xdim, size_t idx,
                                const double* p, size_t pdim, double eps, bool ismean) {
    double* newp = clone_dvector(p, pdim);
    if (!ismean) {
        for (int i = 0; i < pdim; i++) newp[i] = 1.0 / p[i];
    }
    double res =
        gsl_1p_d2logPdx(&gsl_ran_exponential_pdf, xs, xdim, idx, newp, pdim, eps);
    free(newp);
    return res;
}

double gsl_exponential_dlogPdp(const double* xs, size_t xdim, const double* p,
                               size_t pdim, size_t idx, double eps, bool ismean) {
    double* newp = clone_dvector(p, pdim);
    if (!ismean) {
        for (int i = 0; i < pdim; i++) newp[i] = 1.0 / p[i];
    }
    double res = gsl_1p_dlogPdp(&gsl_ran_cauchy_pdf, xs, xdim, newp, pdim, idx, eps);
    free(newp);
    return res;
}

void _inverse_transform_stick2(const double* parameters, double* values, int K) {
    size_t N = K - 1;
    double stick = 1.0;
    int k = 0;
    for (; k < N; k++) {
        values[k] = stick * parameters[k];
        stick -= values[k];
    }
    values[k] = stick;
}

double gsl_dirichlet_dlogPdx(const double* xs, size_t xdim, size_t idx, const double* p,
                             double eps) {
    double* x2 = dvector(xdim - 1);
    double stick = 1;
    for (int i = 0; i < xdim - 1; i++) {
        x2[i] = xs[i] / stick;
        stick -= xs[i];
    }
    double* xm = clone_dvector(x2, xdim - 1);
    xm[idx] -= eps;
    double* xp = clone_dvector(x2, xdim - 1);
    xp[idx] += eps;
    double* newxs = dvector(xdim);

    _inverse_transform_stick2(xm, newxs, xdim);
    double mm = gsl_ran_dirichlet_lnpdf(xdim, p, newxs);

    _inverse_transform_stick2(xp, newxs, xdim);
    double pp = gsl_ran_dirichlet_lnpdf(xdim, p, newxs);

    double result = (pp - mm) / (2 * eps);
    free(newxs);
    free(x2);
    free(xm);
    free(xp);
    return result;
}

char* test_gamma_distribution_dim(Parameters** ps, double* p, size_t pdim,
                                  bool israte) {
    size_t xdim = 2;
    double x[2] = {2, 4};
    Parameters* xs = create_x_dim(x, xdim);
    DistributionModel* dm = NULL;
    if (israte) {
        dm = new_GammaDistributionModel_with_parameters(ps, xs,
                                                        DISTRIBUTION_GAMMA_SHAPE_RATE);
    } else {
        dm = new_GammaDistributionModel_with_parameters(ps, xs,
                                                        DISTRIBUTION_GAMMA_SHAPE_SCALE);
    }
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = gsl_gamma_logP(x, xdim, p, pdim, israte);
    mu_assert(logP == logP2, "logP not matching");

    x[0] = 10;
    Parameters_set_value(xs, 0, x[0]);
    logP = model->logP(model);
    logP2 = gsl_gamma_logP(x, xdim, p, pdim, israte);
    mu_assert(logP == logP2, "logP after reset not matching");

    p[pdim] *= 2;  // change first beta
    Parameters_set_value(ps[1], 0, p[pdim]);
    logP = model->logP(model);
    logP2 = gsl_gamma_logP(x, xdim, p, pdim, israte);
    mu_assert(logP == logP2, "logP after changing p not matching");

    double dlogPdx = model->dlogP(model, Parameters_at(xs, 1));
    double eps = 0.00001;
    double dlogPdx2 = gsl_gamma_dlogPdx(x, xdim, 1, p, pdim, eps, israte);
    mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    model->free(model);
    return NULL;
}

char* test_gamma_distribution_scale_1() {
    size_t pdim = 1;
    double p[2] = {2, 0.5};
    Parameters** ps = create_2_parameters_dim(&p[0], &p[1], pdim);
    char* res = test_gamma_distribution_dim(ps, p, pdim, false);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return res;
}

char* test_gamma_distribution_scale_2() {
    size_t pdim = 2;
    double p[4] = {2, 0.5, 4, 5};  // first 2 are mus...
    Parameters** ps = create_2_parameters_dim(&p[0], &p[pdim], pdim);
    char* res = test_gamma_distribution_dim(ps, p, pdim, false);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return res;
}

char* test_gamma_distribution_rate_1() {
    size_t pdim = 1;
    double p[2] = {2, 0.5};
    Parameters** ps = create_2_parameters_dim(&p[0], &p[1], pdim);
    char* res = test_gamma_distribution_dim(ps, p, pdim, true);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return res;
}

char* test_gamma_distribution_rate_2() {
    size_t pdim = 2;
    double p[4] = {2, 0.5, 4, 5};  // first 2 are mus...
    Parameters** ps = create_2_parameters_dim(&p[0], &p[pdim], pdim);
    char* res = test_gamma_distribution_dim(ps, p, pdim, true);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return res;
}

char* test_lognormal_distribution_dim(Parameters** ps, double* p, size_t pdim) {
    size_t xdim = 2;
    double x[2] = {2, 4};
    Parameters* xs = create_x_dim(x, xdim);

    DistributionModel* dm = new_LogNormalDistributionModel_with_parameters(ps, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = gsl_lognormal_logP(x, xdim, p, pdim);
    mu_assert(logP == logP2, "logP not matching");

    x[0] = 10;
    Parameters_set_value(xs, 0, x[0]);
    logP = model->logP(model);
    logP2 = gsl_lognormal_logP(x, xdim, p, pdim);
    mu_assert(logP == logP2, "logP after changing x not matching");

    p[pdim] *= 2;  // change first sigma
    Parameters_set_value(ps[1], 0, p[pdim]);
    logP = model->logP(model);
    logP2 = gsl_lognormal_logP(x, xdim, p, pdim);
    mu_assert(logP == logP2, "logP after changing p not matching");

    double eps = 0.00001;
    double dlogPdx = model->dlogP(model, Parameters_at(xs, 1));
    double dlogPdx2 = gsl_lognormal_dlogPdx(x, xdim, 1, p, pdim, eps);
    mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    size_t pidx = pdim - 1;  // last p
    double dlogPdmu = model->dlogP(model, Parameters_at(ps[0], pidx));
    double dlogPdmu2 = gsl_lognormal_dlogPdp(x, xdim, p, pdim, 0, pidx, eps);
    mu_assert(fabs(dlogPdmu - dlogPdmu2) < 0.0001, "dlogPdmu not matching");

    double dlogPds = model->dlogP(model, Parameters_at(ps[1], pidx));
    double dlogPds2 = gsl_lognormal_dlogPdp(x, xdim, p, pdim, 1, pidx, eps);
    mu_assert(fabs(dlogPds - dlogPds2) < 0.0001, "dlogPdsigma not matching");

    model->free(model);
    return NULL;
}

char* test_lognormal_distribution_1() {
    size_t pdim = 1;
    double p[2] = {2, 0.5};
    Parameters** ps = create_2_parameters_dim(&p[0], &p[1], pdim);
    char* res = test_lognormal_distribution_dim(ps, p, pdim);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return res;
}

char* test_lognormal_distribution_2() {
    size_t pdim = 2;
    double p[4] = {2, 0.5, 4, 5};  // first 2 are mus...
    Parameters** ps = create_2_parameters_dim(&p[0], &p[pdim], pdim);
    char* res = test_lognormal_distribution_dim(ps, p, pdim);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return res;
}

char* test_cauchy_distribution() {
    size_t pdim = 1;
    double p[2] = {2, 0.5};
    Parameters** ps = create_2_parameters_dim(&p[0], &p[1], 1);

    size_t xdim = 2;
    double x[2] = {2, 4};
    Parameters* xs = create_x_dim(x, xdim);

    DistributionModel* dm = new_CauchyDistributionModel_with_parameters(ps, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = gsl_cauchy_logP(x, xdim, p, 1);
    mu_assert(logP == logP2, "logP not matching");

    x[0] = 10;
    Parameters_set_value(xs, 0, x[0]);
    logP = model->logP(model);
    logP2 = gsl_cauchy_logP(x, xdim, p, 1);
    mu_assert(logP == logP2, "logP after reset not matching");

    double eps = 0.00001;
    double dlogPdx = model->dlogP(model, Parameters_at(xs, 1));
    double dlogPdx2 = gsl_cauchy_dlogPdx(x, xdim, 1, p, pdim, eps);
    mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    double d2logPdx = model->d2logP(model, Parameters_at(xs, 1));
    double d2logPdx2 = gsl_cauchy_d2logPdx(x, xdim, 1, p, pdim, eps);
    mu_assert(fabs(d2logPdx - d2logPdx2) < 0.0001, "d2logPdx not matching");

    model->free(model);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return NULL;
}

char* test_beta_distribution() {
    size_t pdim = 1;
    double p[2] = {2, 0.5};
    Parameters** ps = create_2_parameters_dim(&p[0], &p[1], 1);

    size_t xdim = 2;
    double x[2] = {0.1, 0.2};
    Parameters* xs = create_x_dim(x, xdim);

    DistributionModel* dm = new_BetaDistributionModel_with_parameters(ps, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = gsl_beta_logP(x, xdim, p, 1);
    mu_assert(logP == logP2, "logP not matching");

    x[0] = 0.3;
    Parameters_set_value(xs, 0, x[0]);
    logP = model->logP(model);
    logP2 = gsl_beta_logP(x, xdim, p, 1);
    mu_assert(logP == logP2, "logP after reset not matching");

    // 	double eps = 0.00001;
    // 	double dlogPdx = model->dlogP(model, Parameters_at(xs, 1));
    // 	double dlogPdx2 = gsl_beta_dlogPdx(x, xdim, 1, p, pdim, eps);
    // 	mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    model->free(model);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return NULL;
}

char* test_kumaraswamy_distribution() {
    size_t pdim = 1;
    double p[2] = {2, 0.5};
    Parameters** ps = create_2_parameters_dim(&p[0], &p[1], 1);

    size_t xdim = 2;
    double x[2] = {0.1, 0.2};
    Parameters* xs = create_x_dim(x, xdim);

    DistributionModel* dm = new_KumaraswamyDistributionModel_with_parameters(ps, xs);
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = -3.886587;  // sum(dkumar(c(0.1,0.2), 2,0.5,log=T))
    mu_assert(fabs(logP - logP2) < .000001, "logP not matching");

    x[0] = 0.3;
    Parameters_set_value(xs, 0, x[0]);
    logP = model->logP(model);
    logP2 = -2.745844;  // sum(dkumar(c(0.3,0.2), 2,0.5,log=T))
    mu_assert(fabs(logP - logP2) < .000001, "logP after reset not matching");

    double eps = 0.00001;
    double dlogPdx = model->dlogP(model, Parameters_at(xs, 1));
    double dlogPdx2 = gsl_kumaraswamy_dlogPdx(x, xdim, 1, p, pdim, eps);
    mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    model->free(model);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return NULL;
}

char* test_exponential_distribution_dim(Parameters** ps, double* p, size_t pdim,
                                        bool ismean) {
    size_t xdim = 2;
    double x[2] = {2, 4};
    Parameters* xs = create_x_dim(x, xdim);

    DistributionModel* dm = NULL;
    if (ismean) {
        dm = new_ExponentialDistributionModel_with_parameters(
            ps, xs, DISTRIBUTION_EXPONENTIAL_MEAN);
    } else {
        dm = new_ExponentialDistributionModel_with_parameters(
            ps, xs, DISTRIBUTION_EXPONENTIAL_RATE);
    }
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = gsl_exponential_logP(x, xdim, p, pdim, ismean);
    mu_assert(fabs(logP - logP2) < 0.0001, "logP not matching");

    x[0] = 10;
    Parameters_set_value(xs, 0, x[0]);
    logP = model->logP(model);
    logP2 = gsl_exponential_logP(x, xdim, p, pdim, ismean);
    mu_assert(fabs(logP - logP2) < 0.0001, "logP after reset not matching");

    double eps = 0.00001;
    double dlogPdx = model->dlogP(model, Parameters_at(xs, 1));
    double dlogPdx2 = gsl_exponential_dlogPdx(x, xdim, 1, p, pdim, eps, ismean);
    mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    double d2logPdx = model->d2logP(model, Parameters_at(xs, 1));
    double d2logPdx2 = gsl_exponential_d2logPdx(x, xdim, 1, p, pdim, eps, ismean);
    mu_assert(fabs(d2logPdx - d2logPdx2) < 0.0001, "d2logPdx not matching");

    model->free(model);
    return NULL;
}

char* test_exponential_distribution_mean_1() {
    size_t pdim = 1;
    double p[1] = {2};
    Parameters** ps = create_1_parameters_dim(&p[0], 1);
    char* res = test_exponential_distribution_dim(ps, p, pdim, true);
    free_Parameters(ps[0]);
    free(ps);
    return res;
}

char* test_exponential_distribution_rate_1() {
    size_t pdim = 1;
    double p[1] = {2};
    Parameters** ps = create_1_parameters_dim(&p[0], 1);
    char* res = test_exponential_distribution_dim(ps, p, pdim, false);
    free_Parameters(ps[0]);
    free(ps);
    return res;
}

char* test_exponential_distribution_mean_2() {
    size_t pdim = 2;
    double p[2] = {2, 3};
    Parameters** ps = create_1_parameters_dim(p, pdim);
    char* res = test_exponential_distribution_dim(ps, p, pdim, true);
    free_Parameters(ps[0]);
    free(ps);
    return res;
}

char* test_exponential_distribution_rate_2() {
    size_t pdim = 2;
    double p[2] = {2, 3};
    Parameters** ps = create_1_parameters_dim(p, pdim);
    char* res = test_exponential_distribution_dim(ps, p, pdim, false);
    free_Parameters(ps[0]);
    free(ps);
    return res;
}

char* test_normal_distribution_issigma(bool issigma) {
    size_t pdim = 1;
    double p[2] = {2, 0.5};  // mu and sigma
    Parameters** ps = create_2_parameters_dim(&p[0], &p[1], 1);
    size_t xdim = 2;
    double x[2] = {2, 4};
    Parameters* xs = create_x_dim(x, xdim);

    DistributionModel* dm = NULL;
    if (issigma) {
        dm = new_NormalDistributionModel_with_parameters(
            ps, xs, DISTRIBUTION_NORMAL_MEAN_SIGMA);
    } else {
        dm = new_NormalDistributionModel_with_parameters(ps, xs,
                                                         DISTRIBUTION_NORMAL_MEAN_TAU);
    }
    Model* model = new_DistributionModel2("dist", dm);

    double logP = model->logP(model);
    double logP2 = gsl_normal_logP(x, xdim, p, pdim, issigma);
    mu_assert(logP == logP2, "logP not matching");

    x[0] = 10;
    Parameters_set_value(xs, 0, x[0]);
    logP = model->logP(model);
    logP2 = gsl_normal_logP(x, xdim, p, pdim, issigma);
    mu_assert(logP == logP2, "logP after reset not matching");

    double eps = 0.00001;
    size_t idx = 1;
    double dlogPdx = model->dlogP(model, Parameters_at(xs, idx));
    double dlogPdx2 = gsl_normal_dlogPdx(x, xdim, idx, p, pdim, eps, issigma);
    mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    double d2logPdx = model->d2logP(model, Parameters_at(xs, idx));
    double d2logPdx2 = gsl_normal_d2logPdx(x, xdim, idx, p, pdim, eps, issigma);
    mu_assert(fabs(d2logPdx - d2logPdx2) < 0.0001, "d2logPdx not matching");
    //
    // 	double dlogPdmu = model->dlogP(model, Parameters_at(ps[0], 0));
    // 	double dlogPdmu2 = gsl_normal_1dist_dlogPdp(xs, p, 0, eps, issigma);
    // 	mu_assert(fabs(dlogPdmu - dlogPdmu2) < 0.0001, "dlogPdmu not matching");
    //
    // 	double dlogPds = model->dlogP(model, Parameters_at(ps[1], 0));
    // 	double dlogPds2 = gsl_normal_1dist_dlogPdp(xs, p, 1, eps, issigma);
    // 	mu_assert(fabs(dlogPds - dlogPds2) < 0.0001, "dlogPdsigma not matching");

    model->free(model);
    for (size_t i = 0; i < 2; i++) free_Parameters(ps[i]);
    free(ps);
    return NULL;
}

char* test_normal_distribution_sigma() {
    return test_normal_distribution_issigma(true);
}

char* test_normal_distribution_tau() { return test_normal_distribution_issigma(false); }

char* test_dirichlet_distribution() {
    double p[3] = {3, 0.5, 1};
    Parameters** ps = create_1_parameters_dim(p, 3);

    size_t xdim = 3;
    double x[3] = {0.2, 0.3, 0.5};
    Simplex* simplex = new_Simplex_with_values("simplex", x, xdim);
    Model* msimplex = new_SimplexModel("simplex", simplex);

    DistributionModel* dm = new_DirichletDistributionModel_with_parameters(ps, simplex);
    Model* model = new_DistributionModel3("dist", dm, msimplex);
    msimplex->free(msimplex);

    double logP = model->logP(model);
    double logP2 = gsl_ran_dirichlet_lnpdf(xdim, p, x);
    mu_assert(fabs(logP - logP2) < .000001, "logP not matching");

    // change simplex
    x[0] = 0.3;
    x[1] = 0.2;
    simplex->set_values(simplex, x);
    logP = model->logP(model);
    logP2 = gsl_ran_dirichlet_lnpdf(xdim, p, x);
    mu_assert(fabs(logP - logP2) < .000001, "logP after changing simplex not matching");

    // change parameter of dirichlet
    p[1] = 10;
    Parameters_set_value(dm->parameters[0], 1, p[1]);
    logP = model->logP(model);
    logP2 = gsl_ran_dirichlet_lnpdf(xdim, p, x);
    mu_assert(fabs(logP - logP2) < .000001,
              "logP after changing dirichlet not matching");

    double eps = 0.00001;
    double dlogPdx = model->dlogP(model, Parameters_at(simplex->parameters, 1));
    double dlogPdx2 = gsl_dirichlet_dlogPdx(x, xdim, 1, p, eps);
    mu_assert(fabs(dlogPdx - dlogPdx2) < 0.0001, "dlogPdx not matching");

    model->free(model);
    free_Parameters(ps[0]);
    free(ps);
    return NULL;
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_gamma_distribution_rate_1);
    mu_run_test(test_gamma_distribution_scale_1);
    mu_run_test(test_gamma_distribution_rate_2);
    mu_run_test(test_gamma_distribution_scale_2);
    mu_run_test(test_lognormal_distribution_1);
    mu_run_test(test_lognormal_distribution_2);
    mu_run_test(test_normal_distribution_sigma);
    mu_run_test(test_normal_distribution_tau);
    mu_run_test(test_cauchy_distribution);
    mu_run_test(test_exponential_distribution_mean_1);
    mu_run_test(test_exponential_distribution_rate_1);
    mu_run_test(test_exponential_distribution_mean_2);
    mu_run_test(test_exponential_distribution_rate_2);
    mu_run_test(test_beta_distribution);
    mu_run_test(test_dirichlet_distribution);
    mu_run_test(test_kumaraswamy_distribution);

    return NULL;
}

RUN_TESTS(all_tests);
