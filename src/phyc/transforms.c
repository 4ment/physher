//
//  transforms.c
//  viphy
//
//  Created by Mathieu Fourment on 4/05/2017.
//  Copyright © 2017 University of Technology Sydney. All rights reserved.
//

#include "transforms.h"

#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "hessenberg.h"
#include "matrix.h"
#include "parameters.h"

double log1p_exp(double a) {
    // like log_sum_exp below with b=0.0; prevents underflow
    if (a > 0.0) return a + log1p(exp(-a));
    return log1p(exp(a));
}

double softplus(double x) { return log(1.0 + exp(x)); }

// logit
double logit(const double x) { return log(x / (1.0 - x)); }

// inverse logit
double inverse_logit(const double y) { return 1.0 / (1 + exp(-y)); }

double grad_inverse_logit(const double y) {
    const double exp_y = exp(-y);
    return exp_y / pow(exp_y + 1.0, 2);
}

double grad_log_det_inverse_logit(double y) {
    const double exp_y = exp(y);
    return (1.0 - exp_y) / (exp_y + 1.0);
}

double transform(double x, double lb, double ub) {
    // lower bound
    if (!isinf(lb) && isinf(ub)) {
        return log(x - lb);
    }
    // upper bound
    else if (!isinf(ub) && isinf(lb)) {
        return log(ub - x);
    } else if (!isinf(ub) && !isinf(lb)) {
        return logit((x - lb) / (ub - lb));
    }
    return x;
}

double transform2(double x, double lb, double ub, double* lp) {
    // lower bound
    if (!isinf(lb) && isinf(ub)) {
        *lp += -log(x - lb);
        return log(x - lb);
    }
    // upper bound
    else if (!isinf(ub) && isinf(lb)) {
        *lp += log(-1.0 / (ub - x));
        return log(ub - x);
    } else if (!isinf(ub) && !isinf(lb)) {
        *lp += log((ub - lb) / ((ub - x) * (x - lb)));
        return logit((x - lb) / (ub - lb));
    }
    return x;
}

double inverse_transform(double z, double lb, double ub, double* lp) {
    // lower bound
    /*if(lb == 0 && isinf(ub)){
        *lp += z - log(exp(z) + 1);
        return log(exp(z)+1);
    }*/
    if (!isinf(lb) && isinf(ub)) {
        *lp += z;
        return exp(z) + lb;
    }
    // upper bound
    else if (!isinf(ub) && isinf(lb)) {
        *lp += z;
        return ub - exp(z);
    }
    // no transform
    else if (isinf(ub) && isinf(lb)) {
        return z;
    }

    double inv_logit_x;
    if (z > 0) {
        double exp_minus_x = exp(-z);
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        // printf("%f %f\n", inv_logit_x, x);
        *lp += log(ub - lb) - z - 2 * log1p(exp_minus_x);
        // Prevent x from reaching one unless it really really should.
        if (!isinf(z) && inv_logit_x == 1) inv_logit_x = 1 - 1e-15;
    } else {
        double exp_x = exp(z);
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        *lp += log(ub - lb) + z - 2 * log1p(exp_x);
        // Prevent x from reaching zero unless it really really should.
        if (!isinf(z) && inv_logit_x == 0) inv_logit_x = 1e-15;
    }
    return lb + (ub - lb) * inv_logit_x;
}

double inverse_transform2(double z, double lb, double ub) {
    // lower bound
    /*if(lb == 0 && isinf(ub)){
        return log(exp(z)+1);
    }*/
    if (!isinf(lb) && isinf(ub)) {
        return exp(z) + lb;
    }
    // upper bound
    else if (!isinf(ub) && isinf(lb)) {
        return ub - exp(z);
    }
    // no transform
    else if (isinf(ub) && isinf(lb)) {
        return z;
    }

    double inv_logit_x;
    if (z > 0) {
        double exp_minus_x = exp(-z);
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        // Prevent x from reaching one unless it really really should.
        if (!isinf(z) && inv_logit_x == 1) inv_logit_x = 1 - 1e-15;
    } else {
        double exp_x = exp(z);
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        // printf("%f %f\n", inv_logit_x, x);
        //  Prevent x from reaching zero unless it really really should.
        if (!isinf(z) && inv_logit_x == 0) inv_logit_x = 1e-15;
    }
    return lb + (ub - lb) * inv_logit_x;
}

double grad_inverse_transform(double z, double lb, double ub) {
    // lower bound
    if (!isinf(lb) && isinf(ub)) {
        return exp(z);
    }
    // upper bound
    else if (!isinf(ub) && isinf(lb)) {
        return -exp(z);
    }
    // no transform
    else if (isinf(ub) && isinf(lb)) {
        return 1.0;
    }

    double inv_logit_x;
    if (z > 0) {
        double exp_minus_x = exp(-z);
        inv_logit_x = exp_minus_x / pow((1.0 + exp_minus_x), 2);
    } else {
        double exp_x = exp(z);
        inv_logit_x = exp_x / pow((1.0 + exp_x), 2);
    }
    return (ub - lb) * inv_logit_x;
}

// Calculate gradient of log determinant of Jacobian
// the det should always be positive because all the terms are positive so no sign(det)
// term
double grad_log_det_inverse_transform(double z, double lb, double ub) {
    if (lb == 0.0 && isinf(ub)) {
        return 1.0;
    }
    // no transform
    else if (isinf(ub) && isinf(lb)) {
        return 0;
    }

    double exp_x = exp(z);
    return (1.0 - exp_x) / (1.0 + exp_x);
}

// Simplex

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_simplex2(double* x, const double* y, int K) {
    double stick = 1.0;
    int k = 0;
    for (; k < K - 1; k++) {
        double z = inverse_logit(y[k] - log(K - k - 1));
        x[k] = stick * z;
        stick -= x[k];
    }
    x[k] = stick;
}

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_simplex(double* x, const double* y, int K, double* lp) {
    double stick = 1.0;
    int k = 0;
    for (; k < K - 1; k++) {
        double adj_y_k = y[k] - log(K - k - 1);
        double z = inverse_logit(adj_y_k);
        x[k] = stick * z;
        *lp += log(stick);
        *lp -= log1p_exp(-adj_y_k);
        *lp -= log1p_exp(adj_y_k);
        stick -= x[k];
    }
    x[k] = stick;
}

void grad_inverse_transform_simplex(double* dz, const double* z, int K) {}

/*void grad_log_det_inverse_transform_simplex(const double* x, double* dz, int K){
        //dx[0] = 1.0/(x[0]-1.0) + 1.0/(x[0]+x[1]-1.0);
        //dx[1] = 1.0/(x[0]+x[1]-1.0);
        dz[K-1] = 0;
        for(int k = K-2; k >=0; --k){
                dz[k] = -1.0;
                for(int i = 0; i <= k; --i){
                        dz[k] += x[k];
                }
                dz[k] = 1.0/dz[k] + dz[k+1];
        }
}*/

void grad_log_det_inverse_transform_simplex(const double* z, double* dz, int K) {
    for (int k = 0; k < K - 1; k++) {
        dz[k] = 1.0 / z[k] + 1.0 / (z[k] - 1.0);
    }
}

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_msimplex(double* x, const double* y, int K) {
    double sum = 0;
    int k = 0;
    for (; k < K - 1; k++) {
        x[k] = exp(y[k]);
        sum += x[k];
    }
    x[k] = 1.0 / (1.0 + sum);
    for (int i = 0; i < K - 1; i++) {
        x[i] /= x[k];
    }
}

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_msimplex2(double* x, const double* y, int K, double* lp) {
    double sum = 0;
    double sum_one = 1.0;
    int dim = K - 1;
    int k = 0;
    for (; k < dim; k++) {
        x[k] = exp(y[k]);
        sum_one += x[k];
    }
    x[k] = 1.0 / (sum_one);
    for (int i = 0; i < dim; i++) {
        x[i] /= x[k];
    }

    double inv_sum_one_squared = 1.0 / (sum_one * sum_one);
    double* m = dvector(dim * dim);
    k = 0;
    for (int i = 0; i < dim; i++) {
        double temp = 1;
        int j = 0;
        for (; j < i; j++) temp += x[j];
        j++;
        for (; j < dim; j++) temp += x[j];
        m[dim * i + i] = x[i] * temp / inv_sum_one_squared;

        for (int j = i + 1; j < dim; j++) {
            m[dim * i + j] = -exp(y[i] + y[j]) / inv_sum_one_squared;
            m[dim * j + i] = m[dim * i + j];
        }
    }
    *lp = log(LUDecompose_det2(m, dim));
    free(m);
}

void _transform_simplex_stan(const double* x, double* y, size_t dim, double lower,
                             double upper) {
    double sum = 0;
    for (size_t i = 0; i < dim; i++) {
        double z = x[i] / (1.0 - sum);
        y[i] = logit(z) + log(dim - i);
        sum += x[i];
    }
}

// dim=K-1
void _inverse_transform_simplex_stan(double* x, const double* y, size_t dim,
                                     double lower, double upper) {
    double stick = 1.0;
    for (size_t k = 0; k < dim; k++) {
        x[k] = stick * inverse_logit(y[k] - log(dim - k));
        stick -= x[k];
    }
    x[dim] = stick;
}

// dL/dy += dL/dx * dx/dy
// static void _backward_inverse_transform_exp(double* grad, const double* y, const
// double* ingrad, size_t dim, double lower, double upper){

// dL/dy += dL/dx * dx/dy
// ingrad: dim(ingrad) = dim+1
// grad: dim(grad) = dim x dim+1
static void _backward_inverse_transform_simplex_stan(double* grad, const double* y,
                                                     const double* ingrad, size_t dim,
                                                     double lower, double upper) {
    double stickRemaining = 1.0;
    // size_t j = 0;
    for (size_t index = 0; index < dim; index++) {
        double g = stickRemaining * grad_inverse_logit(y[index] - log(dim - index));
        grad[index] += ingrad[index] * g;
        double cum = g;
        for (size_t i = index + 1; i < dim; i++) {
            g = -cum * inverse_logit(y[i] - log(dim - i));
            grad[index] += ingrad[i] * g;
            cum += g;
        }
        grad[index] += ingrad[dim] * -cum;
        double x = stickRemaining * inverse_logit(y[index] - log(dim - index));
        stickRemaining -= x;
    }
}

// dx/dy
// x: constrained values with dim(x) = dim+1
// y: unconstrained values with dim(y) = dim
// jacobian: dim(jacobian) = dim x dim+1
static void _inverse_transform_simplex_jacobian_stan(double* jacobian, const double* y,
                                                     size_t dim, double lower,
                                                     double upper) {
    memset(jacobian, 0.0, sizeof(double) * dim * (dim + 1));
    double stickRemaining = 1.0;
    size_t j = 0;
    for (size_t index = 0; index < dim; index++, j++) {
        j += index;
        jacobian[j] = stickRemaining * grad_inverse_logit(y[index] - log(dim - index));
        double cum = jacobian[j];
        j++;
        for (size_t i = index + 1; i < dim; i++, j++) {
            jacobian[j] = -cum * inverse_logit(y[i] - log(dim - i));
            cum += jacobian[j];
        }
        jacobian[j] = -cum;
        double x = stickRemaining * inverse_logit(y[index] - log(dim - index));
        stickRemaining -= x;
    }
}

// log(dT(y)/dy) = log(dx/dy) = log(exp(y)) = y
static double _inverse_transform_log_jacobian_simplex_stan(double* res, const double* y,
                                                           size_t dim, double lower,
                                                           double upper) {
    if (res != NULL) {
        error("_inverse_transform_log_jacobian_simplex_stan not implemented");
        memcpy(res, y, sizeof(double) * dim);
    }
    double sumLogDetJac = 0;
    double stickRemaining = 1.0;
    for (size_t i = 0; i < dim; i++) {
        double z = inverse_logit(y[i] - log(dim - i));
        sumLogDetJac += log(z * (1.0 - z) * stickRemaining);
        stickRemaining -= stickRemaining * z;
    }
    return sumLogDetJac;
}

// d log(dT(y)/dy)/dy
static double _inverse_transform_gradient_log_jacobian_simplex_stan(
    double* res, const double* y, size_t dim, double lower, double upper) {
    double lp = 0;
    double stickRemaining = 1.0;
    size_t K = dim + 1;
    double* jacobian = dvector(dim * K);

    _inverse_transform_simplex_jacobian_stan(jacobian, y, dim, lower, upper);
    if (res != NULL) {
        for (size_t i = 0; i < dim; i++) {
            double z = inverse_logit(y[i] - log(dim - i));
            res[i] += z * (1.0 - z) / z - (z * (1.0 - z)) / (1.0 - z);
            for (size_t j = 0; j < i; j++) {
                for (size_t k = 0; k < i; k++) {
                    res[j] += -1.0 / stickRemaining * jacobian[j * K + k];
                }
            }
            stickRemaining -= stickRemaining * z;
        }
    } else {
        for (size_t i = 0; i < dim; i++) {
            lp += 1.0 / y[i] + 1.0 / (y[i] - 1.0) + 1.0 - 2.0 * inverse_logit(y[i]);
        }
    }
    free(jacobian);
    return lp;
}

// Lower bounded

// y = T^{-1}(x)
static void _transform_exp(const double* x, double* y, size_t dim, double lower,
                           double upper) {
    for (size_t i = 0; i < dim; i++) {
        y[i] = log(x[i] - lower);
    }
}

// x = T(y)
static void _inverse_transform_exp(double* x, const double* y, size_t dim, double lower,
                                   double upper) {
    for (size_t i = 0; i < dim; i++) {
        x[i] = exp(y[i]) + lower;
    }
}

// dL/dy += dL/dx * dx/dy
static void _backward_inverse_transform_exp(double* grad, const double* y,
                                            const double* ingrad, size_t dim,
                                            double lower, double upper) {
    for (size_t i = 0; i < dim; i++) {
        grad[i] += ingrad[i] * exp(y[i]);
    }
}

// dx/dy
static void _inverse_transform_jacobian_exp(double* grad, const double* y, size_t dim,
                                            double lower, double upper) {
    for (size_t i = 0; i < dim; i++) {
        grad[i] = exp(y[i]);
    }
}

// log(dT(y)/dy) = log(dx/dy) = log(exp(y)) = y
static double _inverse_transform_log_jacobian_exp(double* res, const double* y,
                                                  size_t dim, double lower,
                                                  double upper) {
    if (res != NULL) {
        memcpy(res, y, sizeof(double) * dim);
    }
    double sumLogJac = 0;
    for (size_t i = 0; i < dim; i++) {
        sumLogJac += y[i];
    }
    return sumLogJac;
}

// d log(dT(y)/dy)/dy = dy/dy = 1
static double _inverse_transform_gradient_log_jacobian_exp(double* res, const double* y,
                                                           size_t dim, double lower,
                                                           double upper) {
    if (res != NULL) {
        for (size_t i = 0; i < dim; i++) {
            res[i] += 1.0;
        }
    }
    return (double)dim;
}

// [0,1]

// y = T^{-1}(x)
static void _transform_logit(const double* x, double* y, size_t dim, double lower,
                             double upper) {
    for (size_t i = 0; i < dim; i++) {
        y[i] = logit(x[i]);
    }
}

// x = T(y)
static void _inverse_transform_logit(double* x, const double* y, size_t dim, double lb,
                                     double ub) {
    for (size_t i = 0; i < dim; i++) {
        x[i] = inverse_logit(y[i]);
    }
}

// dL/dy = dL/dx * dx/dy
static void _backward_inverse_transform_logit(double* grad, const double* y,
                                              const double* ingrad, size_t dim,
                                              double lb, double ub) {
    for (size_t i = 0; i < dim; i++) {
        grad[i] += ingrad[i] * inverse_logit(y[i]) * (1.0 - inverse_logit(y[i]));
    }
}

// dx/y = logit^{-1}(y) (1 - logit^{-1}(y))
static double _inverse_transform_log_jacobian_logit(double* res, const double* y,
                                                    size_t dim, double lb, double ub) {
    if (res != NULL) {
        error("_inverse_transform_log_jacobian_logit not implemented");
    }
    double sumLogJac = 0;
    for (size_t i = 0; i < dim; i++) {
        double negAbsY = -fabs(y[i]);
        sumLogJac += negAbsY - 2.0 * log1p(exp(negAbsY));
    }
    return sumLogJac;
}

// d log(dx/dy)/dy = 1 - 2 * logit^1(y)
static double _inverse_transform_gradient_log_jacobian_logit(double* res,
                                                             const double* y,
                                                             size_t dim, double lower,
                                                             double upper) {
    double lp = 0;
    if (res != NULL) {
        double temp;
        for (size_t i = 0; i < dim; i++) {
            temp = 1.0 - 2.0 * inverse_logit(y[i]);
            res[i] += temp;
            lp += temp;
        }
    } else {
        for (size_t i = 0; i < dim; i++) {
            lp += 1.0 - 2.0 * inverse_logit(y[i]);
        }
    }
    return lp;
}

// Lower and upper bounded

// y = T^{-1}(x)
static void _transform_bounded(const double* x, double* y, size_t dim, double lower,
                               double upper) {
    for (size_t i = 0; i < dim; i++) {
        y[i] = logit((x[i] - lower) / (upper - lower));
    }
}

// x = T(y)
static void _inverse_transform_bounded(double* x, const double* y, size_t dim,
                                       double lb, double ub) {
    double diff = ub - lb;
    for (size_t i = 0; i < dim; i++) {
        double inv_logit_x;
        if (y[i] > 0) {
            double exp_minus_x = exp(-y[i]);
            inv_logit_x = 1.0 / (1.0 + exp_minus_x);
            // Prevent x from reaching one unless it really really should.
            if (!isinf(y[i]) && inv_logit_x == 1) inv_logit_x = 1 - 1e-15;
        } else {
            double exp_x = exp(y[i]);
            inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
            // Prevent x from reaching zero unless it really really should.
            if (!isinf(y[i]) && inv_logit_x == 0) inv_logit_x = 1e-15;
        }
        x[i] = lb + diff * inv_logit_x;
    }
}

// dL/dy = dL/dx * dx/dy
static void _backward_inverse_transform_bounded(double* grad, const double* y,
                                                const double* ingrad, size_t dim,
                                                double lb, double ub) {
    double deltaBound = ub - lb;
    for (size_t i = 0; i < dim; i++) {
        // double inv_logit_x;
        // if (y[i] > 0) {
        //     double exp_minus_x = exp(-y[i]);
        //     inv_logit_x = exp_minus_x / pow((1.0 + exp_minus_x), 2);
        // } else {
        //     double exp_x = exp(y[i]);
        //     inv_logit_x = exp_x / pow((1.0 + exp_x), 2);
        // }
        // grad[i] = ingrad[i] *  deltaBound * inv_logit_x;
        grad[i] +=
            ingrad[i] * deltaBound * inverse_logit(y[i]) * (1.0 - inverse_logit(y[i]));
    }
}

// dx/y = (ub - lb) logit^{-1}(y) (1 - logit^{-1}(y))
// log(dx/dy) = log(ub - lb) -|y| - 2 * log(1 + exp(-|y|))
static double _inverse_transform_log_jacobian_bounded(double* res, const double* y,
                                                      size_t dim, double lb,
                                                      double ub) {
    if (res != NULL) {
        error("_inverse_transform_log_jacobian_bounded not implemented");
    }
    double sumLogJac = 0;
    double logDiff = log(ub - lb);
    for (size_t i = 0; i < dim; i++) {
        double negAbsY = -fabs(y[i]);
        sumLogJac += logDiff + negAbsY - 2.0 * log1p(exp(negAbsY));
    }
    return sumLogJac;
}

// d log(dx/dy)/dy = 1 - 2 * logit^1(y)
static double _inverse_transform_gradient_log_jacobian_bounded(double* res,
                                                               const double* y,
                                                               size_t dim, double lower,
                                                               double upper) {
    double lp = 0;
    if (res != NULL) {
        double temp = 0;
        for (size_t i = 0; i < dim; i++) {
            temp = 1.0 - 2.0 * inverse_logit(y[i]);
            res[i] += temp;
            lp += temp;
        }
    } else {
        for (size_t i = 0; i < dim; i++) {
            lp += 1.0 - 2.0 * inverse_logit(y[i]);
            // double exp_x = exp(y[i]);
            // lp += (1.0-exp_x)/(1.0+exp_x);
        }
    }
    return lp;
}

// Softplus

static void _transform_softplus(const double* x, double* y, size_t dim, double lower,
                                double upper) {
    for (size_t i = 0; i < dim; i++) {
        y[i] = log(exp(x[i]) - 1.0);
    }
}

static void _inverse_transform_softplus(double* x, const double* y, size_t dim,
                                        double lower, double upper) {
    for (size_t i = 0; i < dim; i++) {
        x[i] = log(1.0 + exp(y[i]));
    }
}

static void _get(Transform* self, double* x) {
    // bottom-up update
    if (self->parameter->transform != NULL) {
        double* y = Parameter_values(self->parameter->transform->parameter);
        self->parameter->transform->get(self->parameter->transform, y);
    }
    const double* y = Parameter_values(self->parameter);
    self->inverse_transform(x, y, Parameter_size(self->parameter), self->lower,
                            self->upper);
}

static void _set(Transform* self, const double* x) {
    // top-bottom update
    double* y = Parameter_values(self->parameter);
    self->transform(x, y, Parameter_size(self->parameter), self->lower, self->upper);
    if (self->parameter->transform != NULL) {
        self->parameter->transform->set(self->parameter->transform, y);
    }
}

static double _log_det_jacobian(Transform* self) {
    const double* y = Parameter_values(self->parameter);
    return self->inverse_transform_log_det_jacobian(
        NULL, y, Parameter_size(self->parameter), self->lower, self->upper);
}

static double _gradient_log_det_jacobian(Transform* self) {
    const double* y = Parameter_values(self->parameter);
    return self->inverse_transform_gradient_log_det_jacobian(
        self->parameter->grad, y, Parameter_size(self->parameter), self->lower,
        self->upper);
}

static void _backward(Transform* self, const double* ingrad) {
    const double* y = Parameter_values(self->parameter);
    self->backward_inverse_transform(self->parameter->grad, y, ingrad,
                                     Parameter_size(self->parameter), self->lower,
                                     self->upper);
}

static void _jacobian(Transform* self, double* jacobian) {
    const double* y = Parameter_values(self->parameter);
    self->inverse_transform_jacobian(jacobian, y, Parameter_size(self->parameter),
                                     self->lower, self->upper);
}

Transform* new_Transform_with_parameter(const char* type, double lower, double upper,
                                        Parameter* parameter) {
    Transform* transform = (Transform*)malloc(sizeof(Transform));
    transform->lower = lower;
    transform->upper = upper;
    transform->parameter = parameter;

    transform->get = _get;
    transform->set = _set;
    transform->backward = _backward;
    transform->log_det_jacobian = _log_det_jacobian;
    transform->gradient_log_det_jacobian = _gradient_log_det_jacobian;
    transform->jacobian = _jacobian;

    transform->inverse_transform_jacobian = NULL;

    if (!isinf(lower) && isinf(upper)) {
        transform->transform = _transform_exp;
        transform->inverse_transform = _inverse_transform_exp;
        transform->backward_inverse_transform = _backward_inverse_transform_exp;
        transform->inverse_transform_jacobian = _inverse_transform_jacobian_exp;
        transform->inverse_transform_log_det_jacobian =
            _inverse_transform_log_jacobian_exp;
        transform->inverse_transform_gradient_log_det_jacobian =
            _inverse_transform_gradient_log_jacobian_exp;
        transform->dim = Parameter_size(parameter);
    } else if (lower == 0.0 && upper == 1.0) {
        transform->transform = _transform_logit;
        transform->inverse_transform = _inverse_transform_logit;
        transform->backward_inverse_transform = _backward_inverse_transform_logit;
        transform->inverse_transform_log_det_jacobian =
            _inverse_transform_log_jacobian_logit;
        transform->inverse_transform_gradient_log_det_jacobian =
            _inverse_transform_gradient_log_jacobian_logit;
        transform->dim = Parameter_size(parameter);
    } else if (!isinf(lower) && !isinf(upper)) {
        transform->transform = _transform_bounded;
        transform->inverse_transform = _inverse_transform_bounded;
        transform->backward_inverse_transform = _backward_inverse_transform_bounded;
        transform->inverse_transform_log_det_jacobian =
            _inverse_transform_log_jacobian_bounded;
        transform->inverse_transform_gradient_log_det_jacobian =
            _inverse_transform_gradient_log_jacobian_bounded;
        transform->dim = Parameter_size(parameter);
    } else {
        fprintf(stderr, "Transform %s not implemented\n", transform);
        exit(1);
    }

    return transform;
}

Transform* new_SimplexTransform_with_parameter(const char* type, Parameter* parameter) {
    Transform* transform = (Transform*)malloc(sizeof(Transform));
    transform->lower = 0;
    transform->upper = 1;
    transform->parameter = parameter;

    transform->get = _get;
    transform->set = _set;
    transform->backward = _backward;
    transform->log_det_jacobian = _log_det_jacobian;
    transform->gradient_log_det_jacobian = _gradient_log_det_jacobian;
    transform->jacobian = _jacobian;

    transform->transform = _transform_simplex_stan;
    transform->inverse_transform = _inverse_transform_simplex_stan;
    transform->backward_inverse_transform = _backward_inverse_transform_simplex_stan;
    transform->inverse_transform_log_det_jacobian =
        _inverse_transform_log_jacobian_simplex_stan;
    transform->inverse_transform_gradient_log_det_jacobian =
        _inverse_transform_gradient_log_jacobian_simplex_stan;
    transform->inverse_transform_jacobian = _inverse_transform_simplex_jacobian_stan;
    transform->dim = Parameter_size(parameter) + 1;

    return transform;
}

void free_Transform(Transform* transform) { 
    free_Parameter(transform->parameter);    
    free(transform);
}
