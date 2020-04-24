//
//  transforms.h
//  viphy
//
//  Created by Mathieu Fourment on 4/05/2017.
//  Copyright © 2017 University of Technology Sydney. All rights reserved.
//

#ifndef transforms_h
#define transforms_h

#include <stdio.h>


double softplus(double x);

double logit(const double x);

double inverse_logit(double y);

double grad_inverse_logit(const double y);

double grad_log_det_inverse_logit(double y);



double transform(double x, double lb, double ub);

double transform2(double x, double lb, double ub, double* lp);

double inverse_transform(double x, double lb, double ub, double* lp);

double inverse_transform2(double x, double lb, double ub);

double grad_inverse_transform(double x, double lb, double ub);

// Calculate gradient of log determinant of Jacobian
// the det should always be positive because all the terms are positive so no sign(det) term
double grad_log_det_inverse_transform(double x, double lb, double ub);



void inverse_transform_simplex(double* x, const double* y, int K, double* lp);

void inverse_transform_simplex2(double* x, const double* y, int K);

void grad_log_det_inverse_transform_simplex(const double* z, double* dz, int K);

#endif /* transforms_h */
