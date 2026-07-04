// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _GAUSSIAN_H_
#define _GAUSSIAN_H_

// PDF
double dnorm( const double x, const double mu, const double sigma );

// Log PDF
double dnorml( const double x, const double mu, const double sigma );

// CDF
double pnorm( const double x, const double mu, const double sigma );

// Inverse CDF
double qnorm( const double p, const double mu, const double sigma );

// Generate a Gaussian random variable
double rnorm();
 
#endif
