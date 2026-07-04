// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _GAMMA_H_
#define _GAMMA_H_
 
/*****************************************************
 * x ~ Gamma(alpha,beta)
 *****************************************************/

// Probality density distribution function of Gamma function
double dgamma(const double x, const double alpha, const double beta);

double dloggamma(const double x, const double alpha, const double beta);

// Cumulative distribution function of Gamma function
double pgamma( const double x, const double alpha, const double beta );

// Inverse cumulative distribution function of Gamma function
double qgamma( const double p, const double alpha, const double beta );

double rgamma(double a);

/*****************************************************
 * Incomplete Gamma function
 *****************************************************/

// Cumulative distribution function of incomplete gamma function
double gammp( const double a, const double x );

// Inverse incomplete Gamma function
double invgammp( const double p, const double a);


/*****************************************************
 * Gamma function
 *****************************************************/

// Returns the value ln[Γ(xx)] for xx > 0
double gammln( const double xx );

// Returns the value Γ(xx) for xx > 0
double gamm( const double xx );

#endif
