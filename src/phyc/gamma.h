/*
 *  gamma.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/24/11.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */


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
