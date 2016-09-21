/*
 *  exponential.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 3/7/11.
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


#ifndef _EXPONENTIAL_H_
#define _EXPONENTIAL_H_

double rexp( const double lambda );

double dexp( const double x,  const double lambda);

// Cumulative distribution function
double pexp( const double x, const double lambda);


// Inverse cumulative distribution function
double qexp( const double p, const double lambda);

double exp_mean( const double lambda );

void exponential_discretize( const double lambda, double *bins, const int count );


#endif
