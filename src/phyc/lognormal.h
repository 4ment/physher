/*
 *  lognormal.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 13/8/11.
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

#ifndef _LOGNORMAL_H_
#define _LOGNORMAL_H_

double dlnorm( const double x,  const double mu, const double sigma );

// Cumulative distribution function
double plnorm( const double x, const double mu, const double sigma );


// Inverse cumulative distribution function
double qlnorm( const double p, const double logmu, const double logsigma );

double lognorm_mean( const double logmu, const double logsigma );

void lognormal_discretize( const double mu, const double sigma, double *bins, const int count );

#endif
