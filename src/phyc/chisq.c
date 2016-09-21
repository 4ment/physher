/*
 *  chisq.c
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

#include "chisq.h"

#include "gamma.h"

// Cumulative distribution function
// return P(X <= x)
double pchisq( const double x, const int df ){
	if (x < 0.0 || df < 0.0) return 1.0;
	return gammp( ((double)df)/2, x/2);
}

// Inverse cumulative distribution function
double qchisq( const double p, const int df ){
	return 2.*invgammp(p,0.5*df);
}
