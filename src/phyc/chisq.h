/*
 *  chisq.h
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


#ifndef _CHI_SQUARE_H_
#define _CHI_SQUARE_H_

#include "gamma.h"

// Cumulative distribution function
// return P(X <= x)
double pchisq( const double x, const int df );

// Inverse cumulative distribution function
double qchisq( const double p, const int df );

#endif
