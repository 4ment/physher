/*
 *  dirichlet.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 30/08/2016.
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

#ifndef dirichlet_h
#define dirichlet_h

#include <stdio.h>

#include "gamma.h"

// Density function

double ddirchlet( const double *x, const size_t dim, const double *alphas );

double ddirchletln( const double *x, const size_t dim, const double *alphas );

// Flat dirichlet with alpha == (1,1,..,1)
double ddirchlet_flat( const size_t dim );

void rdirichlet(double*x, const size_t dim, const double* alphas);

#endif /* dirichlet_h */
