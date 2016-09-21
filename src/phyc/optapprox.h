/*
 *  optapprox.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 23/01/2014.
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


#ifndef PhyC_optapprox_h
#define PhyC_optapprox_h

#include "treelikelihood.h"

// Use a Taylor expansion of the likelihood function using a diagonal variance matrix
double loglikelihood_approximation( SingleTreeLikelihood *tlk  );

// Use a Taylor expansion of the likelihood function using a diagonal variance matrix
double loglikelihood_approximation2( SingleTreeLikelihood *tlk  );

// Use a Taylor expansion of the likelihood function using a  variance-covariance matrix
double loglikelihood_approximation3( SingleTreeLikelihood *tlk  );


double loglikelihood_approximation_diag( void *data  );

double loglikelihood_approximation_full( void *data  );

#endif
