/*
 *  dirichlet.c
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

#include "dirichlet.h"

#include <math.h>

double ddirchletln( const double *x, const int dim, const double *alphas ){
    double logp = 0;
    double sum = 0;
    for (int i = 0; i < dim; i++) {
        logp += alphas[i]-1 + log(x[i]);
    }
    for (int i = 0; i < dim; i++) {
        sum += alphas[i];
    }
    logp += gammln(sum);
    
    for (int i = 0; i < dim; i++) {
        logp -= gammln(alphas[i]);
    }
    return logp;
}

double ddirchlet( const double *x, const int dim, const double *alphas ){
    return exp(ddirchletln(x, dim, alphas));
}

double ddirchlet_flat( const int dim ){
    return gamm(dim);
}
