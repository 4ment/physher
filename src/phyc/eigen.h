/*
 *  eigen.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/3/10.
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

#ifndef _EIGEN_H_
#define _EIGEN_H_


//#define LAPACK_ENABLED 1

#include <string.h>

typedef struct EigenDecomposition{
	double **evec;
	double **Invevec;
	double *eval;
	double *evali; // imaginary part of eigenvalues
	size_t dim;
#ifdef LAPACK_ENABLED
    int *isuppz;
    double *M;
#endif
}EigenDecomposition;

EigenDecomposition * new_EigenDecomposition( const size_t dimension );

void EigenDecomposition_decompose( double **a, EigenDecomposition *eigen );

void free_EigenDecomposition( EigenDecomposition *eigencmp );

EigenDecomposition * clone_EigenDecomposition( EigenDecomposition *eigen );



EigenDecomposition * eigen2( double **a, size_t dim );

int hqr2( int N, double **H, double *d, double *e, double **V, int maxIterations );

void normalize( double **a, size_t n );

void jacobi( double **a, int n, double *d, double **v, int *nrot);

#endif
