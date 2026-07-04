// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _EIGEN_H_
#define _EIGEN_H_


//#define LAPACK_ENABLED 1

#include <string.h>
#include <stdbool.h>

typedef struct EigenDecomposition{
	double **evec;
	double **Invevec;
	double *eval;
	double *evali; // imaginary part of eigenvalues
	size_t dim;
	bool failed;
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
