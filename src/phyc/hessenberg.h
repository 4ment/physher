// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _HESSENBERG_H_
#define _HESSENBERG_H_

#include <string.h>

void balance(double **a, int n);

void elmhes(double **a, int *order, const int n );

void orthes ( int n, double **H, double **V, double *ort);

void hqr(double **a, int n, double *wr, double *wi);

void eltran(double **a, double **zz, int *order, int n);

void ludcmp(double **a, int n, int *indx, double *d);
void ludcmp2(double *a, int n, int *indx, double *d);

void lubksb( double **a, int n, int *indx, double *b );

void lubksb2( double *a, int n, int *indx, double *b );

void hqr3(int n, int low, int hgh, double **h, double **zz, double *wr, double *wi);


double LUDecompose_det( double **m, const int dim );
double LUDecompose_det2( double *m, const int dim );

double LUDecompose_det_and_inverse( double *m, const size_t dim );




#endif
