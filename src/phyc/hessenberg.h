/*
 *  hessenberg.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/7/10.
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
