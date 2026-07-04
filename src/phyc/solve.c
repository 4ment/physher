// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "solve.h"

#include "hessenberg.h"
#include "matrix.h"

void inverse(double **a, int dim){
    int i,j,*indx;
    double **y,d,*col;
    
    y = dmatrix(dim,dim);
    indx = ivector(dim);
    col = dvector(dim);
    
    ludcmp(a,dim,indx,&d);
    
    for ( j = 0 ; j < dim; j++ ){
        memset(col, 0, sizeof(double)*dim);
        col[j] = 1.0;
        lubksb(a,dim,indx,col);
        for ( i = 0; i < dim; i++ ) y[i][j] = col[i];
    }
    for ( i = 0; i < dim; i++ )
        for ( j = 0; j < dim; j++ )
            a[i][j] = y[i][j];
    
    free_dmatrix(y,dim);
    free(col);
    free(indx);
}

void inverse2(double *a, int dim){
    int i,j,*indx;
    double **y,d,*col;
    
    y = dmatrix(dim,dim);
    indx = ivector(dim);
    col = dvector(dim);
    
    ludcmp2(a,dim,indx,&d);
    
    for ( j = 0 ; j < dim; j++ ){
        memset(col, 0, sizeof(double)*dim);
        col[j] = 1.0;
        lubksb2(a,dim,indx,col);
        for ( i = 0; i < dim; i++ ) y[i][j] = col[i];
    }
    for ( i = 0; i < dim; i++ )
        for ( j = 0; j < dim; j++ )
            a[i*dim+j] = y[i][j];
    
    free_dmatrix(y,dim);
    free(col);
    free(indx);
}
