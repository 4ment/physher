/*
 *  solve.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 14/12/2015.
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
