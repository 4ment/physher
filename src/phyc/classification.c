/*
 *  classification.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 5/09/13.
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

#include <stdio.h>

#include "classification.h"

#include "matrix.h"


// data should be sorted
int * classification_Jenks_breaks(double *data, int n, int nclass) {
    
    int    **mat1 = imatrix(n + 1, nclass + 1);
    double **mat2 = dmatrix(n + 1, nclass + 1);
    
    for ( int i = 1; i <= nclass; i++ ) {
        mat1[1][i] = 1;
        mat2[1][i] = 0;
        for ( int j = 2; j <= n; j++ ){
            mat2[j][i] = INFINITY;
        }
    }
    
    double variance = 0;
    for (int l = 2; l <= n; l++) {
        double sum = 0;
        double sum_square = 0;
        double w = 0; // number of data points considered so far
        
        for ( int m = 1; m <= l; m++ ) {
            int lower_class_limit = l - m + 1;
            
            double val = data[lower_class_limit-1];
            
            
            sum_square += val * val;
            sum += val;
            
            
            w++;
            variance = sum_square - (sum * sum) / w;
            int i4 = lower_class_limit - 1;
            if ( i4 != 0 ) {
                for ( int j = 2; j <= nclass; j++ ) {
                    if ( mat2[l][j] >= (variance + mat2[i4][j - 1]) ) {
                        mat1[l][j] = lower_class_limit;
                        mat2[l][j] = variance + mat2[i4][j - 1];
                        
                    }
                }
            }
        }
        mat1[l][1] = 1;
        mat2[l][1] = variance;
    }
    
    int k = n;
    int *kclass = ivector(nclass);
    
    kclass[nclass - 1] = n - 1;    
    
    for ( int j = nclass; j >= 2; j-- ) {
        int id =  mat1[k][j] - 2;
        
        kclass[j - 2] = id;
        
        k = mat1[k][j] - 1;
        
        
    }
    
    free_dmatrix(mat2, nclass+1);
    free_imatrix(mat1, nclass+1);
    
    return kclass;
}


