// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "exp.h"

#include <math.h>

#include "matrix.h"
#include "eigen.h"

Matrix *Exp_Taylor_Series( Matrix *m, const int n ){
	size_t i,j;
	int k = 1;
	Matrix *e = new_Matrix( m->nrow, m->ncol );
	Matrix *f = new_Matrix( m->nrow, m->ncol );
	
	for( i = 0; i < m->nrow; i++ ){
		for( j = 0; j < m->ncol; j++ ){
			e->matrix[i][j] = 0;
			f->matrix[i][j] = 0;
		}
		f->matrix[i][i] = 1;
	}
	for ( i = 0; i < n; i++ ) {
		Add(e,f);
		Matrix *temp = Mult(f,m);
        free_Matrix(f);
        f = temp;
		Mult_Scalar(f,1/k);
		k++;
	}
	free_Matrix(f);
	return e;
}

double ** Exp_QR_EigenDecomposition( EigenDecomposition *ed ){
	int i,j,k;
	double expEval;
	double **expp = dmatrix(ed->dim, ed->dim );
	for ( i = 0; i < ed->dim; i++ ){
		expEval = exp(ed->eval[i]);
		for ( j = 0; j < ed->dim; j++ ){
			expp[i][j] = ed->Invevec[i][j] * expEval;
		}
	}
	
	double **exp_matrix = dmatrix(ed->dim, ed->dim );
	for ( i=0; i < ed->dim; i++){
		for ( j=0; j < ed->dim; j++){
			exp_matrix[i][j] = 0;
			for ( k = 0; k < ed->dim; k++ )
				exp_matrix[i][j] += ed->evec[i][k] * expp[k][j];
		}
	}
	free(expp);
	
	return exp_matrix;
}


double ** Exp_QR( double **matrix, const int dim ){
	EigenDecomposition *ed = new_EigenDecomposition( dim );
	EigenDecomposition_decompose( matrix, ed);
	return Exp_QR_EigenDecomposition(ed);
}
