/*
 *  matrix.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 9/22/10.
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


#include "matrix.h"

#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#include "mathconstant.h"
#include "hessenberg.h" // for det

void cholesky(long N, double *A, double *diag){
    long i,j,k;
    for(j = 0; j < N; j++)
        diag[j] = A[N*j+j];
    for(j = 0; j < N; j++){
        for(k = 0; k < j; k++)
            diag[j] -= A[N*k+j]*A[N*k+j];
        diag[j] = sqrt(diag[j]);
        for(i = j+1; i < N;i++){
            for(k = 0; k < j; k++)
                A[N*j+i] -= A[N*k+i]*A[N*k+j];
            A[N*j+i]/=diag[j];
        }
    }
}


#pragma mark -
#pragma mark Vector

struct _Vector{
	size_t count;
	size_t capacity;
	double *vector;
};

Vector * new_Vector( const size_t length ){
	Vector *v = (Vector *)malloc(sizeof(Vector));
    assert(v);
	v->vector = dvector(length);
	v->count	= 0;
	v->capacity = length;
	if (!v) error("allocation failure in vector()"); 
	return v;
}

int Vector_length( Vector *v ){
	return v->count;
}

void Vector_push( Vector *v, double value ){
	if( v->count == v->capacity ){
		v->capacity++;
		v->vector = realloc(v->vector, v->capacity*sizeof(double));
	}
	v->vector[v->count] = value;
	v->count++;
}

void Vector_insert( Vector *v, double value, int index ){
	if( v->count == v->capacity ){
		v->capacity++;
		v->vector = realloc(v->vector, v->capacity*sizeof(double));
	}
	memmove(v->vector+index+1, v->vector+index, v->count*sizeof(double));
	v->vector[index] = value;
}

void Vector_unshift( Vector *v, double value ){
	if( v->count == v->capacity ){
		v->capacity++;
		v->vector = realloc(v->vector, v->capacity*sizeof(double));
	}
	memmove(v->vector+1, v->vector, v->count*sizeof(double));
	v->vector[0] = value;
}

double Vector_at( const Vector *v, int index ){
	return v->vector[index];
}

double Vector_pop( Vector *v ){
	if ( v->count == 0 ) {
		return 0;
	}
	double val = v->vector[v->count-1];
	v->count--;
	return val;
}

double Vector_shift( Vector *v ){
	if ( v->count == 0 ) {
		return 0;
	}	
	double val = v->vector[0];
	Vector_remove(v, 0);
	return val;
}

void Vector_remove_all( Vector *v ){
	v->count = 0;
}

void Vector_remove( Vector *v, int index ){
	memmove( v->vector+index, v->vector+index+1, (v->count-index-1)*sizeof(double));
	v->count--;
}

void Vector_sort_from_iVector( Vector *v, iVector *s ){
	bool done = false;
    int i;
	int size = Vector_length(v);
	while ( !done ) {
		done = true;
		for ( i = 0 ; i < size-1 ; i++ ) {
			if ( iVector_at(s,i) > iVector_at(s,i+1) ) {
				done = false;
				iVector_swap(s, i, i+1);
				Vector_swap(v, i, i+1);
			}
		}
		size--;
	}
}

void Vector_sort_from_ivector( Vector *v, int *s ){
	bool done = false;
    int i;
	int size = Vector_length(v);
	while ( !done ) {
		done = true;
		for ( i = 0 ; i < size-1 ; i++ ) {
			if ( s[i] > s[i+1] ) {
				done = false;
				swap_int( &s[i], &s[i+1] );
				Vector_swap(v, i, i+1);
			}
		}
		size--;
	}
}

// perm should be a vector containing 0,1,2,...,n
void Vector_sort_track( Vector *v, int *perm ){
    dvector_sort_track(v->vector, Vector_length(v), perm);
}

void Vector_swap( Vector *v, int a, int b ){
	double temp = v->vector[a];
	v->vector[a] = v->vector[b];
	v->vector[b] = temp;
}

double * Vector_data( Vector *v ){
	return v->vector;
}


void free_Vector( Vector *v ){
	free(v->vector);
	free(v);
}

#pragma mark -
#pragma mark iVector

struct _iVector{
    size_t count;
    size_t capacity;
    int *vector;
};

iVector * new_iVector( const size_t length ){
    iVector *v = (iVector *)malloc(sizeof(iVector));
    assert(v);
    v->vector = ivector(length);
    v->count	= 0;
    v->capacity = length;
    if (!v) error("allocation failure in vector()");
    return v;
}

int iVector_length( iVector *v ){
    return v->count;
}

void iVector_push( iVector *v, int value ){
    if( v->count == v->capacity ){
        v->capacity++;
        v->vector = realloc(v->vector, v->capacity*sizeof(int));
    }
    v->vector[v->count] = value;
    v->count++;
}

void iVector_insert( iVector *v, int value, int index ){
    if( v->count == v->capacity ){
        v->capacity++;
        v->vector = realloc(v->vector, v->capacity*sizeof(int));
    }
    memcpy(v->vector+1, v->vector, v->count*sizeof(int));
    v->vector[0] = value;
}

int iVector_at( iVector *v, int index ){
    return v->vector[index];
}

int iVector_pop( iVector *v ){
    if ( v->count == 0 ) {
        return 0;
    }
    int val = v->vector[v->count-1];
    v->count--;
    return val;
}

void iVector_remove_all( iVector *v ){
    v->count = 0;
}

void iVector_remove( iVector *v, int index ){
    memcpy( v->vector+index, v->vector+index+1, (v->count-index-1)*sizeof(int));
    v->count--;
}

void iVector_swap( iVector *v, int a, int b ){
    int temp = v->vector[a];
    v->vector[a] = v->vector[b];
    v->vector[b] = temp;
}

void free_iVector( iVector *v ){
    free(v->vector);
    free(v);
}


#pragma mark -
#pragma mark Matrix

Matrix * new_Matrix( const size_t nrow, const size_t ncol ){
	Matrix *m = (Matrix *)malloc(sizeof(Matrix));
    assert(m);
    size_t i;
//	assert(m != NULL);
//	m->matrix = (double **)malloc((size_t) nrow*sizeof(double*));
//	
//	m->nrow = nrow;
//	m->ncol = ncol;
//	
//	if (!m->matrix) error("row allocation failure in vector()");
//	assert(m->matrix != NULL);
//	m->matrix[0] = (double *)calloc(nrow*ncol, sizeof(double));
//	
//	for( size_t i = 1; i<nrow; i++ ) m->matrix[i] = m->matrix[i-1]+ncol;
	
	m->matrix = (double **)malloc((size_t) nrow*sizeof(double*));
    assert(m->matrix);
    for( i = 0; i< nrow; i++ ){
        m->matrix[i] = (double *)calloc(ncol, sizeof(double));
        assert(m->matrix[i]);
    }
	
	m->nrow = nrow;
	m->ncol = ncol;

	return m;
}



Matrix *Transpose( const Matrix *m ){
	Matrix *m2 = new_Matrix(m->ncol, m->nrow);
	long i,j;
	for( i = 0; i < m->nrow; i++)
		for( j = 0; j < m->ncol; j++)
			m2->matrix[i][j] = m->matrix[j][i];
	return m2;
}

void Add( Matrix *m1, const Matrix *m2){
	size_t i,j;
	for( i = 0; i < m1->nrow; i++)
		for( j = 0; j < m1->ncol; j++)
			m1->matrix[i][j] += m2->matrix[i][j];
}

Matrix *Add2( const Matrix *m1, const Matrix *m2){
	Matrix *m = new_Matrix(m1->nrow,m1->ncol);
	size_t i,j;
	for( i = 0; i < m1->nrow; i++)
		for( j = 0; j < m1->ncol; j++)
			m->matrix[i][j] = m1->matrix[i][j] + m2->matrix[i][j];
	return m;
}



Matrix *Mult( const Matrix *m1, const Matrix *m2 ){
	if( m1->ncol != m2->nrow )
		error("fmult: incompatible size\n");
	Matrix *m = new_Matrix(m1->nrow,m2->ncol);
	size_t i,j,k;
	for( i = 0; i < m1->nrow; i++ ){
		for( j = 0; j < m2->ncol; j++ ){
			m->matrix[i][j] = 0;
			for( k=0; k < m1->ncol; k++)
				m->matrix[i][j] += m1->matrix[i][k] * m2->matrix[k][j];
			
		}
	}
	return m;
}

void Mult_Scalar( Matrix *m, double value ){
	size_t i,j;
	for( i = 0; i < m->nrow; i++ )
		for( j = 0; j < m->ncol; j++ )
			m->matrix[i][j] *= value;
}


Matrix *Mult_Scalar2( const Matrix *m, double value ){
	Matrix *m2 = new_Matrix(m->nrow,m->ncol);
	size_t i,j;
	for( i = 0; i < m->nrow; i++ )
		for( j = 0; j < m->ncol; j++ )
			m2->matrix[i][j] = m->matrix[i][j] * value;
	return m2;
}

void Subtract( Matrix *m1, const Matrix *m2 ){
	size_t i,j;
	for( i = 0; i < m1->nrow; i++)
		for( j = 0; j < m1->ncol; j++)
			m1->matrix[i][j] -= m2->matrix[i][j];

}

Matrix *Subtract2( const Matrix *m1, const Matrix *m2 ){
	Matrix *m3 = new_Matrix(m1->nrow,m1->ncol);
	size_t i,j;
	for( i = 0; i < m1->nrow; i++ )
		for( j = 0; j < m1->ncol; j++ )
			m3->matrix[i][j] = m1->matrix[i][j] - m2->matrix[i][j];
	return m3;
}

void free_Matrix( Matrix *m ){
	free_dmatrix( m->matrix, m->nrow );
	free( m );
}


/****************************************************************************/

#pragma mark -
#pragma mark float matrix

float **fmatrix( const size_t nrow, const size_t ncol ){
    size_t i;
    float **m = (float **)malloc(nrow * sizeof(float *) );
    assert(m);
    for( i = 0; i < nrow; i++ ){
        m[i] = (float *)calloc(ncol, sizeof(float) );
        assert(m[i]);
    }
    return m;
}


/****************************************************************************/

#pragma mark -
#pragma mark matrix

double ** dmatrix2( const size_t nrow, const size_t ncol ){
	int i;
	double **matrix = (double **)malloc( nrow*sizeof(double*) );
	assert(matrix);
	matrix[0] = (double *)calloc(nrow*ncol, sizeof(double));
	assert(matrix[0]);
	for(i=1; i<nrow; i++ ) matrix[i] = matrix[i-1]+ncol;
	
	return matrix;
}

void free_dmatrix2(double **m, const long nrow, const long ncol ){
	free(m[0]);
	free(m);
}

double **dmatrix( const size_t nrow, const size_t ncol ){
    size_t i;
    double **m = (double **)malloc(nrow * sizeof(double *) );
	assert(m);
	for( i = 0; i < nrow; i++ ){
		m[i] = (double *)calloc(ncol, sizeof(double) );
		assert(m[i]);
	}
	return m;
}

double **deye( const size_t n ){
    size_t i;
	double ** m = dmatrix( n, n );
	for ( i = 0; i < n; i++) {
		m[i][i] = 1.0;
	}
	return m;
}

double **deye2( double **m, const size_t n ){
	size_t i,j;
	for ( i = 0; i < n; i++) {
		m[i][i] = 1.;
		for ( j = i+1; j < n; j++ ) {
			m[i][j] = m[j][i] = 0.;
		}
	}
	return m;
}

int **imatrix( const size_t nrow, const size_t ncol ){
    size_t i;
	int **m = (int **)calloc(nrow, sizeof(int *) );
	assert(m);
	for( i = 0; i < nrow; i++){
		m[i] = (int *)calloc(ncol, sizeof(int) );
		assert(m[i]);
	}
	return m;
}

long **lmatrix( const size_t nrow, const size_t ncol ){
    size_t i;
	long **m = (long **)calloc(nrow, sizeof(long *) );
	assert(m);
	for( i = 0; i < nrow; i++){
		m[i] = (long *)calloc(ncol, sizeof(long) );
		assert(m[i]);
	}
	return m;
}
unsigned **uimatrix( const size_t nrow, const size_t ncol ){
    size_t i;
    unsigned **m = (unsigned **)calloc(nrow, sizeof(unsigned *) );
	assert(m);
	for( i = 0; i < nrow; i++){
		m[i] = (unsigned *)calloc(ncol, sizeof(unsigned) );
		assert(m[i]);
	}
	return m;
}

char **cmatrix( const size_t nrow, const size_t ncol ){
    size_t i;
	char **m = (char **)malloc(nrow * sizeof(char *) );
	assert(m);
	for( i = 0; i < nrow; i++){
		m[i] = (char *)calloc(ncol, sizeof(char) );
		assert(m[i]);
	}
	return m;
}

bool **bmatrix( const size_t nrow, const size_t ncol ){
    size_t i;
	bool **m = (bool **)malloc(nrow * sizeof(bool *) );
	assert(m);
	for( i = 0; i < nrow; i++){
		m[i] = (bool *)calloc(ncol, sizeof(bool) );
		assert(m[i]);
	}
	return m;
}

uint8_t ** ui8matrix( const size_t nrow, const size_t ncol ){
    size_t i;
	uint8_t **m = (uint8_t **)malloc(nrow * sizeof(uint8_t *) );
	assert(m);
	for( i = 0; i < nrow; i++){
		m[i] = (uint8_t *)calloc(ncol, sizeof(uint8_t) );
		assert(m[i]);
	}
	return m;
}

void free_dmatrix( double **m, const size_t nrow ){
    size_t i;
    for( i =0; i < nrow; i++ ){
        free(m[i]);
        m[i] = NULL;
    }
    free(m);
    m = NULL;
}

void free_fmatrix( float **m, const size_t nrow ){
    size_t i;
    for( i =0; i < nrow; i++ ){
        free(m[i]);
    }
    free(m);
    m = NULL;
}

void free_imatrix( int **m, const size_t nrow ){
    size_t i;
	for( i = 0; i < nrow; i++ ){
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}

void free_lmatrix( long **m, const size_t nrow ){
    size_t i;
	for( i = 0; i < nrow; i++ ){
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}

void free_uimatrix( unsigned **m, const size_t nrow ){
    size_t i;
	for( i = 0; i < nrow; i++ ){
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}

void free_cmatrix( char **m, const size_t nrow ){
	for( size_t i = 0; i < nrow; i++ ){
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}

void free_bmatrix( bool **m, const size_t nrow ){
	for( size_t i = 0; i < nrow; i++ ){
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}


void free_ui8matrix( uint8_t **m, const size_t nrow ){
	for( size_t i = 0; i < nrow; i++ ){
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}


double ** clone_dmatrix(  double **const m, const size_t nrow, const size_t ncol ){
	double **m2 = dmatrix(nrow, ncol);
	memcpy_dmatrix(m2, m, nrow, ncol);
	return m2;
}

int **clone_imatrix( int **const m, const size_t nrow, const size_t ncol ){
	int **m2 = imatrix(nrow, ncol);
	memcpy_imatrix(m2, m, nrow, ncol);
	return m2;
}

void memcpy_dmatrix(  double **dst, double ** const src, const size_t nrow, const size_t ncol){
	for( size_t i = 0; i < nrow; i++){
		memcpy(dst[i], src[i], ncol * sizeof(double) );
	}
}

void memcpy_imatrix(  int **dst, int ** const src, const size_t nrow, const size_t ncol){
	for( size_t i = 0; i < nrow; i++){
		memcpy(dst[i], src[i], ncol * sizeof(int) );
	}
}

void print_dmatrix_to_file(){
    
}
#pragma mark -
#pragma mark 3D matrix

double *** d3matrix( const size_t dim1, const size_t dim2, const size_t dim3 ){
    assert(dim2>0);
	size_t i,j;
	double ***m = calloc( dim1, sizeof(double**));
	assert(m);
	for ( i = 0 ; i < dim1; i++ ) {
		m[i] = malloc( dim2 * sizeof(double*));
		assert(m[i]);
		
		for ( j =0 ; j < dim2; j++ ) {
			m[i][j] = (double*)calloc(  dim3, sizeof(double));
			assert(m[i][j]);
			
		}		
	}
	return m;
}


double *** d3matrix2( const size_t dim1, const size_t dim2 ){
    assert(dim2>0);
	size_t i,j;
	double ***m = (double***)calloc( dim1, sizeof(double**));
	assert(m);
	for ( i = 0 ; i < dim1; i++ ) {
		m[i] = (double**)calloc( dim2, sizeof(double*));
		assert(m[i]);
		
		for ( j =0 ; j < dim2; j++ ) {
			m[i][j] = NULL;
			
		}		
	}
	return m;
}

void free_d3matrix( double ***m, const size_t dim1, const size_t dim2 ){
	size_t i,j;
	for( i = 0; i < dim1; i++ ){
		for( j = 0; j < dim2; j++ ){
			free(m[i][j]);
			m[i][j] = NULL;
		}
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}

char *** c3matrix( const size_t dim1, const size_t dim2, const size_t dim3 ){
    assert(dim2);
	size_t i,j;
	char ***m = (char***)calloc( dim1, sizeof(char**));
	assert(m);
	for ( i = 0 ; i < dim1; i++ ) {
		m[i] = (char**)calloc( dim2, sizeof(char*));
		assert(m[i]);
		
		for ( j = 0 ; j < dim2; j++ ) {
			m[i][j] = (char*)calloc(  dim3, sizeof(char));
			assert(m[i][j]);
			
		}		
	}
	return m;
}

void free_c3matrix( char ***m, const size_t dim1, const size_t dim2, const size_t dim3 ){
	size_t i,j;
	for( i = 0; i < dim1; i++ ){
		for( j = 0; j < dim2; j++ ){
			free(m[i][j]);
			m[i][j] = NULL;
		}
		free(m[i]);
		m[i] = NULL;
	}
	free(m);
	m = NULL;
}

double *** clone_d3matrix( double ***mat, const size_t dim1, const size_t dim2, const size_t dim3 ){
    assert(dim2);
	size_t i,j;
	double ***m = (double***)malloc( dim1 * sizeof(double**));
	assert(m);
	for ( i = 0 ; i < dim1; i++ ) {
		m[i] = (double**)malloc( dim2 * sizeof(double*));
		assert(m[i]);
		for ( j =0 ; j < dim2; j++ ) {
			m[i][j] = (double*)malloc(  dim3 * sizeof(double));
			assert(m[i][j]);
			
			memcpy(m[i][j], mat[i][j], dim3 * sizeof(double));
		}		
	}
	return m;
}

#pragma mark -
#pragma mark vector

double * dvector( const size_t n ){
	double *v = (double *)calloc( n, sizeof(double));
	assert(v);
	return v;
}

unsigned int * uivector( const size_t n ){
	unsigned int *v = (unsigned int *)calloc( n, sizeof(unsigned int));
	assert(v);
	return v;
}

int * ivector( const size_t n ){
	int *v = (int *)calloc( n, sizeof(int));
	assert(v);
	return v;
}

float * fvector( const size_t n ){
	float *v = (float *)calloc( n, sizeof(float));
	assert(v);
	return v;
}

char * cvector( const size_t n ){
	char *v = (char *)calloc( n, sizeof(char));
	assert(v);
	return v;
}

bool * bvector( const size_t n ){
	bool *v = (bool *)calloc( n, sizeof(bool));
	assert(v);
	return v;
}

double *clone_dvector( const double * v, const size_t n ){
	double *clone = dvector( n );
	memcpy(clone, v, n*sizeof(double) );
	return clone;
}

int *clone_ivector( const int * v, const size_t n ){
	int *clone = ivector( n );
	memcpy(clone, v, n*sizeof(int) );
	return clone;	
}

unsigned int *clone_uivector( const unsigned int * v, const size_t n ){
	unsigned int *clone = uivector( n );
	memcpy(clone, v, n*sizeof( int unsigned) );
	return clone;	
}

char *clone_cvector( const char *v, const size_t n ){
	char *clone = cvector( n );
	memcpy(clone, v, n*sizeof(char) );
	return clone;	
}

bool * clone_bvector( const bool *v, const size_t n ){
	bool *clone = bvector( n );
	memcpy(clone, v, n*sizeof(bool) );
	return clone;	
}


void dvector_sort_from_ivector( double *v, int *s, int size ){
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( s[i] > s[i+1] ) {
				done = false;
				swap_int( &s[i], &s[i+1] );
				dswap( &v[i], &v[i+1] );
			}
		}
		size--;
	}
}

// perm should be a vector containing 0,1,2,...,n
// increasing
void dvector_sort_track( double *v, size_t size, int *perm ){
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( v[i] > v[i+1] ) {
				done = false;
				swap_int( &perm[i], &perm[i+1] );
                dswap(&v[i], &v[i+1]);
			}
		}
		size--;
	}
}

#pragma mark -

double * Matrix_mult( const double *A, const double *B, size_t Am, size_t An, size_t Bm, size_t Bn ){
	if( An != Bm ) error("Matrix_mult: incompatible size\n");
	double *m = dvector( Am * Bn );
	
	size_t i,j,k;
	for( i = 0; i < Am; i++ ){
		for( j = 0; j < Bn; j++ ){
			for( k = 0; k < Bm; k++ )
				m[i*Bn+j] += A[i*An+k] * B[k*Bn+j];
		}
	}
	return m;
}

void Matrix_mult2( double* m, const double *A, const double *B, size_t Am, size_t An, size_t Bm, size_t Bn ){
	if( An != Bm ) error("Matrix_mult: incompatible size\n");
	memset(m, 0, sizeof(double)*Am*Bn);
	
	size_t i,j,k;
	for( i = 0; i < Am; i++ ){
		for( j = 0; j < Bn; j++ ){
			for( k = 0; k < Bm; k++ )
				m[i*Bn+j] += A[i*An+k] * B[k*Bn+j];
		}
	}
}

void Matrix_mult3( double* m, const double **A, const double *B, size_t Am, size_t An, size_t Bm, size_t Bn ){
	if( An != Bm ) error("Matrix_mult: incompatible size\n");
	memset(m, 0, sizeof(double)*Am*Bn);
	
	size_t i,j,k;
	for( i = 0; i < Am; i++ ){
		for( j = 0; j < Bn; j++ ){
			for( k = 0; k < Bm; k++ )
				m[i*Bn+j] += A[i][k] * B[k*Bn+j];
		}
	}
}

void Matrix_mult4( double* m, const double *A, const double **B, size_t Am, size_t An, size_t Bm, size_t Bn ){
	if( An != Bm ) error("Matrix_mult: incompatible size\n");
	memset(m, 0, sizeof(double)*Am*Bn);
	
	size_t i,j,k;
	for( i = 0; i < Am; i++ ){
		for( j = 0; j < Bn; j++ ){
			for( k = 0; k < Bm; k++ )
				m[i*Bn+j] += A[i*An+k] * B[k][j];
		}
	}
}
#pragma mark -
#pragma mark print functions

void print_dmatrix( FILE *file, const double **m, const size_t dim1, const size_t dim2, char sep ){
	size_t i,j;
	for( i = 0; i < dim1; i++ ){
		for( j = 0; j < dim2; j++ ){
			fprintf(file, "%f%c",m[i][j], (j == dim2-1 ? '\n' : sep));
		}
	}
}

void print_dvector(  const double *m,  const size_t dim ){
	fprintf(stdout, "\n=-------------------\n");
	for( size_t i = 0; i < dim; i++ ){
		fprintf(stdout, "%e ",m[i]);
	}
	fprintf(stdout, "\n===================\n");
}

void print_uivector( const unsigned *m, const size_t dim ){
	//fprintf(stdout, "\n=-------------------\n");
	for( size_t i = 0; i < dim; i++ ){
		fprintf(stderr, "%u ",m[i]);
	}
	fprintf(stderr, "\n");
	//fprintf(stdout, "\n===================\n");
}

void print_ivector(  const int *m,  const size_t dim ){
	fprintf(stdout, "\n=-------------------\n");
	for( size_t i = 0; i < dim; i++ ){
		fprintf(stdout, "%d ",m[i]);
	}
	fprintf(stdout, "\n===================\n");
}

void print_imatrix( int **mat, const int m, const size_t n ){
	size_t i,j;
	fprintf(stdout, "\n=-------------------\n");
	for( i = 0; i < m; i++ ){
		for( j = 0; j < n; j++ ){
			fprintf(stdout, "%d ",mat[i][j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n===================\n");
}


void compare_dmatrix(  const double **m1,  const double **m2, const size_t nrow, const size_t ncol ){
	for (size_t i = 0; i < nrow; i++) {
		assert(memcmp(m1[i], m2[i], ncol*sizeof(double) ) == 0 );
	}	
}
void compare_d3matrix( const double **m1, const double **m2, const size_t dim1, const size_t dim2, const size_t dim3 ){
	size_t i,j;
	for ( i = 0; i < dim1; i++) {
		for ( j = 0; j < dim2; j++) 
			assert(memcmp(m1[i], m2[i], dim3*sizeof(double) ) == 0 );
	}
}

unsigned int row_index( unsigned int i, unsigned int M ){
    double m = M;
    double row = (-2*m - 1 + sqrt( (4*m*(m+1) - 8*(double)i - 7) )) / -2;
    if( row == (double)(int) row ) row -= 1;
    return (unsigned int) row;
}


unsigned int column_index( unsigned int i, unsigned int M ){
    unsigned int row = row_index( i, M);
    return  i - M * row + row*(row+1) / 2;
}

#pragma mark -

// Approximation of a function at point x near a
double normal_approximation( double *a, double *x, double *cov, size_t dim ){
    assert(dim>0);
	double *diff = dvector(dim);
	size_t i = 0;
	for ( ; i < dim; i++) {
		diff[i] = x[i] - a[i];
	}
	double det = LUDecompose_det_and_inverse(cov, dim); // cov becomes inverse(cov)
	
	double *temp = Matrix_mult(diff, cov, 1, dim, dim, dim);
	
	double d = 0.0;
	for ( i = 0; i < dim; i++ ) {
		d += ( temp[dim] * diff[dim] );
	}
	
	d = dim * LOG_2PI + log(det) + d;
	free(temp);
	free(diff);
	return d;
}

void cov2cor(double* cov, int dim){
	double* d = dvector(dim);
	for(int i = 0; i < dim; i++){
		d[i] = 1.0/sqrt(cov[i*dim+i]);
	}
	for(int i = 0; i < dim; i++){
		for(int j = 0; j < dim; j++){
			cov[j*dim+i] *= d[i] * d[j];
		}
	}
	free(d);
}
