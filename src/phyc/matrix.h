/*
 *  matrix.h
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


#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>


#include "utils.h"

struct _Vector;
typedef struct _Vector Vector;

struct _iVector;
typedef struct _iVector iVector;

struct _dVector;
typedef struct _dVector dVector;

typedef struct __matrix{
	size_t nrow;
	size_t ncol;
	double **matrix;
} Matrix;


void cholesky(long N, double *A, double *diag);

#pragma mark -
// MARK: Vector

Vector *new_Vector( const size_t length );

int Vector_length( Vector *v );

void Vector_push( Vector *v, double value );

double Vector_at( const Vector *v, int index );

void Vector_insert( Vector *v, double value, int index );

void Vector_unshift( Vector *v, double value );

double Vector_pop( Vector *v );

double Vector_shift( Vector *v );

void Vector_remove_all( Vector *v );

void Vector_remove( Vector *v, int index );

void Vector_sort_from_iVector( Vector *v, iVector *s );

void Vector_sort_from_ivector( Vector *v, int *s );

void Vector_swap( Vector *v, int a, int b );

double * Vector_data( Vector *v );

void free_Vector( Vector *v );

void Vector_pack(Vector* v);

#pragma mark -
// MARK: iVector

iVector * new_iVector( const size_t length );

int iVector_length( iVector *v );

void iVector_push( iVector *v, int value );

void iVector_insert( iVector *v, int value, int index );

int iVector_at( iVector *v, int index );

int iVector_pop( iVector *v );

void iVector_remove_all( iVector *v );

void iVector_remove( iVector *v, int index );

void iVector_swap( iVector *v, int a, int b );

void free_iVector( iVector *v );


#pragma mark -
// MARK: Matrix

Matrix *new_Matrix( const size_t nrow, const size_t ncol );


void Add( Matrix *m1, const Matrix *m2);

void Mult_Scalar( Matrix *m, const double value );

void Subtract( Matrix *m1, const Matrix *m2 );


Matrix *Transpose( const Matrix *m );

Matrix *Add2( const Matrix *m1, const Matrix *m2);

Matrix *Mult( const Matrix *m1, const Matrix *m2 );

Matrix *Mult_Scalar2( const Matrix *m, double value );

Matrix *Subtract2( const Matrix *m1, const Matrix *m2 );


void free_Vector( Vector *v );


void free_Matrix( Matrix *m);

/****************************************************************************/

float **fmatrix( const size_t nrow, const size_t ncol );

/************************************************************/

double ** dmatrix( const size_t nrow, const size_t ncol );
int    ** imatrix( const size_t nrow, const size_t ncol );
long   ** lmatrix( const size_t nrow, const size_t ncol );
unsigned ** uimatrix( const size_t nrow, const size_t ncol );
char   ** cmatrix( const size_t nrow, const size_t ncol );
uint8_t ** ui8matrix( const size_t nrow, const size_t ncol );
bool ** bmatrix( const size_t nrow, const size_t ncol );

double **deye( const size_t n );
double **deye2( double **m, const size_t n );

double **clone_dmatrix(  double ** const m, const size_t nrow, const size_t ncol );
int    **clone_imatrix(  int ** const m, const size_t nrow, const size_t ncol );

void free_dmatrix( double **m, const size_t nrow );
void free_fmatrix( float **m, const size_t nrow );
void free_imatrix( int **m, const size_t nrow );
void free_lmatrix( long **m, const size_t nrow );
void free_uimatrix( unsigned **m, const size_t nrow );
void free_cmatrix( char **m, const size_t nrow );
void free_ui8matrix( uint8_t **m, const size_t nrow );
void free_bmatrix( bool **m, const size_t nrow );


void memcpy_dmatrix(  double **dst, double ** const src, const size_t nrow, const size_t ncol);
void memcpy_imatrix(  int **dst, int ** const src, const size_t nrow, const size_t ncol);


void print_dmatrix( FILE *file, const double **m, const size_t dim1, const size_t dim2, char sep );
void print_imatrix( int **mat, const int m, const size_t n );

void compare_dmatrix(  const double **m1, const double **m2, const size_t nrow, const size_t ncol );

/************************************************************/

double *** d3matrix( const size_t dim1, const size_t dim2, const size_t dim3 );
double *** d3matrix2( const size_t dim1, const size_t dim2 );
char *** c3matrix( const size_t dim1, const size_t dim2, const size_t dim3 );

void free_c3matrix( char ***m, const size_t dim1, const size_t dim2, const size_t dim3 );
void free_d3matrix( double ***m, const size_t dim1, const size_t dim2 );

double *** clone_d3matrix( double ***mat, const size_t dim1, const size_t dim2, const size_t dim3 );

void compare_d3matrix( const double **m1, const double **m2, const size_t dim1, const size_t dim2, const size_t dim3 );

/************************************************************/

double  * dvector( const size_t n );

double *clone_dvector( const double * v, const size_t n );

void dvector_sort_from_ivector( double *v, int *s, int size );

void dvector_sort_track( double *v, size_t size, int *perm );

void print_dvector( const double *m, const size_t dim );


unsigned int * uivector( const size_t n );
int     * ivector( const size_t n );
float   * fvector( const size_t n );
char    * cvector( const size_t n );
bool * bvector( const size_t n );

int    *clone_ivector( const int * v, const size_t n );
unsigned int *clone_uivector( const unsigned int * v, const size_t n );
char *clone_cvector( const char *v, const size_t n );
bool * clone_bvector( const bool *v, const size_t n );

void print_uivector( const unsigned *m, const size_t dim );
void print_ivector( const int *m, const size_t dim );


/************************************************************/
unsigned int row_index( unsigned int i, unsigned int M );

unsigned int column_index( unsigned int i, unsigned int M);

double * Matrix_mult( const double *A, const double *B, size_t Am, size_t An, size_t Bm, size_t Bn );

void Matrix_mult2( double* m, const double *A, const double *B, size_t Am, size_t An, size_t Bm, size_t Bn );

void Matrix_mult3( double* m, const double **A, const double *B, size_t Am, size_t An, size_t Bm, size_t Bn );

void Matrix_mult4( double* m, const double *A, const double **B, size_t Am, size_t An, size_t Bm, size_t Bn );

void cov2cor(double* cov, int dim);

#endif
