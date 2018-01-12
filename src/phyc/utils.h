/*
 *  utils.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 9/19/10.
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


#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include <stdbool.h>


#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#if defined (USE_FLOAT)
#define real_t float
#else
#define real_t double
#endif

double logaddexp(double x, double y);

double logFactorial(int n);

double logDoubleFactorial(int n);

#pragma mark -
#pragma mark double

double dmax( double a, double b);

double dmin( double a, double b);

double dmax_vector( const double *array, int n );

double dmin_vector( const double *array, int n );

void dswap( double *a, double *b );

int qsort_desc_dvector( const void *a, const void *b );

int qsort_asc_dvector( const void *a, const void *b );

void randomize_dvector( double *v, double n );


#pragma mark -
#pragma mark int

int imax( int a, int b);

int imin( int a, int b);

int imax_vector( const int *array, int n );

void swap_int( int *a, int *b);

int qsort_desc_ivector( const void *a, const void *b );

int qsort_asc_ivector( const void *a, const void *b );

void randomize_ivector( int *v, int n );

void reverse_ivector( int *v, unsigned int n );


#pragma mark -
#pragma mark size_t

size_t stmax( size_t a, size_t b );

size_t stmin( size_t a, size_t b);

int qsort_desc_ivector( const void *a, const void *b );

int qsort_asc_ivector( const void *a, const void *b );


#pragma mark -
#pragma mark unsigned int

unsigned uimax_vector( const unsigned *array, const int n );

void swap_uint( unsigned int *a, unsigned int *b);

int qsort_desc_uivector( const void *a, const void *b );

int qsort_asc_uivector( const void *a, const void *b );

void randomize_uivector( unsigned int *v, int n );

void uivector_canonical( unsigned *vec, size_t len_vec, int *helper, size_t helper_len );


#pragma mark -
#pragma mark real

real_t rmax( real_t a, real_t b);

real_t rmin( real_t a, real_t b);

real_t rmax_vector( const real_t *array, int n );

real_t rmin_vector( const real_t *array, int n );

void rswap( real_t *a, real_t *b);

int qsort_desc_rvector( const void *a, const void *b );

int qsort_asc_rvector( const void *a, const void *b );

#pragma mark -
#pragma mark sequences

double * linearly_spaced_vector( double lower, double upper, int n );

double * log_spaced_spaced_vector( double lower, double upper, int n );
void log_spaced_spaced_vector2( double *v, double lower, double upper, int n );

double * exp_spaced_spaced_vector( double lower, double upper, int n );
void exp_spaced_spaced_vector2( double *v, double lower, double upper, int n );

#pragma mark -
#pragma mark other

void error( const char message[] );

void shift4( double *a, double *b, double *c, double d );

void rshift4( real_t *a, real_t *b, real_t *c, real_t d );

double SQR( double a );



bool file_exists(const char * filename);

bool isInt( const char *str );

bool isFloat( const char *str );

bool isFloat2( const char *str );

void print_pretty_time( FILE *file, double dseconds );


void * aligned16_malloc( size_t size );

#endif
