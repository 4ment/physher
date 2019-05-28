/*
 *  utils.c
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

#include "utils.h"

#include <assert.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>

#include "random.h"
#include "matrix.h"
#include "mstring.h"

#if defined (USE_FLOAT)
#define real_t float
#else
#define real_t double
#endif


double logaddexp(double x, double y){
	double tmp = x - y;

	if (x == y) return x + M_LN2;
    
	if (tmp > 0)
		return x + log1p(exp(-tmp));
    
	return y + log1p(exp(tmp));
}

double logFactorial(int n){
	double logF = 0;
	for (int i = 2; i <= n; i++) {
		logF += log(i);
	}
	return logF;
}

double logDoubleFactorial(int n){
	double logF = 0;
	for (int i = n; i > 1; i-=2) {
		logF += log(i);
	}
	return logF;
}

#pragma mark -
#pragma mark double

double dmax( double a, double b){
	return ( a > b ? a : b);
}

double dmin( double a, double b){
	return ( a < b ? a : b);
}

double dmin_vector( const double *array, int n ){
	double min = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] < min ) min = array[i];
	}
	return min;
}

double dmax_vector( const double *array, int n ){
	double max = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] > max ) max = array[i];
	}
	return max;
}

size_t which_dmin( const double *array, size_t n ){
	double min = array[0];
	size_t index = 0;
	for (size_t i = 1; i < n; i++) {
		if( array[i] < min ){
			min = array[i];
			index = i;
		}
	}
	return index;
}

inline void dswap( double *a, double *b ){
	double c = *a; *a = *b; *b = c;
}

int qsort_desc_dvector( const void *a, const void *b ){
	const double *aa = a;
	const double *bb = b;
	if( *aa < *bb ) return 1;
	else if( *aa == *bb ) return 0;
	return -1;
}

int qsort_asc_dvector( const void *a, const void *b ){
	const double *aa = a;
	const double *bb = b;
	if( *aa > *bb ) return 1;
	else if( *aa == *bb ) return 0;
	return -1;
}

void randomize_dvector( double *v, double n ){
	int rpos = 0;
	for (int i = 0; i < n; i++) {
		rpos = random_int(n-1);
		dswap( &v[rpos], &v[i]);
	}
}

#pragma mark -
#pragma mark int

int imax( int a, int b){
	return ( a > b ? a : b);
}

int imin( int a, int b){
	return ( a < b ? a : b);
}

int imax_vector( const int *array, int n ){
	int max = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] > max ) max = array[i];
	}
	return max;
}

int imin_vector( const int *array, int n ){
	int min = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] < min ) min = array[i];
	}
	return min;
}

inline void swap_int( int *a, int *b ){
	int c = *a; *a = *b; *b = c;
}

int qsort_desc_ivector( const void *a, const void *b ){
	return (*(int*)b - *(int*)a);
}

int qsort_asc_ivector( const void *a, const void *b ){
	return (*(int*)a - *(int*)b);
}

void randomize_ivector( int *v, int n ){
	int rpos = 0;
	for (int i = 0; i < n; i++) {
		rpos = random_int(n-1);
		swap_int( &v[rpos], &v[i]);
	}	
}

void reverse_ivector( int *v, unsigned int n ){
    int swap;
    for( int a = 0; a < --n; a++ ){
        swap = v[a];
        v[a] = v[n];
        v[n] = swap;
    }
}

#pragma mark -
#pragma mark unsigned int

unsigned uimax_vector( const unsigned *array, int n ){
	unsigned max = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] > max ) max = array[i];
	}
	return max;
}


inline void swap_uint( unsigned int *a, unsigned int *b){
	unsigned int c = *a; *a = *b; *b = c;
}

int qsort_desc_uivector( const void *a, const void *b ){
	return (*(unsigned*)b - *(unsigned*)a);
}

int qsort_asc_uivector( const void *a, const void *b ){
	return (*(unsigned*)a - *(unsigned*)b);
}

void randomize_uivector( unsigned int *v, int n ){
	int rpos = 0;
	for (int i = 0; i < n; i++) {
		rpos = random_int(n-1);
		swap_uint( &v[rpos], &v[i]);
	}	
}

// This function convert an array of unsigned integer representing class membership (e.g. 2012031) to a unique representation (e.g. 0120132)
// For example the arrays 2012031 and 3013021 represent the same partioning into 4 groups.
// This algorithm translate these 2 arrays into the canonical representation 0120132
// @param helper_len is the number of classes
// If there is a missing class in vec then helper will have a -1 at index class-1 
void uivector_canonical( unsigned *vec, size_t len_vec, int *helper, size_t helper_len ){
	int i = 0;
	for ( ; i < helper_len; i++ ) {
		helper[i] = -1;
	}
	int c = 0;
	for ( i = 0; i < len_vec; i++ ) {
		if( helper[ vec[i] ] == -1 ){
			helper[ vec[i] ] = c++;
		}
		vec[i] = helper[ vec[i] ] ;
	}
}

unsigned umax_vector( const unsigned *array, int n ){
	unsigned max = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] > max ) max = array[i];
	}
	return max;
}

#pragma mark -
#pragma mark size_t

size_t stmax( size_t a, size_t b ){
	return ( a > b ? a : b);
}

size_t stmin( size_t a, size_t b ){
	return ( a < b ? a : b);
}

#pragma mark -
#pragma mark real

real_t rmax( real_t a, real_t b){
	return ( a > b ? a : b);
}

real_t rmin( real_t a, real_t b){
	return ( a < b ? a : b);
}

real_t rmin_vector( const real_t *array, int n ){
	real_t min = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] < min ) min = array[i];
	}
	return min;
}

real_t rmax_vector( const real_t *array, int n ){
	real_t max = array[0];
	for (int i = 1; i < n; i++) {
		if( array[i] > max ) max = array[i];
	}
	return max;
}

inline void rswap( real_t *a, real_t *b){
	real_t c = *a; *a = *b; *b = c;
}

int qsort_desc_rvector( const void *a, const void *b ){
	const real_t *aa = a;
	const real_t *bb = b;
	if( *aa < *bb ) return 1;
	else if( *aa == *bb ) return 0;
	return -1;
}

int qsort_asc_rvector( const void *a, const void *b ){
	const real_t *aa = a;
	const real_t *bb = b;
	if( *aa > *bb ) return 1;
	else if( *aa == *bb ) return 0;
	return -1;
}



double * linearly_spaced_vector( double lower, double upper, int n ){
    double *v = dvector(n);
    double space = (upper-lower) / (n-1);
    v[0] = 0;
    for ( int i = 1; i < n; i++ ) {
		v[i] = v[i-1] + space;
	}
    return v;
}

double * log_spaced_spaced_vector( double lower, double upper, int n ){
    double *v = dvector(n);
    double logMin = log(lower);
    double logMax = log(upper);
    double delta = (logMax - logMin) / (n-1);
    
    double accDelta = 0;
    for (int i = 0; i < n; ++i) {
        v[i] = exp(logMin + accDelta);
        accDelta += delta;
    }
    return v;
}

void log_spaced_spaced_vector2( double *v, double lower, double upper, int n ){
    double logMin = log(lower);
    double logMax = log(upper);
    double delta = (logMax - logMin) / (n-1);
    
    double accDelta = 0;
    for (int i = 0; i < n; ++i) {
        v[i] = exp(logMin + accDelta);
        accDelta += delta;
    }
}

double * exp_spaced_spaced_vector( double lower, double upper, int n ){
    double *v = dvector(n);
    double expMin = exp(lower);
    double expMax = exp(upper);
    double delta = (expMax - expMin) / (n-1);
    
    double accDelta = 0;
    for (int i = 0; i < n; ++i) {
        v[i] = log(expMin + accDelta);
        accDelta += delta;
    }
    return v;
}

void exp_spaced_spaced_vector2( double *v, double lower, double upper, int n ){
    double expMin = exp(lower);
    double expMax = exp(upper);
    double delta = (expMax - expMin) / (n-1);
    
    double accDelta = 0;
    for (int i = 0; i < n; ++i) {
        v[i] = log(expMin + accDelta);
        accDelta += delta;
    }
}

#pragma mark -
#pragma mark strings

bool array_of_string_contains(const char *str, const char *array[], int count, bool casesensitive){
	if (casesensitive) {
		for ( int i = 0; i < count; i++) {
			if( strcmp(str, array[i]) == 0 ){
				return true;
			}
		}
	}
	else{
		for ( int i = 0; i < count; i++) {
			if( strcasecmp(str, array[i]) == 0 ){
				return true;
			}
		}
	}
	return false;
}

#pragma mark -
#pragma mark Other

void error( const char message[] ){
	fprintf(stderr, "%s\n", message);
	exit(1);
}

inline void shift4( double *a, double *b, double *c, double d ){
	*a = *b; *b = *c; *c = d;
}

inline void rshift4( real_t *a, real_t *b, real_t *c, real_t d ){
	*a = *b; *b = *c; *c = d;
}

inline double SQR( double a ){
	return a == 0.0 ? 0.0 : a*a; 
}


bool file_exists(const char * filename){
	FILE * file = fopen(filename, "r");
	if ( file != NULL ){
        fclose(file);       
		return true;
    }
    return false;
}

// a number could be 1e10
bool isFloat( const char *str ){
	const char *pch = str;
	int npoint = 0;
	while ( *pch != '\0' ) {
		if ( (*pch >= 48 && *pch <= 57) || *pch == '.' ) {
			if ( *pch == '.' ) npoint++;
			if (npoint > 1 ) {
				return false;
			}
		}
		else return false;
		pch++;
	}
	return true;
}


bool isInt( const char *str ){
	const char *pch = str;
	
	if ( *pch == '\0' ) {
		return false;
	}
	
	if( *pch == '-' ){
		pch++;
		if ( pch == '\0') { // str cannot be "-"
			return false;
		}
	}
	
	while ( *pch >= 48 && *pch <= 57 ) {
		pch++;
	}
	return (*pch == '\0');
}

bool isFloat2( const char *str ){
	const char *pch = str;
	
	if ( *pch == '\0' ) {
		return false;
	}	
	
	if( *pch == '-' ){
		pch++;
		if ( pch == '\0') { // str cannot be "-"
			return false;
		}
	}
	
	int npoint = 0;
	while ( (*pch >= 48 && *pch <= 57) || *pch == '.' ) {
		if ( *pch == '.' ) npoint++;
		if (npoint > 1 ) {
			return false;
		}
		pch++;
	}
	
	if( *pch != '\0' ){
		if ( *pch == 'e' || *pch == 'E' ) {
			pch++;
			return isInt(pch);			
		}
		else return false;
	}
	return true;
}


// this one is not tested but should be better than isFloat2
bool isFloat3( const char *str ){
	const char *pch = str;
	
	if ( *pch == '\0' ) {
		return false;
	}
	
	if( *pch == '-' ){
		pch++;
		if ( pch == '\0') { // str cannot be "-"
			return false;
		}
	}
	
	int npoint = 0;
	while ( (*pch >= 48 && *pch <= 57) || *pch == '.' ) {
		if ( *pch == '.' ) npoint++;
		if (npoint > 1 ) {
			return false;
		}
		pch++;
	}
	
	if( *pch != '\0' ){
		if ( *pch == 'e' || *pch == 'E' ) {
			pch++;
			if( *pch == '\0') return false;
            
            if( *pch == '-' || *pch == '+' ){
                pch++;
            }
            
            if( *pch == '\0') return false;
            
            if( *pch >= 48 && *pch <= 57){
                while ( (*pch >= 48 && *pch <= 57) || *pch == '\0' ) {
                    pch++;
                }
                if( *pch == '\0') return false;
            }
            else {
                return false;
            }
		}
		else return false;
	}
	return true;
}

void print_pretty_time( FILE *file, double dseconds ){
	int seconds = (int)dseconds;
	
	int weeks   = seconds / (604800);
	int days    = (seconds  - (weeks*604800) ) / 86400;
	int hours   = (seconds  - (weeks*604800) - (days*86400) ) /  3600;
	int minutes = (seconds  - (weeks*604800) - (days*86400) -  (hours*3600) ) / 60;
	int secs    = seconds  - (weeks*604800) - (days*86400) -  (hours*3600) - (minutes*60);
	
	if ( weeks != 0 ) {
		fprintf(file, "%d week%s: ", weeks, (weeks==1 ? "" : "s") );
	}
	if ( days != 0 ) {
		fprintf(file, "%d day%s: ", days, (days==1 ? "" : "s") );
	}
	if ( hours != 0 ) {
		fprintf(file, "%d hour%s: ", hours, (hours==1 ? "" : "s") );
	}
	if ( minutes != 0 ) {
		fprintf(file, "%d min%s: ", minutes, (minutes==1 ? "" : "s") );
	}
	
	fprintf(file, "%d sec%s\n", secs, (secs==1 ? "" : "s")  );
}


//void * aligned_malloc( size_t alignment, size_t size ){
//	void *memptr = NULL;
//#if defined (__MACH__) || defined (__WIN32__)
//	memptr = malloc(size);
//	assert(memptr);
//#else
//	int r = posix_memalign(&memptr, alignment,  size);
//	assert(r == 0);
//#endif
//	return memptr;
//}

// works for mac, linux and Windows
void * aligned16_malloc( size_t size ){
	void *memptr = NULL;
#if defined (__MACH__) || defined (__WIN32__)
	memptr = malloc(size);
	assert(memptr);
#else
	int r = posix_memalign(&memptr, 16,  size);
    assert(r == 0);
#endif
	return memptr;
}

int cmp_double_int_pair_asc(const void *a, const void *b){
	const double_int_pair_t *const *aa = a;
	const double_int_pair_t *const *bb = b;
	if ((*aa)->value < (*bb)->value)
		return -1;
	else if ((*aa)->value > (*bb)->value)
		return 1;
	else
		return 0;
}

int cmp_double_int_pair_desc(const void *a, const void *b){
	const double_int_pair_t *const *aa = a;
	const double_int_pair_t *const *bb = b;
	if ((*aa)->value > (*bb)->value)
		return -1;
	else if ((*aa)->value < (*bb)->value)
		return 1;
	else
		return 0;
}

