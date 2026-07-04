// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "random.h"


#define RAND_MAX32 4294967295.0

// [0,n]
int random_int(int n){
	return random_int2(0, n);
}

// [l,n]
int random_int2(int l, int u){
	//return (rand() %(u+1-l)+l);
	return (genrand_int32() % (u+1-l)+l);
}

/*double random_double(double n){
 return random_double2(0, n);
 }
 
 
 double random_double2(double l, double u){
 return (rand() %(u+1.-l)+l);
 }*/



// [0,n]
double random_double2(double n ){
	//return rand() / (((double)RAND_MAX) / n);
	return n * (genrand_int32() / (double)RAND_MAX32);
}

//[0,n)
double random_double3(double n ){
	//return rand() / (((double)RAND_MAX + 1) / n);
	return n * (genrand_int32() / (double)RAND_MAX32 + 1);
}

// [l,u]
double random_double4( double l, double u ){
	return l+(u-l)*(genrand_int32() / ((double)RAND_MAX32));
}


// array has to sum to 1
int roulette_wheel( const double *array, int len ){
    double accum = 0.0;
    double rnum = random_double();
    int i = 0;
    for ( ; i < len; i++ ) {
        accum += array[i];
        if( accum >= rnum ) break;
    }
    return i;
}

int roulette_wheel2( const double *array, int len, double tot ){
    double accum = 0.0;
    double rnum = random_double2(tot);
    int i = 0;
    for ( ; i < len; i++ ) {
        accum += array[i];
        if( accum >= rnum ) break;
    }
    return i;
}

