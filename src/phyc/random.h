// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _RANDOM_H_
#define _RANDOM_H_

#define nseed init_genrand

// [0,1]
#define random_double genrand_real1

// prototypes for Mersenne Twister pseudorandom generator functions
// file: mt19937ar.c

void init_genrand(unsigned long s);

unsigned long genrand_int32(void);



// [0,u]
int random_int(int u);

// [l,u]
int random_int2(int l, int u);




// [0,1]
double genrand_real1(void);

// [0,1)
double genrand_real2(void);


// [0,n]
double random_double2(double n );

// [0,n)
double random_double3(double n );

// [l,u]
double random_double4( double l, double u );

// sum of array must be equal to 1
int roulette_wheel( const double *array, int len );

int roulette_wheel2( const double *array, int len, double tot );

#endif
