/*
 *  random.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/23/11.
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
