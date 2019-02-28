//
//  utilsgsl.c
//  physher
//
//  Created by Mathieu Fourment on 28/2/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "utilsgsl.h"

// array has to sum to 1
size_t roulette_wheel_gsl(gsl_rng* rng, const double *array, size_t len ){
	double accum = 0.0;
	double rnum = gsl_rng_uniform(rng);
	size_t i = 0;
	for ( ; i < len; i++ ) {
		accum += array[i];
		if( accum >= rnum ) break;
	}
	return i;
}
