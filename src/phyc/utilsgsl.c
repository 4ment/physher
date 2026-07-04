// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
