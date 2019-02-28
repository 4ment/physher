//
//  utilsgsl.h
//  physher
//
//  Created by Mathieu Fourment on 28/2/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef utilsgsl_h
#define utilsgsl_h

#include <gsl/gsl_rng.h>

size_t roulette_wheel_gsl(gsl_rng* rng, const double *array, size_t len );

#endif /* utilsgsl_h */
