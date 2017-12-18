//
//  marginal.h
//  physher
//
//  Created by Mathieu Fourment on 18/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef marginal_h
#define marginal_h

#include <stdio.h>

#include "matrix.h"

double log_harmonic_mean(const Vector** values);

double log_marginal_stepping_stone(const Vector** values, size_t temp_count, const double* temperatures, double* lrssk);

#endif /* marginal_h */
