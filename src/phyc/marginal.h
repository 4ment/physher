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

double log_marginal_path_sampling(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk);

double log_marginal_path_sampling_modified(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk);

#endif /* marginal_h */
