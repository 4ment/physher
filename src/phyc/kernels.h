//
//  kernels.h
//  physher
//
//  Created by Mathieu Fourment on 26/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef kernels_h
#define kernels_h

#include <stdio.h>

double gaussian_kernel_density(const double* values, size_t length, double x, double h);

#endif /* kernels_h */
