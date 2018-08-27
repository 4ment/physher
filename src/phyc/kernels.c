//
//  kernels.c
//  physher
//
//  Created by Mathieu Fourment on 26/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "kernels.h"

#include <math.h>

#include "mathconstant.h"

// h: bandwith
double gaussian_kernel_density(const double* values, size_t length, double x, double h){
	double prob = 1.0/length/SQRT_2PI;
	for (size_t i = 0; i < length; i++) {
		prob += exp(-0.5*pow((x-values[i])/h, 2));
	}
	return prob;
}
