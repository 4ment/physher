//
//  utilssse.h
//  physher
//
//  Created by Mathieu Fourment on 29/8/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef utilssse_h
#define utilssse_h

#include <stdio.h>

void mult_vector_vector_inplace(double* vec1, const double* vec2, size_t length, size_t chunk);

void mult_vector_vector(double* out, const double* vec1, const double* vec2, size_t length, size_t chunk);

void add_vector_vector_mult_vector(double* out, const double* vec1, const double* vec2, const double* vec3, size_t length, size_t chunk);

#endif /* utilssse_h */
