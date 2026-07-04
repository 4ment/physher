// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef utilssse_h
#define utilssse_h

#include <stdio.h>

void mult_vector_vector_inplace(double* vec1, const double* vec2, size_t length, size_t chunk);

void mult_vector_vector(double* out, const double* vec1, const double* vec2, size_t length, size_t chunk);

void add_vector_vector_mult_vector(double* out, const double* vec1, const double* vec2, const double* vec3, size_t length, size_t chunk);

#endif /* utilssse_h */
