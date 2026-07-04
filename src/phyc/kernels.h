// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef kernels_h
#define kernels_h

#include <stdio.h>

double gaussian_kernel_density(const double* values, size_t length, double x, double h);

#endif /* kernels_h */
