// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef dirichlet_h
#define dirichlet_h

#include <stdio.h>

#include "gamma.h"

// Density function

double ddirchlet( const double *x, const size_t dim, const double *alphas );

double ddirchletln( const double *x, const size_t dim, const double *alphas );

// Flat dirichlet with alpha == (1,1,..,1)
double ddirchlet_flat( const size_t dim );

void rdirichlet(double*x, const size_t dim, const double* alphas);

#endif /* dirichlet_h */
