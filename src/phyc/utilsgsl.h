// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef utilsgsl_h
#define utilsgsl_h

#include <gsl/gsl_rng.h>

size_t roulette_wheel_gsl(gsl_rng* rng, const double *array, size_t len );

#endif /* utilsgsl_h */
