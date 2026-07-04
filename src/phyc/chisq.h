// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _CHI_SQUARE_H_
#define _CHI_SQUARE_H_

#include "gamma.h"

// Cumulative distribution function
// return P(X <= x)
double pchisq( const double x, const int df );

// Inverse cumulative distribution function
double qchisq( const double p, const int df );

#endif
