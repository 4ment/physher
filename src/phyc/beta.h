// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef beta_h
#define beta_h

#include <stdio.h>

double invbetai(double p, double a, double b);

double dbetaprime(double x, double alpha, double beta);

double dlogbetaprime(double x, double alpha, double beta);

#endif /* beta_h */
