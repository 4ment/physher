// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef discreteoperator_h
#define discreteoperator_h

#include <stdio.h>

#include "operator.h"

#include <gsl/gsl_rng.h>

bool operator_discrete_bitflip(Operator* op, double* logHR);

bool operator_discrete_exchange(Operator* op, double* logHR);

#endif /* discreteoperator_h */
