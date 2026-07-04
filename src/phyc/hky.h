// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef hky_h
#define hky_h

#include <stdio.h>

#include "substmodel.h"

struct SubstitutionModel* new_HKY(Parameter* freqs);

struct SubstitutionModel* new_HKY_with_values(const double* freqs, const double kappa);

struct SubstitutionModel* new_HKY_with_parameters(Parameter* freqs, Parameter* kappa);

#endif /* hky_h */
