// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef K80_h
#define K80_h

#include <stdio.h>

#include "substmodel.h"

struct SubstitutionModel* new_K80_with_values(const double kappa);

struct SubstitutionModel* new_K80_with_parameters(Parameter* kappa);

#endif /* K80_h */
