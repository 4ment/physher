// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef gy94_h
#define gy94_h

#include <stdio.h>

#include "substmodel.h"

SubstitutionModel* new_GY94(Parameter* freqs, unsigned gen_code);

SubstitutionModel* new_GY94_with_values(Parameter* freqs, const double omega,
                                        const double kappa, unsigned gen_code);

#endif /* gy94_h */
