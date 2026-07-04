// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef mg94_h
#define mg94_h

#include <stdio.h>

#include "substmodel.h"

SubstitutionModel* new_MG94(Parameter* freqs, unsigned gen_code);

SubstitutionModel* new_MG94_with_values(Parameter* freqs, const double alpha,
                                        const double beta, const double kappa,
                                        unsigned gen_code);

SubstitutionModel* new_MG94_with_parameters(Parameter* freqs, Parameter* alpha,
                                            Parameter* beta, Parameter* kappa,
                                            unsigned gen_code);

#endif /* mg94_h */
