// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef gensubst_h
#define gensubst_h

#include <stdio.h>

#include "discreteparameter.h"
#include "substmodel.h"

struct SubstitutionModel* new_GeneralJC69Model_with_parameters(DataType* datatype,
                                                               Parameter* freqs_simplex,
                                                               bool normalize);

struct SubstitutionModel* new_GeneralModel_with_parameters(
    DataType* datatype, DiscreteParameter* model, const Parameters* rates,
    Parameter* freqs, int relativeTo, bool normalize);

void general_dQdp(struct SubstitutionModel* m);

#endif /* gensubst_h */
