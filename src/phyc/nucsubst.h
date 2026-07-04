// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef nucsubst_h
#define nucsubst_h

#include "substmodel.h"

struct SubstitutionModel* new_ReversibleNucleotideModel_with_parameters(
    const char* model, Parameter* freqs, const Parameters* rates);

#endif /* nucsubst_h */
