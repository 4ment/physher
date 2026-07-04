// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef nonstat_h
#define nonstat_h

#include "substmodel.h"

struct SubstitutionModel* new_NONSTATNucleotideModel();

struct SubstitutionModel* new_NONSTATNucleotideModel_with_parameters(
    Parameter* freqs, const Parameters* rates);

#endif /* nonstat_h */
