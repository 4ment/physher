// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef unrest_h
#define unrest_h

#include "substmodel.h"

struct SubstitutionModel * new_UnrestrictedNucleotideModel();

struct SubstitutionModel* new_UnrestrictedNucleotideModel_with_parameters(
    Parameter* freqs, const Parameters* rates);

#endif /* unrest_h */
