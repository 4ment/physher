// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef Math_neutralitytest_h
#define Math_neutralitytest_h

#include "sequence.h"

double Watterson_theta_estimator( const Sequences *sequences );

double Tajima_D( const Sequences *sequences );

double FuLi_Dstar( const Sequences *sequences );

double FuLi_Fstar( const Sequences *sequences );


#endif
