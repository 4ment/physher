// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef __PhyC__phyresampling__
#define __PhyC__phyresampling__

#include <stdio.h>

#include "sequence.h"
#include "sitepattern.h"

typedef enum resampling_scheme{RESAMPLING_BOOTSTRAP, RESAMPLING_JACKKNIFE, RESAMPLING_JACKKNIFE_PROPORTION}resampling_scheme;

Sequences * Sequences_bootstrap( const Sequences *sequences );

Sequences * Sequences_jackknife( const Sequences *sequences, int index );

Sequences * Sequences_jackknife_n( const Sequences *sequences, int n );



SitePattern * SitePattern_bootstrap( const SitePattern *sitepattern );

SitePattern * SitePattern_jackknife( const SitePattern *original, int index );

SitePattern * SitePattern_jackknife_n( SitePattern *original, int n);

SitePattern * SitePattern_reweight( const SitePattern *original, const double *weights );

#endif /* defined(__PhyC__phyresampling__) */
