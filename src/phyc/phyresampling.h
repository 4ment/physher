/*
 *  phyclustering.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/01/2015.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

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
