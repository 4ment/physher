/*
 *  localclock.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/12/11.
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

#ifndef _LOCAL_CLOCK_H_
#define _LOCAL_CLOCK_H_

#include "treelikelihood.h"
#include "heterotachy.h"


#pragma mark -
// MARK: Heuristic search

ClockSearch * new_LocalClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads );



#pragma mark -
// MARK: GA

ClockSearch * new_GALocalClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize );

ClockSearch * new_GALocalClockSearch2( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize );

#endif
