/*
 *  discreteclock.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/18/11.
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

#ifndef _DISCRETE_CLOCK_H_
#define _DISCRETE_CLOCK_H_

#include "treelikelihood.h"
#include "heterotachy.h"


ClockSearch * new_GADiscreteClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize );

void DiscreteClock_greedy( SingleTreeLikelihood *tlk );

void DiscreteClock_greedy2( SingleTreeLikelihood *tlk );
#endif
