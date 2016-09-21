/*
 *  neutralitytest.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 7/2/12.
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


#ifndef Math_neutralitytest_h
#define Math_neutralitytest_h

#include "sequence.h"

double Watterson_theta_estimator( const Sequences *sequences );

double Tajima_D( const Sequences *sequences );

double FuLi_Dstar( const Sequences *sequences );

double FuLi_Fstar( const Sequences *sequences );


#endif
