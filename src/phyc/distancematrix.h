/*
 *  distancematrix.h
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

#ifndef __PhyC__distancematrix__
#define __PhyC__distancematrix__

#include <stdio.h>

#include "sequence.h"
#include "sitepattern.h"

typedef enum distancematrix_model{
    DISTANCE_MATRIX_UNCORRECTED,
    DISTANCE_MATRIX_JC69,
    DISTANCE_MATRIX_K2P,
    
    DISTANCE_MATRIX_KIMURA
}distancematrix_model;

double ** Sequences_distance( const Sequences *sequences, distancematrix_model model );

float ** Sequences_distance_float( const Sequences *sequences, distancematrix_model model );

double ** SitePattern_distance( const SitePattern *patterns, distancematrix_model model );

#endif /* defined(__PhyC__distancematrix__) */
