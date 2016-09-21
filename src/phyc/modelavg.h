/*
 *  modelavg.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 27/6/12.
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

#ifndef Math_modelavg_h
#define Math_modelavg_h

#include <stdio.h>

#include "branchmodel.h"

typedef struct ModelToAverage{
	double IC;
	double weight;
	double *params;
	int n;
} ModelToAverage;

typedef struct ModelAveraged{
	double *mean;
	double *min;
	double *max;
	double n;
} ModelAveraged;


ModelToAverage * new_ModelToAverage( const double IC, const int n );

void free_ModelToAverage( ModelToAverage *m );


Tree *Model_average_from_log( const char *filename, double p, char ***orderedNames, int start, int end );

void free_ModelAveraged( ModelAveraged *ma );



#endif
