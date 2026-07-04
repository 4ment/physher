// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
