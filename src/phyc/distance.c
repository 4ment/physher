// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "distance.h"

#include "sitepattern.h"
#include "matrix.h"


double ** calculate_distance_matrix( SitePattern *sp ){
	int j = 0;
	int k = 0;
	double **matrix = dmatrix(sp->size, sp->size);
	for ( int i = 0; i < sp->size; i++ ) {
		for ( j = i+1; j < sp->size; j++ ) {
			for ( k = 0; j < sp->count; k++ ) {
				// FIXME: if there is a gap or state unknown (X) it should be different
				if ( sp->patterns[i][k] != sp->patterns[j][k] ) {
					matrix[i][j] += sp->weights[k];
				}
			}
			matrix[i][j] /= sp->nsites;
			matrix[i][j] = matrix[j][i];
		}
	}
	return matrix;
}
