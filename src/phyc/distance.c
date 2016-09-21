/*
 *  distance.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 5/27/11.
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
