/*
 *  exp.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/1/10.
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

#ifndef _MATRIX_EXP_H_
#define _MATRIX_EXP_H_

#include "matrix.h"
#include "eigen.h"

Matrix *Exp_Taylor_Series( Matrix *m, const int n ); 

double ** Exp_QR_EigenDecomposition( EigenDecomposition *ed );

double ** Exp_QR( double **matrix, const int dim );

#endif
