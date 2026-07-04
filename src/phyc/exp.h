// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _MATRIX_EXP_H_
#define _MATRIX_EXP_H_

#include "matrix.h"
#include "eigen.h"

Matrix *Exp_Taylor_Series( Matrix *m, const int n ); 

double ** Exp_QR_EigenDecomposition( EigenDecomposition *ed );

double ** Exp_QR( double **matrix, const int dim );

#endif
