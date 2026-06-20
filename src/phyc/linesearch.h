/*
 *  linesearch.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/11/11.
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


#ifndef _LINESEARCH_H_
#define _LINESEARCH_H_

#include "parameters.h"
#include "optimizer.h"

typedef struct {
    double alpha;
    double fx;
    int nfev;
    int status;
} LineSearchResult;

LineSearchResult strong_wolfe_line_search(Parameters* parameters, opt_func fun,
                                          opt_grad_func grad_f, void* data, const double *x,
                                          const double *p, const double *g, double f0,
                                          double c1, double c2, double alpha0, double amax);

opt_result lnsrch(Parameters *parameters, double* x,  opt_func func, void *data, double fold, double *g, double *p, double *fmin, double stpmax, double alam);

#endif
