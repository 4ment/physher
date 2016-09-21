/*
 *  lm.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/29/10.
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

#ifndef _LINEAR_MODEL_H_
#define _LINEAR_MODEL_H_


void regression( double *x, double *y, int n, double *slope, double *intercept );

void residuals( const double *x, const double *y, int n, double s, double intcp, double *res );
double residual( double x, double y, double s, double intcp );
double sse( double *x, double *y, int n, double s, double intcp );
double mse( double *x, double *y, int n, double s, double intcp);
double var( double *x, double *y, int n, double s, double intcp );

void CI_xIntercept(double *x, double *y, int n, double slope, double intercept, double *lower, double *upper, double alpha );
double CI_slope( double *x, double *y, int n, double slope, double intercept, double prob );
double CI_intercept( double *x, double *y, int n, double slope, double intercept, double prob );

#endif
