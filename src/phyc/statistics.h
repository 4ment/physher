/*
 *  statistics.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 15/8/12.
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


#ifndef PhyC_statistics_h
#define PhyC_statistics_h

double correlation( double *x, double *y, int dim );

double covariance( const double *x, const double *y, double meanX, double meanY, int dim );

double mean( const double *x, int dim );

double variance( const double *x, int dim, double mean );

double standard_deviation( const double *x, int dim, double mean );

#endif
