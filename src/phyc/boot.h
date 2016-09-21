/*
 *  boot.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 6/12/2013.
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

#ifndef PhyC_boot_h
#define PhyC_boot_h


double * boot_ci_norm( const double *thetas, int n, double theta_hat, double alpha );

double boot_ci_norm2( const double *thetas, int n, double theta_hat, double q );

double bootci_BCa(double *thetai, int n, const double *thetab, int b, double theta_hat, double q);

double bootci_BCa_weighted(double *thetai, int n, const double *weights, const double *thetab, int b, double theta_hat, double q);

double bootci_BCa_weighted_debug(double *thetai, int n, const double *weights, const double *thetab, int b, double theta_hat, double q);

#endif
