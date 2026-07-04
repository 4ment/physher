// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_boot_h
#define PhyC_boot_h


double * boot_ci_norm( const double *thetas, int n, double theta_hat, double alpha );

double boot_ci_norm2( const double *thetas, int n, double theta_hat, double q );

double bootci_BCa(double *thetai, int n, const double *thetab, int b, double theta_hat, double q);

double bootci_BCa_weighted(double *thetai, int n, const double *weights, const double *thetab, int b, double theta_hat, double q);

double bootci_BCa_weighted_debug(double *thetai, int n, const double *weights, const double *thetab, int b, double theta_hat, double q);

#endif
