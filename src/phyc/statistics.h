// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_statistics_h
#define PhyC_statistics_h

double correlation( double *x, double *y, int dim );

double covariance( const double *x, const double *y, double meanX, double meanY, int dim );

double mean( const double *x, int dim );

double variance( const double *x, int dim, double mean );

double standard_deviation( const double *x, int dim, double mean );

#endif
