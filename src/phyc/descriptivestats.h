// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_descriptivestats_h
#define PhyC_descriptivestats_h

double dmean( const double *v, int len );

double dweighted_average( const double *v, const double *weights, int len );

double dmedian( const double *v, int len );

double dmedian_ordered( const double *v, int len );

double dpercentile( const double *v, int len, double percentile );

double dpercentile_ordered( const double *v, int len, double percentile );


int imedian( const int *v, int len );

int imedian_ordered( const int *v, int len );

int ipercentile( const int *v, int len, double percentile );

int ipercentile_ordered( const int *v, int len, double percentile );

#endif
