// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_treestat_h
#define PhyC_treestat_h

#include "tree.h"

double TreeStat_rate_correlation( const Tree *tree );

double TreeStat_mean_rate_scaled( const Tree *tree );

double TreeStat_mean_rate_tips_scaled( const Tree *tree );

double TreeStat_mean_rate_internal_scaled( const Tree *tree );

double TreeStat_mean_rate( const Tree *tree, double *min, double *max );

#endif
