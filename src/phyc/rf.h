// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later


#ifndef Math_rf_h
#define Math_rf_h


double rf_distance( bool **s1, bool **s2, int splitCount, int length );

double rf_norm_distance( bool **s1, bool **s2, int splitCount, int length );

double Branch_score( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length );

double K_tree_score( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length );

#endif
