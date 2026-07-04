// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef Math_roottotip_h
#define Math_roottotip_h

double * lm_tree( Tree *tree, bool forward, bool use_correlation );

double ** max_lm_tree( Tree *tree, bool forward, bool use_correlation );

double * lm_tree_cluster( Tree *tree, bool forward, int k );

#endif
