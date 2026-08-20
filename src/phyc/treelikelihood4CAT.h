// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_treelikelihood4CAT_h
#define PhyC_treelikelihood4CAT_h

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood

void update_partials_4_cat(SingleTreeLikelihood *tlk, int partialsIndex,
                           int partialsIndex1, int matrixIndex1, int partialsIndex2,
                           int matrixIndex2);

#pragma mark -
#pragma mark Upper Likelihood

void calculate_branch_partials_4_cat(SingleTreeLikelihood *tlk, double *rootPartials,
                                     int upperPartialsIndex, int partialsIndex,
                                     int matrixIndex);

#pragma mark -
#pragma mark SSE

#ifdef SSE3_ENABLED

void update_partials_4_SSE_cat(SingleTreeLikelihood *tlk, int partialsIndex,
                               int partialsIndex1, int matrixIndex1, int partialsIndex2,
                               int matrixIndex2);

void calculate_branch_partials_4_SSE_cat(SingleTreeLikelihood *tlk, double *rootPartials,
                                         int upperPartialsIndex, int partialsIndex,
                                         int matrixIndex);

#endif

#endif
