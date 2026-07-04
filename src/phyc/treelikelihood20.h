// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_treelikelihood20_h
#define PhyC_treelikelihood20_h

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood SSE

#ifdef SSE3_ENABLED

void update_partials_20_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void partials_states_and_states_20_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_20_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, double *partials );

void partials_states_and_undefined_20_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_20_SSE( const SingleTreeLikelihood *tlk, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined_20_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

#endif

#pragma mark -
#pragma mark Upper Likelihood

void calculate_branch_partials_20(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);

#pragma mark -
#pragma mark Upper Likelihood SSE

#ifdef SSE3_ENABLED
void calculate_branch_partials_20_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);
#endif

#endif
