// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_treelikelihood4_h
#define PhyC_treelikelihood4_h

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood

void update_partials_flexible_4( SingleTreeLikelihood *tlk, double* partials, int partialsIndex1, double* partials1, double* matrix1, int partialsIndex2, double* partials2, double* matrix2 );

void update_partials_4( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

//void update_partials_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void integrate_partials_4( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_4( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods);

void partials_undefined_and_undefined_4( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_states_and_undefined_4( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_states_and_states_4( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );


void partials_undefined_4( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, double *partials3 );

void partials_states_4( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, double *partials );


#pragma mark -
#pragma mark OpenMP

#ifdef _OPENMP
void update_partials_4_openmp( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void partials_states_and_states_4_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_4_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3);

void partials_undefined_and_undefined_4_openmp( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3);
#endif

#pragma mark -
#pragma mark SSE

#ifdef SSE3_ENABLED

void update_partials_flexible_4_SSE( SingleTreeLikelihood *tlk, double* partials, int partialsIndex1, double* partials1, double* matrix1, int partialsIndex2, double* partials2, double* matrix2 );

void update_partials_4_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void integrate_partials_4_SSE( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_4_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

void calculate_branch_partials_4_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);
#endif


#pragma mark -
#pragma mark AVX

#ifdef AVX_ENABLED
void update_partials_4_AVX( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void integrate_partials_4_AVX( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_4_AVX( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );
#endif


#pragma mark -
#pragma mark Upper Likelihood

void calculate_branch_partials_4(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);

void update_partials_upper_4( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_4( const SingleTreeLikelihood *tlk, Node *node );

#pragma mark -
#pragma mark Upper Likelihood SSE

#ifdef SSE3_ENABLED
void update_partials_upper_sse_4( SingleTreeLikelihood *tlk, Node *node );
#endif

#endif
