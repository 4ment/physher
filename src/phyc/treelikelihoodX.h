// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_treelikelihoodX_h
#define PhyC_treelikelihoodX_h

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood

void update_partials_general( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void integrate_partials_general( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_general( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

void partials_states_and_states( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined( const SingleTreeLikelihood *tlk, const double *partials2, const double *matrices2, double *partials3 );
void partials_states( const SingleTreeLikelihood *tlk, int idx1, const double * matrices1, double * partials );

void update_partials_general_transpose( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );
void partials_states_and_states_transpose( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );
void partials_states_and_undefined_transpose( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );


void partials_undefined_and_undefined_omp( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );


#pragma mark -
#pragma mark Lower Likelihood SSE

void update_partials_general_even_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void integrate_partials_general_even_SSE( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_general_even_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

#pragma mark -
#pragma mark Upper Likelihood

void calculate_branch_likelihood(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);

void calculate_branch_partials(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);

void node_log_likelihoods_upper_general( const SingleTreeLikelihood *tlk, Node *node );

#pragma mark -
#pragma mark Upper Likelihood SSE

void calculate_branch_likelihood_even_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);

void calculate_branch_partials_even_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);

#endif
