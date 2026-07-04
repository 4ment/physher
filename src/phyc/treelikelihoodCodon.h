// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_treelikelihoodCodon_h
#define PhyC_treelikelihoodCodon_h

#include "treelikelihood.h"


void node_log_likelihoods_codon( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

void integrate_partials_codon( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void update_partials_codon( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void partials_states_and_states_codon( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_codon( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined_codon( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

#pragma mark -
#pragma mark OpenMP

#ifdef _OPENMP
void update_partials_codon_openmp( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void partials_states_and_states_codon_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_codon_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined_codon_openmp( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );
#endif

#pragma mark -
#pragma mark SSE3


#ifdef SSE3_ENABLED
void node_log_likelihoods_codon_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

void integrate_partials_codon_SSE( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );


void update_partials_codon_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void partials_states_and_states_codon_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_codon_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined_codon_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );


void update_partials_codon_odd_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void partials_states_and_states_codon_odd_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_codon_odd_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined_codon_odd_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

#endif


#pragma mark -
#pragma mark Upper Likelihood

void update_partials_upper_codon( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_codon( const SingleTreeLikelihood *tlk, Node *node );

#pragma mark -
#pragma mark Upper Likelihood SSE

#ifdef SSE3_ENABLED
void update_partials_upper_sse_codon( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_sse_codon( const SingleTreeLikelihood *tlk, Node *node );
#endif

#endif
