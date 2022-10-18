/*
 *  treelikelihoodCodon.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 27/9/12.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */


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
void update_partials_codon_openmp( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

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
