/*
 *  treelikelihoodX.h
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


#ifndef PhyC_treelikelihoodX_h
#define PhyC_treelikelihoodX_h

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood

void update_partials_general( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void integrate_partials_general( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_general( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

void node_likelihoods_general( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

void partials_states_and_states( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );


void update_partials_general_transpose( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );
void partials_states_and_states_transpose( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );
void partials_states_and_undefined_transpose( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );


void partials_undefined_and_undefined_omp( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );


#pragma mark -
#pragma mark Lower Likelihood SSE

void update_partials_general_even_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 );

void integrate_partials_general_even_SSE( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_general_even_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

void node_likelihoods_general_even_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );

#pragma mark -
#pragma mark Upper Likelihood

void calculate_branch_likelihood(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex);

void update_partials_upper_general( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_general( const SingleTreeLikelihood *tlk, Node *node );


#pragma mark -
#pragma mark ASR

void update_partials_general_ancestral( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );


void update_partials_general_ancestral_rel( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );
void integrate_partials_general_rel( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void update_partials_general_ancestral_relb( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );


#pragma mark -
#pragma mark Lower Likelihood REL

void update_partials_general_rel( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

#endif
