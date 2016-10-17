/*
 *  treelikelihood4.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 24/9/12.
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

#ifndef PhyC_treelikelihood4_h
#define PhyC_treelikelihood4_h

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood

void update_partials_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void update_partials_uncertainty_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void integrate_partials_4( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_4( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods);


void partials_undefined_and_undefined_4( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_states_and_undefined_4( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_states_and_states_4( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_uncertainty_4( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_states_and_states_uncertainty_4( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );


void update_partials_4_ancestral( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );




void update_partials_noexp_integrate_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void partials_undefined_and_undefined_noexp_integrate_4( const SingleTreeLikelihood *tlk, const double *partials1, const double *exps1, const double *partials2, const double *exps2, double *partials3);

void partials_states_and_undefined_noexp_integrate_4( const SingleTreeLikelihood *tlk, int idx1, const double *exps1, const double *partials2, const double *exps2, double *partials3 );

void partials_states_and_states_noexp_integrate_4( const SingleTreeLikelihood *tlk, int idx1, const double *exps1, int idx2, const double *exps2, double *partials );


#pragma mark -
#pragma mark OpenMP

void update_partials_4_openmp( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void partials_states_and_states_4_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_4_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3);

void partials_undefined_and_undefined_4_openmp( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3);


#pragma mark -
#pragma mark SSE

#ifdef SSE3_ENABLED
void update_partials_4_SSE( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void update_partials_uncertainty_4_SSE( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void integrate_partials_4_SSE( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_4_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );
#endif


#pragma mark -
#pragma mark AVX

#ifdef AVX_ENABLED
void update_partials_4_AVX( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );
void update_partials_uncertainty_4_AVX( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void integrate_partials_4_AVX( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials );

void node_log_likelihoods_4_AVX( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods );
#endif


#pragma mark -
#pragma mark Upper Likelihood

void update_partials_upper_4( SingleTreeLikelihood *tlk, Node *node );

void update_partials_uncertainty_upper_4( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_4( const SingleTreeLikelihood *tlk, Node *node );

#pragma mark -
#pragma mark Upper Likelihood SSE

void update_partials_upper_sse_4( SingleTreeLikelihood *tlk, Node *node );

void update_partials_uncertainty_upper_sse_4( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_sse_4( const SingleTreeLikelihood *tlk, Node *node );


#endif
