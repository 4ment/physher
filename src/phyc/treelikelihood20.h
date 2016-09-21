/*
 *  treelikelihood20.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 17/09/2014.
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

#ifndef PhyC_treelikelihood20_h
#define PhyC_treelikelihood20_h

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood SSE

#ifdef SSE3_ENABLED

void update_partials_20_SSE( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 );

void partials_states_and_states_20_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials );

void partials_states_and_undefined_20_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

void partials_undefined_and_undefined_20_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 );

#endif

#pragma mark -
#pragma mark Upper Likelihood

void update_partials_upper_20( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_20( const SingleTreeLikelihood *tlk, Node *node );

#pragma mark -
#pragma mark Upper Likelihood SSE

#ifdef SSE3_ENABLED
void update_partials_upper_sse_20( SingleTreeLikelihood *tlk, Node *node );

void node_log_likelihoods_upper_sse_20( const SingleTreeLikelihood *tlk, Node *node );
#endif

#endif
