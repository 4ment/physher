/*
 *  treelikelihood4.c
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

#include "treelikelihood4.h"

#include <stdio.h>
#include <math.h>


#ifdef SSE3_ENABLED
#if defined(__aarch64__)
#include "neon2sse.h"
#else
#include <xmmintrin.h> // SSE
#include <pmmintrin.h> // SSE3
//#include <tmmintrin.h> // SSSE3
#endif
#endif

#if 0
#define RESTRICT restrict
#else
#define RESTRICT
#endif

#include "treelikelihood.h"

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
#define GNUC47
#endif

#pragma mark -
#pragma mark Lower Likelihood


void update_partials_flexible_4( SingleTreeLikelihood *tlk, double* partials, int partialsIndex1, double* partials1, double* matrix1, int partialsIndex2, double* partials2, double* matrix2 ) {
	
	if(partialsIndex2 < 0){
		if(  tlk->partials[0][partialsIndex1] != NULL ){
			partials_undefined_4(tlk,
								 partials1,
								 matrix1,
								 partials);
		}
		else{
			
			partials_states_4(tlk,
							  tlk->mapping[partialsIndex1],
							  matrix1,
							  partials);
		}
	}
	else if( tlk->partials[0][partialsIndex1] != NULL ){
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_undefined_and_undefined_4(tlk,
											   partials1,
											   matrix1,
											   partials2,
											   matrix2,
											   partials);
		}
		else {
			partials_states_and_undefined_4(tlk,
											tlk->mapping[partialsIndex2],
											matrix2,
											partials1,
											matrix1,
											partials);
		}
	}
	else{
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_states_and_undefined_4(tlk,
											tlk->mapping[partialsIndex1],
											matrix1,
											partials2,
											matrix2,
											partials);
		}
		else{
			partials_states_and_states_4(tlk,
										 tlk->mapping[partialsIndex1],
										 matrix1,
										 tlk->mapping[partialsIndex2],
										 matrix2,
										 partials);
		}
	}
	
//	if ( tlk->scale ) {
//		SingleTreeLikelihood_scalePartials( tlk, partialsIndex);
//	}
}

void update_partials_4( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {

	if(partialsIndex2 < 0){
		if(  tlk->partials[0][partialsIndex1] != NULL ){
			partials_undefined_4(tlk,
							   tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
							   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
							   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		else{
			
			partials_states_4(tlk,
							 tlk->mapping[partialsIndex1],
							 tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
							 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
	}
	else if( tlk->partials[0][partialsIndex1] != NULL ){
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_undefined_and_undefined_4(tlk,
											   tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											   tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		else {
			partials_states_and_undefined_4(tlk,
											tlk->mapping[partialsIndex2],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_states_and_undefined_4(tlk,
											tlk->mapping[partialsIndex1],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
			
		}
		else{
			partials_states_and_states_4(tlk,
										 tlk->mapping[partialsIndex1],
										 tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
										 tlk->mapping[partialsIndex2],
										 tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
										 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
	}
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex, partialsIndex1, partialsIndex2);
	}
}

void integrate_partials_4( const SingleTreeLikelihood * tlk, const double * RESTRICT inPartials, const double * RESTRICT proportions, double * RESTRICT outPartials ){
	int k;
	double *pPartials = outPartials;
	const double *pInPartials = inPartials;
    double prop = proportions[0];
	
    // should be faster than multiplying by 1
    if( tlk->cat_count == 1 ){
        memcpy(outPartials, inPartials, tlk->pattern_count*4*sizeof(double));
    }
    else {
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            *pPartials++ = *pInPartials++ * prop;
            *pPartials++ = *pInPartials++ * prop;
            *pPartials++ = *pInPartials++ * prop;
            *pPartials++ = *pInPartials++ * prop;
        }
        
        
        for ( int l = 1; l < tlk->cat_count; l++ ) {
            pPartials = outPartials;
            prop = proportions[l];
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                *pPartials += *pInPartials++ * prop;
                pPartials++;
                
                *pPartials += *pInPartials++ * prop;
                pPartials++;
                
                *pPartials += *pInPartials++ * prop;
                pPartials++;
                
                *pPartials += *pInPartials++ * prop;
                pPartials++;
            }
        }
    }
}

void node_log_likelihoods_4( const SingleTreeLikelihood *tlk, const double * RESTRICT partials, const double * RESTRICT frequencies, double * RESTRICT outLogLikelihoods ){
	int v = 0;
	for ( int k = 0; k < tlk->pattern_count; k++ ) {
		
		outLogLikelihoods[k] = 	frequencies[0] * partials[v]; v++;
		outLogLikelihoods[k] += frequencies[1] * partials[v]; v++;
		outLogLikelihoods[k] += frequencies[2] * partials[v]; v++;
		outLogLikelihoods[k] += frequencies[3] * partials[v]; v++;
		
		outLogLikelihoods[k] = log(outLogLikelihoods[k]);
		
		if ( tlk->scale ) {
			outLogLikelihoods[k] += getLogScalingFactor( tlk, k);
		}
	}
}

void partials_states_and_states_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT matrices2, double * RESTRICT partials ){
	int k,w;
	int u = 0;
	int state1, state2;
	double *pPartials = partials;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
        for ( k = 0; k < tlk->pattern_count; k++ ) {
			state1 = tlk->sp->patterns[idx1][k];
			state2 = tlk->sp->patterns[idx2][k];
			
			//w = l * tlk->matrix_size;
			w = u;
			
			if (state1 < 4 && state2 < 4) {
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2];
				
			}
			else if (state1 < 4 ) {
				// child 1 has a gap or unknown state so treat it as unknown
				*pPartials++ = matrices1[w + state1]; w += 4;
				*pPartials++ = matrices1[w + state1]; w += 4;
				*pPartials++ = matrices1[w + state1]; w += 4;
				*pPartials++ = matrices1[w + state1];
			}
			else if (state2 < 4 ) {
				// child 2 has a gap or unknown state so treat it as unknown
				*pPartials++ = matrices2[w + state2]; w += 4;
				*pPartials++ = matrices2[w + state2]; w += 4;
				*pPartials++ = matrices2[w + state2]; w += 4;
				*pPartials++ = matrices2[w + state2];
				
			}
			else {
				// both children have a gap or unknown state so set partials to 1
				memcpy(pPartials, TWENTY_DOUBLE_ONES, sizeof(double)*4);
				pPartials += 4;
			}
		}
		u += 16;
	}
}

// compute e^{tQ}.partials
void partials_states_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, double * RESTRICT partials ){
	int k,w;
	int u = 0;
	int state1;
	double *pPartials = partials;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			state1 = tlk->sp->patterns[idx1][k];
			
			w = u;
			
			if (state1 < 4 ) {
				*pPartials++ = matrices1[w + state1]; w += 4;
				*pPartials++ = matrices1[w + state1]; w += 4;
				*pPartials++ = matrices1[w + state1]; w += 4;
				*pPartials++ = matrices1[w + state1];
			}
			else {
				// P.[1 1 1 1]^T = [1 1 1 1] when P is probability matrix
				//memcpy(pPartials, TWENTY_DOUBLE_ONES, sizeof(double)*4);
				// When derivatives are calculated P is NOT a probability matrix so we do the calculation
				*pPartials++ = matrices1[0] + matrices1[1] + matrices1[2] + matrices1[3];
				*pPartials++ = matrices1[4] + matrices1[5] + matrices1[6] + matrices1[7];
				*pPartials++ = matrices1[8] + matrices1[9] + matrices1[10] + matrices1[11];
				*pPartials++ = matrices1[12] + matrices1[13] + matrices1[14] + matrices1[15];
			}
		}
		u += 16;
	}
}

// Auto-vectorization
#ifdef GNUC47
void partials_states_and_undefined_4( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3){
	const double *partials2 = __builtin_assume_aligned(apartials2, 16);
	const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
	const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#else
void partials_states_and_undefined_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, const double * RESTRICT partials2, const double * RESTRICT matrices2, double * RESTRICT partials3){
#endif
	
	double sum;
	int v = 0;
	int k;
	int w = 0;
	int state1;
	
	double *pPartials = partials3;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[idx1][k];
			
			if ( state1 < 4) {
				
				sum  = matrices2[w]     * partials2[v];
				sum += matrices2[w + 1] * partials2[v + 1];
				sum += matrices2[w + 2] * partials2[v + 2];
				sum += matrices2[w + 3] * partials2[v + 3];
				
				*pPartials++ = matrices1[w + state1] * sum;
				
				
				sum  = matrices2[w + 4] * partials2[v];
				sum += matrices2[w + 5] * partials2[v + 1];
				sum += matrices2[w + 6] * partials2[v + 2];
				sum += matrices2[w + 7] * partials2[v + 3];
				
				*pPartials++ = matrices1[w + 4 + state1] * sum;
				
				
				sum  = matrices2[w + 8]  * partials2[v];
				sum += matrices2[w + 9]  * partials2[v + 1];
				sum += matrices2[w + 10] * partials2[v + 2];
				sum += matrices2[w + 11] * partials2[v + 3];
				
				*pPartials++ = matrices1[w + 8 + state1] * sum;
				
				
				sum  = matrices2[w + 12] * partials2[v];
				sum += matrices2[w + 13] * partials2[v + 1];
				sum += matrices2[w + 14] * partials2[v + 2];
				sum += matrices2[w + 15] * partials2[v + 3];
				
				*pPartials++ = matrices1[w + 12 + state1] * sum;
				
			}
			else {
				// Child 1 has a gap or unknown state so don't use it
				
				*pPartials  = matrices2[w]     * partials2[v];
				*pPartials += matrices2[w + 1] * partials2[v + 1];
				*pPartials += matrices2[w + 2] * partials2[v + 2];
				*pPartials += matrices2[w + 3] * partials2[v + 3];
				
				pPartials++;
				
				*pPartials  = matrices2[w + 4] * partials2[v];
				*pPartials += matrices2[w + 5] * partials2[v + 1];
				*pPartials += matrices2[w + 6] * partials2[v + 2];
				*pPartials += matrices2[w + 7] * partials2[v + 3];
				
				pPartials++;
				
				*pPartials  = matrices2[w + 8]  * partials2[v];
				*pPartials += matrices2[w + 9]  * partials2[v + 1];
				*pPartials += matrices2[w + 10] * partials2[v + 2];
				*pPartials += matrices2[w + 11] * partials2[v + 3];
				
				pPartials++;
				
				*pPartials  = matrices2[w + 12] * partials2[v];
				*pPartials += matrices2[w + 13] * partials2[v + 1];
				*pPartials += matrices2[w + 14] * partials2[v + 2];
				*pPartials += matrices2[w + 15] * partials2[v + 3];
				
				pPartials++;
				
			}
			v += 4;
		}
		w += 16;
	}
}

#ifdef GNUC47
void partials_undefined_and_undefined_4( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double * restrict partials3 ){
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#else
void partials_undefined_and_undefined_4( const SingleTreeLikelihood *tlk, const double * RESTRICT partials1, const double * RESTRICT matrices1, const double * RESTRICT partials2, const double * RESTRICT matrices2, double * RESTRICT partials3){
#endif
    double sum1, sum2;
    int v = 0;
    int k;
    int w = 0;
	
    double *pPartials = partials3;
    int cat_count = tlk->cat_count;
    int sp_count = tlk->pattern_count;
    
    for ( int l = 0; l < cat_count; l++ ) {
        
        for ( k = 0; k < sp_count; k++ ) {
            
            sum1  = matrices1[w]   * partials1[v];
            sum2  = matrices2[w]   * partials2[v];
            
            sum1 += matrices1[w+1] * partials1[v + 1];
            sum2 += matrices2[w+1] * partials2[v + 1];
            
            sum1 += matrices1[w+2] * partials1[v + 2];
            sum2 += matrices2[w+2] * partials2[v + 2];
            
            sum1 += matrices1[w+3] * partials1[v + 3];
            sum2 += matrices2[w+3] * partials2[v + 3];
            
            *pPartials++ = sum1 * sum2;
            
            
            sum1  = matrices1[w+4] * partials1[v];
            sum2  = matrices2[w+4] * partials2[v];
            
            sum1 += matrices1[w+5] * partials1[v + 1];
            sum2 += matrices2[w+5] * partials2[v + 1];
            
            sum1 += matrices1[w+6] * partials1[v + 2];
            sum2 += matrices2[w+6] * partials2[v + 2];
            
            sum1 += matrices1[w+7] * partials1[v + 3];
            sum2 += matrices2[w+7] * partials2[v + 3];
            
            *pPartials++ = sum1 * sum2;
            
            
            sum1  = matrices1[w+8]  * partials1[v];
            sum2  = matrices2[w+8]  * partials2[v];
            
            sum1 += matrices1[w+9]  * partials1[v + 1];
            sum2 += matrices2[w+9]  * partials2[v + 1];
            
            sum1 += matrices1[w+10] * partials1[v + 2];
            sum2 += matrices2[w+10] * partials2[v + 2];
            
            sum1 += matrices1[w+11] * partials1[v + 3];
            sum2 += matrices2[w+11] * partials2[v + 3];
            
            *pPartials++ = sum1 * sum2;
            
            
            sum1  = matrices1[w+12] * partials1[v];
            sum2  = matrices2[w+12] * partials2[v];
            
            sum1 += matrices1[w+13] * partials1[v + 1];
            sum2 += matrices2[w+13] * partials2[v + 1];
            
            sum1 += matrices1[w+14] * partials1[v + 2];
            sum2 += matrices2[w+14] * partials2[v + 2];
            
            sum1 += matrices1[w+15] * partials1[v + 3];
            sum2 += matrices2[w+15] * partials2[v + 3];
            
            *pPartials++ = sum1 * sum2;
            
            v += 4;
        }
        w += 16;
    }
    
}
	
#ifdef GNUC47
void partials_undefined_4( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, double * restrict partials3 ){
	const double *partials1 = __builtin_assume_aligned(apartials1, 16);
	const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
#else
void partials_undefined_4( const SingleTreeLikelihood *tlk, const double * RESTRICT partials1, const double * RESTRICT matrices1, double * RESTRICT partials3){
#endif
	double sum1;
	int v = 0;
	int k;
	int w = 0;
	
	double *pPartials = partials3;
	int cat_count = tlk->cat_count;
	int sp_count = tlk->pattern_count;
	
	for ( int l = 0; l < cat_count; l++ ) {
		
		for ( k = 0; k < sp_count; k++ ) {
			
			sum1  = matrices1[w]   * partials1[v];
			sum1 += matrices1[w+1] * partials1[v + 1];
			sum1 += matrices1[w+2] * partials1[v + 2];
			sum1 += matrices1[w+3] * partials1[v + 3];
			
			*pPartials++ = sum1;
			
			
			sum1  = matrices1[w+4] * partials1[v];
			sum1 += matrices1[w+5] * partials1[v + 1];
			sum1 += matrices1[w+6] * partials1[v + 2];
			sum1 += matrices1[w+7] * partials1[v + 3];
			
			*pPartials++ = sum1;
			
			
			sum1  = matrices1[w+8]  * partials1[v];
			sum1 += matrices1[w+9]  * partials1[v + 1];
			sum1 += matrices1[w+10] * partials1[v + 2];
			sum1 += matrices1[w+11] * partials1[v + 3];
			
			*pPartials++ = sum1;
			
			
			sum1  = matrices1[w+12] * partials1[v];
			sum1 += matrices1[w+13] * partials1[v + 1];
			sum1 += matrices1[w+14] * partials1[v + 2];
			sum1 += matrices1[w+15] * partials1[v + 3];
			
			*pPartials++ = sum1;
			
			v += 4;
		}
		w += 16;
	}
	
}

#pragma mark -
#pragma mark OpenMP

#ifdef _OPENMP
void update_partials_4_openmp( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
	
	if( tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1] != NULL ){
		if(  tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2] != NULL ){
			partials_undefined_and_undefined_4_openmp(tlk,
												   tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		else {
			partials_states_and_undefined_4_openmp(tlk,
												tlk->mapping[partialsIndex2],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
												tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2] != NULL ){
			partials_states_and_undefined_4_openmp(tlk,
												tlk->mapping[partialsIndex1],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
												tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
			
		}
		else{
			partials_states_and_states_4_openmp(tlk,
											 tlk->mapping[partialsIndex1],
											 tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											 tlk->mapping[partialsIndex2],
											 tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
	}
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex, partialsIndex1, partialsIndex2);
	}
}

void partials_states_and_states_4_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
	
    int nThreads = tlk->nthreads;
    
    #pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for ( int lk = 0; lk < tlk->pattern_count*tlk->cat_count; lk++ ) {
        int l = lk / tlk->pattern_count;
        int k = lk % tlk->pattern_count;
        
        int state1 = tlk->sp->patterns[idx1][k];
        int state2 = tlk->sp->patterns[idx2][k];
        
        int w = l * tlk->matrix_size;
        
        double *pPartials = partials + (l*tlk->pattern_count + k)*4;
        
        if (state1 < 4 && state2 < 4) {
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
            
        }
        else if (state1 < 4 ) {
            // child 1 has a gap or unknown state so treat it as unknown
            *pPartials++ = matrices1[w + state1]; w += 4;
            *pPartials++ = matrices1[w + state1]; w += 4;
            *pPartials++ = matrices1[w + state1]; w += 4;
            *pPartials++ = matrices1[w + state1]; w += 4;
        }
        else if (state2 < 4 ) {
            // child 2 has a gap or unknown state so treat it as unknown
            *pPartials++ = matrices2[w + state2]; w += 4;
            *pPartials++ = matrices2[w + state2]; w += 4;
            *pPartials++ = matrices2[w + state2]; w += 4;
            *pPartials++ = matrices2[w + state2]; w += 4;
            
        }
        else {
            // both children have a gap or unknown state so set partials to 1
            *pPartials++ = 1.0;
            *pPartials++ = 1.0;
            *pPartials++ = 1.0;
            *pPartials++ = 1.0;
        }
    }
}


void partials_states_and_undefined_4_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3){
	
    int nThreads = tlk->nthreads;
    
    #pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for ( int lk = 0; lk < tlk->pattern_count*tlk->cat_count; lk++ ) {
        int l = lk / tlk->pattern_count;
        int k = lk % tlk->pattern_count;
        
        int w = l * tlk->matrix_size;
        int v = (l*tlk->pattern_count + k) * 4;
        
        int state1 = tlk->sp->patterns[idx1][k];
        
        double sum;
    
        double *pPartials = partials3+v;
        
        if ( state1 < 4) {
            
            sum  = matrices2[w]     * partials2[v];
            sum += matrices2[w + 1] * partials2[v + 1];
            sum += matrices2[w + 2] * partials2[v + 2];
            sum += matrices2[w + 3] * partials2[v + 3];
            
            *pPartials++ = matrices1[w + state1] * sum;
            
            
            sum  = matrices2[w + 4] * partials2[v];
            sum += matrices2[w + 5] * partials2[v + 1];
            sum += matrices2[w + 6] * partials2[v + 2];
            sum += matrices2[w + 7] * partials2[v + 3];
            
            *pPartials++ = matrices1[w + 4 + state1] * sum;
            
            
            sum  = matrices2[w + 8]  * partials2[v];
            sum += matrices2[w + 9]  * partials2[v + 1];
            sum += matrices2[w + 10] * partials2[v + 2];
            sum += matrices2[w + 11] * partials2[v + 3];
            
            *pPartials++ = matrices1[w + 8 + state1] * sum;
            
            
            sum  = matrices2[w + 12] * partials2[v];
            sum += matrices2[w + 13] * partials2[v + 1];
            sum += matrices2[w + 14] * partials2[v + 2];
            sum += matrices2[w + 15] * partials2[v + 3];
            
            *pPartials++ = matrices1[w + 12 + state1] * sum;
            
        }
        else {
            // Child 1 has a gap or unknown state so don't use it
            
            *pPartials  = matrices2[w]     * partials2[v];
            *pPartials += matrices2[w + 1] * partials2[v + 1];
            *pPartials += matrices2[w + 2] * partials2[v + 2];
            *pPartials += matrices2[w + 3] * partials2[v + 3];
            
            pPartials++;
            
            *pPartials  = matrices2[w + 4] * partials2[v];
            *pPartials += matrices2[w + 5] * partials2[v + 1];
            *pPartials += matrices2[w + 6] * partials2[v + 2];
            *pPartials += matrices2[w + 7] * partials2[v + 3];
            
            pPartials++;
            
            *pPartials  = matrices2[w + 8]  * partials2[v];
            *pPartials += matrices2[w + 9]  * partials2[v + 1];
            *pPartials += matrices2[w + 10] * partials2[v + 2];
            *pPartials += matrices2[w + 11] * partials2[v + 3];
            
            pPartials++;
            
            *pPartials  = matrices2[w + 12] * partials2[v];
            *pPartials += matrices2[w + 13] * partials2[v + 1];
            *pPartials += matrices2[w + 14] * partials2[v + 2];
            *pPartials += matrices2[w + 15] * partials2[v + 3];
            
            pPartials++;
            
        }
    }
    
}


void partials_undefined_and_undefined_4_openmp( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3){
    
    int nThreads = tlk->nthreads;
    
    #pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for ( int lk = 0; lk < tlk->pattern_count*tlk->cat_count; lk++ ) {
        int l = lk / tlk->pattern_count;
        int k = lk % tlk->pattern_count;
        
        int w = l * tlk->matrix_size;
        int v = (l*tlk->pattern_count + k) * 4;
        
        double sum1, sum2;
        
        double *pPartials = partials3 + v;
        
        sum1  = matrices1[w]   * partials1[v];
        sum2  = matrices2[w]   * partials2[v];
        
        sum1 += matrices1[w+1] * partials1[v + 1];
        sum2 += matrices2[w+1] * partials2[v + 1];
        
        sum1 += matrices1[w+2] * partials1[v + 2];
        sum2 += matrices2[w+2] * partials2[v + 2];
        
        sum1 += matrices1[w+3] * partials1[v + 3];
        sum2 += matrices2[w+3] * partials2[v + 3];
        
        *pPartials++ = sum1 * sum2;
        
        
        sum1  = matrices1[w+4] * partials1[v];
        sum2  = matrices2[w+4] * partials2[v];
        
        sum1 += matrices1[w+5] * partials1[v + 1];
        sum2 += matrices2[w+5] * partials2[v + 1];
        
        sum1 += matrices1[w+6] * partials1[v + 2];
        sum2 += matrices2[w+6] * partials2[v + 2];
        
        sum1 += matrices1[w+7] * partials1[v + 3];
        sum2 += matrices2[w+7] * partials2[v + 3];
        
        *pPartials++ = sum1 * sum2;
        
        
        sum1  = matrices1[w+8]  * partials1[v];
        sum2  = matrices2[w+8]  * partials2[v];
        
        sum1 += matrices1[w+9]  * partials1[v + 1];
        sum2 += matrices2[w+9]  * partials2[v + 1];
        
        sum1 += matrices1[w+10] * partials1[v + 2];
        sum2 += matrices2[w+10] * partials2[v + 2];
        
        sum1 += matrices1[w+11] * partials1[v + 3];
        sum2 += matrices2[w+11] * partials2[v + 3];
        
        *pPartials++ = sum1 * sum2;
        
        
        sum1  = matrices1[w+12] * partials1[v];
        sum2  = matrices2[w+12] * partials2[v];
        
        sum1 += matrices1[w+13] * partials1[v + 1];
        sum2 += matrices2[w+13] * partials2[v + 1];
        
        sum1 += matrices1[w+14] * partials1[v + 2];
        sum2 += matrices2[w+14] * partials2[v + 2];
        
        sum1 += matrices1[w+15] * partials1[v + 3];
        sum2 += matrices2[w+15] * partials2[v + 3];
        
        *pPartials++ = sum1 * sum2;
        
    }
}
#endif

#pragma mark -
#pragma mark Lower Likelihood SSE

#ifdef SSE3_ENABLED

void integrate_partials_4_SSE( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials ){
	int k;
	double *pPartials = outPartials;
	const double *pInPartials = inPartials;
	
	double prop[2] __attribute__ ((aligned (16)));
	__m128d in, pn;
	__m128d *pSSE = (__m128d*)outPartials;
    
	prop[0] = prop[1] = proportions[0];
	pn = _mm_load_pd(prop);
    
    // should be faster than multiplying by 1
    if( tlk->cat_count == 1 ){
        memcpy(outPartials, inPartials, tlk->pattern_count*4*sizeof(double));
    }
    else {
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            in = _mm_load_pd(pInPartials);
            _mm_store_pd(pPartials, _mm_mul_pd(in, pn));
            pPartials   += 2;
            pInPartials += 2;
            
            in = _mm_load_pd(pInPartials);
            _mm_store_pd(pPartials, _mm_mul_pd(in, pn));
            pPartials   += 2;
            pInPartials += 2;
        }
    
	
	
        for ( int l = 1; l < tlk->cat_count; l++ ) {
            prop[0] = prop[1] = proportions[l];
            pn = _mm_load_pd(prop);
            pPartials = outPartials;
            pSSE = (__m128d*)outPartials;
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                
                in = _mm_load_pd(pInPartials);
                in = _mm_mul_pd(in, pn);
                _mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
                pSSE++;
                pInPartials += 2;
                pPartials   += 2;
                
                
                in = _mm_load_pd(pInPartials);
                in = _mm_mul_pd(in, pn);
                _mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
                pSSE++;
                pInPartials += 2;
                pPartials   += 2;
                
            }
        }
    }
}

void node_log_likelihoods_4_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods ){
    
	const double *pInPartials = partials;
	double *pOutPartials = outLogLikelihoods;
	
	double temp[2] __attribute__ ((aligned (16)));
    
	__m128d f0 = _mm_load_pd(&frequencies[0]);
	__m128d f2 = _mm_load_pd(&frequencies[2]);
	__m128d in0, in2;
	
	for ( int k = 0; k < tlk->pattern_count; k++ ) {
		
		in0 = _mm_load_pd(pInPartials);
		pInPartials += 2;
		in2 = _mm_load_pd(pInPartials);
		pInPartials += 2;
		
		in0 = _mm_mul_pd(in0, f0);
		in2 = _mm_mul_pd(in2, f2);
        
		in0 = _mm_add_pd(in0, in2);
		
		_mm_store_pd(temp, in0);
		
		*pOutPartials = log(temp[0]+temp[1]);
		
		if ( tlk->scale ) {
			*pOutPartials += getLogScalingFactor( tlk, k);
		}
        
        pOutPartials++;
	}
}
    
#ifdef GNUC47
void partials_states_and_states_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, int idx2, const double * restrict amatrices2, double *partials ){
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#else
void partials_states_and_states_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
#endif
    
	int k;
	int state1, state2;
	int u = 0;
	int w;
    
	const double *m1 = matrices1;
	const double *m2 = matrices2;
    double *pPartials = partials;
    
	__m128d m1v0, m1v2,m2v0, m2v2;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[idx1][k];
			state2 = tlk->sp->patterns[idx2][k];
			
            // w = tlk->sm->get_site_category(tlk->sm, k)*16 + u;
			w = u;
			
			if (state1 < 4 && state2 < 4) {
				
				m1 = &matrices1[w+4*state1];
				m2 = &matrices2[w+4*state2];
				
				m1v0 = _mm_load_pd(&m1[0]);
				m1v2 = _mm_load_pd(&m1[2]);
				
				m2v0 = _mm_load_pd(&m2[0]);
				m2v2 = _mm_load_pd(&m2[2]);
                
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				_mm_store_pd(pPartials, _mm_mul_pd(m1v2, m2v2));
				pPartials += 2;
                
			}
			else if (state1 < 4) {
				// child 1 has a gap or unknown state so treat it as unknown
				m1 = &matrices1[w+4*state1];
                
				*pPartials++ = *m1++;
				*pPartials++ = *m1++;
				*pPartials++ = *m1++;
				*pPartials++ = *m1;
			}
			else if (state2 < 4 ) {
				// child 2 has a gap or unknown state so treat it as unknown
				m2 = &matrices2[w+4*state2];
                
				*pPartials++ = *m2++;
				*pPartials++ = *m2++;
				*pPartials++ = *m2++;
				*pPartials++ = *m2;
				
			}
			else {
				// both children have a gap or unknown state so set partials to 1
				*pPartials++ = 1.0;
				*pPartials++ = 1.0;
				*pPartials++ = 1.0;
				*pPartials++ = 1.0;
			}
		}
		u += 16;
	}
}

	
#ifdef GNUC47
void partials_states_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, double *partials ){
	const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
#else
void partials_states_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, double *partials ){
#endif
	
	int k;
	int state1;
	int u = 0;
	
	double *pPartials = partials;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[idx1][k];
			// u = tlk->sm->get_site_category(tlk->sm, k)*16 + l*16;
			
			if (state1 < 4) {
				memcpy(pPartials, &matrices1[u+4*state1], 4*sizeof(double));
				pPartials += 4;
			}
			else {
				// P.[1 1 1 1]^T = [1 1 1 1] when P is probability matrix
				// When derivatives are calculated P is NOT a probability matrix so we do the calculation
				*pPartials++ = matrices1[0] + matrices1[4] + matrices1[8] + matrices1[12];
				*pPartials++ = matrices1[1] + matrices1[5] + matrices1[9] + matrices1[13];
				*pPartials++ = matrices1[2] + matrices1[6] + matrices1[10] + matrices1[14];
				*pPartials++ = matrices1[3] + matrices1[7] + matrices1[11] + matrices1[15];
			}
		}
		u += 16;
	}
}
        
#ifdef GNUC47
void partials_states_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3 ){
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#else
void partials_states_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
#endif
	
	int v = 0;
	int k;
	int w = 0;
	int state1;
	
	__m128d p2v0, p2v2, m2v0, m2v2, *m1, temp;
	__m128d* pPartials = (__m128d*)partials3;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[idx1][k];
			// w = tlk->sm->get_site_category(tlk->sm, k)*16 + l*16;
			
			
			p2v0 = _mm_load_pd(&partials2[v]);
			p2v2 = _mm_load_pd(&partials2[v+2]);
			
			if ( state1 < 4) {
				m1 = (__m128d*)&matrices1[w + state1*4];
				
				
				m2v0 = _mm_load_pd(&matrices2[w]);
				m2v2 = _mm_load_pd(&matrices2[w+2]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				temp = _mm_add_pd(m2v0, m2v2);

				
				m2v0 = _mm_load_pd(&matrices2[w+4]);
				m2v2 = _mm_load_pd(&matrices2[w+6]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				m2v2 = _mm_unpacklo_pd(temp,m2v0);
				temp = _mm_unpackhi_pd(temp,m2v0);
				*pPartials++ = _mm_mul_pd(*m1, _mm_add_pd(m2v2, temp));
				m1++;
				
				
				
				m2v0 = _mm_load_pd(&matrices2[w+8]);
				m2v2 = _mm_load_pd(&matrices2[w+10]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				temp = _mm_add_pd(m2v0, m2v2);

				
				m2v0 = _mm_load_pd(&matrices2[w+12]);
				m2v2 = _mm_load_pd(&matrices2[w+14]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				m2v2 = _mm_unpacklo_pd(temp,m2v0);
				temp = _mm_unpackhi_pd(temp,m2v0);
				*pPartials++ = _mm_mul_pd(*m1, _mm_add_pd(m2v2, temp));
			}
			else {
				// Child 1 has a gap or unknown state so don't use it
				
				m2v0 = _mm_load_pd(&matrices2[w]);
				m2v2 = _mm_load_pd(&matrices2[w+2]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				temp = _mm_add_pd(m2v0, m2v2);
				
				
				m2v0 = _mm_load_pd(&matrices2[w+4]);
				m2v2 = _mm_load_pd(&matrices2[w+6]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);

				m2v2 = _mm_unpacklo_pd(temp,m2v0);
				temp = _mm_unpackhi_pd(temp,m2v0);
				*pPartials++ = _mm_add_pd(m2v2, temp);
				
				
				m2v0 = _mm_load_pd(&matrices2[w+8]);
				m2v2 = _mm_load_pd(&matrices2[w+10]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				temp = _mm_add_pd(m2v0, m2v2);
				
				
				m2v0 = _mm_load_pd(&matrices2[w+12]);
				m2v2 = _mm_load_pd(&matrices2[w+14]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				m2v2 = _mm_unpacklo_pd(temp,m2v0);
				temp = _mm_unpackhi_pd(temp,m2v0);
				*pPartials++ = _mm_add_pd(m2v2, temp);
				
			}
			v += 4;
		}
		w += 16;
	}
}

            // Auto-vectorization
#ifdef GNUC47
void partials_undefined_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3 ){
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#else
void partials_undefined_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
#endif
	
	int v = 0;
	int k;
	int w = 0;
	
	__m128d* pPartials = (__m128d*)partials3;
	__m128d m1v0, m1v2, m2v0, m2v2, p1v0, p1v2, p2v0, p2v2, temp;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			// w = tlk->sm->get_site_category(tlk->sm, k)*16 + l*16;
			
			p1v0 = _mm_load_pd(&partials1[v]);
			p1v2 = _mm_load_pd(&partials1[v+2]);
			
			p2v0 = _mm_load_pd(&partials2[v]);
			p2v2 = _mm_load_pd(&partials2[v+2]);
			
			
			
			m1v0 = _mm_load_pd(&matrices1[w]);
			m1v2 = _mm_load_pd(&matrices1[w+2]);
			
			m2v0 = _mm_load_pd(&matrices2[w]);
			m2v2 = _mm_load_pd(&matrices2[w+2]);
			
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			
			m1v0 = _mm_add_pd(m1v0, m1v2);
			
			m2v0 = _mm_mul_pd(m2v0, p2v0);
			m2v2 = _mm_mul_pd(m2v2, p2v2);
			
			m2v0 = _mm_add_pd(m2v0, m2v2);
			
			temp = _mm_hadd_pd(m1v0, m2v0);
			
			
			m1v0 = _mm_load_pd(&matrices1[w+4]);
			m1v2 = _mm_load_pd(&matrices1[w+6]);
			
			m2v0 = _mm_load_pd(&matrices2[w+4]);
			m2v2 = _mm_load_pd(&matrices2[w+6]);
			
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			
			m1v0 = _mm_add_pd(m1v0, m1v2);
			
			m2v0 = _mm_mul_pd(m2v0, p2v0);
			m2v2 = _mm_mul_pd(m2v2, p2v2);
			
			m2v0 = _mm_add_pd(m2v0, m2v2);
			
			m1v0 = _mm_hadd_pd(m1v0, m2v0);
			
			m2v0 = _mm_unpacklo_pd(temp,m1v0);
			temp = _mm_unpackhi_pd(temp,m1v0);
			*pPartials++ = _mm_mul_pd(m2v0, temp);
			
			
			m1v0 = _mm_load_pd(&matrices1[w+8]);
			m1v2 = _mm_load_pd(&matrices1[w+10]);
			
			m2v0 = _mm_load_pd(&matrices2[w+8]);
			m2v2 = _mm_load_pd(&matrices2[w+10]);
			
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			
			m1v0 = _mm_add_pd(m1v0, m1v2);
			
			m2v0 = _mm_mul_pd(m2v0, p2v0);
			m2v2 = _mm_mul_pd(m2v2, p2v2);
			
			m2v0 = _mm_add_pd(m2v0, m2v2);
			
			temp = _mm_hadd_pd(m1v0, m2v0);
			
			
			m1v0 = _mm_load_pd(&matrices1[w+12]);
			m1v2 = _mm_load_pd(&matrices1[w+14]);
			
			m2v0 = _mm_load_pd(&matrices2[w+12]);
			m2v2 = _mm_load_pd(&matrices2[w+14]);
			
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			
			m1v0 = _mm_add_pd(m1v0, m1v2);
			
			m2v0 = _mm_mul_pd(m2v0, p2v0);
			m2v2 = _mm_mul_pd(m2v2, p2v2);
			
			m2v0 = _mm_add_pd(m2v0, m2v2);
			
			m1v0 = _mm_hadd_pd(m1v0, m2v0);
			
			m2v0 = _mm_unpacklo_pd(temp,m1v0);
			temp = _mm_unpackhi_pd(temp,m1v0);
			*pPartials++ = _mm_mul_pd(m2v0, temp);
			
			v += 4;
		}
		w += 16;
	}
}

#ifdef GNUC47
void partials_undefined_4_SSE( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, double *partials3 ){
	const double *partials1 = __builtin_assume_aligned(apartials1, 16);
	const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
#else
void partials_undefined_4_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, double *partials3 ){
#endif
	
	int v = 0;
	int k;
	int w = 0;
	
	__m128d *pPartials = (__m128d*)partials3;
	__m128d m1v0, m1v2, p1v0, p1v2, temp;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			// w = tlk->sm->get_site_category(tlk->sm, k)*16 + l*16;
			
			p1v0 = _mm_load_pd(&partials1[v]);
			p1v2 = _mm_load_pd(&partials1[v+2]);
			
			
			m1v0 = _mm_load_pd(&matrices1[w]);
			m1v2 = _mm_load_pd(&matrices1[w+2]);
			
			
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			temp = _mm_add_pd(m1v0, m1v2);
			
			
			m1v0 = _mm_load_pd(&matrices1[w+4]);
			m1v2 = _mm_load_pd(&matrices1[w+6]);
			
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			m1v0 = _mm_add_pd(m1v0, m1v2);
			
			m1v2 = _mm_unpacklo_pd(temp,m1v0);
			temp = _mm_unpackhi_pd(temp,m1v0);
			*pPartials++ = _mm_add_pd(m1v2, temp);
			
			
			m1v0 = _mm_load_pd(&matrices1[w+8]);
			m1v2 = _mm_load_pd(&matrices1[w+10]);
			
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			temp = _mm_add_pd(m1v0, m1v2);
			
			
			m1v0 = _mm_load_pd(&matrices1[w+12]);
			m1v2 = _mm_load_pd(&matrices1[w+14]);
			m1v0 = _mm_mul_pd(m1v0, p1v0);
			m1v2 = _mm_mul_pd(m1v2, p1v2);
			
			m1v0 = _mm_add_pd(m1v0, m1v2);
			
			m1v2 = _mm_unpacklo_pd(temp,m1v0);
			temp = _mm_unpackhi_pd(temp,m1v0);
			*pPartials++ = _mm_add_pd(m1v2, temp);
			
			v += 4;
		}
		w += 16;
	}
}

void update_partials_flexible_4_SSE( SingleTreeLikelihood *tlk, double* partials, int partialsIndex1, double* partials1, double* matrix1, int partialsIndex2, double* partials2, double* matrix2 ) {
	
	if(partialsIndex2 < 0){
		if(  tlk->partials[0][partialsIndex1] != NULL ){
			partials_undefined_4_SSE(tlk,
								 partials1,
								 matrix1,
								 partials);
		}
		else{
			
			partials_states_4_SSE(tlk,
							  tlk->mapping[partialsIndex1],
							  matrix1,
							  partials);
		}
	}
	else if( tlk->partials[0][partialsIndex1] != NULL ){
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_undefined_and_undefined_4_SSE(tlk,
											   partials1,
											   matrix1,
											   partials2,
											   matrix2,
											   partials);
		}
		else {
			partials_states_and_undefined_4_SSE(tlk,
											tlk->mapping[partialsIndex2],
											matrix2,
											partials1,
											matrix1,
											partials);
		}
	}
	else{
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_states_and_undefined_4_SSE(tlk,
											tlk->mapping[partialsIndex1],
											matrix1,
											partials2,
											matrix2,
											partials);
		}
		else{
			partials_states_and_states_4_SSE(tlk,
										 tlk->mapping[partialsIndex1],
										 matrix1,
										 tlk->mapping[partialsIndex2],
										 matrix2,
										 partials);
		}
	}
	
	//	if ( tlk->scale ) {
	//		SingleTreeLikelihood_scalePartials( tlk, partialsIndex);
	//	}
}
	
void update_partials_4_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
	
	if(partialsIndex2 < 0){
		if(  tlk->partials[0][partialsIndex1] != NULL ){
			partials_undefined_4_SSE(tlk,
												   tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		else{
			partials_states_4_SSE(tlk,
											 tlk->mapping[partialsIndex1],
											 tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
	}
	else if( tlk->partials[0][partialsIndex1] != NULL ){
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_undefined_and_undefined_4_SSE(tlk,
											   tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											   tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		else {
			partials_states_and_undefined_4_SSE(tlk,
											tlk->mapping[partialsIndex2],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[0][partialsIndex2] != NULL ){
			partials_states_and_undefined_4_SSE(tlk,
											tlk->mapping[partialsIndex1],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
			
		}
		else{
			partials_states_and_states_4_SSE(tlk,
										 tlk->mapping[partialsIndex1],
										 tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
										 tlk->mapping[partialsIndex2],
										 tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
										 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
	}
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex, partialsIndex1, partialsIndex2);
	}
}
	
static void _calculate_branch_likelihood_undefined_SSE(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count);
	int v = 0;
	__m128d* m;
	__m128d p1, p2, temp;
	
	double t[2] __attribute__ ((aligned (16)));
	
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			m = (__m128d*)&matrices[l*16];
			p1 = _mm_load_pd(partials+v);
			p2 = _mm_load_pd(partials+v+2);
			
			temp = _mm_mul_pd(*m, p1); m++;
			temp = _mm_add_pd(temp, _mm_mul_pd(*m, p2)); m++;
			_mm_store_pd(t, temp);
			rootPartials[u++] += (t[0]+t[1]) * upperPartials[v++] * weight;
			
			temp = _mm_mul_pd(*m, p1); m++;
			temp = _mm_add_pd(temp, _mm_mul_pd(*m, p2)); m++;
			_mm_store_pd(t, temp);
			rootPartials[u++] += (t[0]+t[1]) * upperPartials[v++] * weight;
			
			temp = _mm_mul_pd(*m, p1); m++;
			temp = _mm_add_pd(temp, _mm_mul_pd(*m, p2)); m++;
			_mm_store_pd(t, temp);
			rootPartials[u++] += (t[0]+t[1]) * upperPartials[v++] * weight;
			
			temp = _mm_mul_pd(*m, p1); m++;
			temp = _mm_add_pd(temp, _mm_mul_pd(*m, p2)); m++;
			_mm_store_pd(t, temp);
			rootPartials[u++] += (t[0]+t[1]) * upperPartials[v++] * weight;
		}
	}
}

static void _calculate_branch_likelihood_upper_undefined_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count);
	
	double t[2] __attribute__ ((aligned (16)));
	__m128d temp;
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			__m128d p1 = _mm_load_pd(partials+v);
			__m128d p2 = _mm_load_pd(partials+v+2);
			int state = tlk->sp->patterns[upperPartialsIndex][k];
			
			if(state < 4){
				__m128d* m = (__m128d*)&matrices[l * 16 + state*4];
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2));
				_mm_store_pd(t,temp);
				rootPartials[u+state] += (t[0]+t[1]) * weight;
				u+=4;
			}
			else{
				__m128d* m = (__m128d*)&matrices[l*16];
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
			}
			v+=4;
		}
	}
}

// matrices are transposed
static void _calculate_branch_likelihood_state_SSE(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, int partialsIndex, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*4*tlk->pattern_count);
	int v = 0;
	__m128d* up = (__m128d*)upperPartials;
	__m128d weight;
	__m128d* rp, *m1, *m2;
	for(int l = 0; l < tlk->cat_count; l++) {
		rp = (__m128d*)rootPartials;
		double ww = tlk->sm->get_proportion(tlk->sm, l);
		weight = _mm_set1_pd(ww);
		const double* mat = &matrices[l*16];
		int u = 0;
		for(int k = 0; k < tlk->pattern_count; k++) {
			const int state = tlk->sp->patterns[partialsIndex][k];
			if(state < 4){
				m1 = (__m128d*)&mat[state*4];
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; up++;
				u+=4;v+=4;
			}
			else{
				const double* transMatrixPtr = &matrices[l*16];
				rootPartials[u++] += (transMatrixPtr[0] + transMatrixPtr[4] + transMatrixPtr[8] + transMatrixPtr[12]) * upperPartials[v++] * ww;
				rootPartials[u++] += (transMatrixPtr[1] + transMatrixPtr[5] + transMatrixPtr[9] + transMatrixPtr[13]) * upperPartials[v++] * ww;
				rootPartials[u++] += (transMatrixPtr[2] + transMatrixPtr[6] + transMatrixPtr[10] + transMatrixPtr[14]) * upperPartials[v++] * ww;
				rootPartials[u++] += (transMatrixPtr[3] + transMatrixPtr[7] + transMatrixPtr[11] + transMatrixPtr[15]) * upperPartials[v++] * ww;
				
//				m1 = (__m128d*)mat;
//				m2 = (__m128d*)&mat[4];
//				__m128 temp = _mm_add_pd(*m1,*m2); m1+=2; m2+=2;
//				temp = _mm_add_pd(temp, _mm_add_pd(*m1,*m2));
//				temp = _mm_mul_pd(*up, _mm_mul_pd(temp, weight)); up++;
//				*rp = _mm_add_pd(*rp, temp); rp++;
//				
//				m1 = (__m128d*)&mat[2];
//				m2 = (__m128d*)&mat[6];
//				temp = _mm_add_pd(*m1,*m2); m1+=2; m2+=2;
//				temp = _mm_add_pd(temp, _mm_add_pd(*m1,*m2));
//				temp = _mm_mul_pd(*up, _mm_mul_pd(temp, weight)); up++;
//				*rp = _mm_add_pd(*rp, temp); rp++;
				
				rp += 2;
				up += 2;
			}
		}
	}
	
}

void calculate_branch_likelihood_4_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex){
	// partialIndex is a taxon so upperPartialsIndex has to be internal
	// matrices are transposed
	if( tlk->partials[0][partialsIndex] == NULL ){
		_calculate_branch_likelihood_state_SSE(tlk, rootPartials,
											   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
											   tlk->mapping[partialsIndex],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	// upperPartialsIndex is a taxon
	// possible for the child of the root with a leaf sibling
	else if( tlk->partials[0][upperPartialsIndex] == NULL ){
		_calculate_branch_likelihood_upper_undefined_SSE(tlk, rootPartials,
														 tlk->mapping[upperPartialsIndex],
														 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
														 tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	else{
		_calculate_branch_likelihood_undefined_SSE(tlk, rootPartials,
												   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
}
	
static void _calculate_branch_partials_undefined_SSE(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*4*tlk->pattern_count*tlk->cat_count);
	int v = 0;
	__m128d* m;
	__m128d p1, p2, temp, temp2;
	
	for(int l = 0; l < tlk->cat_count; l++) {
		for(int k = 0; k < tlk->pattern_count; k++) {
			m = (__m128d*)&matrices[tlk->sm->get_site_category(tlk->sm, k)*16 + l*16];
			p1 = _mm_load_pd(&partials[v]);
			p2 = _mm_load_pd(&partials[v+2]);
			
			temp = _mm_mul_pd(*m, p1); m++;
			temp = _mm_add_pd(temp, _mm_mul_pd(*m, p2)); m++;
			
			temp2 = _mm_mul_pd(*m, p1); m++;
			temp2 = _mm_add_pd(temp2, _mm_mul_pd(*m, p2)); m++;
			
			temp = _mm_hadd_pd(temp, temp2);
			temp2 = _mm_load_pd(&upperPartials[v]);
			_mm_store_pd(&rootPartials[v], _mm_mul_pd(temp, temp2));
			v+=2;
			
			temp = _mm_mul_pd(*m, p1); m++;
			temp = _mm_add_pd(temp, _mm_mul_pd(*m, p2)); m++;
			
			temp2 = _mm_mul_pd(*m, p1); m++;
			temp2 = _mm_add_pd(temp2, _mm_mul_pd(*m, p2)); m++;
			
			temp = _mm_hadd_pd(temp, temp2);
			temp2 = _mm_load_pd(&upperPartials[v]);
			_mm_store_pd(&rootPartials[v], _mm_mul_pd(temp, temp2));
			v+=2;
		}
	}
}

static void _calculate_branch_partials_upper_undefined_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*4*tlk->pattern_count*tlk->cat_count);
	
	double t[2] __attribute__ ((aligned (16)));
	__m128d temp;
	int v = 0;
	int u = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		for(int k = 0; k < tlk->pattern_count; k++) {
			__m128d p1 = _mm_load_pd(&partials[v]);
			__m128d p2 = _mm_load_pd(&partials[v+2]);
			int state = tlk->sp->patterns[upperPartialsIndex][k];
			
			if(state < 4){
				__m128d* m = (__m128d*)&matrices[tlk->sm->get_site_category(tlk->sm, k)*16 + l * 16 + state*4];
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2));
				_mm_store_pd(t,temp);
				rootPartials[u+state] = t[0] + t[1];
				u+=4;
			}
			else{
				__m128d* m = (__m128d*)&matrices[tlk->sm->get_site_category(tlk->sm, k)*16 + l*16];
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] = t[0] + t[1];
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] = t[0] + t[1];
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] = t[0] + t[1];
				
				temp = _mm_mul_pd(*m, p1); m++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,p2)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] = t[0] + t[1];
			}
			v+=4;
		}
	}
}

// matrices are transposed
static void _calculate_branch_partials_state_SSE(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, int partialsIndex, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*4*tlk->pattern_count*tlk->cat_count);
	__m128d* up = (__m128d*)upperPartials;
	__m128d* rp = (__m128d*)rootPartials;
	__m128d* m1, *m2;
	for(int l = 0; l < tlk->cat_count; l++) {
		for(int k = 0; k < tlk->pattern_count; k++) {
			const double* transMatrixPtr = &matrices[tlk->sm->get_site_category(tlk->sm, k)*16 + l*16];
			const int state = tlk->sp->patterns[partialsIndex][k];
			if(state < 4){
				m1 = (__m128d*)&transMatrixPtr[state*4];
				*rp = _mm_mul_pd(*m1, *up); rp++; m1++; up++;
				*rp = _mm_mul_pd(*m1, *up); rp++; up++;
			}
			else{
				m1 = (__m128d*)transMatrixPtr;
				m2 = (__m128d*)&transMatrixPtr[4];
				*rp = _mm_add_pd(*m1,*m2); m1+=4; m2+=4;
				*rp = _mm_add_pd(*rp, _mm_add_pd(*m1,*m2));
				*rp = _mm_mul_pd(*up, *rp); up++; rp++;

				m1 = (__m128d*)&transMatrixPtr[2];
				m2 = (__m128d*)&transMatrixPtr[6];
				*rp = _mm_add_pd(*m1,*m2); m1+=4; m2+=4;
				*rp = _mm_add_pd(*rp, _mm_add_pd(*m1,*m2));
				*rp = _mm_mul_pd(*up, *rp); up++; rp++;
			}
		}
	}
	
}

void calculate_branch_partials_4_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex){
	// partialIndex is a taxon so upperPartialsIndex has to be internal
	// matrices are transposed
	if( tlk->partials[0][partialsIndex] == NULL ){
		_calculate_branch_partials_state_SSE(tlk, rootPartials,
											   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
											   tlk->mapping[partialsIndex],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	// upperPartialsIndex is a taxon
	// possible for the child of the root with a leaf sibling
	else if( tlk->partials[0][upperPartialsIndex] == NULL ){
		_calculate_branch_partials_upper_undefined_SSE(tlk, rootPartials,
														 tlk->mapping[upperPartialsIndex],
														 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
														 tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	else{
		_calculate_branch_partials_undefined_SSE(tlk, rootPartials,
												   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
}
#endif

#pragma mark -
#pragma mark AVX

//#define AVX_ENABLED

#ifdef AVX_ENABLED
#include <immintrin.h>

void integrate_partials_4_AVX( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials ){
    int k;
    double *pPartials = outPartials;
    const double *pInPartials = inPartials;
    
    __m256d in, pn;
    //__m256d *pSSE = NULL;
    __m256d p;
    
    pn = _mm256_set1_pd(proportions[0]);
    
    for ( k = 0; k < tlk->pattern_count; k++ ) {
        in = _mm256_load_pd(pInPartials);
        _mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
        
        pPartials   += 4;
        pInPartials += 4;
    }
    
    
    
    for ( int l = 1; l < tlk->cat_count; l++ ) {
        pn = _mm256_set1_pd(proportions[l]);
        
        pPartials = outPartials;
        //pSSE = (__m256d*)outPartials;
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            in = _mm256_load_pd(pInPartials);
            in = _mm256_mul_pd(in, pn);
            p = _mm256_load_pd(pPartials);
            _mm256_store_pd( pPartials, _mm256_add_pd(in, p) );
            //_mm256_store_pd( pPartials, _mm256_add_pd(in, *pSSE) );
            //pSSE++;
            pInPartials += 4;
            pPartials   += 4;
        }
    }
}


    // frequencies is not aligned,using set_pd instead of load_pd
    void node_log_likelihoods_4_AVX2( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods ){
        
        const double *pInPartials = partials;
        double *pOutPartials = outLogLikelihoods;
        
        double temp[4] __attribute__ ((aligned (32)));
        //memcpy(temp,frequencies,4*sizeof(double));
        //__m256d f0 = _mm256_load_pd(temp);
        __m256d f0 = _mm256_set_pd(frequencies[3],frequencies[2],frequencies[1],frequencies[0]);
        __m256d in0,in1;
        __m128d sum;
        
        for ( int k = 0; k < tlk->pattern_count; k+=2 ) {
            
            in0 = _mm256_load_pd(pInPartials);
            pInPartials += 4;
            
            in0 = _mm256_mul_pd(in0, f0);
            
            
            in1 = _mm256_load_pd(pInPartials);
            pInPartials += 4;
            
            in1 = _mm256_mul_pd(in1, f0);
            
            in0 = _mm256_hadd_pd(in0, in1);
			
            sum = _mm_add_pd(_mm256_extractf128_pd(in0, 0), _mm256_extractf128_pd(in0, 1));
            
            _mm_store_pd(pOutPartials, sum);
            
            pOutPartials[0] = log(pOutPartials[0]);
            pOutPartials[1] = log(pOutPartials[1]);
            
            if ( tlk->scale ) {
                pOutPartials[0] += getLogScalingFactor( tlk, k);
                pOutPartials[1] += getLogScalingFactor( tlk, k+1);
            }
            pOutPartials+=2;
        }
        
        if( tlk->pattern_count & 1 ){
            in0 = _mm256_load_pd(pInPartials);
            
            in0 = _mm256_mul_pd(in0, f0);
            
            _mm256_store_pd(temp, in0);
            
            *pOutPartials = log(temp[0]+temp[1]+temp[2]+temp[3]);
            
            if ( tlk->scale ) {
                *pOutPartials += getLogScalingFactor( tlk, tlk->pattern_count-1);
            }
        }
    }
    
    // frequencies is not aligned so I am using set_pd instead of load_pd
    void node_log_likelihoods_4_AVX( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods ){
        
        const double *pInPartials = partials;
        double *pOutPartials = outLogLikelihoods;
        
        double temp[4] __attribute__ ((aligned (32)));
        //memcpy(temp,frequencies,4*sizeof(double));
        //__m256d f0 = _mm256_load_pd(temp);
        __m256d f0 = _mm256_set_pd(frequencies[3],frequencies[2],frequencies[1],frequencies[0]);
        __m256d in0;
        
        for ( int k = 0; k < tlk->pattern_count; k++ ) {
            
            in0 = _mm256_load_pd(pInPartials);
            pInPartials += 4;
            
            in0 = _mm256_mul_pd(in0, f0);
            
            _mm256_store_pd(temp, in0);
            
            *pOutPartials = log(temp[0]+temp[1]+temp[2]+temp[3]);
            
            if ( tlk->scale ) {
                *pOutPartials += getLogScalingFactor( tlk, k);
            }
            pOutPartials++;
        }
    }
    
    void partials_states_and_states_4_AVX( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
        int k;
        int state1, state2;
        int u = 0;
        int w;
        
        double *pPartials = partials;
        const double *m1 = NULL;
        const double *m2 = NULL;
        
        __m256d m1v0, m2v0;
        
        for ( int l = 0; l < tlk->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                
                state1 = tlk->sp->patterns[idx1][k];
                state2 = tlk->sp->patterns[idx2][k];
                
                w = u;
                
                if (state1 < tlk->m->nstate && state2 < tlk->m->nstate) {
                    
                    m1v0 = _mm256_load_pd(&matrices1[w+4*state1]);
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+4*state2]);
                    
                    _mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
                    pPartials += 4;
                    
                }
                else if (state1 < tlk->m->nstate) {
                    // child 1 has a gap or unknown state so treat it as unknown
                    m1 = &matrices1[w+4*state1];
                    
                    *pPartials++ = *m1++;
                    *pPartials++ = *m1++;
                    *pPartials++ = *m1++;
                    *pPartials++ = *m1;
                }
                else if (state2 < tlk->m->nstate ) {
                    // child 2 has a gap or unknown state so treat it as unknown
                    m2 = &matrices2[w+4*state2];
                    
                    *pPartials++ = *m2++;
                    *pPartials++ = *m2++;
                    *pPartials++ = *m2++;
                    *pPartials++ = *m2;
                    
                }
                else {
                    // both children have a gap or unknown state so set partials to 1
                    *pPartials++ = 1.0;
                    *pPartials++ = 1.0;
                    *pPartials++ = 1.0;
                    *pPartials++ = 1.0;
                }
            }
            u += tlk->matrix_size;
        }
    }
    
    void partials_states_and_undefined_4_AVX( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
        
        int v = 0;
        int k;
        int w = 0;
        int state1;
        
        double *pPartials = partials3;
        //const double *m1 = matrices1;
        
        __m256d p2v0, m2v0;
        __m256d tt;
        __m128d sum, m1v0;
        
        //double temp[4] __attribute__ ((aligned (32)));
        
        for ( int l = 0; l < tlk->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                
                state1 = tlk->sp->patterns[idx1][k];
                
                
                p2v0 = _mm256_load_pd(&partials2[v]);
                
                if ( state1 < 4) {
                    
                    //m1 = &matrices1[w + state1*4];
                    m1v0 = _mm_load_pd(&matrices1[w + state1*4]);
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w]);
                    
                    tt = _mm256_mul_pd(m2v0, p2v0);
                    //m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = *m1 * (temp[0]+temp[1]+temp[2]+temp[3]);
                    //m1++;
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+4]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = *m1 * (temp[0]+temp[1]+temp[2]+temp[3]);
                    //m1++;
                    
                    tt = _mm256_hadd_pd(tt, m2v0);
                    
                    sum = _mm_add_pd(_mm256_extractf128_pd(tt, 0), _mm256_extractf128_pd(tt, 1));
                    sum = _mm_mul_pd(m1v0,sum);
                    _mm_store_pd(pPartials, sum);
                    pPartials += 2;
                    
                    m1v0 = _mm_load_pd(&matrices1[w + state1*4+2]);
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+8]);
                    
                    tt = _mm256_mul_pd(m2v0, p2v0);
                    //m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = *m1 *  (temp[0]+temp[1]+temp[2]+temp[3]);
                    //m1++;
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+12]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = *m1 * (temp[0]+temp[1]+temp[2]+temp[3]);
                    tt = _mm256_hadd_pd(tt, m2v0);
                    
                    sum = _mm_add_pd(_mm256_extractf128_pd(tt, 0), _mm256_extractf128_pd(tt, 1));
                    sum = _mm_mul_pd(m1v0,sum);
                    _mm_store_pd(pPartials, sum);
                    pPartials += 2;
                }
                else {
                    // Child 1 has a gap or unknown state so don't use it
                    
                    m2v0 = _mm256_load_pd(&matrices2[w]);
                    
                    tt = _mm256_mul_pd(m2v0, p2v0);
                    //m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+4]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                    tt = _mm256_hadd_pd(tt, m2v0);
                    
                    sum = _mm_add_pd(_mm256_extractf128_pd(tt, 0), _mm256_extractf128_pd(tt, 1));
                    _mm_store_pd(pPartials, sum);
                    pPartials += 2;
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+8]);
                    
                    tt = _mm256_mul_pd(m2v0, p2v0);
                    //m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+12]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    //_mm256_store_pd(temp, m2v0);
                    
                    //*pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                    tt = _mm256_hadd_pd(tt, m2v0);
                    
                    sum = _mm_add_pd(_mm256_extractf128_pd(tt, 0), _mm256_extractf128_pd(tt, 1));
                    _mm_store_pd(pPartials, sum);
                    pPartials += 2;
                    
                }
                v += 4;
            }
            w += tlk->matrix_size;
        }
    }
    
    void partials_states_and_undefined_4_AVX2( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
        
        int v = 0;
        int k;
        int w = 0;
        int state1;
        
        double *pPartials = partials3;
        const double *m1 = matrices1;
        
        __m256d p2v0, m2v0;
        double temp[4] __attribute__ ((aligned (32)));
        
        for ( int l = 0; l < tlk->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                
                state1 = tlk->sp->patterns[idx1][k];
                
                
                p2v0 = _mm256_load_pd(&partials2[v]);
                
                if ( state1 < 4) {
                    
                    m1 = &matrices1[w + state1*4];
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = *m1 * (temp[0]+temp[1]+temp[2]+temp[3]);
                    m1++;
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+4]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = *m1 * (temp[0]+temp[1]+temp[2]+temp[3]);
                    m1++;
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+8]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = *m1 *  (temp[0]+temp[1]+temp[2]+temp[3]);
                    m1++;
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+12]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = *m1 * (temp[0]+temp[1]+temp[2]+temp[3]);
                    
                }
                else {
                    // Child 1 has a gap or unknown state so don't use it
                    
                    m2v0 = _mm256_load_pd(&matrices2[w]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+4]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+8]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                    
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+12]);
                    
                    m2v0 = _mm256_mul_pd(m2v0, p2v0);
                    
                    _mm256_store_pd(temp, m2v0);
                    
                    *pPartials++ = temp[0]+temp[1]+temp[2]+temp[3];
                    
                }
                v += 4;
            }
            w += tlk->matrix_size;
        }
    }
    
    /*
     __m256d x1, x2;
     // calculate 4 two-element horizontal sums:
     // lower 64 bits contain x1[0] + x1[1]
     // next 64 bits contain x2[0] + x2[1]
     // next 64 bits contain x1[2] + x1[3]
     // next 64 bits contain x2[2] + x2[3]
     __m256d sum = _mm256_hadd_pd(x1, x2);
     // extract upper 128 bits of result
     __m128d sum_high = _mm256_extractf128_pd(sum1, 1);
     // add upper 128 bits of sum to its lower 128 bits
     __m128d result = _mm_add_pd(sum_high, (__m128d) sum);
     // lower 64 bits of result contain the sum of x1[0], x1[1], x1[2], x1[3]
     // upper 64 bits of result contain the sum of x2[0], x2[1], x2[2], x2[3]
     */
    void partials_undefined_and_undefined_4_AVX2( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
        
        int v = 0;
        int k;
        int w = 0;
        
        double *pPartials = partials3;
        
        __m256d m1v0, m2v0, p1, p2;
        __m256d m1v2, m2v2;
        //__m128d sum1,sum2,sum;
        double temp[4] __attribute__ ((aligned (32)));
        for ( int l = 0; l < tlk->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                
                p1 = _mm256_load_pd(&partials1[v]);
                
                p2 = _mm256_load_pd(&partials2[v]);
                
                
                m1v0 = _mm256_load_pd(&matrices1[w]);
                
                m2v0 = _mm256_load_pd(&matrices2[w]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1);
                
                m2v0 = _mm256_mul_pd(m2v0, p2);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                //sum1 = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));
                
                
                
                m1v2 = _mm256_load_pd(&matrices1[w+4]);
                
                m2v2 = _mm256_load_pd(&matrices2[w+4]);
                
                
                m1v2 = _mm256_mul_pd(m1v2, p1);
                
                m2v2 = _mm256_mul_pd(m2v2, p2);
                
                m1v2 = _mm256_hadd_pd(m1v2, m2v2);
                
                
                m2v0 = _mm256_blend_pd(m1v0, m1v2, 0b1100);
                m2v2 = _mm256_permute2f128_pd(m1v0,m1v2, 0x21);
                m2v0 = _mm256_add_pd(m2v2,m2v0);
                _mm256_store_pd(temp, m2v0);
                *pPartials++ = temp[0] * temp[1];
                *pPartials++ = temp[2] * temp[3];
                
                
                m1v0 = _mm256_load_pd(&matrices1[w+8]);
                
                m2v0 = _mm256_load_pd(&matrices2[w+8]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1);
                
                m2v0 = _mm256_mul_pd(m2v0, p2);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                //sum1 = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));
                
                
                
                m1v2 = _mm256_load_pd(&matrices1[w+12]);
                
                m2v2 = _mm256_load_pd(&matrices2[w+12]);
                
                
                m1v2 = _mm256_mul_pd(m1v2, p1);
                
                m2v2 = _mm256_mul_pd(m2v2, p2);
                
                m1v2 = _mm256_hadd_pd(m1v2, m2v2);
                
                //sum2 = _mm_add_pd(_mm256_extractf128_pd(m1v2, 0), _mm256_extractf128_pd(m1v2, 1));;
                
                
                m2v0 = _mm256_blend_pd(m1v0, m1v2, 0b1100);
                m2v2 = _mm256_permute2f128_pd(m1v0,m1v2, 0x21);
                m2v0 = _mm256_add_pd(m2v2,m2v0);
                _mm256_store_pd(temp, m2v0);
                *pPartials++ = temp[0] * temp[1];
                *pPartials++ = temp[2] * temp[3];
                
                v += 4;
            }
            w += tlk->matrix_size;
        }
    }
    void partials_undefined_and_undefined_4_AVX( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
        
        int v = 0;
        int k;
        int w = 0;
        
        double *pPartials = partials3;
        
        __m256d m1v0, m2v0, p1, p2;
        __m256d m1v2, m2v2;
        __m128d sum1,sum2,sum;
        
        for ( int l = 0; l < tlk->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                
                p1 = _mm256_load_pd(&partials1[v]);
                
                p2 = _mm256_load_pd(&partials2[v]);
                
                
                m1v0 = _mm256_load_pd(&matrices1[w]);
                
                m2v0 = _mm256_load_pd(&matrices2[w]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1);
                
                m2v0 = _mm256_mul_pd(m2v0, p2);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                sum1 = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));
                
                //_mm_store_pd(temp, sum);
                
                //*pPartials++ = temp[0]*temp[1];
                
                
                m1v2 = _mm256_load_pd(&matrices1[w+4]);
                
                m2v2 = _mm256_load_pd(&matrices2[w+4]);
                
                
                m1v2 = _mm256_mul_pd(m1v2, p1);
                
                m2v2 = _mm256_mul_pd(m2v2, p2);
                
                m1v2 = _mm256_hadd_pd(m1v2, m2v2);
                
                sum2 = _mm_add_pd(_mm256_extractf128_pd(m1v2, 0), _mm256_extractf128_pd(m1v2, 1));
                
                //_mm_store_pd(temp, sum);
                
                //*pPartials++ = temp[0]*temp[1];
                
                sum  = _mm_shuffle_pd(sum1,sum2, _MM_SHUFFLE2(1,0));
                sum1 = _mm_shuffle_pd(sum1,sum2, _MM_SHUFFLE2(0,1));
                _mm_store_pd(pPartials, _mm_mul_pd(sum,sum1));
                pPartials += 2;
                
                
                m1v0 = _mm256_load_pd(&matrices1[w+8]);
                
                m2v0 = _mm256_load_pd(&matrices2[w+8]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1);
                
                m2v0 = _mm256_mul_pd(m2v0, p2);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                sum1 = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));
                
                //_mm_store_pd(temp, sum);
                
                //*pPartials++ = temp[0]*temp[1];
                
                
                m1v2 = _mm256_load_pd(&matrices1[w+12]);
                
                m2v2 = _mm256_load_pd(&matrices2[w+12]);
                
                
                m1v2 = _mm256_mul_pd(m1v2, p1);
                
                m2v2 = _mm256_mul_pd(m2v2, p2);
                
                m1v2 = _mm256_hadd_pd(m1v2, m2v2);
                
                sum2 = _mm_add_pd(_mm256_extractf128_pd(m1v2, 0), _mm256_extractf128_pd(m1v2, 1));;
                
                //_mm_store_pd(temp, sum);
                
                //*pPartials++ = temp[0]*temp[1];
                
                sum  = _mm_shuffle_pd(sum1,sum2, _MM_SHUFFLE2(1,0));
                sum1 = _mm_shuffle_pd(sum1,sum2, _MM_SHUFFLE2(0,1));
                _mm_store_pd(pPartials, _mm_mul_pd(sum,sum1));
                pPartials += 2;
                
                v += 4;
            }
            w += tlk->matrix_size;
        }
    }
    
    void partials_undefined_and_undefined_4_AVX3( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
        
        int v = 0;
        int k;
        int w = 0;
        
        double *pPartials = partials3;
        
        __m256d m1v0, m2v0, p1v0, p2v0;
        __m128d sum;
        double temp[2] __attribute__ ((aligned (16)));
        
        for ( int l = 0; l < tlk->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->pattern_count; k++ ) {
                
                p1v0 = _mm256_load_pd(&partials1[v]);
                
                p2v0 = _mm256_load_pd(&partials2[v]);
                
                
                m1v0 = _mm256_load_pd(&matrices1[w]);
                
                m2v0 = _mm256_load_pd(&matrices2[w]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1v0);
                
                m2v0 = _mm256_mul_pd(m2v0, p2v0);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                sum = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));
                
                _mm_store_pd(temp, sum);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                m1v0 = _mm256_load_pd(&matrices1[w+4]);
                
                m2v0 = _mm256_load_pd(&matrices2[w+4]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1v0);
                
                m2v0 = _mm256_mul_pd(m2v0, p2v0);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                sum = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));
                
                _mm_store_pd(temp, sum);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                m1v0 = _mm256_load_pd(&matrices1[w+8]);
                
                m2v0 = _mm256_load_pd(&matrices2[w+8]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1v0);
                
                m2v0 = _mm256_mul_pd(m2v0, p2v0);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                sum = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));
                
                _mm_store_pd(temp, sum);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                m1v0 = _mm256_load_pd(&matrices1[w+12]);
                
                m2v0 = _mm256_load_pd(&matrices2[w+12]);
                
                
                m1v0 = _mm256_mul_pd(m1v0, p1v0);
                
                m2v0 = _mm256_mul_pd(m2v0, p2v0);
                
                m1v0 = _mm256_hadd_pd(m1v0, m2v0);
                
                sum = _mm_add_pd(_mm256_extractf128_pd(m1v0, 0), _mm256_extractf128_pd(m1v0, 1));;
                
                _mm_store_pd(temp, sum);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                v += 4;
            }
            w += tlk->matrix_size;
        }
    }

void update_partials_4_AVX( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
	
	if( tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1] != NULL ){
		if(  tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2] != NULL ){
			partials_undefined_and_undefined_4_AVX(tlk,
												   tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		else {
			partials_states_and_undefined_4_AVX(tlk,
												tlk->mapping[partialsIndex2],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
												tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2] != NULL ){
			partials_states_and_undefined_4_AVX(tlk,
												tlk->mapping[partialsIndex1],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
												tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
												tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
			
		}
		else{
			partials_states_and_states_4_AVX(tlk,
											 tlk->mapping[partialsIndex1],
											 tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
											 tlk->mapping[partialsIndex2],
											 tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
											 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
	}
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex, partialsIndex1, partialsIndex2);
	}
}

#endif


#pragma mark -
#pragma mark Upper likelihood Reversible

static void _calculate_branch_likelihood_undefined2(SingleTreeLikelihood *tlk, double* rootPartials, int idx, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			const double* partialsChildPtr = partials+v;
			int state = tlk->sp->patterns[idx][k];
			
			if(state < 4){
				const double* transMatrixPtr = &matrices[l*16+state*tlk->m->nstate];
				double sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[3];
				rootPartials[u+state] += sum * weight;
				u+=4;
			}
			else{
				const double* transMatrixPtr = &matrices[l*16];
				double sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
                rootPartials[u] += sum * weight;u++;
				
				sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
				rootPartials[u] += sum * weight;u++;
				
				sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
				rootPartials[u] += sum * weight;u++;
				
				sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
				sum += *transMatrixPtr * partialsChildPtr[3];
				rootPartials[u] += sum * weight;u++;
			}
			v+=4;
		}
	}
}
	
static void _calculate_branch_likelihood_undefined(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			const double* partialsChildPtr = partials+v;
			const double* transMatrixPtr = &matrices[l*16];
			
			double sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum  = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;

			sum  = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum  = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3];
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
		}
	}
}
	
static void _calculate_branch_likelihood_state(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, int partialsIndex, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
            int state = tlk->sp->patterns[partialsIndex][k];
            if(state < 4){
                int w =  l * 16 + state;

                rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 4;u++;v++;
                rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 4;u++;v++;
                rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 4;u++;v++;
                rootPartials[u] += matrices[w] * upperPartials[v] * weight;u++;v++;
            }
            else{
                const double* transMatrixPtr = &matrices[l*16];
                
                rootPartials[u] += (transMatrixPtr[0] + transMatrixPtr[1] + transMatrixPtr[2] + transMatrixPtr[3]) * upperPartials[v] * weight;u++;v++;
                rootPartials[u] += (transMatrixPtr[4] + transMatrixPtr[5] + transMatrixPtr[6] + transMatrixPtr[7]) * upperPartials[v] * weight;u++;v++;
                rootPartials[u] += (transMatrixPtr[8] + transMatrixPtr[9] + transMatrixPtr[10] + transMatrixPtr[11]) * upperPartials[v] * weight;u++;v++;
                rootPartials[u] += (transMatrixPtr[12] + transMatrixPtr[13] + transMatrixPtr[14] + transMatrixPtr[15]) * upperPartials[v] * weight;u++;v++;
            }
		}
	}
	
}

    //TODO: implement this one for sse and the other treelikleihood files
static void _calculate_branch_likelihood_state2(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, const double* partials, const double* matrices){
    memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count);
    int v = 0;
    for(int l = 0; l < tlk->cat_count; l++) {
        int u = 0;
        const double weight = tlk->sm->get_proportion(tlk->sm, l);
        for(int k = 0; k < tlk->pattern_count; k++) {
            const double* partialsChildPtr = partials+v;
            int state = tlk->sp->patterns[upperPartialsIndex][k];
            const double* transMatrixPtr = NULL;
            double sum;
            
            if(state < 4){
                transMatrixPtr = &matrices[l * 16 + state*4];
                
                sum = transMatrixPtr[0] * partialsChildPtr[0];
                sum += transMatrixPtr[1] * partialsChildPtr[1];
                sum += transMatrixPtr[2] * partialsChildPtr[2];
                sum += transMatrixPtr[3] * partialsChildPtr[3];
                rootPartials[u+state] += sum * weight;
                u+=4;
            }
            else{
                transMatrixPtr = &matrices[l*16];
                
                sum = transMatrixPtr[0] * partialsChildPtr[0];
                sum += transMatrixPtr[1] * partialsChildPtr[1];
                sum += transMatrixPtr[2] * partialsChildPtr[2];
                sum += transMatrixPtr[3] * partialsChildPtr[3];
                rootPartials[u] += sum * weight;u++;
                
                sum  = transMatrixPtr[4] * partialsChildPtr[0];
                sum += transMatrixPtr[5] * partialsChildPtr[1];
                sum += transMatrixPtr[6] * partialsChildPtr[2];
                sum += transMatrixPtr[7] * partialsChildPtr[3];
                rootPartials[u] += sum * weight;u++;
                
                sum  = transMatrixPtr[8] * partialsChildPtr[0];
                sum += transMatrixPtr[9] * partialsChildPtr[1];
                sum += transMatrixPtr[10] * partialsChildPtr[2];
                sum += transMatrixPtr[11] * partialsChildPtr[3];
                rootPartials[u] += sum * weight;u++;
                
                sum  = transMatrixPtr[12] * partialsChildPtr[0];
                sum += transMatrixPtr[13] * partialsChildPtr[1];
                sum += transMatrixPtr[14] * partialsChildPtr[2];
                sum += transMatrixPtr[15] * partialsChildPtr[3];
                rootPartials[u] += sum * weight;u++;
            }
            v+=4;
        }
    }
    
}

void calculate_branch_likelihood_4(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex){
    // Sibling is a leaf
	if( tlk->partials[0][partialsIndex] == NULL ){
		_calculate_branch_likelihood_state(tlk, rootPartials,
										   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
										   tlk->mapping[partialsIndex],
										   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
    }
    // upper is coming from the other side of the tree (right node) and right node is a leaf
    // right node should not be a leaf
//    else if(upperPartialsIndex < tlk->sp->size){
//        exit(1);
//    }
    // upper is coming from the other side of the tree (right node) and right node is an internal node
    /*else if(upperPartialsIndex < Tree_node_count(tlk->tree)){
        _calculate_branch_likelihood_undefined(tlk, rootPartials, tlk->partials[upperPartialsIndex], tlk->partials[partialsIndex], tlk->matrices[matrixIndex]);
    }*/
    else if( tlk->partials[0][upperPartialsIndex] == NULL ){
		_calculate_branch_likelihood_state2(tlk, rootPartials,
											tlk->mapping[upperPartialsIndex],
											tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	else{
		_calculate_branch_likelihood_undefined(tlk, rootPartials,
											   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
											   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
}

static void _calculate_branch_partials_undefined(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count*tlk->cat_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		for(int k = 0; k < tlk->pattern_count; k++) {
			const double* partialsChildPtr = partials+v;
			const double* transMatrixPtr = &matrices[tlk->sm->get_site_category(tlk->sm, k) +l*16];
			
			double sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			rootPartials[v] = sum * upperPartials[v]; v++;
			
			sum  = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			rootPartials[v] = sum * upperPartials[v]; v++;
			
			sum  = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			rootPartials[v] = sum * upperPartials[v]; v++;
			
			sum  = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3];
			rootPartials[v] = sum * upperPartials[v]; v++;
		}
	}
}

static void _calculate_branch_partials_state(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, int partialsIndex, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		for(int k = 0; k < tlk->pattern_count; k++) {
			int state = tlk->sp->patterns[partialsIndex][k];
			if(state < 4){
				int w =  tlk->sm->get_site_category(tlk->sm, k) + l * 16 + state;
				
				rootPartials[v] = matrices[w] * upperPartials[v]; w += 4; v++;
				rootPartials[v] = matrices[w] * upperPartials[v]; w += 4; v++;
				rootPartials[v] = matrices[w] * upperPartials[v]; w += 4; v++;
				rootPartials[v] = matrices[w] * upperPartials[v]; v++;
			}
			else{
				const double* transMatrixPtr = &matrices[tlk->sm->get_site_category(tlk->sm, k) + l*16];
				
				rootPartials[v] = (transMatrixPtr[0] + transMatrixPtr[1] + transMatrixPtr[2] + transMatrixPtr[3]) * upperPartials[v]; v++;
				rootPartials[v] = (transMatrixPtr[4] + transMatrixPtr[5] + transMatrixPtr[6] + transMatrixPtr[7]) * upperPartials[v]; v++;
				rootPartials[v] = (transMatrixPtr[8] + transMatrixPtr[9] + transMatrixPtr[10] + transMatrixPtr[11]) * upperPartials[v]; v++;
				rootPartials[v] = (transMatrixPtr[12] + transMatrixPtr[13] + transMatrixPtr[14] + transMatrixPtr[15]) * upperPartials[v]; v++;
			}
		}
	}
	
}

//TODO: implement this one for sse and the other treelikleihood files
static void _calculate_branch_partials_state2(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->m->nstate*tlk->pattern_count*tlk->cat_count);
	int v = 0;
	int u = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		for(int k = 0; k < tlk->pattern_count; k++) {
			const double* partialsChildPtr = partials+v;
			int state = tlk->sp->patterns[upperPartialsIndex][k];
			const double* transMatrixPtr = NULL;
			double sum;
			
			if(state < 4){
				transMatrixPtr = &matrices[tlk->sm->get_site_category(tlk->sm, k) + l * 16 + state*4];
				
				sum = transMatrixPtr[0] * partialsChildPtr[0];
				sum += transMatrixPtr[1] * partialsChildPtr[1];
				sum += transMatrixPtr[2] * partialsChildPtr[2];
				sum += transMatrixPtr[3] * partialsChildPtr[3];
				rootPartials[u+state] = sum;
				u+=4;
			}
			else{
				transMatrixPtr = &matrices[tlk->sm->get_site_category(tlk->sm, k) + l*16];
				
				sum = transMatrixPtr[0] * partialsChildPtr[0];
				sum += transMatrixPtr[1] * partialsChildPtr[1];
				sum += transMatrixPtr[2] * partialsChildPtr[2];
				sum += transMatrixPtr[3] * partialsChildPtr[3];
				rootPartials[u++] = sum;
				
				sum  = transMatrixPtr[4] * partialsChildPtr[0];
				sum += transMatrixPtr[5] * partialsChildPtr[1];
				sum += transMatrixPtr[6] * partialsChildPtr[2];
				sum += transMatrixPtr[7] * partialsChildPtr[3];
				rootPartials[u++] = sum;
				
				sum  = transMatrixPtr[8] * partialsChildPtr[0];
				sum += transMatrixPtr[9] * partialsChildPtr[1];
				sum += transMatrixPtr[10] * partialsChildPtr[2];
				sum += transMatrixPtr[11] * partialsChildPtr[3];
				rootPartials[u++] = sum;
				
				sum  = transMatrixPtr[12] * partialsChildPtr[0];
				sum += transMatrixPtr[13] * partialsChildPtr[1];
				sum += transMatrixPtr[14] * partialsChildPtr[2];
				sum += transMatrixPtr[15] * partialsChildPtr[3];
				rootPartials[u++] = sum;
			}
			v+=4;
		}
	}
	
}

	void calculate_branch_partials_4(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex){
	// Sibling is a leaf
	if( tlk->partials[0][partialsIndex] == NULL ){
		_calculate_branch_partials_state(tlk, rootPartials,
										   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
										   tlk->mapping[partialsIndex],
										   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	else if( tlk->partials[0][upperPartialsIndex] == NULL ){
		_calculate_branch_partials_state2(tlk, rootPartials,
											tlk->mapping[upperPartialsIndex],
											tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
											tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	else{
		_calculate_branch_partials_undefined(tlk, rootPartials,
											   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
											   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
}


#pragma mark -
#pragma mark Upper likelihood Non reversible
	
// Called by a node whose parent is NOT the root and the node's sibling is a leaf
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _update_upper_partials_state( SingleTreeLikelihood *tlk, const double * restrict amatrix_upper, const double * restrict apartials_upper, const double * restrict amatrix_lower, int sibling_index, double *partials ){
	const double *matrix_upper   = __builtin_assume_aligned(amatrix_upper, 16);
	const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
	const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
#else
static void _update_upper_partials_state( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, int sibling_index, double *partials ){
#endif
    int w,k;
    //int w2;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[ sibling_index ][k];
            
//			w2 = w = l * 16;
            
//            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
//            
//            
//            if( state < 4){
//                *pPartials *= matrix_lower[w+state];
//            }
//            pPartials++;
//            w += 4;
//            
//            
//            w2 = l * 16 + 1;
//            
//            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
//            
//            
//            if( state < 4){
//                *pPartials *= matrix_lower[w+state];
//            }
//            pPartials++;
//            w += 4;
//            
//            
//            w2 = l * 16 + 2;
//            
//            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
//            
//            
//            if( state < 4){
//                *pPartials *= matrix_lower[w+state];
//            }
//            pPartials++;
//            w += 4;
//            
//            
//            w2 = l * 16 + 3;
//            
//            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
//            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
//            
//            
//            if( state < 4){
//                *pPartials *= matrix_lower[w+state];
//            }
//            pPartials++;
            
            w = l * 16;
            pPartials[0]  = matrix_upper[w]   * partials_upper[v];
            pPartials[1]  = matrix_upper[w+1] * partials_upper[v];
            pPartials[2]  = matrix_upper[w+2] * partials_upper[v];
            pPartials[3]  = matrix_upper[w+3] * partials_upper[v];

            pPartials[0]  += matrix_upper[w+4] * partials_upper[v+1];
            pPartials[1]  += matrix_upper[w+5] * partials_upper[v+1];
            pPartials[2]  += matrix_upper[w+6] * partials_upper[v+1];
            pPartials[3]  += matrix_upper[w+7] * partials_upper[v+1];

            pPartials[0]  += matrix_upper[w+8]  * partials_upper[v+2];
            pPartials[1]  += matrix_upper[w+9]  * partials_upper[v+2];
            pPartials[2]  += matrix_upper[w+10] * partials_upper[v+2];
            pPartials[3]  += matrix_upper[w+11] * partials_upper[v+2];

            pPartials[0]  += matrix_upper[w+12] * partials_upper[v+3];
            pPartials[1]  += matrix_upper[w+13] * partials_upper[v+3];
            pPartials[2]  += matrix_upper[w+14] * partials_upper[v+3];
            pPartials[3]  += matrix_upper[w+15] * partials_upper[v+3];


            if( state < 4){
                pPartials[0] *= matrix_lower[state+w];
                pPartials[1] *= matrix_lower[state+w+4];
                pPartials[2] *= matrix_lower[state+w+8];
                pPartials[3] *= matrix_lower[state+w+12];
            }
            pPartials+=4;
			v += 4;
		}
    }
}

// Called by a node whose parent is NOT the root and the node's sibling is NOT a leaf
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _update_upper_partials_undefined( SingleTreeLikelihood *tlk, const double * restrict amatrix_upper, const double * restrict apartials_upper, const double * restrict amatrix_lower, const double * restrict apartials_lower, double *partials ){
#else
static void _update_upper_partials_undefined( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ){
#endif
    int w,w2,k;
    int v = 0;
    double sum1,sum2;
    double *pPartials = partials;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *matrix_upper   = __builtin_assume_aligned(amatrix_upper, 16);
	const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
    const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
    const double *partials_lower = __builtin_assume_aligned(apartials_lower, 16);
#endif
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			w2 = w = l * 16;
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            sum2  = matrix_lower[w] * partials_lower[v]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 3]; w++;
            
            *pPartials++ = sum1 * sum2 ;
            
            
            w2 = l * 16 + 1;
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            sum2  = matrix_lower[w] * partials_lower[v]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 3]; w++;
            
            
            *pPartials++ = sum1 * sum2 ;
            
            
            w2 = l * 16 + 2;
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            sum2  = matrix_lower[w] * partials_lower[v]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 3]; w++;
            
            *pPartials++ = sum1 * sum2 ;
            
            
            w2 = l * 16 + 3;
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            sum2  = matrix_lower[w] * partials_lower[v]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum2 += matrix_lower[w] * partials_lower[v + 3];
            
            *pPartials++ = sum1 * sum2;
            
			
			v += 4;
		}
    }
}

// Called by a child of the root but not if the child's sibling is NOT leaf
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _update_upper_partials_root_and_undefined( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, const double * restrict frequencies, double *partials_upper ){
#else
static void _update_upper_partials_root_and_undefined( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper ){
#endif
    int w,k;
    int v = 0;
	double *pPartials = partials_upper;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *partials1 = __builtin_assume_aligned(apartials1, 16);
	const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
#endif
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			w = l * 16;
            
            *pPartials  = matrices1[w] * partials1[v]; w++;
            *pPartials += matrices1[w] * partials1[v + 1]; w++;
            *pPartials += matrices1[w] * partials1[v + 2]; w++;
            *pPartials += matrices1[w] * partials1[v + 3]; w++;
            
            *pPartials++ *= frequencies[0];
            
            
            *pPartials  = matrices1[w] * partials1[v]; w++;
            *pPartials += matrices1[w] * partials1[v + 1]; w++;
            *pPartials += matrices1[w] * partials1[v + 2]; w++;
            *pPartials += matrices1[w] * partials1[v + 3]; w++;
            
            *pPartials++ *= frequencies[1];
            
            
            *pPartials  = matrices1[w] * partials1[v]; w++;
            *pPartials += matrices1[w] * partials1[v + 1]; w++;
            *pPartials += matrices1[w] * partials1[v + 2]; w++;
            *pPartials += matrices1[w] * partials1[v + 3]; w++;
            
            *pPartials++ *= frequencies[2];
            
            
            *pPartials  = matrices1[w] * partials1[v]; w++;
            *pPartials += matrices1[w] * partials1[v + 1]; w++;
            *pPartials += matrices1[w] * partials1[v + 2]; w++;
            *pPartials += matrices1[w] * partials1[v + 3]; w++;
            
            *pPartials++ *= frequencies[3];
            
			
			v += 4;
		}
    }
}

    
// Called by a child of the root and the child's sibling is a leaf
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _update_upper_partials_root_and_state( const SingleTreeLikelihood *tlk, const double * restrict amatrices1, int idx1, const double * restrict frequencies, double *partials_upper ){
#else
static void _update_upper_partials_root_and_state( const SingleTreeLikelihood *tlk, const double *matrices1, int idx1, const double *frequencies, double *partials_upper ){
#endif
    int w;
	double *pPartials = partials_upper;
    int state1;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
#endif
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( int k = 0; k < tlk->pattern_count; k++ ) {

            state1 = tlk->sp->patterns[idx1][k];
            
			w = l * 16;
            
            if( state1 < 4 ){
                *pPartials++ = matrices1[w+state1] * frequencies[0]; w+=4;
                *pPartials++ = matrices1[w+state1] * frequencies[1]; w+=4;
                *pPartials++ = matrices1[w+state1] * frequencies[2]; w+=4;
                *pPartials++ = matrices1[w+state1] * frequencies[3];
            }
            else {
                memcpy(pPartials, frequencies, sizeof(double)*4);
                pPartials += 4;

            }
            
		}
    }

}
    
void update_partials_upper_4( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);
	const double* freqs = tlk->get_root_frequencies(tlk);

    if( Node_isroot(parent) ){
        if( Node_isleaf(sibling) ){
            _update_upper_partials_root_and_state(tlk,
												  tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
												  tlk->mapping[Node_id(sibling)],
												  freqs,
												  tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]]);
        }
        else {
            _update_upper_partials_root_and_undefined(tlk,
													  tlk->partials[tlk->current_partials_indexes[Node_id(sibling)]][Node_id(sibling)],
													  tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
													  freqs,
													  tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]]);
        }
    }
    else if( Node_isleaf(sibling) ){
        _update_upper_partials_state(tlk,
									 tlk->matrices[tlk->current_matrices_indexes[Node_id(parent)]][ Node_id(parent)],
									 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(parent)]]][tlk->upper_partial_indexes[Node_id(parent)]],
									 tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
									 tlk->mapping[Node_id(sibling)],
									 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]]);
    }
    else {
        _update_upper_partials_undefined(tlk,
										 tlk->matrices[tlk->current_matrices_indexes[Node_id(parent)]][ Node_id(parent)],
										 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(parent)]]][tlk->upper_partial_indexes[Node_id(parent)]],
										 tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
										 tlk->partials[tlk->current_partials_indexes[Node_id(sibling)]][Node_id(sibling)],
										 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]]);
    }
}

    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _partial_lower_upper( const SingleTreeLikelihood *tlk, const double * restrict apartials_upper, const double * restrict apartials_lower, const double * restrict amatrix_lower, const double * restrict proportions, double *pattern_lk ){
#else
static void _partial_lower_upper( const SingleTreeLikelihood *tlk, const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ){
#endif
    int w,k;
    int v = 0;
    double p,sum;
    
    memset(pattern_lk, 0, tlk->pattern_count*sizeof(double));
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
	const double *partials_lower = __builtin_assume_aligned(apartials_lower, 16);
	const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
#endif
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			w = l * 16;
            
            sum  = matrix_lower[w] * partials_lower[v + 0]; w++;
            sum += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum += matrix_lower[w] * partials_lower[v + 3]; w++;
            
            p = sum * partials_upper[v];
            
            
            sum  = matrix_lower[w] * partials_lower[v + 0]; w++;
            sum += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum += matrix_lower[w] * partials_lower[v + 3]; w++;
            
            p += sum * partials_upper[v + 1];
            
            
            sum  = matrix_lower[w] * partials_lower[v + 0]; w++;
            sum += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum += matrix_lower[w] * partials_lower[v + 3]; w++;
            
            p += sum * partials_upper[v + 2];
            
            
            sum  = matrix_lower[w] * partials_lower[v + 0]; w++;
            sum += matrix_lower[w] * partials_lower[v + 1]; w++;
            sum += matrix_lower[w] * partials_lower[v + 2]; w++;
            sum += matrix_lower[w] * partials_lower[v + 3]; w++;
            
            p += sum * partials_upper[v + 3];
			
            pattern_lk[k] += p * proportions[l];
			v += 4;
		}
    }
}
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _partial_lower_upper_leaf( const SingleTreeLikelihood *tlk, const double * restrict apartials_upper, int idx, const double * restrict amatrix_lower, const double * restrict proportions, double * restrict pattern_lk ){
#else
static void _partial_lower_upper_leaf( const SingleTreeLikelihood *tlk, const double *partials_upper, int idx, const double *matrix_lower, const double *proportions, double *pattern_lk ){
#endif
    int w,k;
    int v = 0;
    double p = 0;
    int state;
    
    memset(pattern_lk, 0, tlk->pattern_count*sizeof(double));
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
	const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
#endif
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			state = tlk->sp->patterns[idx][k];
            
			w = l * 16;
            if( state < 4 ){
                p  = matrix_lower[w+state] * partials_upper[v]; w += 4;
                p += matrix_lower[w+state] * partials_upper[v + 1]; w += 4;
                p += matrix_lower[w+state] * partials_upper[v + 2]; w += 4;
                p += matrix_lower[w+state] * partials_upper[v + 3];
            }
            else {
                p  = partials_upper[v];
                p += partials_upper[v + 1];
                p += partials_upper[v + 2];
                p += partials_upper[v + 3];
            }
            pattern_lk[k] += p * proportions[l];
			v += 4;
		}
    }
}

void node_log_likelihoods_upper_4( const SingleTreeLikelihood *tlk, Node *node ){
	int node_index = Node_id(node);
	if ( !Node_isleaf(node) ) {
		_partial_lower_upper(tlk,
							 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[node_index]]][tlk->upper_partial_indexes[node_index]],
							 tlk->partials[tlk->current_partials_indexes[node_index]][node_index],
							 tlk->matrices[tlk->current_matrices_indexes[node_index]][node_index],
							 tlk->sm->get_proportions(tlk->sm),
							 tlk->pattern_lk+tlk->sp->count );
	}
	else {
		_partial_lower_upper_leaf(tlk,
								  tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[node_index]]][tlk->upper_partial_indexes[node_index]],
								  tlk->mapping[node_index],
								  tlk->matrices[tlk->current_matrices_indexes[node_index]][node_index],
								  tlk->sm->get_proportions(tlk->sm),
								  tlk->pattern_lk+tlk->sp->count );
	}
}

#pragma mark -
#pragma mark Upper Likelihood Non reversible SSE

#ifdef SSE3_ENABLED
	
// Called by a node whose parent is NOT the root and the node's sibling is a leaf
// matrix_lower is transposed
static void _update_upper_partials_state_sse2( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, int sibling_index, double *partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[ sibling_index ][k];
            
			w2 = w = l * 16;
            
            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
            
            
            if( state < 4){
                *pPartials *= matrix_lower[w+4*state];
            }
            pPartials++;
            w++;
            
            
            w2 = l * 16 + 1;
            
            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
            
            
            if( state < 4){
                *pPartials *= matrix_lower[w+4*state];
            }
            pPartials++;
            w++;
            
            
            w2 = l * 16 + 2;
            
            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
            
            
            if( state < 4){
                *pPartials *= matrix_lower[w+4*state];
            }
            pPartials++;
            w++;
            
            
            w2 = l * 16 + 3;
            
            *pPartials  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            *pPartials += matrix_upper[w2] * partials_upper[v + 3];
            
            
            if( state < 4){
                *pPartials *= matrix_lower[w+4*state];
            }
            pPartials++;
            
			
			v += 4;
		}
    }
}

// Called by a node whose parent is NOT the root and the node's sibling is a leaf
// matrix_lower is transposed
static void _update_upper_partials_state_sse( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, int sibling_index, double *partials ){
    int w,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[ sibling_index ][k];
            
            w = l * 16;
            __m128d pu = _mm_set1_pd(partials_upper[v]);
            __m128d *mu = (__m128d*)&matrix_upper[w];
            __m128d p1;
            __m128d p2;
            
            p1 = _mm_mul_pd(*mu, pu); mu++;
            p2 = _mm_mul_pd(*mu, pu); mu++;
            
            pu = _mm_set1_pd(partials_upper[v+1]);
            p1 = _mm_add_pd(p1, _mm_mul_pd(*mu, pu)); mu++;
            p2 = _mm_add_pd(p2, _mm_mul_pd(*mu, pu)); mu++;
            
            pu = _mm_set1_pd(partials_upper[v+2]);
            p1 = _mm_add_pd(p1, _mm_mul_pd(*mu, pu)); mu++;
            p2 = _mm_add_pd(p2, _mm_mul_pd(*mu, pu)); mu++;
            
            pu = _mm_set1_pd(partials_upper[v+3]);
            p1 = _mm_add_pd(p1, _mm_mul_pd(*mu, pu)); mu++;
            p2 = _mm_add_pd(p2, _mm_mul_pd(*mu, pu));
            
            if( state < 4){
                p1 = _mm_mul_pd(p1, _mm_load_pd(&matrix_lower[4*state+w]));
                p2 = _mm_mul_pd(p2, _mm_load_pd(&matrix_lower[4*state+w+2]));
            }
            _mm_store_pd(&pPartials[0], p1);
            _mm_store_pd(&pPartials[2], p2);

            pPartials+=4;
           
            v += 4;
        }
    }
}
    
// Called by a node whose parent is NOT the root and the node's sibling is NOT a leaf
// All matrices are NOT transposed
static void _update_upper_partials_undefined_sse2( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ){
    int w,w2,k;
    int v = 0;
    double sum1;
    double *pPartials = partials;
    
    __m128d *pl;
    __m128d *ml;
    __m128d temp;
    double t[2] __attribute__ ((aligned (16)));
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			w2 = w = l * 16;
            
            ml = (__m128d*)&matrix_lower[w];
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            pl = (__m128d*)&partials_lower[v];
            temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
            _mm_store_pd(t,temp);
            
            *pPartials++ = sum1 * (t[0]+t[1]) ;
            
            
            w2 = l * 16 + 1;
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            pl = (__m128d*)&partials_lower[v];
            temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
            _mm_store_pd(t,temp);
            
            *pPartials++ = sum1 * (t[0]+t[1]) ;
            
            
            w2 = l * 16 + 2;
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            pl = (__m128d*)&partials_lower[v];
            temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
            _mm_store_pd(t,temp);
            
            *pPartials++ = sum1 * (t[0]+t[1]) ;
            
            
            w2 = l * 16 + 3;
            
            sum1  = matrix_upper[w2] * partials_upper[v]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 1]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 2]; w2 += 4;
            sum1 += matrix_upper[w2] * partials_upper[v + 3];
            
            pl = (__m128d*)&partials_lower[v];
            temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl));
            _mm_store_pd(t,temp);
            
            *pPartials++ = sum1 * (t[0]+t[1]) ;
            
			
			v += 4;
		}
    }
}

// Called by a node whose parent is NOT the root and the node's sibling is NOT a leaf
// All matrices are NOT transposed
static void _update_upper_partials_undefined_sse( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ){
    int w,k;
    int v = 0;
    double *pPartials = partials;
    
    __m128d *pl;
    __m128d *ml;
    __m128d temp,temp2;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            w = l * 16;            
            
            __m128d pu = _mm_set1_pd(partials_upper[v]);
            __m128d *mu = (__m128d*)&matrix_upper[w];
            __m128d p1;
            __m128d p2;
            
            p1 = _mm_mul_pd(*mu, pu); mu++;
            p2 = _mm_mul_pd(*mu, pu); mu++;
            
            pu = _mm_set1_pd(partials_upper[v+1]);
            p1 = _mm_add_pd(p1, _mm_mul_pd(*mu, pu)); mu++;
            p2 = _mm_add_pd(p2, _mm_mul_pd(*mu, pu)); mu++;
            
            pu = _mm_set1_pd(partials_upper[v+2]);
            p1 = _mm_add_pd(p1, _mm_mul_pd(*mu, pu)); mu++;
            p2 = _mm_add_pd(p2, _mm_mul_pd(*mu, pu)); mu++;
            
            pu = _mm_set1_pd(partials_upper[v+3]);
            p1 = _mm_add_pd(p1, _mm_mul_pd(*mu, pu)); mu++;
            p2 = _mm_add_pd(p2, _mm_mul_pd(*mu, pu));
            
            
            ml = (__m128d*)&matrix_lower[w];

            pl = (__m128d*)&partials_lower[v];
            temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
            
            pl = (__m128d*)&partials_lower[v];
            temp2 = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp2 = _mm_add_pd(temp2, _mm_mul_pd(*ml,*pl)); ml++;
            
            _mm_store_pd(&pPartials[0], _mm_mul_pd(p1, _mm_add_pd(_mm_unpackhi_pd(temp,temp2),_mm_unpacklo_pd(temp,temp2))));
            
            
            
            pl = (__m128d*)&partials_lower[v];
            temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
            
            pl = (__m128d*)&partials_lower[v];
            temp2 = _mm_mul_pd(*ml,*pl); ml++; pl++;
            temp2 = _mm_add_pd(temp2, _mm_mul_pd(*ml,*pl));
            
            _mm_store_pd(&pPartials[2], _mm_mul_pd(p2, _mm_add_pd(_mm_unpackhi_pd(temp,temp2),_mm_unpacklo_pd(temp,temp2))));
            
            pPartials+=4;
            
            v += 4;
        }
    }
}
    
    
// Called by a child of the root but not if the child's sibling is NOT leaf
// All matrices are NOT transposed
static void _update_upper_partials_root_and_undefined_sse( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper ){
    int w,k;
    int v = 0;
	double *pPartials = partials_upper;
    
    __m128d *m, *p, temp;
    double t[2] __attribute__ ((aligned (16)));
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			w = l * 16;
            
            m = (__m128d*)&matrices1[w];
            
            p = (__m128d*)&partials1[v];
            
            temp = _mm_mul_pd(*p, *m); p++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); m++;
            
            _mm_store_pd(t, temp);
            
            *pPartials++ = frequencies[0] * (t[0]+t[1]);
            
            
            p = (__m128d*)&partials1[v];
            
            temp = _mm_mul_pd(*p, *m); p++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); m++;
            
            _mm_store_pd(t, temp);
            
            *pPartials++ = frequencies[1] * (t[0]+t[1]);
            
            
            p = (__m128d*)&partials1[v];
            
            temp = _mm_mul_pd(*p, *m); p++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); m++;
            
            _mm_store_pd(t, temp);
            
            *pPartials++ = frequencies[2] * (t[0]+t[1]);
            
            
            p = (__m128d*)&partials1[v];
            
            temp = _mm_mul_pd(*p, *m); p++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m));
            
            _mm_store_pd(t, temp);
            
            *pPartials++ = frequencies[3] * (t[0]+t[1]);
			
			v += 4;
		}
    }
}

// Called by a child of the root and the child's sibling is a leaf
// matrices1 is transposed
static void _update_upper_partials_root_and_state_sse( const SingleTreeLikelihood *tlk, const double *matrices1, int idx1, const double *frequencies, double *partials_upper ){
    int w;
	double *pPartials = partials_upper;
    int state1;
    __m128d *m;
    __m128d f0 = _mm_load_pd(&frequencies[0]);
	__m128d f2 = _mm_load_pd(&frequencies[2]);
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( int k = 0; k < tlk->pattern_count; k++ ) {
			
            state1 = tlk->sp->patterns[idx1][k];
            
			w = l * 16;
            
            if( state1 < 4 ){
                m = (__m128d*)&matrices1[w+4*state1];
                
                _mm_store_pd(pPartials, _mm_mul_pd(*m,f0)); pPartials += 2; m++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,f2)); pPartials += 2;
            }
            else {
                memcpy(pPartials, frequencies, sizeof(double)*4);
                pPartials += 4;
            }
            
		}
    }
}

void update_partials_upper_sse_4( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);
	const double* freqs = tlk->get_root_frequencies(tlk);
	
    if( Node_isroot(parent) ){
        // The matrix of the sibling is transposed
        if( Node_isleaf(sibling) ){
            _update_upper_partials_root_and_state_sse(tlk,
													  tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
													  tlk->mapping[Node_id(sibling)],
													  freqs,
													  tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]] );
        }
        else {
            _update_upper_partials_root_and_undefined_sse(tlk,
														  tlk->partials[tlk->current_partials_indexes[Node_id(sibling)]][Node_id(sibling)],
														  tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
														  freqs,
														  tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]] );
        }
    }
    // The matrix of the sibling is transposed
    // The pparent node cannot be leaf
    else if( Node_isleaf(sibling) ){
        _update_upper_partials_state_sse(tlk,
										 tlk->matrices[tlk->current_matrices_indexes[Node_id(parent)]][Node_id(parent)],
										 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(parent)]]][tlk->upper_partial_indexes[Node_id(parent)]],
										 tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
										 tlk->mapping[Node_id(sibling)],
										 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]]);
    }
    else {
        _update_upper_partials_undefined_sse(tlk,
											 tlk->matrices[tlk->current_matrices_indexes[Node_id(parent)]][Node_id(parent)],
											 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(parent)]]][tlk->upper_partial_indexes[Node_id(parent)]],
											 tlk->matrices[tlk->current_matrices_indexes[Node_id(sibling)]][Node_id(sibling)],
											 tlk->partials[tlk->current_partials_indexes[Node_id(sibling)]][Node_id(sibling)],
											 tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[Node_id(node)]]][tlk->upper_partial_indexes[Node_id(node)]]);
    }
}

static void _partial_lower_upper_sse( const SingleTreeLikelihood *tlk, const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int v = 0;
    double p;
    int cat_count = tlk->cat_count;
    int sp_count = tlk->pattern_count;
    
    __m128d *m, *pl, temp;
    
    double t[2] __attribute__ ((aligned (16)));
    
    memset(pattern_lk, 0, sp_count*sizeof(double));
    
	for ( int l = 0; l < cat_count; l++ ) {
		
		for ( k = 0; k < sp_count; k++ ) {
			
			w = l * 16;
            
            m = (__m128d*)&matrix_lower[w];
            
            pl = (__m128d*)&partials_lower[v];
            
            temp = _mm_mul_pd(*pl, *m); pl++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); m++;
            _mm_store_pd(t, temp);
            p = (t[0]+t[1]) * partials_upper[v];
            
            pl = (__m128d*)&partials_lower[v];
            
            temp = _mm_mul_pd(*pl, *m); pl++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); m++;
            _mm_store_pd(t, temp);
            p += (t[0]+t[1]) * partials_upper[v + 1];
            
            pl = (__m128d*)&partials_lower[v];
            
            temp = _mm_mul_pd(*pl, *m); pl++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); m++;
            _mm_store_pd(t, temp);
            p += (t[0]+t[1]) * partials_upper[v + 2];
            
            pl = (__m128d*)&partials_lower[v];
            
            temp = _mm_mul_pd(*pl, *m); pl++; m++;
            temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m));
            _mm_store_pd(t, temp);
            p += (t[0]+t[1]) * partials_upper[v + 3];
			
            pattern_lk[k] += p * proportions[l];
            
			v += 4;
		}
    }
}

// matrix_lower is transposed
static void _partial_lower_upper_leaf_sse( const SingleTreeLikelihood *tlk, const double *partials_upper, int idx, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int state;
    double p;
    
    int cat_count = tlk->cat_count;
    int sp_count = tlk->pattern_count;
    
    __m128d *m;
    __m128d temp;
    __m128d *pu = (__m128d*)partials_upper;
    
    double t[2] __attribute__ ((aligned (16)));
    
    memset(pattern_lk, 0, sp_count*sizeof(double));
    
	for ( int l = 0; l < cat_count; l++ ) {
		p = proportions[l];
        
		for ( k = 0; k < sp_count; k++ ) {
			state = tlk->sp->patterns[idx][k];
            
			w = l * 16;
            
            if( state < 4 ){
                m = (__m128d*)&matrix_lower[w+4*state];
                
                temp = _mm_mul_pd(*m,*pu); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); pu++;
            }
            else {
                temp = *pu; pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
            }
            
            _mm_store_pd(t, temp);
            
            pattern_lk[k] += (t[0]+t[1]) * p;
		}
    }
}

#endif
