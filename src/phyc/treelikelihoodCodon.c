/*
 *  treelikelihoodCodon.c
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


#include "treelikelihoodCodon.h"

#include <stdio.h>
#include <math.h>

#ifdef SSE3_ENABLED
#include <xmmintrin.h>
#include <pmmintrin.h> // SSE3
/*#include <tmmintrin.h>*/
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "matrix.h"
#include "treelikelihood.h"
#include "treelikelihoodX.h"


#define _TRANSPOSE_ 2
#pragma mark -
#pragma mark codon


void node_log_likelihoods_codon( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods ){

	const double *f = frequencies;
	const double *p = partials;
	double *out = outLogLikelihoods;

	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	
	for ( int k = 0; k < tlk->pattern_count; k++, out++ ) {
		f = frequencies;
		
		*out  = *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		*out += *f++ * *p++;
		
			
		// finish the state, number of codons - 60 with a loop
		for ( int i = 0; i < extra; i++ ) {
			*out += *f++ * *p++;
		}						
				
		*out = log(*out);
		
		if ( tlk->scale ) {
			*out += getLogScalingFactor( tlk, k);
		}
	}
	
}

void node_likelihoods_codon( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods ){
    
    const double *f = frequencies;
    const double *p = partials;
    double *out = outLogLikelihoods;
    
    const int nstate = tlk->sm->nstate;
    const int extra = nstate - 60;
    
    for ( int k = 0; k < tlk->pattern_count; k++, out++ ) {
        f = frequencies;
        
        *out  = *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        *out += *f++ * *p++;
        
        
        // finish the state, number of codons - 60 with a loop
        for ( int i = 0; i < extra; i++ ) {
            *out += *f++ * *p++;
        }
        
        if ( tlk->scale ) {
            *out += getLogScalingFactor( tlk, k);
        }
    }
    
}

void integrate_partials_codon( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials ){	
	int k;
	double *pPartials = outPartials;
	const double *pInPartials = inPartials;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	
	for ( k = 0; k < tlk->pattern_count; k++ ) {
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		*pPartials++ = *pInPartials++ * proportions[0];
		
		// finish the state, number of codons - 60 with a loop
		for ( int i = 0; i < extra; i++ ) {
			*pPartials++ = *pInPartials++ * proportions[0];
		}
	}
	
	
	for ( int l = 1; l < tlk->cat_count; l++ ) {
		pPartials = outPartials;
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			
			// finish the state, number of codons - 60 with a loop
			for ( int i = 0; i < extra; i++ ) {
				*pPartials += *pInPartials++ * proportions[l]; pPartials++;
			}
			
		}
	}
	
}

void update_partials_codon( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
	
	if( tlk->partials[partialsIndex1] != NULL ){
		if(  tlk->partials[partialsIndex2] != NULL ){
			partials_undefined_and_undefined_codon(tlk,
											 tlk->partials[partialsIndex1],
											 tlk->matrices[matrixIndex1],
											 tlk->partials[partialsIndex2],
											 tlk->matrices[matrixIndex2],
											 tlk->partials[partialsIndex]);
		}
		else {
			partials_states_and_undefined_codon(tlk,
										  tlk->mapping[partialsIndex2],
										  tlk->matrices[matrixIndex2],
										  tlk->partials[partialsIndex1],
										  tlk->matrices[matrixIndex1],
										  tlk->partials[partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[partialsIndex2] != NULL ){
			partials_states_and_undefined_codon(tlk,
										  tlk->mapping[partialsIndex1],
										  tlk->matrices[matrixIndex1],
										  tlk->partials[partialsIndex2],
										  tlk->matrices[matrixIndex2],
										  tlk->partials[partialsIndex]);
			
		}
		else{
			partials_states_and_states_codon(tlk,
									   tlk->mapping[partialsIndex1],
									   tlk->matrices[matrixIndex1],
									   tlk->mapping[partialsIndex2],
									   tlk->matrices[matrixIndex2],
									   tlk->partials[partialsIndex]);
		}
	}
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex);
	}
}

// Auto-vectorization
void partials_states_and_undefined_codon( const SingleTreeLikelihood *tlk, int idx1, const double * restrict matrices1, const double * restrict partials2, const double * restrict matrices2, double * restrict partials3){
    
	int v = 0;
	int k;
	int w = 0;
	int state1;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
        const double *p2 = __builtin_assume_aligned(partials2, 16);
        
        const double *m1 = __builtin_assume_aligned(matrices1, 16);
        const double *m2 = __builtin_assume_aligned(matrices2, 16);
#else
        const double *p2 = partials2;
        
        const double *m1 = matrices1;
        const double *m2 = matrices2;
#endif
        
	double *pPartials = partials3;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {

		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			m2 = matrices2 + w;
			
			state1 = tlk->sp->patterns[k][idx1];
			
			if ( state1 < nstate ) {
#ifndef _TRANSPOSE_
				m1 = matrices1 + w;				
#else
				m1 = matrices1 + w + nstate*state1;
#endif
				
				for ( int j = 0; j < nstate; j++ ) {
					
					p2 = partials2 + v;
					
					*pPartials  = *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					// finish the state, number of codons - 60 with a loop
					for ( int i = 0; i < extra; i++ ) {
						*pPartials += *m2 * *p2; m2++; p2++;
					}
				
					
#ifndef _TRANSPOSE_				
					*pPartials *= m1[state1];
					m1 += nstate;
#else	
					*pPartials *= m1[j];
#endif
					pPartials++;
				}
				
                
			}
			else {
				// Child 1 has a gap or unknown state so don't use it
				
				for ( int j = 0; j < nstate; j++ ) {
					
					p2 = partials2 + v;
					
					*pPartials  = *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					*pPartials += *m2 * *p2; m2++; p2++;
					
					// finish the state, number of codons - 60 with a loop
					for ( int i = 0; i < extra; i++ ) {
						*pPartials += *m2 * *p2; m2++; p2++;
					}
					pPartials++;
					
				}
                
			}
			v += nstate;
		}
		w += tlk->matrix_size;
	}
	
}


void partials_undefined_and_undefined_codon( const SingleTreeLikelihood *tlk, const double * restrict partials1, const double * restrict matrices1, const double * restrict partials2, const double * restrict matrices2, double * restrict partials3 ){
	double sum1, sum2;
	
	int v = 0;
	int k;
	int w = 0;
	
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *m1 = __builtin_assume_aligned(matrices1, 16);
    const double *m2 = __builtin_assume_aligned(matrices2, 16);
    
    const double *p1 = __builtin_assume_aligned(partials1, 16);
    const double *p2 = __builtin_assume_aligned(partials2, 16);
#else
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	const double *p1 = partials1;
	const double *p2 = partials2;
#endif
    

	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	
	double *pPartials = partials3;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
        
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			m1 = matrices1 + w;
			m2 = matrices2 + w;
			
			for ( int j = 0; j < nstate; j++ ) {
				
				p1 = partials1 + v;
				p2 = partials2 + v;
				
				sum1   = *m1 * *p1; m1++; p1++;
				sum2   = *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				sum1  += *m1 * *p1; m1++; p1++;
				sum2  += *m2 * *p2; m2++; p2++;
				
				// finish the state, number of codons - 60 with a loop
				for ( int i = 0; i < extra; i++ ) {
					sum1 += *m1 * *p1; m1++; p1++;
					sum2 += *m2 * *p2; m2++; p2++;
				}
				
				*pPartials++ = sum1 * sum2;
			}
			
			v += nstate;
		}
		w += tlk->matrix_size;
	}
    
}


void partials_states_and_states_codon( const SingleTreeLikelihood *tlk, int idx1, const double * restrict matrices1, int idx2, const double * restrict matrices2, double * restrict partials ){
	int k,w,i;
	int u = 0;
	int state1, state2;
	
	double *pPartials = partials;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			state2 = tlk->sp->patterns[k][idx2];
			
			w = u;
			
#ifndef _TRANSPOSE_
			if (state1 < nstate && state2 < nstate) {
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
                
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
                
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
                
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
                
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				
				// finish the state, number of codons - 60 with a loop
				for ( i = 0; i < extra; i++ ) {
					*pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
				}
			}
			else if (state1 < nstate ) {
				// child 1 has a gap or unknown state so treat it as unknown
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				*pPartials++ = matrices1[w + state1]; w += nstate;
				
				// finish the state, number of codons - 60 with a loop
				for ( i = 0; i < extra; i++ ) {
					*pPartials++ = matrices1[w + state1]; w += nstate;
				}
				
			}
			else if (state2 < nstate ) {
				// child 2 has a gap or unknown state so treat it as unknown
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				*pPartials++ = matrices2[w + state2]; w += nstate;
				
				// finish the state, number of codons - 60 with a loop
				for ( i = 0; i < extra; i++ ) {
					*pPartials++ = matrices2[w + state2]; w += nstate;
				}
				
			}
#else
			if (state1 < nstate && state2 < nstate) {
				int w1 = w + state1*nstate;
				int w2 = w + state2*nstate;
				
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
                
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
                
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
                
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
                
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				
				// finish the state, number of codons - 60 with a loop
				for ( i = 0; i < extra; i++ ) {
					*pPartials++ = matrices1[w1] * matrices2[w2]; w1++; w2++;
				}
			}
			else if (state1 < nstate ) {
				// child 2 has a gap or unknown state so treat it as unknown
				int w1 = w + state1*nstate;
				memcpy(pPartials, &matrices1[w1], nstate*sizeof(double));
				pPartials+=nstate;
			}
			else if (state2 < nstate ) {
				// child 1 has a gap or unknown state so treat it as unknown
				int w2 = w + state2*nstate;
				memcpy(pPartials, &matrices2[w2], nstate*sizeof(double));
				pPartials+=nstate;
				
			}
#endif
			
			else {
				// both children have a gap or unknown state so set partials to 1
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
				
				// finish the state, number of codons - 60 with a loop
				for ( i = 0; i < extra; i++ ) {
					*pPartials++ = 1.0;
				}
			}
		}
		u += tlk->matrix_size;
	}
    
}

#pragma mark -
#pragma mark OpenMP

#ifdef _OPENMP

void update_partials_codon_openmp( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
	
    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_codon_openmp(tlk,
                                                          tlk->partials[nodeIndex1],
                                                          tlk->matrices[nodeIndex1],
                                                          tlk->partials[nodeIndex2],
                                                          tlk->matrices[nodeIndex2],
                                                          tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_codon_openmp(tlk,
                                                       tlk->mapping[nodeIndex2],
                                                       tlk->matrices[nodeIndex2],
                                                       tlk->partials[nodeIndex1],
                                                       tlk->matrices[nodeIndex1],
                                                       tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_codon_openmp(tlk,
                                                       tlk->mapping[nodeIndex1],
                                                       tlk->matrices[nodeIndex1],
                                                       tlk->partials[nodeIndex2],
                                                       tlk->matrices[nodeIndex2],
                                                       tlk->partials[nodeIndex3]);
            
        }
        else{
            partials_states_and_states_codon_openmp(tlk,
                                                    tlk->mapping[nodeIndex1],
                                                    tlk->matrices[nodeIndex1],
                                                    tlk->mapping[nodeIndex2],
                                                    tlk->matrices[nodeIndex2],
                                                    tlk->partials[nodeIndex3]);
        }
    }
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
	}
}

void partials_states_and_states_codon_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
	
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	int nThreads = tlk->nthreads;
    
    #pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for ( int lk = 0; lk < tlk->pattern_count*tlk->cat_count; lk++ ) {
        int l = lk / tlk->pattern_count;
        int k = lk % tlk->pattern_count;
        
        int state1 = tlk->sp->patterns[k][idx1];
        int state2 = tlk->sp->patterns[k][idx2];
        
        int w = l * tlk->matrix_size;
        
        double *pPartials = partials + (l*tlk->pattern_count + k)*nstate;
        
        if (state1 < nstate && state2 < nstate) {
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            
            // finish the state, number of codons - 60 with a loop
            for ( int i = 0; i < extra; i++ ) {
                *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += nstate;
            }
            
        }
        else if (state1 < nstate ) {
            // child 1 has a gap or unknown state so treat it as unknown
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            *pPartials++ = matrices1[w + state1]; w += nstate;
            
            // finish the state, number of codons - 60 with a loop
            for ( int i = 0; i < extra; i++ ) {
                *pPartials++ = matrices1[w + state1]; w += nstate;
            }
            
        }
        else if (state2 < nstate ) {
            // child 2 has a gap or unknown state so treat it as unknown
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            *pPartials++ = matrices2[w + state2]; w += nstate;
            
            // finish the state, number of codons - 60 with a loop
            for ( int i = 0; i < extra; i++ ) {
                *pPartials++ = matrices2[w + state2]; w += nstate;
            }
            
        }
        else {
            // both children have a gap or unknown state so set partials to 1
            *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
            *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
            *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
            *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
            *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
            *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
            
            // finish the state, number of codons - 60 with a loop
            for ( int i = 0; i < extra; i++ ) {
                *pPartials++ = 1.0;
            }
        }
    }
	
    
}


void partials_states_and_undefined_codon_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
    
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	int nThreads = tlk->nthreads;
    
#pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for ( int lk = 0; lk < tlk->pattern_count*tlk->cat_count; lk++ ) {
        int l = lk / tlk->pattern_count;
        int k = lk % tlk->pattern_count;
        
        int w = l * tlk->matrix_size;
        int v = (l*tlk->pattern_count + k) * nstate;
        const double *m2 = matrices2 + w;
        
        int state1 = tlk->sp->patterns[k][idx1];
        
        const double *p2 = NULL;
        double *pPartials = partials3+v;
        
        if ( state1 < nstate ) {
            
            const double *m1 = matrices1 + w;
            
            for ( int j = 0; j < nstate; j++ ) {
                
                p2 = partials2 + v;
                
                *pPartials  = *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                // finish the state, number of codons - 60 with a loop
                for ( int i = 0; i < extra; i++ ) {
                    *pPartials += *m2 * *p2; m2++; p2++;
                }
                
                *pPartials *= m1[state1];
                pPartials++;
                m1 += nstate;
            }
            
            
        }
        else {
            // Child 1 has a gap or unknown state so don't use it
            
            for ( int j = 0; j < nstate; j++ ) {
                
                p2 = partials2 + v;
                
                *pPartials  = *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                *pPartials += *m2 * *p2; m2++; p2++;
                
                // finish the state, number of codons - 60 with a loop
                for ( int i = 0; i < extra; i++ ) {
                    *pPartials += *m2 * *p2; m2++; p2++;
                }
                pPartials++;
                
            }
            
        }
    }
	
	
}

/*
 With openmp 3.0
 #pragma omp parallel for schedule(dynamic,1) collapse(2) num_threads(2)
 for ( int l = 0; l < tlk->cat_count; l++ ) {
 
 for ( int k = 0; k < tlk->pattern_count; k++ ) {}}
 */

void partials_undefined_and_undefined_codon_openmp( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
    
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	int nThreads = tlk->nthreads;
    
#pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for ( int lk = 0; lk < tlk->pattern_count*tlk->cat_count; lk++ ) {
        int l = lk / tlk->pattern_count;
        int k = lk % tlk->pattern_count;
        
        int w = l * tlk->matrix_size;
        int v = (l*tlk->pattern_count + k) * nstate;
        
        const double *m1 = matrices1 + w;
        const double *m2 = matrices2 + w;
        
        double *pPartials = partials3 + v;
        
        double sum1, sum2;
        const double *p1 = NULL;
        const double *p2 = NULL;
        
        for ( int j = 0; j < nstate; j++ ) {
            
            p1 = partials1 + v;
            p2 = partials2 + v;
            
            sum1   = *m1 * *p1; m1++; p1++;
            sum2   = *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            sum1  += *m1 * *p1; m1++; p1++;
            sum2  += *m2 * *p2; m2++; p2++;
            
            // finish the state, number of codons - 60 with a loop
            for ( int i = 0; i < extra; i++ ) {
                sum1 += *m1 * *p1; m1++; p1++;
                sum2 += *m2 * *p2; m2++; p2++;
            }
            
            *pPartials++ = sum1 * sum2;
        }
        
        
	}
    
}

#endif



#pragma mark -
#pragma mark SSE

#ifdef SSE3_ENABLED

//FIXME: should replace _mm_hadd_pd by _mm_add_pd if possible
void node_log_likelihoods_codon_SSE( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods ){
    
	const double *pInPartials = partials;
	double *pOutPartials = outLogLikelihoods;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 2;
	
	double temp[2] __attribute__ ((aligned (16)));
	
	__m128d *f = (__m128d*)frequencies;
	__m128d in, sum;
	
	
	for ( int k = 0; k < tlk->pattern_count; k++, pOutPartials++ ) {
		
		f = (__m128d*)frequencies;
		
		sum = _mm_setzero_pd();
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		in  = _mm_load_pd(pInPartials);
		in  = _mm_mul_pd(in, *f++);
		sum = _mm_add_pd(sum, in);
		pInPartials += 2;
		
		if ( extra != 0 ) {
			
			// finish the state, number of codons - 60 with a loop
			
			for ( int i = 0; i < extra_pairs; i++ ) {
				in  = _mm_load_pd(pInPartials);
				in  = _mm_mul_pd(in, *f++);
				sum = _mm_add_pd(sum, in);
				pInPartials += 2;
			}
			
			_mm_store_pd(temp, sum);
			
			// odd number so one more
			if ( extra & 1 ) {
				*pOutPartials = log( temp[0] + temp[1] + *pInPartials * frequencies[nstate-1] );
				pInPartials += 2;
			}
			else {
				*pOutPartials = log(temp[0]+temp[1]);
			}
			
		}
		else {
			
			_mm_store_pd(temp, sum);
			
			*pOutPartials = log(temp[0]+temp[1]);
		}
		
		
		if ( tlk->scale ) {
			*pOutPartials += getLogScalingFactor( tlk, k);
		}
        
        if ( extra & 1 ) {
            pOutPartials++;
        }
	}
	
    //	for ( int i = 0; i < tlk->pattern_count; i ++ ) {
    //		nodelog += outLogLikelihoods[i];
    //		//fprintf(stderr, " %f", outLogLikelihoods[i]);
    //	}
    //	fprintf(stderr, "\nlognode %f\n", nodelog);
}


void integrate_partials_codon_SSE( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials ){
	int k;
	double *pPartials = outPartials;
	const double *pInPartials = inPartials;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 2;
	
	double prop[2] __attribute__ ((aligned (16)));
	__m128d in, pn;
	__m128d *pSSE;
	
	prop[0] = prop[1] = proportions[0];
	pn = _mm_load_pd(prop);
	
	for ( k = 0; k < tlk->pattern_count; k++ ) {
		
		//pSSE = (__m128d*)pPartials;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		in = _mm_load_pd(pInPartials);
		_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
		pPartials   += 2;
		pInPartials += 2;
		
		if ( extra != 0 ) {
			
			// finish the state, number of codons - 60 with a loop
			for ( int i = 0; i < extra_pairs; i++ ) {
				in = _mm_load_pd(pInPartials);
				_mm_store_pd(pPartials, _mm_mul_pd(in, pn));
				pPartials   += 2;
				pInPartials += 2;
			}
			
			// odd number so one more
			if ( extra & 1 ) {
				*pPartials = *pInPartials * proportions[0];
                pPartials   += 2;
				pInPartials += 2;
			}
			
		}
		
	}
	
	
	for ( int l = 1; l < tlk->cat_count; l++ ) {
		prop[0] = prop[1] = proportions[l];
		pn = _mm_load_pd(prop);
		pPartials = outPartials;
		//pSSE = (__m128d*)outPartials;
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			pSSE = (__m128d*)pPartials;
			
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
			
			
			if ( extra != 0 ) {
				
				// finish the state, number of codons - 60 with a loop
				
				for ( int i = 0; i < extra_pairs; i++ ) {
					in = _mm_load_pd(pInPartials);
					in = _mm_mul_pd(in, pn);
					_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
					pSSE++;
					pInPartials += 2;
					pPartials   += 2;
				}
				
				// odd number so one more
				if ( extra & 1 ) {
					*pPartials += *pInPartials * proportions[l];
                    pSSE++;
                    pInPartials += 2;
                    pPartials   += 2;
				}
				
			}
			
		}
	}
	
    //	for ( int i = 0; i < tlk->partials_size; i ++ ) {
    //		integrate += outPartials[i];
    //	}
    //	fprintf(stderr, "integrate %f\n", integrate);
}

void update_partials_codon_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
	
	if( tlk->partials[partialsIndex1] != NULL ){
		if(  tlk->partials[partialsIndex2] != NULL ){
			partials_undefined_and_undefined_codon_SSE(tlk,
												   tlk->partials[partialsIndex1],
												   tlk->matrices[matrixIndex1],
												   tlk->partials[partialsIndex2],
												   tlk->matrices[matrixIndex2],
												   tlk->partials[partialsIndex]);
		}
		else {
			partials_states_and_undefined_codon_SSE(tlk,
										  tlk->mapping[partialsIndex2],
										  tlk->matrices[matrixIndex2],
										  tlk->partials[partialsIndex1],
										  tlk->matrices[matrixIndex1],
										  tlk->partials[partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[partialsIndex2] != NULL ){
			partials_states_and_undefined_codon_SSE(tlk,
										  tlk->mapping[partialsIndex1],
										  tlk->matrices[matrixIndex1],
										  tlk->partials[partialsIndex2],
										  tlk->matrices[matrixIndex2],
										  tlk->partials[partialsIndex]);
			
		}
		else{
			partials_states_and_states_codon_SSE(tlk,
									   tlk->mapping[partialsIndex1],
									   tlk->matrices[matrixIndex1],
									   tlk->mapping[partialsIndex2],
									   tlk->matrices[matrixIndex2],
									   tlk->partials[partialsIndex]);
		}
	}
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex);
	}
}

void update_partials_codon_odd_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
	
	if( tlk->partials[partialsIndex1] != NULL ){
		if(  tlk->partials[partialsIndex2] != NULL ){
			partials_undefined_and_undefined_codon_odd_SSE(tlk,
													   tlk->partials[partialsIndex1],
													   tlk->matrices[matrixIndex1],
													   tlk->partials[partialsIndex2],
													   tlk->matrices[matrixIndex2],
													   tlk->partials[partialsIndex]);
		}
		else {
			partials_states_and_undefined_codon_odd_SSE(tlk,
													tlk->mapping[partialsIndex2],
													tlk->matrices[matrixIndex2],
													tlk->partials[partialsIndex1],
													tlk->matrices[matrixIndex1],
													tlk->partials[partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[partialsIndex2] != NULL ){
			partials_states_and_undefined_codon_odd_SSE(tlk,
													tlk->mapping[partialsIndex1],
													tlk->matrices[matrixIndex1],
													tlk->partials[partialsIndex2],
													tlk->matrices[matrixIndex2],
													tlk->partials[partialsIndex]);
			
		}
		else{
			partials_states_and_states_codon_odd_SSE(tlk,
												 tlk->mapping[partialsIndex1],
												 tlk->matrices[matrixIndex1],
												 tlk->mapping[partialsIndex2],
												 tlk->matrices[matrixIndex2],
												 tlk->partials[partialsIndex]);
		}
	}
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex);
	}
}

void partials_states_and_states_codon_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
	int k;
	int state1, state2;
	int u = 0;
	int w;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 2;
	
	double *pPartials = partials;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
    
	__m128d m1v0, m2v0;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			state2 = tlk->sp->patterns[k][idx2];
			
			w = u;
			
			if (state1 < nstate && state2 < nstate) {
				
				m1 = &matrices1[w + nstate*state1];
				m2 = &matrices2[w + nstate*state2];
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				if ( extra != 0 ) {
                    
					// finish the state, number of codons - 60 with a loop
					
					for ( int i = 0; i < extra_pairs; i++ ) {
						m1v0 = _mm_load_pd(m1);
						m2v0 = _mm_load_pd(m2);
						_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
						pPartials += 2;
						m1+=2; m2+=2;
					}
					
					// odd number so one more
					if ( extra & 1 ) {
						*pPartials++ = *m1 * *m2;
					}
					
				}
				
			}
			else if (state1 < nstate) {
				// child 1 has a gap or unknown state so treat it as unknown
				m1 = &matrices1[w + nstate*state1];
				
				memcpy(pPartials, m1, nstate*sizeof(double));
				
				pPartials += nstate;
				
			}
			else if (state2 < nstate ) {
				// child 2 has a gap or unknown state so treat it as unknown
				m2 = &matrices2[w + nstate*state2];
				
				memcpy(pPartials, m2, nstate*sizeof(double));
				
				pPartials += nstate;
				
			}
			else {
				// both children have a gap or unknown state so set partials to 1
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				
				for ( int i = 0; i < extra; i++ ) {
					*pPartials++ = 1.0;
				}
			}
		}
		u += tlk->matrix_size;
	}
	
    //	for ( int i = 0; i < tlk->partials_size; i ++ ) {
    //		tot_state2 += partials[i];
    //	}
    //	fprintf(stderr, "state2 %f\n", tot_state2);
}

void partials_states_and_states_codon_odd_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
	int k;
	int state1, state2;
	int u = 0;
	int w;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 2;
	
	double *pPartials = partials;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
    
	__m128d m1v0, m2v0;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			state2 = tlk->sp->patterns[k][idx2];
			
			w = u;
			
			if (state1 < nstate && state2 < nstate) {
				
				m1 = &matrices1[w + (nstate+1)*state1];
				m2 = &matrices2[w + (nstate+1)*state2];
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_load_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
				
				m1v0 = _mm_load_pd(m1);
				m2v0 = _mm_loadu_pd(m2);
				_mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
				pPartials += 2;
				m1+=2; m2+=2;
		
                    
                // finish the state, number of codons - 60 with a loop
                
                for ( int i = 0; i < extra_pairs; i++ ) {
                    m1v0 = _mm_load_pd(m1);
                    m2v0 = _mm_load_pd(m2);
                    _mm_store_pd(pPartials, _mm_mul_pd(m1v0, m2v0));
                    pPartials += 2;
                    m1+=2; m2+=2;
                }
                
                *pPartials++ = *m1 * *m2;
					
				
			}
			else if (state1 < nstate) {
				// child 1 has a gap or unknown state so treat it as unknown
				m1 = &matrices1[w + (nstate+1)*state1];
				// we don't need to copy the last double (padding)
				memcpy(pPartials, m1, nstate*sizeof(double));
				
				pPartials += nstate;
				
			}
			else if (state2 < nstate ) {
				// child 2 has a gap or unknown state so treat it as unknown
				m2 = &matrices2[w + (nstate+1)*state2];
				
				memcpy(pPartials, m2, nstate*sizeof(double));
				
				pPartials += nstate;
				
			}
			else {
				// both children have a gap or unknown state so set partials to 1
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				
				for ( int i = 0; i < extra; i++ ) {
					*pPartials++ = 1.0;
				}
			}
            
            // there is an odd number of states, skip 1 so it is aligned for _mm_store_pd
            pPartials++;
		}
		u += tlk->matrix_size;
	}
	
}



void partials_states_and_undefined_codon_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
	
	int v = 0;
	int k;
	int w = 0;
	int state1;
	
	double *pPartials = partials3;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	
	__m128d m2v0, sum;
	double temp[2] __attribute__ ((aligned (16)));
	__m128d *p2 = (__m128d *)partials2;
	
	const int nstate = tlk->sm->nstate;
	int extra = nstate - 60;
	int extra_pairs = extra / 2;
    
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
            
			m2 = matrices2 + w;
			
			if ( state1 < nstate ) {
				
				m1 = &matrices1[w + state1*nstate];
				
				for ( int j = 0; j < nstate; j++ ) {
					
					sum = _mm_setzero_pd();
					
					p2 = (__m128d*)&partials2[v];
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
                    
					
					if ( extra != 0 ) {
                        
						// finish the last states, number of codons - 60 with a loop
						
						for ( int i = 0; i < extra_pairs; i++ ) {
							
							m2v0 = _mm_load_pd(m2);
							m2v0 = _mm_mul_pd(m2v0, *p2++);
							sum  = _mm_add_pd(sum, m2v0);
							m2 += 2;
							
						}
						
						_mm_store_pd(temp, sum);
						
						// odd number so one more
						if ( extra & 1 ) {
							*pPartials++ = *m1 * ( temp[0] + temp[1] + (*m2 * partials2[v+nstate-1]) );
							m2++;
						}
						else {
							*pPartials++ = *m1 * ( temp[0] + temp[1] );
						}
                        
					}
					else {
						_mm_store_pd(temp, sum);
						*pPartials++ = *m1 * ( temp[0] + temp[1] );
					}
					
					m1++;
					
				}
				
				
			}
			else {
				// Child 1 has a gap or unknown state so don't use it
				
				for ( int j = 0; j < nstate; j++ ) {
					
					sum = _mm_setzero_pd();
					
					p2 = (__m128d*)&partials2[v];
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
                    
					
					if ( extra != 0 ) {
						
						// finish the last states, number of codons - 60 with a loop
						
						for ( int i = 0; i < extra_pairs; i++ ) {
							
							m2v0 = _mm_load_pd(m2);
							m2v0 = _mm_mul_pd(m2v0, *p2++);
							sum  = _mm_add_pd(sum, m2v0);
							m2 += 2;
							
						}
						
						_mm_store_pd(temp, sum);
						
						// odd number so one more
						if ( extra & 1 ) {
							*pPartials++ = temp[0] + temp[1] + ( *m2 * partials2[v+nstate-1] );
							m2++;
						}
						else {
							*pPartials++ = temp[0] + temp[1];
						}
						
					}
					else {
						_mm_store_pd(temp, sum);
						*pPartials++ = temp[0] + temp[1];
					}
					
					m1++;
					
				}
				
			}
			v += nstate;
		}
		w += tlk->matrix_size;
	}
	
    //	for ( int i = 0; i < tlk->partials_size; i ++ ) {
    //		tot_state_undef += partials3[i];
    //	}
    //	fprintf(stderr, "state undef %f\n", tot_state_undef);
}

void partials_states_and_undefined_codon_odd_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
	
	int v = 0;
	int k;
	int w = 0;
	int state1;
	
	double *pPartials = partials3;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	
	__m128d m2v0, sum;
	double temp[2] __attribute__ ((aligned (16)));
	__m128d *p2 = (__m128d *)partials2;
	
	const int nstate = tlk->sm->nstate;
	int extra = nstate - 60;
	int extra_pairs = extra / 2;
    
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
            
			m2 = matrices2 + w; // matrices2 is padded
			
			if ( state1 < nstate ) {
				
				m1 = &matrices1[w + state1*(nstate+1)]; //matrices1 s padded
				
				for ( int j = 0; j < nstate; j++ ) {
					
					sum = _mm_setzero_pd();
					
					p2 = (__m128d*)&partials2[v]; // p2 always aligned because it is padded
					
                    //0
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
                    //10
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
                    // 20
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
                    // 29
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
                    
                    // finish the last states, number of codons - 60 with a loop
                    
                    for ( int i = 0; i < extra_pairs; i++ ) {
                        m2v0 = _mm_load_pd(m2);
                        m2v0 = _mm_mul_pd(m2v0, *p2++);
                        sum  = _mm_add_pd(sum, m2v0);
                        m2 += 2;
                    }
                    
                    _mm_store_pd(temp, sum);
                    
                    
                    *pPartials++ = *m1 * ( temp[0] + temp[1] + (*m2 * partials2[v+nstate-1]) );
                    m2 += 2; // +2 because of padding
					m1 += 2;
				}
				
				
			}
			else {
				// Child 1 has a gap or unknown state so don't use it
				
				for ( int j = 0; j < nstate; j++ ) {
					
					sum = _mm_setzero_pd();
					
					p2 = (__m128d*)&partials2[v];
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
					m2v0 = _mm_load_pd(m2);
					m2v0 = _mm_mul_pd(m2v0, *p2++);
					sum  = _mm_add_pd(sum, m2v0);
					m2 += 2;
					
                    
					
                    // finish the last states, number of codons - 60 with a loop
                    
                    for ( int i = 0; i < extra_pairs; i++ ) {
                        m2v0 = _mm_load_pd(m2);
                        m2v0 = _mm_mul_pd(m2v0, *p2++);
                        sum  = _mm_add_pd(sum, m2v0);
                        m2 += 2;
                    }
                    
                    _mm_store_pd(temp, sum);
                    
                    *pPartials++ = temp[0] + temp[1] + ( *m2 * partials2[v+nstate-1] );
                    m2 += 2;
					m1 += 2;
					
				}
				
			}
			v += nstate+1;
		}
		w += tlk->matrix_size;
	}
	
}

void partials_undefined_and_undefined_codon_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
	
	int v = 0;
	int k;
	int w = 0;
	
	double *pPartials = partials3;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 2;
	
	__m128d *p1 = (__m128d*)partials1;
	__m128d *p2 = (__m128d*)partials2;
	
	__m128d m1v0;
	
	double temp[2] __attribute__ ((aligned (16)));
	
	__m128d m1s, m2s;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			m1 = matrices1 + w;
			m2 = matrices2 + w;
			
			for ( int j = 0; j < nstate; j++ ) {
				
				m1s = _mm_setzero_pd();
				m2s = _mm_setzero_pd();
				
				p1 = (__m128d*)&partials1[v];
				p2 = (__m128d*)&partials2[v];
                
				
				// matrices 1
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				
				// matrices 2
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				
				
				if ( extra != 0 ) {
					
					// finish the state, number of codons - 60 with a loop
					
					for ( int i = 0; i < extra_pairs; i++ ) {
						
						m1v0 = _mm_load_pd(m1); m1 += 2;
						m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
						m1s = _mm_add_pd(m1s, m1v0);
						
						m1v0 = _mm_load_pd(m2); m2 += 2;
						m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
						m2s = _mm_add_pd(m2s, m1v0);
						
					}
					
					m1s = _mm_hadd_pd(m1s, m2s);
					
					_mm_store_pd(temp, m1s);
					
					// odd number so one more
					if ( extra & 1 ) {
						*pPartials++ = ( temp[0] + *m1 * partials1[v+nstate-1] ) * ( temp[1] + *m2 * partials2[v+nstate-1] );
						m1++;
						m2++;
					}
					else {
						
						*pPartials++ = temp[0] * temp[1];
					}
					
				}
				else {
					
					m1s = _mm_hadd_pd(m1s, m2s);
					
					_mm_store_pd(temp, m1s);
					
					*pPartials++ = temp[0] * temp[1];
				}
				
				
			}
			
			v += nstate;
			
		}
		w += tlk->matrix_size;
	}
	
}


void partials_undefined_and_undefined_codon_odd_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
	
	int v = 0;
	int k;
	int w = 0;
	
	double *pPartials = partials3;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 2;
	
	__m128d *p1 = (__m128d*)partials1;
	__m128d *p2 = (__m128d*)partials2;
	
	__m128d m1v0;
	
	double temp[2] __attribute__ ((aligned (16)));
	
	__m128d m1s, m2s;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			m1 = matrices1 + w;
			m2 = matrices2 + w;
			
			for ( int j = 0; j < nstate; j++ ) {
				
				m1s = _mm_setzero_pd();
				m2s = _mm_setzero_pd();
				
				p1 = (__m128d*)&partials1[v]; // aligned
				p2 = (__m128d*)&partials2[v]; // aligned
                
				
				// matrices 1
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				m1v0 = _mm_load_pd(m1); m1 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
				m1s = _mm_add_pd(m1s, m1v0);
				
				
				// matrices 2
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
				m1v0 = _mm_load_pd(m2); m2 += 2;
				m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
				m2s = _mm_add_pd(m2s, m1v0);
				
					
                // finish the state, number of codons - 60 with a loop
                
                for ( int i = 0; i < extra_pairs; i++ ) {
                    
                    m1v0 = _mm_load_pd(m1); m1 += 2;
                    m1v0 = _mm_mul_pd(m1v0, *p1); p1++;
                    m1s = _mm_add_pd(m1s, m1v0);
                    
                    m1v0 = _mm_load_pd(m2); m2 += 2;
                    m1v0 = _mm_mul_pd(m1v0, *p2); p2++;
                    m2s = _mm_add_pd(m2s, m1v0);
                    
                }
                
                m1s = _mm_hadd_pd(m1s, m2s);
                
                _mm_store_pd(temp, m1s);
                
                // Last one
                *pPartials++ = ( temp[0] + *m1 * partials1[v+nstate-1] ) * ( temp[1] + *m2 * partials2[v+nstate-1] );
                
                m1 += 2;
                m2 += 2;
			}
        
             // padding for alignment
			pPartials++;
			v += nstate+1;
			
		}
		w += tlk->matrix_size;
	}
	
}

#endif

#pragma mark -
#pragma mark AVX

#ifdef AVX_ENABLED
#include <immintrin.h>

void node_log_likelihoods_codon_AVX( const SingleTreeLikelihood *tlk, const double *partials, const double *frequencies, double *outLogLikelihoods ){
    
	const double *pInPartials = partials;
	double *pOutPartials = outLogLikelihoods;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 4;
	
	double temp[4] __attribute__ ((aligned (16)));
	
	__m256d *f = (__m256d*)frequencies;
	__m256d in, sum;
	
	
	for ( int k = 0; k < tlk->pattern_count; k++, pOutPartials++ ) {
		
		f = (__m256d*)frequencies;
		
		sum = _mm256_setzero_pd();
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		in  = _mm256_load_pd(pInPartials);
		in  = _mm256_mul_pd(in, *f++);
		sum = _mm256_add_pd(sum, in);
		pInPartials += 4;
		
		
		if ( extra != 0 ) {
			
			// finish the state, number of codons - 60 with a loop
			
			for ( int i = 0; i < extra_pairs; i++ ) {
				in  = _mm256_load_pd(pInPartials);
				in  = _mm256_mul_pd(in, *f++);
				sum = _mm256_add_pd(sum, in);
				pInPartials += 4;
			}
			
			_mm_store_pd(temp, sum);
			
			// odd number so one more
			if ( extra & 1 ) {
				*pOutPartials = log( temp[0] + temp[1] + temp[2] + temp[3] + *pInPartials * frequencies[nstate-1] );
				pInPartials++;
			}
			else {
				*pOutPartials = log(temp[0] + temp[1] + temp[2] + temp[3]);
			}
			
		}
		else {
			
			_mm_store_pd(temp, sum);
			
			*pOutPartials = log(temp[0] + temp[1] + temp[2] + temp[3]);
		}
		
		
		if ( tlk->scale ) {
			*pOutPartials += getLogScalingFactor( tlk, k);
		}
	}
	
    //	for ( int i = 0; i < tlk->pattern_count; i ++ ) {
    //		nodelog += outLogLikelihoods[i];
    //		//fprintf(stderr, " %f", outLogLikelihoods[i]);
    //	}
    //	fprintf(stderr, "\nlognode %f\n", nodelog);
}

void integrate_partials_codon_AVX( const SingleTreeLikelihood *tlk, const double *inPartials, const double *proportions, double *outPartials ){
	int k;
	double *pPartials = outPartials;
	const double *pInPartials = inPartials;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 4;
	
	__m128d in, pn;
	__m128d *pSSE = (__m128d*)outPartials;
	
	pn = _mm256_set_pd(proportions[0],proportions[0],proportions[0],proportions[0]);
	
	for ( k = 0; k < tlk->pattern_count; k++ ) {
		
		pSSE = (__m128d*)pPartials;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		in = _mm256_load_pd(pInPartials);
		_mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
		pPartials   += 4;
		pInPartials += 4;
		
		
		if ( extra != 0 ) {
			
			// finish the state, number of codons - 60 with a loop
			
			for ( int i = 0; i < extra_pairs; i++ ) {                
                in = _mm256_load_pd(pInPartials);
                _mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
                pPartials   += 4;
                pInPartials += 4;
			}
			
			// odd number so one more
			if ( extra & 1 ) {
				*pPartials++ = *pInPartials++ * proportions[0];
			}
			
		}
		
	}
	
	
	for ( int l = 1; l < tlk->cat_count; l++ ) {
		pn = _mm256_set_pd(proportions[l],proportions[l],proportions[l],proportions[l]);

		pPartials = outPartials;
		//pSSE = (__m128d*)outPartials;
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			pSSE = (__m256d*)pPartials;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;
			
			in = _mm_load_pd(pInPartials);
			in = _mm_mul_pd(in, pn);
			_mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
			pSSE++;
			pInPartials += 4;
			pPartials   += 4;			
						
			
			if ( extra != 0 ) {
				
				// finish the state, number of codons - 60 with a loop
				
				for ( int i = 0; i < extra_pairs; i++ ) {                    
                    in = _mm_load_pd(pInPartials);
                    in = _mm_mul_pd(in, pn);
                    _mm_store_pd( pPartials, _mm_add_pd(in, *pSSE) );
                    pSSE++;
                    pInPartials += 4;
                    pPartials   += 4;
				}
				
				// odd number so one more
				if ( extra & 1 ) {
					*pPartials += *pInPartials * proportions[l];
					pPartials++;
					pInPartials++;
				}
				
			}
			
		}
	}
	
    //	for ( int i = 0; i < tlk->partials_size; i ++ ) {
    //		integrate += outPartials[i];
    //	}
    //	fprintf(stderr, "integrate %f\n", integrate);
}

void update_partials_codon_AVX( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
	if( tlk->integrate_cat ){
		if( tlk->mapping[nodeIndex1] == -1 ){
			if(  tlk->mapping[nodeIndex2] == -1 ){
				partials_undefined_and_undefined_codon_AVX(tlk,
														   tlk->partials[nodeIndex1],
														   tlk->matrices[nodeIndex1],
														   tlk->partials[nodeIndex2],
														   tlk->matrices[nodeIndex2],
														   tlk->partials[nodeIndex3]);
			}
			else {
				partials_states_and_undefined_codon_AVX(tlk,
														tlk->mapping[nodeIndex2],
														tlk->matrices[nodeIndex2],
														tlk->partials[nodeIndex1],
														tlk->matrices[nodeIndex1],
														tlk->partials[nodeIndex3]);
			}
			
		}
		else{
			if(  tlk->mapping[nodeIndex2] == -1 ){
				partials_states_and_undefined_codon_AVX(tlk,
														tlk->mapping[nodeIndex1],
														tlk->matrices[nodeIndex1],
														tlk->partials[nodeIndex2],
														tlk->matrices[nodeIndex2],
														tlk->partials[nodeIndex3]);
				
			}
			else{
				partials_states_and_states_codon_AVX(tlk,
													 tlk->mapping[nodeIndex1],
													 tlk->matrices[nodeIndex1],
													 tlk->mapping[nodeIndex2],
													 tlk->matrices[nodeIndex2],
													 tlk->partials[nodeIndex3]);
			}
		}
	}
	else{
		
	}
	
	
	
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
	}
}

void partials_states_and_states_codon_AVX( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
	int k;
	int state1, state2;
	int u = 0;
	int w;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 4;
	
	double *pPartials = partials;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	
	__m256d m1v0, m2v0;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			state2 = tlk->sp->patterns[k][idx2];
			
			w = u;
			
			if (state1 < nstate && state2 < nstate) {
				
				m1 = &matrices1[w + nstate*state1];
				m2 = &matrices2[w + nstate*state2];
				
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
                
				m1v0 = _mm256_load_pd(m1);
				m2v0 = _mm256_load_pd(m2);
				_mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
				pPartials += 4;
				m1+=4; m2+=4;
				
				
				if ( extra != 0 ) {
                    
					// finish the state, number of codons - 60 with a loop
					
					for ( int i = 0; i < extra_pairs; i++ ) {                        
                        m1v0 = _mm256_load_pd(m1);
                        m2v0 = _mm256_load_pd(m2);
                        _mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
                        pPartials += 4;
                        m1+=4; m2+=4;
					}
					
					// odd number so one more
					if ( extra & 1 ) {
						*pPartials++ = *m1 * *m2;
					}
					
				}
				
			}
			else if (state1 < nstate) {
				// child 1 has a gap or unknown state so treat it as unknown
				m1 = &matrices1[w + nstate*state1];
				
				memcpy(pPartials, m1, nstate*sizeof(double));
				
				pPartials += nstate;
				
			}
			else if (state2 < nstate ) {
				// child 2 has a gap or unknown state so treat it as unknown
				m2 = &matrices2[w + nstate*state2];
				
				memcpy(pPartials, m2, nstate*sizeof(double));
				
				pPartials += nstate;
				
			}
			else {
				// both children have a gap or unknown state so set partials to 1
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0; *pPartials++ = 1.0;	*pPartials++ = 1.0;
				
				for ( int i = 0; i < extra; i++ ) {
					*pPartials++ = 1.0;
				}
			}
		}
		u += tlk->matrix_size;
	}
	
    //	for ( int i = 0; i < tlk->partials_size; i ++ ) {
    //		tot_state2 += partials[i];
    //	}
    //	fprintf(stderr, "state2 %f\n", tot_state2);
}



void partials_states_and_undefined_codon_AVX( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
	
	int v = 0;
	int k;
	int w = 0;
	int state1;
	
	double *pPartials = partials3;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	
	__m256d m2v0, sum;
	double temp[4] __attribute__ ((aligned (16)));
	__m256d *p2 = (__m256d *)partials2;
	
	const int nstate = tlk->sm->nstate;
	int extra = nstate - 60;
	int extra_pairs = extra / 4;
    
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
            
			m2 = matrices2 + w;
			
			if ( state1 < nstate ) {
				
				m1 = &matrices1[w + state1*nstate];
				
				for ( int j = 0; j < nstate; j++ ) {
					
					sum = _mm256_setzero_pd();
					
					p2 = (__m256d*)&partials2[v];
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					
                    
					
					if ( extra != 0 ) {
                        
						// finish the last states, number of codons - 60 with a loop
						
						for ( int i = 0; i < extra_pairs; i++ ) {
                            
                            m2v0 = _mm256_load_pd(m2);
                            m2v0 = _mm256_mul_pd(m2v0, *p2++);
                            sum  = _mm256_add_pd(sum, m2v0);
                            m2 += 4;
							
						}
						
						_mm_store_pd(temp, sum);
						
						// odd number so one more
						if ( extra & 1 ) {
							*pPartials++ = *m1 * ( temp[0] + temp[1] + temp[2] + temp[3] + (*m2 * partials2[v+nstate-1]) );
							m2++;
						}
						else {
							*pPartials++ = *m1 * ( temp[0] + temp[1] + temp[2] + temp[3] );
						}
                        
					}
					else {
						_mm_store_pd(temp, sum);
						*pPartials++ = *m1 * ( temp[0] + temp[1] + temp[2] + temp[3] );
					}
					
					m1++;
					
				}
				
				
			}
			else {
				// Child 1 has a gap or unknown state so don't use it
				
				for ( int j = 0; j < nstate; j++ ) {
					
					sum = _mm256_setzero_pd();
					
					p2 = (__m256d*)&partials2[v];
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;
					
					m2v0 = _mm256_load_pd(m2);
					m2v0 = _mm256_mul_pd(m2v0, *p2++);
					sum  = _mm256_add_pd(sum, m2v0);
					m2 += 4;					
					                    
					
					if ( extra != 0 ) {
						
						// finish the last states, number of codons - 60 with a loop
						
						for ( int i = 0; i < extra_pairs; i++ ) {
                            
                            m2v0 = _mm256_load_pd(m2);
                            m2v0 = _mm256_mul_pd(m2v0, *p2++);
                            sum  = _mm256_add_pd(sum, m2v0);
                            m2 += 4;
							
						}
						
						_mm_store_pd(temp, sum);
						
						// odd number so one more
						if ( extra & 1 ) {
							*pPartials++ = temp[0] + temp[1] + temp[2] + temp[3] + ( *m2 * partials2[v+nstate-1] );
							m2++;
						}
						else {
							*pPartials++ = temp[0] + temp[1] + temp[2] + temp[3];
						}
						
					}
					else {
						_mm_store_pd(temp, sum);
						*pPartials++ = temp[0] + temp[1] + temp[2] + temp[3];
					}
					
					m1++;
					
				}
				
			}
			v += nstate;
		}
		w += tlk->matrix_size;
	}
	
    //	for ( int i = 0; i < tlk->partials_size; i ++ ) {
    //		tot_state_undef += partials3[i];
    //	}
    //	fprintf(stderr, "state undef %f\n", tot_state_undef);
}

void partials_undefined_and_undefined_codon_AVX( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
	
	int v = 0;
	int k;
	int w = 0;
	
	double *pPartials = partials3;
	const double *m1 = matrices1;
	const double *m2 = matrices2;
	
	const int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
	const int extra_pairs = extra / 4;
	
	__m256d *p1 = (__m256d*)partials1;
	__m256d *p2 = (__m256d*)partials2;
	
	__m256d m1v0;
	
	double temp[4] __attribute__ ((aligned (16)));
	
	__m256d m1s, m2s;
	
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			m1 = matrices1 + w;
			m2 = matrices2 + w;
			
			for ( int j = 0; j < nstate; j++ ) {
				
				m1s = _mm_setzero_pd();
				m2s = _mm_setzero_pd();
				
				p1 = (__m256d*)&partials1[v];
				p2 = (__m256d*)&partials2[v];
                
				
				// matrices 1
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				m1v0 = _mm256_load_pd(m1); m1 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
				m1s = _mm256_add_pd(m1s, m1v0);
				
				
				
				
				// matrices 2
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
								
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				m1v0 = _mm256_load_pd(m2); m2 += 4;
				m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
				m2s = _mm256_add_pd(m2s, m1v0);
				
				
				if ( extra != 0 ) {
					
					// finish the state, number of codons - 60 with a loop
					
					if ( extra_pairs != 0 ){
						
                        m1v0 = _mm256_load_pd(m1); m1 += 4;
                        m1v0 = _mm256_mul_pd(m1v0, *p1); p1++;
                        m1s = _mm256_add_pd(m1s, m1v0);
                        
                        m1v0 = _mm256_load_pd(m2); m2 += 4;
                        m1v0 = _mm256_mul_pd(m1v0, *p2); p2++;
                        m2s = _mm256_add_pd(m2s, m1v0);
						
					}
					
					m1s = _mm256_hadd_pd(m1s, m2s);
					
					_mm256_store_pd(temp, m1s);
					
					// odd number so 1-3 more!!!!
                    // there should be a loop here, the other function
                    // the odd test is worng too in this case
					if ( extra & 1 ) {
                        
                        for ( int i = v + 60 + extra_pairs; i < v+nstate; i++ ) {
                            *pPartials = ( temp[0] + temp[1] + *m1 * partials1[i] ) * ( temp[2] + temp[3] + *m2 * partials2[i] );
                            m1++;
                            m2++;
                            pPartials++;
                        }
                        
					}
					else {
						
						*pPartials++ = (temp[0] + temp[1]) * (temp[2] * temp[3]);
					}
					
				}
				else {
					
					m1s = _mm256_hadd_pd(m1s, m2s);
					
					_mm256_store_pd(temp, m1s);
					
					*pPartials++ = temp[0] * temp[1] * temp[2] * temp[3];
				}
				
				
			}
			
			v += nstate;
			
		}
		w += tlk->matrix_size;
	}
	
}

#endif

#pragma mark -
#pragma mark Upper Likelihood

// Called by a node whose parent is NOT the root and the node's sibling is a leaf
static void _update_upper_partials_state( SingleTreeLikelihood *tlk, const double * restrict matrix_upper, const double * restrict partials_upper, const double * restrict matrix_lower, int sibling_index, double * restrict partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    int nstate = tlk->sm->nstate;
	
	memset(partials, 0, tlk->partials_size*sizeof(double));
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
            w2 = l * tlk->matrix_size;
            
            for ( int i = 0; i < nstate; i++ ) {
                
                pPartials[0]   += matrix_upper[w2] * partials_upper[v];
                pPartials[1]   += matrix_upper[w2+1] * partials_upper[v];
                pPartials[2]   += matrix_upper[w2+2] * partials_upper[v];
                pPartials[3]   += matrix_upper[w2+3] * partials_upper[v];
                pPartials[4]   += matrix_upper[w2+4] * partials_upper[v];
                pPartials[5]   += matrix_upper[w2+5] * partials_upper[v];
                pPartials[6]   += matrix_upper[w2+6] * partials_upper[v];
                pPartials[7]   += matrix_upper[w2+7] * partials_upper[v];
                pPartials[8]   += matrix_upper[w2+8] * partials_upper[v];
                pPartials[9]   += matrix_upper[w2+9] * partials_upper[v];
                
                pPartials[10]   += matrix_upper[w2+10] * partials_upper[v];
                pPartials[11]   += matrix_upper[w2+11] * partials_upper[v];
                pPartials[12]   += matrix_upper[w2+12] * partials_upper[v];
                pPartials[13]   += matrix_upper[w2+13] * partials_upper[v];
                pPartials[14]   += matrix_upper[w2+14] * partials_upper[v];
                pPartials[15]   += matrix_upper[w2+15] * partials_upper[v];
                pPartials[16]   += matrix_upper[w2+16] * partials_upper[v];
                pPartials[17]   += matrix_upper[w2+17] * partials_upper[v];
                pPartials[18]   += matrix_upper[w2+18] * partials_upper[v];
                pPartials[19]   += matrix_upper[w2+19] * partials_upper[v];
                
                pPartials[20]   += matrix_upper[w2+20] * partials_upper[v];
                pPartials[21]   += matrix_upper[w2+21] * partials_upper[v];
                pPartials[22]   += matrix_upper[w2+22] * partials_upper[v];
                pPartials[23]   += matrix_upper[w2+23] * partials_upper[v];
                pPartials[24]   += matrix_upper[w2+24] * partials_upper[v];
                pPartials[25]   += matrix_upper[w2+25] * partials_upper[v];
                pPartials[26]   += matrix_upper[w2+26] * partials_upper[v];
                pPartials[27]   += matrix_upper[w2+27] * partials_upper[v];
                pPartials[28]   += matrix_upper[w2+28] * partials_upper[v];
                pPartials[29]   += matrix_upper[w2+29] * partials_upper[v];                
                
                pPartials[30]   += matrix_upper[w2+30] * partials_upper[v];
                pPartials[31]   += matrix_upper[w2+31] * partials_upper[v];
                pPartials[32]   += matrix_upper[w2+32] * partials_upper[v];
                pPartials[33]   += matrix_upper[w2+33] * partials_upper[v];
                pPartials[34]   += matrix_upper[w2+34] * partials_upper[v];
                pPartials[35]   += matrix_upper[w2+35] * partials_upper[v];
                pPartials[36]   += matrix_upper[w2+36] * partials_upper[v];
                pPartials[37]   += matrix_upper[w2+37] * partials_upper[v];
                pPartials[38]   += matrix_upper[w2+38] * partials_upper[v];
                pPartials[39]   += matrix_upper[w2+39] * partials_upper[v];
                
                pPartials[40]   += matrix_upper[w2+40] * partials_upper[v];
                pPartials[41]   += matrix_upper[w2+41] * partials_upper[v];
                pPartials[42]   += matrix_upper[w2+42] * partials_upper[v];
                pPartials[43]   += matrix_upper[w2+43] * partials_upper[v];
                pPartials[44]   += matrix_upper[w2+44] * partials_upper[v];
                pPartials[45]   += matrix_upper[w2+45] * partials_upper[v];
                pPartials[46]   += matrix_upper[w2+46] * partials_upper[v];
                pPartials[47]   += matrix_upper[w2+47] * partials_upper[v];
                pPartials[48]   += matrix_upper[w2+48] * partials_upper[v];
                pPartials[49]   += matrix_upper[w2+49] * partials_upper[v];
                
                pPartials[50]   += matrix_upper[w2+50] * partials_upper[v];
                pPartials[51]   += matrix_upper[w2+51] * partials_upper[v];
                pPartials[52]   += matrix_upper[w2+52] * partials_upper[v];
                pPartials[53]   += matrix_upper[w2+53] * partials_upper[v];
                pPartials[54]   += matrix_upper[w2+54] * partials_upper[v];
                pPartials[55]   += matrix_upper[w2+55] * partials_upper[v];
                pPartials[56]   += matrix_upper[w2+56] * partials_upper[v];
                pPartials[57]   += matrix_upper[w2+57] * partials_upper[v];
                pPartials[58]   += matrix_upper[w2+58] * partials_upper[v];
                pPartials[59]   += matrix_upper[w2+59] * partials_upper[v];
                
                for ( int j = 60; j < nstate; j++ ) {
                    pPartials[j]  += matrix_upper[w2+j] * partials_upper[v];
				}
                v++;
                w2 += nstate;
            }
            
            if( state < nstate){        
	            w = l * tlk->matrix_size;
            	for ( int j = 0; j < nstate; j++ ) {
                    pPartials[j] *= matrix_lower[w+state];
                    w += nstate;
				}
            }
        	pPartials+=nstate;
		}
    }
}

static void _update_upper_partials_state2( SingleTreeLikelihood *tlk, const double * restrict matrix_upper, const double * restrict partials_upper, const double * restrict matrix_lower, int sibling_index, double * restrict partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
            w = l * tlk->matrix_size;
            
            for ( int i = 0; i < nstate; i++ ) {
                w2 = l * tlk->matrix_size + i;
                const double * pPartials_u = &partials_upper[v];
                
                *pPartials   = matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                for ( int j = 0; j < extra; j++ ) {
                    *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
				}
                
                if( state < nstate){
                    *pPartials *= matrix_lower[w+state];
                }
                pPartials++;
                w += nstate;
            }
			
			v += nstate;
		}
    }
}

static void _update_upper_partials_state3( SingleTreeLikelihood *tlk, const double * restrict matrix_upper, const double * restrict partials_upper, const double * restrict matrix_lower, int sibling_index, double * restrict partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
    double *mat = dvector(nstate*nstate);
    for(int i = 0; i < nstate; i++ ){
    for(int k = 0; k < nstate; k++ ){
    	mat[i*nstate+k] = matrix_upper[k*nstate+i];
    }
    }
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
            w2 = w = l * tlk->matrix_size;

            for ( int i = 0; i < nstate; i++ ) {
                const double * pPartials_u = &partials_upper[v];
                
                *pPartials   = mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                *pPartials  += mat[w2] * *pPartials_u++; w2++;
                
                for ( int j = 0; j < extra; j++ ) {
                    *pPartials  += mat[w2] * *pPartials_u++; w2++;
				}
                
                if( state < nstate){
                    *pPartials *= matrix_lower[w+state];
                }
                pPartials++;
                w += nstate;
            }
			
			v += nstate;
		}
    }
    free(mat);
}

// Called by a node whose parent is NOT the root and the node's sibling is NOT a leaf
static void _update_upper_partials_undefined( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ){
    int w,w2,k;
    int v = 0;
    double sum1,sum2;
    double *pPartials = partials;
    
    int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            w = l * tlk->matrix_size;
		            
            for ( int i = 0; i < nstate; i++ ) {
                w2 = l * tlk->matrix_size + i;
                const double * pPartials_u = &partials_upper[v];
                const double * pPartials_l = &partials_lower[v];
                
                sum1   = matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum2   = matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                for ( int j = 0; j < extra; j++ ) {
                    sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                    sum2  += matrix_lower[w] * *pPartials_l++;  w++;
				}
                
                *pPartials++ = sum1 * sum2 ;
            }
			
			v += nstate;
		}
    }
}

// Called by a child of the root but not if the child's sibling is NOT leaf
static void _update_upper_partials_root_and_undefined( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper ){
    int w,k;
    int v = 0;
	double *pPartials = partials_upper;
	const double *pPartials1 = partials1;
    int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			w = l * tlk->matrix_size;
            
            const double *m = &matrices1[w];
            
            for ( int i = 0; i < nstate; i++ ) {
                pPartials1 = &partials1[v];
                
                *pPartials   = *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                for ( int j = 0; j < extra; j++ ) {
					*pPartials  += *m++ * *pPartials1++;
				}
                
                *pPartials++ *= frequencies[i];
            }            
			
			v += nstate;
		}
    }
}

// Called by a child of the root and the child's sibling is a leaf
static void _update_upper_partials_root_and_state( const SingleTreeLikelihood *tlk, const double *matrices1, int idx1, const double *frequencies, double *partials_upper ){
    int w;
	double *pPartials = partials_upper;
    int state1;
    int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( int k = 0; k < tlk->pattern_count; k++ ) {
			
            state1 = tlk->sp->patterns[k][idx1];
            
			w = l * tlk->matrix_size;
            
            const double *ff = frequencies;
            
            if( state1 < nstate ){
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                *pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
                
                for ( int i = 0; i < extra; i++ ) {
					*pPartials++ = matrices1[w+state1] * *ff++; w+=nstate;
				}
            }
            else {
                memcpy(pPartials, frequencies, sizeof(double)*nstate);
                pPartials += nstate;
//				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
//				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
//				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
//				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
//				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
//				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
//				
//				// finish the state, number of codons - 60 with a loop
//				for ( int i = 0; i < extra; i++ ) {
//					*pPartials++ = 1.0;
//				}
            }
            
		}
    }
}

void update_partials_upper_codon( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);
    const double* freqs = tlk->get_root_frequencies(tlk);
	
    if( Node_isroot(parent) ){
        if( Node_isleaf(sibling) ){
            _update_upper_partials_root_and_state(tlk, tlk->matrices[ Node_id(sibling) ], tlk->mapping[ Node_id(sibling) ], freqs, tlk->partials[tlk->upper_partial_indexes[Node_id(node)]] );
        }
        else {
            _update_upper_partials_root_and_undefined(tlk, tlk->partials[ Node_id(sibling) ],  tlk->matrices[ Node_id(sibling) ], freqs, tlk->partials[tlk->upper_partial_indexes[Node_id(node)]] );
        }
    }
    else if( Node_isleaf(sibling) ){
        _update_upper_partials_state(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials[tlk->upper_partial_indexes[Node_id(parent)]], tlk->matrices[ Node_id(sibling) ], tlk->mapping[ Node_id(sibling) ], tlk->partials[tlk->upper_partial_indexes[Node_id(node)]]);
    }
    else {
        _update_upper_partials_undefined(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials[tlk->upper_partial_indexes[Node_id(parent)]], tlk->matrices[ Node_id(sibling) ], tlk->partials[ Node_id(sibling) ], tlk->partials[tlk->upper_partial_indexes[Node_id(node)]]);
    }
}

static void _partial_lower_upper( const SingleTreeLikelihood *tlk, const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int v = 0;
    double sum;
    double p = 0;
    int nstate = tlk->sm->nstate;
	const int extra = nstate - 60;
    
    memset(pattern_lk, 0, tlk->pattern_count*sizeof(double));
    
	for ( int l = 0; l < tlk->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->pattern_count; k++ ) {
			
			w = l * tlk->matrix_size;
            p = 0;
            
            for ( int i = 0; i < nstate; i++ ) {
                const double *pPartials_l = &partials_lower[v];
                
                sum   = matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                for ( int j = 0; j < extra; j++ ) {
                    sum  += matrix_lower[w] * *pPartials_l++; w++;                    
                }
                p += sum * partials_upper[v + i];
            }
			
            pattern_lk[k] += p * proportions[l];
			v += nstate;
		}
    }
}

static void _partial_lower_upper_leaf( const SingleTreeLikelihood *tlk, const double *partials_upper, int idx, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int v = 0;
    double p = 0;
    int state;
    int nstate    = tlk->sm->nstate;
    int sp_count  = tlk->pattern_count;
    int cat_count = tlk->cat_count;
    int matrix_size = tlk->matrix_size;
	const int extra = nstate - 60;
    
    memset(pattern_lk, 0, sp_count*sizeof(double));
    
	for ( int l = 0; l < cat_count; l++ ) {
		
		for ( k = 0; k < sp_count; k++ ) {
			state = tlk->sp->patterns[k][idx];
            
			w = l * matrix_size;
            const double *pPartials_u = &partials_upper[v];
            
            if( state < nstate ){
                p   = matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;
                
                for ( int j = 0; j < extra; j++ ) {
                    p  += matrix_lower[w+state] * *pPartials_u++; w += nstate;                    
                }
                
            }
            else {
                p  = *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                for ( int j = 0; j < extra; j++ ) {
                    p += *pPartials_u++;
                }
            }
            pattern_lk[k] += p * proportions[l];
			v += nstate;
		}
    }
}

void node_log_likelihoods_upper_codon( const SingleTreeLikelihood *tlk, Node *node ){
	int node_index = Node_id(node);
    
    if ( !Node_isleaf(node) ) {
        _partial_lower_upper(tlk, tlk->partials[tlk->upper_partial_indexes[node_index]], tlk->partials[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->pattern_lk+tlk->sp->count );
    }
    else {
        _partial_lower_upper_leaf(tlk, tlk->partials[tlk->upper_partial_indexes[node_index]], tlk->mapping[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->pattern_lk+tlk->sp->count );
    }
}

#pragma mark -
#pragma mark Upper Likelihood SSE

#ifdef SSE3_ENABLED
// Called by a node whose parent is NOT the root and the node's sibling is a leaf
// matrix_lower is transposed
static void _update_upper_partials_state_sse_codon( SingleTreeLikelihood *tlk, const double * restrict matrix_upper, const double * restrict partials_upper, const double * restrict matrix_lower, int sibling_index, double * restrict partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    int nstate = tlk->sm->nstate;
    
    memset(partials, 0, tlk->partials_size*sizeof(double));
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
            w2 = l * tlk->matrix_size;
            
            for ( int i = 0; i < nstate; i++ ) {
            
                pPartials[0]   += matrix_upper[w2] * partials_upper[v];
                pPartials[1]   += matrix_upper[w2+1] * partials_upper[v];
                pPartials[2]   += matrix_upper[w2+2] * partials_upper[v];
                pPartials[3]   += matrix_upper[w2+3] * partials_upper[v];
                pPartials[4]   += matrix_upper[w2+4] * partials_upper[v];
                pPartials[5]   += matrix_upper[w2+5] * partials_upper[v];
                pPartials[6]   += matrix_upper[w2+6] * partials_upper[v];
                pPartials[7]   += matrix_upper[w2+7] * partials_upper[v];
                pPartials[8]   += matrix_upper[w2+8] * partials_upper[v];
                pPartials[9]   += matrix_upper[w2+9] * partials_upper[v];
                
                pPartials[10]   += matrix_upper[w2+10] * partials_upper[v];
                pPartials[11]   += matrix_upper[w2+11] * partials_upper[v];
                pPartials[12]   += matrix_upper[w2+12] * partials_upper[v];
                pPartials[13]   += matrix_upper[w2+13] * partials_upper[v];
                pPartials[14]   += matrix_upper[w2+14] * partials_upper[v];
                pPartials[15]   += matrix_upper[w2+15] * partials_upper[v];
                pPartials[16]   += matrix_upper[w2+16] * partials_upper[v];
                pPartials[17]   += matrix_upper[w2+17] * partials_upper[v];
                pPartials[18]   += matrix_upper[w2+18] * partials_upper[v];
                pPartials[19]   += matrix_upper[w2+19] * partials_upper[v];
                
                pPartials[20]   += matrix_upper[w2+20] * partials_upper[v];
                pPartials[21]   += matrix_upper[w2+21] * partials_upper[v];
                pPartials[22]   += matrix_upper[w2+22] * partials_upper[v];
                pPartials[23]   += matrix_upper[w2+23] * partials_upper[v];
                pPartials[24]   += matrix_upper[w2+24] * partials_upper[v];
                pPartials[25]   += matrix_upper[w2+25] * partials_upper[v];
                pPartials[26]   += matrix_upper[w2+26] * partials_upper[v];
                pPartials[27]   += matrix_upper[w2+27] * partials_upper[v];
                pPartials[28]   += matrix_upper[w2+28] * partials_upper[v];
                pPartials[29]   += matrix_upper[w2+29] * partials_upper[v];                
                
                pPartials[30]   += matrix_upper[w2+30] * partials_upper[v];
                pPartials[31]   += matrix_upper[w2+31] * partials_upper[v];
                pPartials[32]   += matrix_upper[w2+32] * partials_upper[v];
                pPartials[33]   += matrix_upper[w2+33] * partials_upper[v];
                pPartials[34]   += matrix_upper[w2+34] * partials_upper[v];
                pPartials[35]   += matrix_upper[w2+35] * partials_upper[v];
                pPartials[36]   += matrix_upper[w2+36] * partials_upper[v];
                pPartials[37]   += matrix_upper[w2+37] * partials_upper[v];
                pPartials[38]   += matrix_upper[w2+38] * partials_upper[v];
                pPartials[39]   += matrix_upper[w2+39] * partials_upper[v];
                
                pPartials[40]   += matrix_upper[w2+40] * partials_upper[v];
                pPartials[41]   += matrix_upper[w2+41] * partials_upper[v];
                pPartials[42]   += matrix_upper[w2+42] * partials_upper[v];
                pPartials[43]   += matrix_upper[w2+43] * partials_upper[v];
                pPartials[44]   += matrix_upper[w2+44] * partials_upper[v];
                pPartials[45]   += matrix_upper[w2+45] * partials_upper[v];
                pPartials[46]   += matrix_upper[w2+46] * partials_upper[v];
                pPartials[47]   += matrix_upper[w2+47] * partials_upper[v];
                pPartials[48]   += matrix_upper[w2+48] * partials_upper[v];
                pPartials[49]   += matrix_upper[w2+49] * partials_upper[v];
                
                pPartials[50]   += matrix_upper[w2+50] * partials_upper[v];
                pPartials[51]   += matrix_upper[w2+51] * partials_upper[v];
                pPartials[52]   += matrix_upper[w2+52] * partials_upper[v];
                pPartials[53]   += matrix_upper[w2+53] * partials_upper[v];
                pPartials[54]   += matrix_upper[w2+54] * partials_upper[v];
                pPartials[55]   += matrix_upper[w2+55] * partials_upper[v];
                pPartials[56]   += matrix_upper[w2+56] * partials_upper[v];
                pPartials[57]   += matrix_upper[w2+57] * partials_upper[v];
                pPartials[58]   += matrix_upper[w2+58] * partials_upper[v];
                pPartials[59]   += matrix_upper[w2+59] * partials_upper[v];
                
                for ( int j = 60; j < nstate; j++ ) {
                    pPartials[j]  += matrix_upper[w2+j] * partials_upper[v];
				}
                v++;
                w2 += nstate;
            }
            
            if( state < nstate){        
	            w = l * tlk->matrix_size + state*nstate;
            	for ( int j = 0; j < nstate; j++ ) {
                    pPartials[j] *= matrix_lower[w];
                    w++;
				}
            }
        	pPartials+=nstate;
        }
    }
}

// Called by a node whose parent is NOT the root and the node's sibling is a leaf
// matrix_lower is transposed
static void _update_upper_partials_state_sse_codon2( SingleTreeLikelihood *tlk, const double * restrict matrix_upper, const double * restrict partials_upper, const double * restrict matrix_lower, int sibling_index, double * restrict partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    int nstate = tlk->sm->nstate;
    const int extra = nstate - 60;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
            w = l * tlk->matrix_size + nstate * state;
            
            for ( int i = 0; i < nstate; i++ ) {
                w2 = l * tlk->matrix_size + i;
                
                const double * pPartials_u = &partials_upper[v];
                
                *pPartials   = matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                for ( int j = 0; j < extra; j++ ) {
                    *pPartials  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                }
                
                if( state < nstate){
                    *pPartials *= matrix_lower[w];
                }
                pPartials++;
                w++;
            }
            
            v += nstate;
        }
    }
}

// Called by a node whose parent is NOT the root and the node's sibling is NOT a leaf
// All matrices are NOT transposed
static void _update_upper_partials_undefined_sse_codon( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ){
    int w,w2,k;
    int v = 0;
    double sum1,sum2;
    double *pPartials = partials;
    
    int nstate = tlk->sm->nstate;
    const int extra = nstate - 60;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            w = l * tlk->matrix_size;
            
            for ( int i = 0; i < nstate; i++ ) {
                w2 = l * tlk->matrix_size + i;
                const double * pPartials_u = &partials_upper[v];
                const double * pPartials_l = &partials_lower[v];
                
                sum1   = matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                
                sum2   = matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                sum2  += matrix_lower[w] * *pPartials_l++; w++;
                
                for ( int j = 0; j < extra; j++ ) {
                    sum1  += matrix_upper[w2] * *pPartials_u++; w2 += nstate;
                    sum2  += matrix_lower[w] * *pPartials_l++;  w++;
                }
                
                *pPartials++ = sum1 * sum2 ;
            }
            
            v += nstate;
        }
    }
}

// Called by a child of the root but not if the child's sibling is NOT leaf
// All matrices are NOT transposed
static void _update_upper_partials_root_and_undefined_sse_codon( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper ){
    int w,k;
    int v = 0;
    double *pPartials = partials_upper;
    const double *pPartials1 = partials1;
    int nstate = tlk->sm->nstate;
    const int extra = nstate - 60;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            w = l * tlk->matrix_size;
            
            const double *m = &matrices1[w];
            
            for ( int i = 0; i < nstate; i++ ) {
                pPartials1 = &partials1[v];
                
                *pPartials   = *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                *pPartials  += *m++ * *pPartials1++;
                
                for ( int j = 0; j < extra; j++ ) {
                    *pPartials  += *m++ * *pPartials1++;
                }
                
                *pPartials++ *= frequencies[i];
            }
            
            v += nstate;
        }
    }
}

// Called by a child of the root and the child's sibling is a leaf
// matrices1 is transposed
static void _update_upper_partials_root_and_state_sse_codon( const SingleTreeLikelihood *tlk, const double *matrices1, int idx1, const double *frequencies, double *partials_upper ){
    int w;
    double *pPartials = partials_upper;
    int state1;
    int nstate = tlk->sm->nstate;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( int k = 0; k < tlk->pattern_count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            
            w = l * tlk->matrix_size + state1 * nstate;
            
            const double *m1 = &matrices1[w];
            
            if( state1 < nstate ){
                *pPartials++ = m1[0] * frequencies[0];
                *pPartials++ = m1[1] * frequencies[1];
                *pPartials++ = m1[2] * frequencies[2];
                *pPartials++ = m1[3] * frequencies[3];
                *pPartials++ = m1[4] * frequencies[4];
                *pPartials++ = m1[5] * frequencies[5];
                *pPartials++ = m1[6] * frequencies[6];
                *pPartials++ = m1[7] * frequencies[7];
                *pPartials++ = m1[8] * frequencies[8];
                *pPartials++ = m1[9] * frequencies[9];
                
                *pPartials++ = m1[10] * frequencies[10];
                *pPartials++ = m1[11] * frequencies[11];
                *pPartials++ = m1[12] * frequencies[12];
                *pPartials++ = m1[13] * frequencies[13];
                *pPartials++ = m1[14] * frequencies[14];
                *pPartials++ = m1[15] * frequencies[15];
                *pPartials++ = m1[16] * frequencies[16];
                *pPartials++ = m1[17] * frequencies[17];
                *pPartials++ = m1[18] * frequencies[18];
                *pPartials++ = m1[19] * frequencies[19];
                
                *pPartials++ = m1[20] * frequencies[20];
                *pPartials++ = m1[21] * frequencies[21];
                *pPartials++ = m1[22] * frequencies[22];
                *pPartials++ = m1[23] * frequencies[23];
                *pPartials++ = m1[24] * frequencies[24];
                *pPartials++ = m1[25] * frequencies[25];
                *pPartials++ = m1[26] * frequencies[26];
                *pPartials++ = m1[27] * frequencies[27];
                *pPartials++ = m1[28] * frequencies[28];
                *pPartials++ = m1[29] * frequencies[29];
                
                *pPartials++ = m1[30] * frequencies[30];
                *pPartials++ = m1[31] * frequencies[31];
                *pPartials++ = m1[32] * frequencies[32];
                *pPartials++ = m1[33] * frequencies[33];
                *pPartials++ = m1[34] * frequencies[34];
                *pPartials++ = m1[35] * frequencies[35];
                *pPartials++ = m1[36] * frequencies[36];
                *pPartials++ = m1[37] * frequencies[37];
                *pPartials++ = m1[38] * frequencies[38];
                *pPartials++ = m1[39] * frequencies[39];
                
                *pPartials++ = m1[40] * frequencies[40];
                *pPartials++ = m1[41] * frequencies[41];
                *pPartials++ = m1[42] * frequencies[42];
                *pPartials++ = m1[43] * frequencies[43];
                *pPartials++ = m1[44] * frequencies[44];
                *pPartials++ = m1[45] * frequencies[45];
                *pPartials++ = m1[46] * frequencies[46];
                *pPartials++ = m1[47] * frequencies[47];
                *pPartials++ = m1[48] * frequencies[48];
                *pPartials++ = m1[49] * frequencies[49];
                
                *pPartials++ = m1[50] * frequencies[50];
                *pPartials++ = m1[51] * frequencies[51];
                *pPartials++ = m1[52] * frequencies[52];
                *pPartials++ = m1[53] * frequencies[53];
                *pPartials++ = m1[54] * frequencies[54];
                *pPartials++ = m1[55] * frequencies[55];
                *pPartials++ = m1[56] * frequencies[56];
                *pPartials++ = m1[57] * frequencies[57];
                *pPartials++ = m1[58] * frequencies[58];
                *pPartials++ = m1[59] * frequencies[59];
                
                for ( int i = 60 ; i < nstate; i++ ) {
                    *pPartials++ = m1[i] * frequencies[i];
                }
            }
            else {
                memcpy(pPartials, frequencies, sizeof(double)*nstate);
                pPartials += nstate;
                //				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                //				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                //				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                //				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                //				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                //				*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;*pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                //
                //				// finish the state, number of codons - 60 with a loop
                //				for ( int i = 0; i < extra; i++ ) {
                //					*pPartials++ = 1.0;
                //				}
            }
            
        }
    }
}

void update_partials_upper_sse_codon( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);
    const double* freqs = tlk->get_root_frequencies(tlk);
	
    if( Node_isroot(parent) ){
        // The matrix of the sibling is transposed
        if( Node_isleaf(sibling) ){
            _update_upper_partials_root_and_state_sse_codon(tlk, tlk->matrices[ Node_id(sibling) ], tlk->mapping[ Node_id(sibling) ], freqs, tlk->partials[tlk->upper_partial_indexes[Node_id(node)]] );
        }
        else {
            _update_upper_partials_root_and_undefined_sse_codon(tlk, tlk->partials[ Node_id(sibling) ],  tlk->matrices[ Node_id(sibling) ],  freqs, tlk->partials[tlk->upper_partial_indexes[Node_id(node)]] );
        }
    }
    // The matrix of the sibling is transposed
    // The pparent node cannot be leaf
    else if( Node_isleaf(sibling) ){
        _update_upper_partials_state_sse_codon(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials[tlk->upper_partial_indexes[Node_id(parent)]], tlk->matrices[ Node_id(sibling) ], tlk->mapping[ Node_id(sibling) ], tlk->partials[tlk->upper_partial_indexes[Node_id(node)]]);
    }
    else {
        _update_upper_partials_undefined_sse_codon(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials[tlk->upper_partial_indexes[Node_id(parent)]], tlk->matrices[ Node_id(sibling) ], tlk->partials[ Node_id(sibling) ], tlk->partials[tlk->upper_partial_indexes[Node_id(node)]]);
    }
}

static void _partial_lower_upper_sse_codon( const SingleTreeLikelihood *tlk, const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int v = 0;
    double sum;
    double p = 0;
    int nstate = tlk->sm->nstate;
    const int extra = nstate - 60;
    
    memset(pattern_lk, 0, tlk->pattern_count*sizeof(double));
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            w = l * tlk->matrix_size;
            p = 0;
            
            for ( int i = 0; i < nstate; i++ ) {
                const double *pPartials_l = &partials_lower[v];
                
                sum   = matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                sum  += matrix_lower[w] * *pPartials_l++; w++;
                
                for ( int j = 0; j < extra; j++ ) {
                    sum  += matrix_lower[w] * *pPartials_l++; w++;
                }
                p += sum * partials_upper[v + i];
            }
            
            pattern_lk[k] += p * proportions[l];
            v += nstate;
        }
    }
}

// matrix_lower is transposed
static void _partial_lower_upper_leaf_sse_codon2( const SingleTreeLikelihood *tlk, const double *partials_upper, int idx, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int v = 0;
    double p;
    int state;
    int nstate    = tlk->sm->nstate;
    int sp_count  = tlk->pattern_count;
    int cat_count = tlk->cat_count;
    int matrix_size = tlk->matrix_size;
    
    __m128d *m;
    __m128d temp;
    __m128d *pu;
    
    double t[2] __attribute__ ((aligned (16)));
    
    memset(pattern_lk, 0, sp_count*sizeof(double));
    
    for ( int l = 0; l < cat_count; l++ ) {
        
        for ( k = 0; k < sp_count; k++ ) {
            state = tlk->sp->patterns[k][idx];
            p = 0;
            w = l * matrix_size + nstate*state;
            int j = 60;
            
            if( state < nstate ){
            
            	if( nstate & 1 ){
            		if( v & 1){
            		
            		}
            		
					if( k & 1 ){
						p  = matrix_lower[w] * partials_upper[v];
						w++;
						v++;
						j++;
					}
            	}
            	
                m = (__m128d*)&matrix_lower[w];
                pu = (__m128d*)&partials_upper[v];
                
                
                temp = _mm_mul_pd(*m,*pu); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); pu++;
                
                for ( ; j < nstate; j++ ) {
                    p  += matrix_lower[w+j] * partials_upper[v+j];
                }
                
            }
            else {
            
            	if( nstate & 1 ){
					if( k & 1 ){
						p  = partials_upper[v];
						v++;
						j++;
					}
            	}
            	
                pu = (__m128d*)&partials_upper[v];
                
                temp = *pu; pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                temp = _mm_add_pd(temp, *pu); pu++;
                
                for ( ; j < nstate; j++ ) {
                    p += partials_upper[v+j];
                }
            }
            _mm_store_pd(t, temp);
            
            pattern_lk[k] += (t[0]+t[1]+p) * proportions[l];
            v += nstate;
        }
    }
}

// matrix_lower is transposed
static void _partial_lower_upper_leaf_sse_codon( const SingleTreeLikelihood *tlk, const double *partials_upper, int idx, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int v = 0;
    double p = 0;
    int state;
    int nstate    = tlk->sm->nstate;
    int sp_count  = tlk->pattern_count;
    int cat_count = tlk->cat_count;
    int matrix_size = tlk->matrix_size;
    const int extra = nstate - 60;
    
    memset(pattern_lk, 0, sp_count*sizeof(double));
    
    for ( int l = 0; l < cat_count; l++ ) {
        
        for ( k = 0; k < sp_count; k++ ) {
            state = tlk->sp->patterns[k][idx];
            
            w = l * matrix_size + state * nstate;
            const double *pPartials_u = &partials_upper[v];
            
            if( state < nstate ){
                p   = matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                p  += matrix_lower[w] * *pPartials_u++; w++;
                
                for ( int j = 0; j < extra; j++ ) {
                    p  += matrix_lower[w] * *pPartials_u++; w++;
                }
                
            }
            else {
                p  = *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                p += *pPartials_u++;
                
                for ( int j = 0; j < extra; j++ ) {
                    p += *pPartials_u++;
                }
            }
            pattern_lk[k] += p * proportions[l];
            v += nstate;
        }
    }
}

void node_log_likelihoods_upper_sse_codon( const SingleTreeLikelihood *tlk, Node *node ){
    int node_index = Node_id(node);
    
    if ( !Node_isleaf(node) ) {
        _partial_lower_upper_sse_codon(tlk, tlk->partials[tlk->upper_partial_indexes[node_index]], tlk->partials[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->pattern_lk+tlk->sp->count );
    }
    else {
        _partial_lower_upper_leaf_sse_codon(tlk, tlk->partials[tlk->upper_partial_indexes[node_index]], tlk->mapping[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->pattern_lk+tlk->sp->count );
    }
}
#endif
