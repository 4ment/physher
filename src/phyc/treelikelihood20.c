/*
 *  treelikelihood4.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/09/2014.
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

#include "treelikelihood20.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "matrix.h"

#ifdef SSE3_ENABLED
//#include <xmmintrin.h> // SSE
//#include <emmintrin.h> // SSE2
#include <pmmintrin.h> // SSE3
//#include <tmmintrin.h> // SSSE3
#endif

#if 1
#define RESTRICT restrict
#else
#define RESTRICT
#endif

#pragma mark -
#pragma mark Lower Likelihood SSE

#include "treelikelihoodX.h"

#ifdef SSE3_ENABLED

void update_partials_20_SSE( SingleTreeLikelihood *tlk, int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
	
	if( tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1] != NULL ){
		if(  tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2] != NULL ){
			partials_undefined_and_undefined_20_SSE(tlk,
												   tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		else {
			partials_states_and_undefined_20_SSE(tlk,
										  tlk->mapping[partialsIndex2],
										  tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
										  tlk->partials[tlk->current_partials_indexes[partialsIndex1]][partialsIndex1],
										  tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
										  tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
		
	}
	else{
		if(  tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2] != NULL ){
			partials_states_and_undefined_20_SSE(tlk,
										  tlk->mapping[partialsIndex1],
										  tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
										  tlk->partials[tlk->current_partials_indexes[partialsIndex2]][partialsIndex2],
										  tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
										  tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
			
		}
		else{
			partials_states_and_states_20_SSE(tlk,
									   tlk->mapping[partialsIndex1],
									   tlk->matrices[tlk->current_matrices_indexes[matrixIndex1]][matrixIndex1],
									   tlk->mapping[partialsIndex2],
									   tlk->matrices[tlk->current_matrices_indexes[matrixIndex2]][matrixIndex2],
									   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex]);
		}
	}

	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, partialsIndex);
	}
}


void partials_states_and_states_20_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT matrices2, double *partials ){

    int k;
    int state1, state2;
    int u = 0;
    int w;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    matrices1 = __builtin_assume_aligned(matrices1, 16);
    matrices2 = __builtin_assume_aligned(matrices2, 16);
#endif
    
    double *pPartials = partials;
    
    __m128d *m1,*m2;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            
            w = u;
            
            if (state1 < 20 && state2 < 20) {
                
                m1 = (__m128d*)&matrices1[w+20*state1];
                m2 = (__m128d*)&matrices2[w+20*state2];
                
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2)); m1++; m2++; pPartials+=2;
                _mm_store_pd(pPartials, _mm_mul_pd(*m1, *m2));  pPartials+=2;
            }
            else if (state1 < 20) {
                // child 2 has a gap or unknown state so treat it as unknown
                memcpy(pPartials, &matrices1[w+20*state1], sizeof(double)*20);
                pPartials += 20;
            }
            else if (state2 < 20 ) {
                // child 1 has a gap or unknown state so treat it as unknown
                memcpy(pPartials, &matrices2[w+20*state2], sizeof(double)*20);
                pPartials += 20;
            }
            else {
                // both children have a gap or unknown state so set partials to 1
                *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0; *pPartials++ = 1.0;
                
            }
        }
        u += 400;
    }
}
    

void partials_states_and_undefined_20_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, const double * RESTRICT partials2, const double * RESTRICT matrices2, double *partials3 ){
    
    int v = 0;
    int k;
    int w = 0;
    int state1;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    partials2 = __builtin_assume_aligned(partials2, 16);
    matrices1 = __builtin_assume_aligned(matrices1, 16);
    matrices2 = __builtin_assume_aligned(matrices2, 16);
#endif
    
    double *pPartials = partials3;
    const double *m1 = matrices1;
    
    __m128d *m2,*p2,res;
    
    double temp[2] __attribute__ ((aligned (16)));
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            
            m2 = (__m128d*)&matrices2[w];
            
            if ( state1 < 20) {
                
                m1 = &matrices1[w + state1*20];
                
                for ( int j = 0; j < 5; j++ ) {
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = *m1 * (temp[0]+temp[1]);
                    m1++;
                    
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = *m1 * (temp[0]+temp[1]);
                    m1++;
                    
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = *m1 * (temp[0]+temp[1]);
                    m1++;
                    
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = *m1 * (temp[0]+temp[1]);
                    m1++;
                }
            
            }
            else {
                // Child 1 has a gap or unknown state so don't use it
                
                for ( int j = 0; j < 5; j++ ) {
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = temp[0]+temp[1];
                    
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = temp[0]+temp[1];
                    
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = temp[0]+temp[1];
                    
                    
                    p2 = (__m128d*)&partials2[v];
                    res = _mm_mul_pd(*m2, *p2); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++; p2++;
                    res = _mm_add_pd(res, _mm_mul_pd(*m2, *p2)); m2++;
                    
                    _mm_store_pd(temp, res);
                    
                    *pPartials++ = temp[0]+temp[1];
                }
                
            }
            v += 20;
        }
        w += 400;
    }
}
    

void partials_undefined_and_undefined_20_SSE( const SingleTreeLikelihood *tlk, const double * RESTRICT partials1, const double * RESTRICT matrices1, const double * RESTRICT partials2, const double * RESTRICT matrices2, double *partials3 ){
    
    int v = 0;
    int k;
    int w = 0;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    partials1 = __builtin_assume_aligned(partials1, 16);
    partials2 = __builtin_assume_aligned(partials2, 16);
    
    matrices1 = __builtin_assume_aligned(matrices1, 16);
    matrices2 = __builtin_assume_aligned(matrices2, 16);
#endif
    
    double *pPartials = partials3;
    
    double temp[2] __attribute__ ((aligned (16)));
    
    __m128d *m1,*m2,*p,res,res2;
    
    for ( int l = 0; l < tlk->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->pattern_count; k++ ) {
            
            m1 = (__m128d*)&matrices1[w];
            m2 = (__m128d*)&matrices2[w];
            
            for ( int j = 0; j < 5; j++ ) {
                
                p = (__m128d*)&partials1[v];
                
                res = _mm_mul_pd(*m1, *p); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++;
                
                p = (__m128d*)&partials2[v];
                
                res2 = _mm_mul_pd(*m2, *p); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++;
                
                _mm_store_pd(temp, _mm_hadd_pd(res, res2));
                
                *pPartials++ = temp[0]*temp[1];
                
                
                p = (__m128d*)&partials1[v];
                
                res = _mm_mul_pd(*m1, *p); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++;
                
                p = (__m128d*)&partials2[v];
                
                res2 = _mm_mul_pd(*m2, *p); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++;
                
                _mm_store_pd(temp, _mm_hadd_pd(res, res2));
                
                *pPartials++ = temp[0]*temp[1];
                
                
                p = (__m128d*)&partials1[v];
                
                res = _mm_mul_pd(*m1, *p); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++;
                
                p = (__m128d*)&partials2[v];
                
                res2 = _mm_mul_pd(*m2, *p); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++;
                
                _mm_store_pd(temp, _mm_hadd_pd(res, res2));
                
                *pPartials++ = temp[0]*temp[1];
                
                
                p = (__m128d*)&partials1[v];
                
                res = _mm_mul_pd(*m1, *p); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++; p++;
                res = _mm_add_pd(res, _mm_mul_pd(*m1, *p)); m1++;
                
                p = (__m128d*)&partials2[v];
                
                res2 = _mm_mul_pd(*m2, *p); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++; p++;
                res2 = _mm_add_pd(res2, _mm_mul_pd(*m2, *p)); m2++;
                
                _mm_store_pd(temp, _mm_hadd_pd(res, res2));
                
                *pPartials++ = temp[0]*temp[1];
                
            }
            v += 20;
        }
        w += 400;
    }
}
#endif
    
#pragma mark -
#pragma mark Upper Likelihood

static void _calculate_branch_likelihood_undefined(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*20*tlk->pattern_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			const double* partialsChildPtr = partials+v;
			const double* transMatrixPtr = &matrices[l*400];
			
			double sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
			sum = *transMatrixPtr * partialsChildPtr[0]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[1]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[2]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[3]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[4]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[5]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[6]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[7]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[8]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[9]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[10]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[11]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[12]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[13]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[14]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[15]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[16]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[17]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[18]; transMatrixPtr++;
			sum += *transMatrixPtr * partialsChildPtr[19]; transMatrixPtr++;
			rootPartials[u] += sum * upperPartials[v] * weight;u++;v++;
			
		}
	}
}

static void _calculate_branch_likelihood_state(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, int partialsIndex, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->sm->nstate*tlk->pattern_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			int state = tlk->sp->patterns[k][partialsIndex];
			if(state < 20){
				int w =  l * 400 + state;
				
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight; w += 20;u++;v++;
				rootPartials[u] += matrices[w] * upperPartials[v] * weight;u++;v++;
			}
			else{
				const double* transMatrixPtr = &matrices[l*400];
				for(int i = 0; i < 20; i++){
					double sum = 0;
					for(int j = 0; j < 20; j++){
						sum += transMatrixPtr[i*20+j];
					}
					rootPartials[u++] += sum * weight * upperPartials[v++];
				}
			}
		}
	}
	
}

static void _calculate_branch_likelihood_state2(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->sm->nstate*tlk->pattern_count);
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			const double* partialsChildPtr = partials+v;
			int state = tlk->sp->patterns[k][upperPartialsIndex];
			const double* transMatrixPtr = NULL;
			double sum;
			
			if(state < 20){
				transMatrixPtr = &matrices[l * 400 + state*20];
				
				sum = transMatrixPtr[0] * partialsChildPtr[0];
				sum += transMatrixPtr[1] * partialsChildPtr[1];
				sum += transMatrixPtr[2] * partialsChildPtr[2];
				sum += transMatrixPtr[3] * partialsChildPtr[3];
				sum += transMatrixPtr[4] * partialsChildPtr[4];
				sum += transMatrixPtr[5] * partialsChildPtr[5];
				sum += transMatrixPtr[6] * partialsChildPtr[6];
				sum += transMatrixPtr[7] * partialsChildPtr[7];
				sum += transMatrixPtr[8] * partialsChildPtr[8];
				sum += transMatrixPtr[9] * partialsChildPtr[9];
				sum += transMatrixPtr[10] * partialsChildPtr[10];
				sum += transMatrixPtr[11] * partialsChildPtr[11];
				sum += transMatrixPtr[12] * partialsChildPtr[12];
				sum += transMatrixPtr[13] * partialsChildPtr[13];
				sum += transMatrixPtr[14] * partialsChildPtr[14];
				sum += transMatrixPtr[15] * partialsChildPtr[15];
				sum += transMatrixPtr[16] * partialsChildPtr[16];
				sum += transMatrixPtr[17] * partialsChildPtr[17];
				sum += transMatrixPtr[18] * partialsChildPtr[18];
				sum += transMatrixPtr[19] * partialsChildPtr[19];
				rootPartials[u+state] += sum * weight;
				u+=20;
			}
			else{
				transMatrixPtr = &matrices[l*400];
				for(int i = 0; i < 20; i++){
					double sum = 0;
					for(int j = 0; j < 20; j++){
						sum += transMatrixPtr[i*20+j] * partialsChildPtr[j];
					}
					rootPartials[u] += sum * weight; u++;
				}
			}
			v+=20;
		}
	}
	
}

void calculate_branch_likelihood_20(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex){
	// Sibling is a leaf
	if( tlk->partials[0][partialsIndex] == NULL ){
		_calculate_branch_likelihood_state(tlk, rootPartials,
										   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
										   tlk->mapping[partialsIndex],
										   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	// upper is coming from the other side of the tree (right node) and right node is a leaf
	// right node should not be a leaf
	else if(upperPartialsIndex < tlk->sp->size){
		exit(1);
	}
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
	


#pragma mark -
#pragma mark Upper Likelihood SSE

#ifdef SSE3_ENABLED
static void _calculate_branch_likelihood_undefined_20_SSE(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*20*tlk->pattern_count);
	int v = 0;
	__m128d* m, *p;
	__m128d temp;

	double t[2] __attribute__ ((aligned (16)));

	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			m = (__m128d*)&matrices[l*400];
			for(int j = 0; j < 20; j++){
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * upperPartials[v+j] * weight;
			}
			v += 20;

		}
	}
}

static void _calculate_branch_likelihood_upper_undefined_20_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, const double* partials, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*tlk->sm->nstate*tlk->pattern_count);
	
	double t[2] __attribute__ ((aligned (16)));
	__m128d temp;
	int v = 0;
	for(int l = 0; l < tlk->cat_count; l++) {
		int u = 0;
		const double weight = tlk->sm->get_proportion(tlk->sm, l);
		for(int k = 0; k < tlk->pattern_count; k++) {
			__m128d* p = (__m128d*)(partials+v);
			int state = tlk->sp->patterns[k][upperPartialsIndex];
			
			if(state < 20){
				__m128d* m = (__m128d*)&matrices[l * 400 + state*20];
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p));
				
				_mm_store_pd(t,temp);
				rootPartials[u+state] += (t[0]+t[1]) * weight;
				u+=20;
			}
			else{
				__m128d* m = (__m128d*)&matrices[l*400];
				
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
				p = (__m128d*)(partials+v);
				temp = _mm_mul_pd(*m, *p); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++; p++;
				temp = _mm_add_pd(temp, _mm_mul_pd(*m,*p)); m++;
				_mm_store_pd(t,temp);
				rootPartials[u++] += (t[0]+t[1]) * weight;
				
			}
			v+=20;
		}
	}
}

// matrices are transposed
static void _calculate_branch_likelihood_state_20_SSE(SingleTreeLikelihood *tlk, double* rootPartials, const double* upperPartials, int partialsIndex, const double* matrices){
	memset(rootPartials, 0, sizeof(double)*20*tlk->pattern_count);
	int v = 0;
	__m128d* up = (__m128d*)upperPartials;
	__m128d weight;
	__m128d* rp, *m1;
	for(int l = 0; l < tlk->cat_count; l++) {
		rp = (__m128d*)rootPartials;
		double ww = tlk->sm->get_proportion(tlk->sm, l);
		weight = _mm_set1_pd(ww);
		const double* mat = &matrices[l*400];
		int u = 0;
		for(int k = 0; k < tlk->pattern_count; k++) {
			const int state = tlk->sp->patterns[k][partialsIndex];
			if(state < 20){
				m1 = (__m128d*)&mat[state*20];
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; m1++; up++;
				*rp = _mm_add_pd(*rp, _mm_mul_pd(weight, _mm_mul_pd(*m1, *up))); rp++; up++;
				u+=20;v+=20;
			}
			else{
				const double* transMatrixPtr = &matrices[l*400];
				for(int i = 0; i < 20; i++){
					double sum = 0;
					for(int j = 0; j < 20; j++){
						sum += transMatrixPtr[i+j*20];
					}
					rootPartials[u++] += sum * ww * upperPartials[v];
					v++;
				}
				rp += 10;
				up += 10;
			}
		}
	}
	
}
	
void calculate_branch_likelihood_20_SSE(SingleTreeLikelihood *tlk, double* rootPartials, int upperPartialsIndex, int partialsIndex, int matrixIndex){
	// partialIndex is a taxon so upperPartialsIndex has to be internal
	// matrices are transposed
	if( tlk->partials[0][partialsIndex] == NULL ){
		_calculate_branch_likelihood_state_20_SSE(tlk, rootPartials,
											   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
											   tlk->mapping[partialsIndex],
											   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	// upperPartialsIndex is a taxon
	// possible for the child of the root with a leaf sibling
	else if( tlk->partials[0][upperPartialsIndex] == NULL ){
		_calculate_branch_likelihood_upper_undefined_20_SSE(tlk, rootPartials,
														 tlk->mapping[upperPartialsIndex],
														 tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
														 tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
	else{
		_calculate_branch_likelihood_undefined_20_SSE(tlk, rootPartials,
												   tlk->partials[tlk->current_partials_indexes[upperPartialsIndex]][upperPartialsIndex],
												   tlk->partials[tlk->current_partials_indexes[partialsIndex]][partialsIndex],
												   tlk->matrices[tlk->current_matrices_indexes[matrixIndex]][matrixIndex]);
	}
}
#endif
