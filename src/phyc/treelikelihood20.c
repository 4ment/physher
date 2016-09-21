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
void update_partials_20_SSE( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_20_SSE(tlk,
                                                   tlk->partials[nodeIndex1],
                                                   tlk->matrices[nodeIndex1],
                                                   tlk->partials[nodeIndex2],
                                                   tlk->matrices[nodeIndex2],
                                                   tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_20_SSE(tlk,
                                                tlk->mapping[nodeIndex2],
                                                tlk->matrices[nodeIndex2],
                                                tlk->partials[nodeIndex1],
                                                tlk->matrices[nodeIndex1],
                                                tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_20_SSE(tlk,
                                                tlk->mapping[nodeIndex1],
                                                tlk->matrices[nodeIndex1],
                                                tlk->partials[nodeIndex2],
                                                tlk->matrices[nodeIndex2],
                                                tlk->partials[nodeIndex3]);
            
        }
        else{
            partials_states_and_states_20_SSE(tlk,
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
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

// Called by a node whose parent is NOT the root and the node's sibling is a leaf
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _update_upper_partials_state( SingleTreeLikelihood *tlk, const double * restrict amatrix_upper, const double * restrict apartials_upper, const double * restrict amatrix_lower, int sibling_index, double *partials ){
#else
static void _update_upper_partials_state( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, int sibling_index, double *partials ){
#endif
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *matrix_upper   = __builtin_assume_aligned(amatrix_upper,   16);
    const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
    const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower,   16);
#endif
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
            w = l * 400;
            
            for ( int j = 0; j < 5; j++ ) {
            
                w2 = l * 400 + j*4;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+state];
                }
                pPartials++;
                w += 20;
                
                
                w2 = l * 400 + j*4 + 1;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+state];
                }
                pPartials++;
                w += 20;
                
                
                w2 = l * 400 + j*4 + 2;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+state];
                }
                pPartials++;
                w += 20;
                
                
                w2 = l * 400 + j*4 + 3;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+state];
                }
                pPartials++;
                w += 20;
            }
            
            v += 20;
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
    const double *matrix_upper   = __builtin_assume_aligned(amatrix_upper,   16);
    const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
    const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower,   16);
    const double *partials_lower = __builtin_assume_aligned(apartials_lower, 16);
#endif
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            w = l*400;
            
            for( int j = 0; j < 5; j++ ){
                
                w2 = l*400 + j*4;
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                sum2  = matrix_lower[w] * partials_lower[v];      w++;
                sum2 += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                *pPartials++ = sum1 * sum2 ;
                
                
                w2 = l*400 + j*4 + 1;
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                sum2  = matrix_lower[w] * partials_lower[v];      w++;
                sum2 += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                
                *pPartials++ = sum1 * sum2 ;
                
                
                w2 = l*400 + j*4 + 2;
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                sum2  = matrix_lower[w] * partials_lower[v];      w++;
                sum2 += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                *pPartials++ = sum1 * sum2 ;
                
                
                w2 = l*400 + j*4 + 3;
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                sum2  = matrix_lower[w] * partials_lower[v];      w++;
                sum2 += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum2 += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum2 += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                *pPartials++ = sum1 * sum2;
                
            }
            v += 20;
        }
    }
}
        
        // Called by a child of the root but not if the child's sibling is NOT leaf
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _update_upper_partials_root_and_undefined( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, const double * restrict frequencies, double *partials_upper ){
#else
static void _update_upper_partials_root_and_undefined( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper ){
#endif
    int w,k,j;
    int v = 0;
    double *pPartials = partials_upper;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
#endif
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            w = l * 400;
            
            for( j = 0; j < 5; j++ ){
                
                *pPartials  = matrices1[w] * partials1[v];      w++;
                *pPartials += matrices1[w] * partials1[v + 1];  w++;
                *pPartials += matrices1[w] * partials1[v + 2];  w++;
                *pPartials += matrices1[w] * partials1[v + 3];  w++;
                *pPartials += matrices1[w] * partials1[v + 4];  w++;
                *pPartials += matrices1[w] * partials1[v + 5];  w++;
                *pPartials += matrices1[w] * partials1[v + 6];  w++;
                *pPartials += matrices1[w] * partials1[v + 7];  w++;
                *pPartials += matrices1[w] * partials1[v + 8];  w++;
                *pPartials += matrices1[w] * partials1[v + 9];  w++;
                *pPartials += matrices1[w] * partials1[v + 10]; w++;
                *pPartials += matrices1[w] * partials1[v + 11]; w++;
                *pPartials += matrices1[w] * partials1[v + 12]; w++;
                *pPartials += matrices1[w] * partials1[v + 13]; w++;
                *pPartials += matrices1[w] * partials1[v + 14]; w++;
                *pPartials += matrices1[w] * partials1[v + 15]; w++;
                *pPartials += matrices1[w] * partials1[v + 16]; w++;
                *pPartials += matrices1[w] * partials1[v + 17]; w++;
                *pPartials += matrices1[w] * partials1[v + 18]; w++;
                *pPartials += matrices1[w] * partials1[v + 19]; w++;
                
                *pPartials++ *= frequencies[j*4];
                
                
                *pPartials  = matrices1[w] * partials1[v];      w++;
                *pPartials += matrices1[w] * partials1[v + 1];  w++;
                *pPartials += matrices1[w] * partials1[v + 2];  w++;
                *pPartials += matrices1[w] * partials1[v + 3];  w++;
                *pPartials += matrices1[w] * partials1[v + 4];  w++;
                *pPartials += matrices1[w] * partials1[v + 5];  w++;
                *pPartials += matrices1[w] * partials1[v + 6];  w++;
                *pPartials += matrices1[w] * partials1[v + 7];  w++;
                *pPartials += matrices1[w] * partials1[v + 8];  w++;
                *pPartials += matrices1[w] * partials1[v + 9];  w++;
                *pPartials += matrices1[w] * partials1[v + 10]; w++;
                *pPartials += matrices1[w] * partials1[v + 11]; w++;
                *pPartials += matrices1[w] * partials1[v + 12]; w++;
                *pPartials += matrices1[w] * partials1[v + 13]; w++;
                *pPartials += matrices1[w] * partials1[v + 14]; w++;
                *pPartials += matrices1[w] * partials1[v + 15]; w++;
                *pPartials += matrices1[w] * partials1[v + 16]; w++;
                *pPartials += matrices1[w] * partials1[v + 17]; w++;
                *pPartials += matrices1[w] * partials1[v + 18]; w++;
                *pPartials += matrices1[w] * partials1[v + 19]; w++;
                
                *pPartials++ *= frequencies[j*4+1];
                
                
                *pPartials  = matrices1[w] * partials1[v];      w++;
                *pPartials += matrices1[w] * partials1[v + 1];  w++;
                *pPartials += matrices1[w] * partials1[v + 2];  w++;
                *pPartials += matrices1[w] * partials1[v + 3];  w++;
                *pPartials += matrices1[w] * partials1[v + 4];  w++;
                *pPartials += matrices1[w] * partials1[v + 5];  w++;
                *pPartials += matrices1[w] * partials1[v + 6];  w++;
                *pPartials += matrices1[w] * partials1[v + 7];  w++;
                *pPartials += matrices1[w] * partials1[v + 8];  w++;
                *pPartials += matrices1[w] * partials1[v + 9];  w++;
                *pPartials += matrices1[w] * partials1[v + 10]; w++;
                *pPartials += matrices1[w] * partials1[v + 11]; w++;
                *pPartials += matrices1[w] * partials1[v + 12]; w++;
                *pPartials += matrices1[w] * partials1[v + 13]; w++;
                *pPartials += matrices1[w] * partials1[v + 14]; w++;
                *pPartials += matrices1[w] * partials1[v + 15]; w++;
                *pPartials += matrices1[w] * partials1[v + 16]; w++;
                *pPartials += matrices1[w] * partials1[v + 17]; w++;
                *pPartials += matrices1[w] * partials1[v + 18]; w++;
                *pPartials += matrices1[w] * partials1[v + 19]; w++;
                
                *pPartials++ *= frequencies[j*4+2];
                
                
                *pPartials  = matrices1[w] * partials1[v];      w++;
                *pPartials += matrices1[w] * partials1[v + 1];  w++;
                *pPartials += matrices1[w] * partials1[v + 2];  w++;
                *pPartials += matrices1[w] * partials1[v + 3];  w++;
                *pPartials += matrices1[w] * partials1[v + 4];  w++;
                *pPartials += matrices1[w] * partials1[v + 5];  w++;
                *pPartials += matrices1[w] * partials1[v + 6];  w++;
                *pPartials += matrices1[w] * partials1[v + 7];  w++;
                *pPartials += matrices1[w] * partials1[v + 8];  w++;
                *pPartials += matrices1[w] * partials1[v + 9];  w++;
                *pPartials += matrices1[w] * partials1[v + 10]; w++;
                *pPartials += matrices1[w] * partials1[v + 11]; w++;
                *pPartials += matrices1[w] * partials1[v + 12]; w++;
                *pPartials += matrices1[w] * partials1[v + 13]; w++;
                *pPartials += matrices1[w] * partials1[v + 14]; w++;
                *pPartials += matrices1[w] * partials1[v + 15]; w++;
                *pPartials += matrices1[w] * partials1[v + 16]; w++;
                *pPartials += matrices1[w] * partials1[v + 17]; w++;
                *pPartials += matrices1[w] * partials1[v + 18]; w++;
                *pPartials += matrices1[w] * partials1[v + 19]; w++;
                
                *pPartials++ *= frequencies[j*4+3];
                
            }
            v += 20;
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( int k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            
            w = l * 400;
            
            if( state1 < 20 ){
                *pPartials++ = matrices1[w+state1] * frequencies[0];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[1];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[2];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[3];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[4];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[5];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[6];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[7];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[8];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[9];  w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[10]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[11]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[12]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[13]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[14]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[15]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[16]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[17]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[18]; w+=20;
                *pPartials++ = matrices1[w+state1] * frequencies[19];
            }
            else {
                memcpy(pPartials, frequencies, sizeof(double)*20);
                pPartials += 20;
                
            }
            
        }
    }
    
}
    
void update_partials_upper_20( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);
    
    if( Node_isroot(parent) ){
        if( Node_isleaf(sibling) ){
            _update_upper_partials_root_and_state(tlk, tlk->matrices[ Node_id(sibling) ], tlk->mapping[ Node_id(sibling) ], tlk->sm->m->_freqs, tlk->partials_upper[ Node_id(node) ] );
        }
        else {
            _update_upper_partials_root_and_undefined(tlk, tlk->partials[ Node_id(sibling) ],  tlk->matrices[ Node_id(sibling) ],  tlk->sm->m->_freqs, tlk->partials_upper[ Node_id(node) ] );
        }
    }
    else if( Node_isleaf(sibling) ){
        _update_upper_partials_state(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials_upper[ Node_id(parent) ], tlk->matrices[ Node_id(sibling) ], tlk->mapping[ Node_id(sibling) ], tlk->partials_upper[ Node_id(node) ]);
    }
    else {
        _update_upper_partials_undefined(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials_upper[ Node_id(parent) ], tlk->matrices[ Node_id(sibling) ], tlk->partials[ Node_id(sibling) ], tlk->partials_upper[ Node_id(node) ]);
    }
}
                
                
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _partial_lower_upper( const SingleTreeLikelihood *tlk, const double * restrict apartials_upper, const double * restrict apartials_lower, const double * restrict amatrix_lower, const double * restrict proportions, double *pattern_lk ){
#else
static void _partial_lower_upper( const SingleTreeLikelihood *tlk, const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ){
#endif
    int w,k,j,l;
    int v = 0;
    double p,sum;
    
    memset(pattern_lk, 0, tlk->sp->count*sizeof(double));
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
    const double *partials_lower = __builtin_assume_aligned(apartials_lower, 16);
    const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
#endif
    
    for ( l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            w = l * 400;
            p = 0;
            
            for( j = 0; j < 5; j++ ){
                
                sum  = matrix_lower[w] * partials_lower[v];      w++;
                sum += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                p += sum * partials_upper[v + j*4];
                
                
                sum  = matrix_lower[w] * partials_lower[v];      w++;
                sum += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                p += sum * partials_upper[v + j*4 + 1];
                
                
                sum  = matrix_lower[w] * partials_lower[v];      w++;
                sum += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                p += sum * partials_upper[v + j*4 + 2];
                
                
                sum  = matrix_lower[w] * partials_lower[v];      w++;
                sum += matrix_lower[w] * partials_lower[v + 1];  w++;
                sum += matrix_lower[w] * partials_lower[v + 2];  w++;
                sum += matrix_lower[w] * partials_lower[v + 3];  w++;
                sum += matrix_lower[w] * partials_lower[v + 4];  w++;
                sum += matrix_lower[w] * partials_lower[v + 5];  w++;
                sum += matrix_lower[w] * partials_lower[v + 6];  w++;
                sum += matrix_lower[w] * partials_lower[v + 7];  w++;
                sum += matrix_lower[w] * partials_lower[v + 8];  w++;
                sum += matrix_lower[w] * partials_lower[v + 9];  w++;
                sum += matrix_lower[w] * partials_lower[v + 10]; w++;
                sum += matrix_lower[w] * partials_lower[v + 11]; w++;
                sum += matrix_lower[w] * partials_lower[v + 12]; w++;
                sum += matrix_lower[w] * partials_lower[v + 13]; w++;
                sum += matrix_lower[w] * partials_lower[v + 14]; w++;
                sum += matrix_lower[w] * partials_lower[v + 15]; w++;
                sum += matrix_lower[w] * partials_lower[v + 16]; w++;
                sum += matrix_lower[w] * partials_lower[v + 17]; w++;
                sum += matrix_lower[w] * partials_lower[v + 18]; w++;
                sum += matrix_lower[w] * partials_lower[v + 19]; w++;
                
                p += sum * partials_upper[v + j*4 + 3];
            }
            
            v += 20;
            
            pattern_lk[k] += p * proportions[l];
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
    
    memset(pattern_lk, 0, tlk->sp->count*sizeof(double));
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
    const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower,   16);
#endif
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            state = tlk->sp->patterns[k][idx];
            
            w = l * 400;
            if( state < 20 ){
                p  = matrix_lower[w+state] * partials_upper[v];      w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 1];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 2];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 3];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 4];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 5];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 6];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 7];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 8];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 9];  w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 10]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 11]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 12]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 13]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 14]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 15]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 16]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 17]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 18]; w += 20;
                p += matrix_lower[w+state] * partials_upper[v + 19];
            }
            else {
                p  = partials_upper[v];
                p += partials_upper[v + 1];
                p += partials_upper[v + 2];
                p += partials_upper[v + 3];
                p += partials_upper[v + 4];
                p += partials_upper[v + 5];
                p += partials_upper[v + 6];
                p += partials_upper[v + 7];
                p += partials_upper[v + 8];
                p += partials_upper[v + 9];
                p += partials_upper[v + 10];
                p += partials_upper[v + 11];
                p += partials_upper[v + 12];
                p += partials_upper[v + 13];
                p += partials_upper[v + 14];
                p += partials_upper[v + 15];
                p += partials_upper[v + 16];
                p += partials_upper[v + 17];
                p += partials_upper[v + 18];
                p += partials_upper[v + 19];
            }
            pattern_lk[k] += p * proportions[l];
            v += 20;
        }
    }
}

void node_log_likelihoods_upper_20( const SingleTreeLikelihood *tlk, Node *node ){
    int node_index = Node_id(node);
    
    if ( !Node_isleaf(node) ) {
        _partial_lower_upper(tlk, tlk->partials_upper[node_index], tlk->partials[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
    else {
        _partial_lower_upper_leaf(tlk, tlk->partials_upper[node_index], tlk->mapping[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
}
                    


#pragma mark -
#pragma mark Upper Likelihood SSE

//#define SSE3_ENABLED 1

#ifdef SSE3_ENABLED
                                                                    
// Called by a node whose parent is NOT the root and the node's sibling is a leaf
// matrix_lower is transposed
static void _update_upper_partials_state_sse( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, int sibling_index, double *partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
            w = l * 400;
            for( int j = 0; j < 5; j++ ){
                
                w2 = l * 400 + j*4;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+20*state];
                }
                pPartials++;
                w++;
                
                
                w2 = l * 400 + j*4 + 1;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+20*state];
                }
                pPartials++;
                w++;
                
                
                w2 = l * 400 + j*4 + 2;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+20*state];
                }
                pPartials++;
                w++;
                
                
                w2 = l * 400 + j*4 + 3;
                
                *pPartials  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                *pPartials += matrix_upper[w2] * partials_upper[v + 19];
                
                
                if( state < 20){
                    *pPartials *= matrix_lower[w+20*state];
                }
                pPartials++;
                w++;
             
            }
            v += 20;
        }
    }
}

    
// Called by a node whose parent is NOT the root and the node's sibling is NOT a leaf
// All matrices are NOT transposed
static void _update_upper_partials_undefined_sse( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ){
    int w,w2,k;
    int v = 0;
    double sum1;
    double *pPartials = partials;
    
    __m128d *pl;
    __m128d *ml;
    __m128d temp;
    double t[2] __attribute__ ((aligned (16)));
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            w = l*400;
            
            ml = (__m128d*)&matrix_lower[w];
            
            for( int j = 0; j < 5; j++ ){
                
                w2 = l*400 + j*4;
                
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                pl = (__m128d*)&partials_lower[v];
                temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
                _mm_store_pd(t,temp);
                
                *pPartials++ = sum1 * (t[0]+t[1]) ;
                
                
                w2 = l*400 + j*4 + 1;
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                pl = (__m128d*)&partials_lower[v];
                temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
                _mm_store_pd(t,temp);
                
                *pPartials++ = sum1 * (t[0]+t[1]) ;
                
                
                w2 = l*400 + j*4 + 2;
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                pl = (__m128d*)&partials_lower[v];
                temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;
                
                _mm_store_pd(t,temp);
                
                *pPartials++ = sum1 * (t[0]+t[1]) ;
                
                
                w2 = l*400 + j*4 + 3;
                
                sum1  = matrix_upper[w2] * partials_upper[v];      w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 1];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 2];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 3];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 4];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 5];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 6];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 7];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 8];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 9];  w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 10]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 11]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 12]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 13]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 14]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 15]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 16]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 17]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 18]; w2 += 20;
                sum1 += matrix_upper[w2] * partials_upper[v + 19];
                
                pl = (__m128d*)&partials_lower[v];
                temp = _mm_mul_pd(*ml,*pl); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++; pl++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*ml,*pl)); ml++;

                _mm_store_pd(t,temp);
                
                *pPartials++ = sum1 * (t[0]+t[1]) ;
                
            }
            v += 20;
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            w = l * 400;
            
            m = (__m128d*)&matrices1[w];
            
            for( int j = 0; j < 5; j++ ){
                
                p = (__m128d*)&partials1[v];
                
                temp = _mm_mul_pd(*p, *m); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); m++;
                
                _mm_store_pd(t, temp);
                
                *pPartials++ = frequencies[j*4] * (t[0]+t[1]);
                
                
                p = (__m128d*)&partials1[v];
                
                temp = _mm_mul_pd(*p, *m); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); m++;
                
                _mm_store_pd(t, temp);
                
                *pPartials++ = frequencies[j*4+1] * (t[0]+t[1]);
                
                
                p = (__m128d*)&partials1[v];
                
                temp = _mm_mul_pd(*p, *m); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); m++;
                
                _mm_store_pd(t, temp);
                
                *pPartials++ = frequencies[j*4+2] * (t[0]+t[1]);
                
                
                p = (__m128d*)&partials1[v];
                
                temp = _mm_mul_pd(*p, *m); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); p++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*p, *m)); m++;
                
                _mm_store_pd(t, temp);
                
                *pPartials++ = frequencies[j*4+3] * (t[0]+t[1]);
                
            }
            v += 20;
        }
    }
}

// Called by a child of the root and the child's sibling is a leaf
// matrices1 is transposed
static void _update_upper_partials_root_and_state_sse( const SingleTreeLikelihood *tlk, const double *matrices1, int idx1, const double *frequencies, double *partials_upper ){
    int l,k,w;
    double *pPartials = partials_upper;
    int state1;
    __m128d *m;
    __m128d *f;
    
    for ( l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            
            w = l * 400;
            
            if( state1 < 20 ){
                m = (__m128d*)&matrices1[w+20*state1];
                f = (__m128d*)frequencies;
                
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2; m++; f++;
                _mm_store_pd(pPartials, _mm_mul_pd(*m,*f)); pPartials += 2;
            }
            else {
                memcpy(pPartials, frequencies, sizeof(double)*20);
                pPartials += 20;
            }
            
        }
    }
}

void update_partials_upper_sse_20( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);

    if( Node_isroot(parent) ){
        // The matrix of the sibling is transposed
        if( Node_isleaf(sibling) ){
            _update_upper_partials_root_and_state_sse(tlk, tlk->matrices[Node_id(sibling) ], tlk->mapping[Node_id(sibling) ], tlk->sm->m->_freqs, tlk->partials_upper[Node_id(node) ] );
        }
        else {
            _update_upper_partials_root_and_undefined_sse(tlk, tlk->partials[Node_id(sibling) ],  tlk->matrices[Node_id(sibling) ],  tlk->sm->m->_freqs, tlk->partials_upper[Node_id(node) ] );
        }
    }
    // The matrix of the sibling is transposed
    // The pparent node cannot be leaf
    else if( Node_isleaf(sibling) ){
        _update_upper_partials_state_sse(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials_upper[ Node_id(parent) ], tlk->matrices[ Node_id(sibling) ], tlk->mapping[ Node_id(sibling) ], tlk->partials_upper[ Node_id(node) ]);
    }
    else {
        _update_upper_partials_undefined_sse(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials_upper[ Node_id(parent) ], tlk->matrices[ Node_id(sibling) ], tlk->partials[ Node_id(sibling) ], tlk->partials_upper[ Node_id(node) ]);
    }
}

static void _partial_lower_upper_sse( const SingleTreeLikelihood *tlk, const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int j,k;
    int v = 0;
    double p;
    int cat_count = tlk->sm->cat_count;
    int sp_count = tlk->sp->count;
    
    __m128d *m, *pl, temp;
    
    double t[2] __attribute__ ((aligned (16)));
    
    memset(pattern_lk, 0, sp_count*sizeof(double));
    
    for ( int l = 0; l < cat_count; l++ ) {
        
        for ( k = 0; k < sp_count; k++ ) {
            
            
            p = 0;
            
            m = (__m128d*)&matrix_lower[l*400];
            
            for( j = 0; j < 5; j++ ){
                
                pl = (__m128d*)&partials_lower[v];
                
                temp = _mm_mul_pd(*pl, *m); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); m++;
                
                _mm_store_pd(t, temp);
                p += (t[0]+t[1]) * partials_upper[v + j*4];
                
                pl = (__m128d*)&partials_lower[v];
                
                temp = _mm_mul_pd(*pl, *m); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); m++;
                
                _mm_store_pd(t, temp);
                p += (t[0]+t[1]) * partials_upper[v + j*4 + 1];
                
                pl = (__m128d*)&partials_lower[v];
                
                temp = _mm_mul_pd(*pl, *m); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); m++;
                
                _mm_store_pd(t, temp);
                p += (t[0]+t[1]) * partials_upper[v + j*4 + 2];
                
                pl = (__m128d*)&partials_lower[v];
                
                temp = _mm_mul_pd(*pl, *m); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); pl++; m++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*pl, *m)); m++;
                
                _mm_store_pd(t, temp);
                p += (t[0]+t[1]) * partials_upper[v + j*4 + 3];
                
            }
            v += 20;
            
            pattern_lk[k] += p * proportions[l];
        }
    }
}

    
// matrix_lower is transposed
static void _partial_lower_upper_leaf_sse( const SingleTreeLikelihood *tlk, const double *partials_upper, int idx, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int state;
    double p;
    
    int cat_count = tlk->sm->cat_count;
    int sp_count = tlk->sp->count;
    
    __m128d *m;
    __m128d temp;
    __m128d *pu = (__m128d*)partials_upper;
    
    double t[2] __attribute__ ((aligned (16)));
    
    memset(pattern_lk, 0, sp_count*sizeof(double));
    
    for ( int l = 0; l < cat_count; l++ ) {
        p = proportions[l];
        
        for ( k = 0; k < sp_count; k++ ) {
            state = tlk->sp->patterns[k][idx];
            
            w = l * 400;
            
            if( state < 20 ){
                m = (__m128d*)&matrix_lower[w+20*state];
                
                temp = _mm_mul_pd(*m,*pu); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); m++; pu++;
                temp = _mm_add_pd(temp, _mm_mul_pd(*m,*pu)); pu++;
            }
            else {
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
            }
            
            _mm_store_pd(t, temp);
            
            pattern_lk[k] += (t[0]+t[1]) * p;
        }
    }
}

void node_log_likelihoods_upper_sse_20( const SingleTreeLikelihood *tlk, Node *node ){
    int node_index = Node_id(node);
    
    if ( !Node_isleaf(node) ) {
        _partial_lower_upper_sse(tlk, tlk->partials_upper[node_index], tlk->partials[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
    else {
        _partial_lower_upper_leaf_sse(tlk, tlk->partials_upper[node_index], tlk->mapping[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
}
#endif
