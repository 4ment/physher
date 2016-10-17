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
#include <xmmintrin.h> // SSE
#include <pmmintrin.h> // SSE3
//#include <tmmintrin.h> // SSSE3
#endif

#if 0
#define RESTRICT restrict
#else
#define RESTRICT
#endif

#include "treelikelihood.h"

#pragma mark -
#pragma mark Lower Likelihood

void update_partials_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_4(tlk,
                                               tlk->partials[nodeIndex1],
                                               tlk->matrices[nodeIndex1],
                                               tlk->partials[nodeIndex2],
                                               tlk->matrices[nodeIndex2],
                                               tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_4(tlk,
                                            tlk->mapping[nodeIndex2],
                                            tlk->matrices[nodeIndex2],
                                            tlk->partials[nodeIndex1],
                                            tlk->matrices[nodeIndex1],
                                            tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_4(tlk,
                                            tlk->mapping[nodeIndex1],
                                            tlk->matrices[nodeIndex1],
                                            tlk->partials[nodeIndex2],
                                            tlk->matrices[nodeIndex2],
                                            tlk->partials[nodeIndex3]);
            
        }
        else{
            partials_states_and_states_4(tlk,
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

void update_partials_uncertainty_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_4(tlk,
                                               tlk->partials[nodeIndex1],
                                               tlk->matrices[nodeIndex1],
                                               tlk->partials[nodeIndex2],
                                               tlk->matrices[nodeIndex2],
                                               tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_uncertainty_4(tlk,
                                            tlk->mapping[nodeIndex2],
                                            tlk->matrices[nodeIndex2],
                                            tlk->partials[nodeIndex1],
                                            tlk->matrices[nodeIndex1],
                                            tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_uncertainty_4(tlk,
                                            tlk->mapping[nodeIndex1],
                                            tlk->matrices[nodeIndex1],
                                            tlk->partials[nodeIndex2],
                                            tlk->matrices[nodeIndex2],
                                            tlk->partials[nodeIndex3]);
            
        }
        else{
            partials_states_and_states_uncertainty_4(tlk,
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

void update_partials_uncertainty2_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    partials_undefined_and_undefined_4(tlk,
                                       tlk->partials[nodeIndex1],
                                       tlk->matrices[nodeIndex1],
                                       tlk->partials[nodeIndex2],
                                       tlk->matrices[nodeIndex2],
                                       tlk->partials[nodeIndex3]);
    if ( tlk->scale ) {
        SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
    }
}

void update_partials_noexp_integrate_4( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_noexp_integrate_4(tlk,
                                                               tlk->partials[nodeIndex1],
                                                               tlk->matrices[nodeIndex1],
                                                               tlk->partials[nodeIndex2],
                                                               tlk->matrices[nodeIndex2],
                                                               tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_noexp_integrate_4(tlk,
                                                            tlk->mapping[nodeIndex2],
                                                            tlk->matrices[nodeIndex2],
                                                            tlk->partials[nodeIndex1],
                                                            tlk->matrices[nodeIndex1],
                                                            tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_noexp_integrate_4(tlk,
                                                            tlk->mapping[nodeIndex1],
                                                            tlk->matrices[nodeIndex1],
                                                            tlk->partials[nodeIndex2],
                                                            tlk->matrices[nodeIndex2],
                                                            tlk->partials[nodeIndex3]);
            
        }
        else{
            partials_states_and_states_noexp_integrate_4(tlk,
                                                         tlk->mapping[nodeIndex1],
                                                         tlk->matrices[nodeIndex1],
                                                         tlk->mapping[nodeIndex2],
                                                         tlk->matrices[nodeIndex2],
                                                         tlk->partials[nodeIndex3]);
        }
    }

    if( nodeIndex3 != Tree_root(tlk->tree)->id ){
        int count = tlk->sm->cat_count*tlk->sp->count;
        double temp[4];
        double *pPartials = tlk->partials[nodeIndex3];
        double **invevec = tlk->sm->m->eigendcmp->Invevec;
        
        for ( int l = 0; l < count; l++ ) {
            
            memcpy(temp, pPartials, sizeof(double)*4);
            
            for ( int i = 0; i < 4 ; i++ ) {
                *pPartials = 0;
                for ( int x = 0; x < 4 ; x++ ) {
                    *pPartials += temp[x] * invevec[i][x];
                }
                pPartials++;
            }
        }
	}
    
	if ( tlk->scale ) {
		SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
	}
}


void integrate_partials_4( const SingleTreeLikelihood * tlk, const double * RESTRICT inPartials, const double * RESTRICT proportions, double * RESTRICT outPartials ){
	int k;
	double *pPartials = outPartials;
	const double *pInPartials = inPartials;
    double prop = proportions[0];
	
    // should be faster than multiplying by 1
    if( tlk->sm->cat_count == 1 ){
        memcpy(outPartials, inPartials, tlk->sp->count*4*sizeof(double));
    }
    else {
        for ( k = 0; k < tlk->sp->count; k++ ) {
            *pPartials++ = *pInPartials++ * prop;
            *pPartials++ = *pInPartials++ * prop;
            *pPartials++ = *pInPartials++ * prop;
            *pPartials++ = *pInPartials++ * prop;
        }
        
        
        for ( int l = 1; l < tlk->sm->cat_count; l++ ) {
            pPartials = outPartials;
            prop = proportions[l];
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
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
	for ( int k = 0; k < tlk->sp->count; k++ ) {
		
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            
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
                *pPartials++ = 1.0;
                *pPartials++ = 1.0;
                *pPartials++ = 1.0;
                *pPartials++ = 1.0;
            }
        }
        u += 16;
    }
}

void partials_states_and_states_uncertainty_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT matrices2, double * RESTRICT partials ){
    int k,w;
    int u = 0;
    int state1, state2;
    double *pPartials = partials;
    double sum;
    
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            w = u;
            
            if (state1 < 4 && state2 < 4) {
                *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
                *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
                *pPartials++ = matrices1[w + state1] * matrices2[w + state2]; w += 4;
                *pPartials++ = matrices1[w + state1] * matrices2[w + state2];
                
            }
            else if (state1 < 4 && state2 < 15 ) {
                const double *partials2 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state2];
                sum  = matrices2[w]     * partials2[0];
                sum += matrices2[w + 1] * partials2[1];
                sum += matrices2[w + 2] * partials2[2];
                sum += matrices2[w + 3] * partials2[3];
                
                *pPartials++ = matrices1[w + state1] * sum;
                
                
                sum  = matrices2[w + 4] * partials2[0];
                sum += matrices2[w + 5] * partials2[1];
                sum += matrices2[w + 6] * partials2[2];
                sum += matrices2[w + 7] * partials2[3];
                
                *pPartials++ = matrices1[w + 4 + state1] * sum;
                
                
                sum  = matrices2[w + 8]  * partials2[0];
                sum += matrices2[w + 9]  * partials2[1];
                sum += matrices2[w + 10] * partials2[2];
                sum += matrices2[w + 11] * partials2[3];
                
                *pPartials++ = matrices1[w + 8 + state1] * sum;
                
                
                sum  = matrices2[w + 12] * partials2[0];
                sum += matrices2[w + 13] * partials2[1];
                sum += matrices2[w + 14] * partials2[2];
                sum += matrices2[w + 15] * partials2[3];
                
                *pPartials++ = matrices1[w + 12 + state1] * sum;
            }
            else if (state2 < 4 && state1 < 15 ) {
                const double *partials1 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state1];
                sum  = matrices1[w]     * partials1[0];
                sum += matrices1[w + 1] * partials1[1];
                sum += matrices1[w + 2] * partials1[2];
                sum += matrices1[w + 3] * partials1[3];
                
                *pPartials++ = matrices2[w + state2] * sum;
                
                
                sum  = matrices1[w + 4] * partials1[0];
                sum += matrices1[w + 5] * partials1[1];
                sum += matrices1[w + 6] * partials1[2];
                sum += matrices1[w + 7] * partials1[3];
                
                *pPartials++ = matrices2[w + 4 + state2] * sum;
                
                
                sum  = matrices1[w + 8]  * partials1[0];
                sum += matrices1[w + 9]  * partials1[1];
                sum += matrices1[w + 10] * partials1[2];
                sum += matrices1[w + 11] * partials1[3];
                
                *pPartials++ = matrices2[w + 8 + state2] * sum;
                
                
                sum  = matrices1[w + 12] * partials1[0];
                sum += matrices1[w + 13] * partials1[1];
                sum += matrices1[w + 14] * partials1[2];
                sum += matrices1[w + 15] * partials1[3];
                
                *pPartials++ = matrices2[w + 12 + state2] * sum;
                
            }
            else if (state1 < 15 && state2 < 15 ) {
                const double *partials1 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state1];
                const double *partials2 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state2];
                double sum1 = 0;
                double sum2 = 0;
                sum1  = matrices1[w]   * partials1[0];
                sum2  = matrices2[w]   * partials2[0];
                
                sum1 += matrices1[w+1] * partials1[1];
                sum2 += matrices2[w+1] * partials2[1];
                
                sum1 += matrices1[w+2] * partials1[2];
                sum2 += matrices2[w+2] * partials2[2];
                
                sum1 += matrices1[w+3] * partials1[3];
                sum2 += matrices2[w+3] * partials2[3];
                
                *pPartials++ = sum1 * sum2;
                
                
                sum1  = matrices1[w+4] * partials1[0];
                sum2  = matrices2[w+4] * partials2[0];
                
                sum1 += matrices1[w+5] * partials1[1];
                sum2 += matrices2[w+5] * partials2[1];
                
                sum1 += matrices1[w+6] * partials1[2];
                sum2 += matrices2[w+6] * partials2[2];
                
                sum1 += matrices1[w+7] * partials1[3];
                sum2 += matrices2[w+7] * partials2[3];
                
                *pPartials++ = sum1 * sum2;
                
                
                sum1  = matrices1[w+8]  * partials1[0];
                sum2  = matrices2[w+8]  * partials2[0];
                
                sum1 += matrices1[w+9]  * partials1[1];
                sum2 += matrices2[w+9]  * partials2[1];
                
                sum1 += matrices1[w+10] * partials1[2];
                sum2 += matrices2[w+10] * partials2[2];
                
                sum1 += matrices1[w+11] * partials1[3];
                sum2 += matrices2[w+11] * partials2[3];
                
                *pPartials++ = sum1 * sum2;
                
                
                sum1  = matrices1[w+12] * partials1[0];
                sum2  = matrices2[w+12] * partials2[0];
                
                sum1 += matrices1[w+13] * partials1[1];
                sum2 += matrices2[w+13] * partials2[1];
                
                sum1 += matrices1[w+14] * partials1[2];
                sum2 += matrices2[w+14] * partials2[2];
                
                sum1 += matrices1[w+15] * partials1[3];
                sum2 += matrices2[w+15] * partials2[3];
                
                *pPartials++ = sum1 * sum2;
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

void partials_states_and_states_noexp_integrate_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT exps1, int idx2, const double * RESTRICT exps2, double * RESTRICT partials ){
	int k;
	int u = 0;
	int state1, state2;
	double sum1,sum2;

    double **evec = tlk->sm->m->eigendcmp->evec;
    double **invevec = tlk->sm->m->eigendcmp->Invevec;
    double *pPartials = partials;
    
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			state2 = tlk->sp->patterns[k][idx2];
    
			if (state1 < 4 && state2 < 4) {
                
                for ( int x = 0; x < 4 ; x++ ) {
                    sum1 = sum2 = 0;
                    for ( int i = 0; i < 4 ; i++ ) {
                        sum1 += evec[x][i] * exps1[l*4+i] * invevec[i][state1];
                        sum2 += evec[x][i] * exps2[l*4+i] * invevec[i][state2];
                    }
                    *pPartials++ = sum1 * sum2;
                }
				
			}
			else if (state1 < 4 ) {
				// child 2 has a gap or unknown state so treat it as unknown
                
                for ( int x = 0; x < 4 ; x++ ) {
                    *pPartials = 0;
                    for ( int i = 0; i < 4 ; i++ ) {
                        *pPartials += evec[x][i] * exps1[l*4+i] * invevec[i][state1];
                    }
                    pPartials++;
                }
                
			}
			else if (state2 < 4 ) {
				// child 1 has a gap or unknown state so treat it as unknown
                
                for ( int x = 0; x < 4 ; x++ ) {
                    *pPartials = 0;
                    for ( int i = 0; i < 4 ; i++ ) {
                        *pPartials += evec[x][i] * exps2[l*4+i] * invevec[i][state2];
                    }
                    pPartials++;
                }
                
				
			}
			else {
                printf("partials_states_and_states_noexp_integrate_4 not done yet!\n");
                exit(1);
				// both children have a gap or unknown state so set partials to 1
				partials[u]   = 1.0;
				partials[u+1] = 1.0;
				partials[u+2] = 1.0;
				partials[u+3] = 1.0;
			}
		}
		u += 4;
	}
}



static void _partials_tip_and_tip_4_ancestral( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT matrices2, int idx3, double * RESTRICT partials ){
	int k;
	int u = 0;
	int state1, state2, state3;
	double *pPartials = partials;
    const int nstate   = tlk->sm->nstate;
	
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			state2 = tlk->sp->patterns[k][idx2];
			state3 = tlk->sp->patterns[k][idx3];
						
            int p = u + nstate * state3;
            
			if (state1 < 4 && state2 < 4) {
				*pPartials++ = matrices1[p + state1] * matrices2[p + state2];
				
			}
			else if (state1 < 4 ) {
				// child 1 has a gap or unknown state so treat it as unknown
				*pPartials++ = matrices1[p + state1];
			}
			else if (state2 < 4 ) {
				// child 2 has a gap or unknown state so treat it as unknown
				*pPartials++ = matrices2[p + state2];
				
			}
			else {
				// both children have a gap or unknown state so set partials to 1
				*pPartials++ = 1.0;
			}
		}
		u += 16;
	}
}

// Auto-vectorization
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_states_and_undefined_4( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3){
#else
void partials_states_and_undefined_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, const double * RESTRICT partials2, const double * RESTRICT matrices2, double * RESTRICT partials3){
#endif
    double sum;
	int v = 0;
	int k;
	int w = 0;
	int state1;
	

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *partials2 = __builtin_assume_aligned(apartials2, 16);
	
	const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    
	double *pPartials = partials3;
	
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			
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

    // Auto-vectorization
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    void partials_states_and_undefined_uncertainty_4( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3){
#else
    void partials_states_and_undefined_uncertainty_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, const double * RESTRICT partials2, const double * RESTRICT matrices2, double * RESTRICT partials3){
#endif
        double sum1,sum2;
        int v = 0;
        int k;
        int w = 0;
        int state1;
        
        
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
        const double *partials2 = __builtin_assume_aligned(apartials2, 16);
        
        const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
        const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
        
        double *pPartials = partials3;
        
        for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
                state1 = tlk->sp->patterns[k][idx1];
                
                if ( state1 < 4) {
                    
                    sum1  = matrices2[w]     * partials2[v];
                    sum1 += matrices2[w + 1] * partials2[v + 1];
                    sum1 += matrices2[w + 2] * partials2[v + 2];
                    sum1 += matrices2[w + 3] * partials2[v + 3];
                    
                    *pPartials++ = matrices1[w + state1] * sum1;
                    
                    
                    sum1  = matrices2[w + 4] * partials2[v];
                    sum1 += matrices2[w + 5] * partials2[v + 1];
                    sum1 += matrices2[w + 6] * partials2[v + 2];
                    sum1 += matrices2[w + 7] * partials2[v + 3];
                    
                    *pPartials++ = matrices1[w + 4 + state1] * sum1;
                    
                    
                    sum1  = matrices2[w + 8]  * partials2[v];
                    sum1 += matrices2[w + 9]  * partials2[v + 1];
                    sum1 += matrices2[w + 10] * partials2[v + 2];
                    sum1 += matrices2[w + 11] * partials2[v + 3];
                    
                    *pPartials++ = matrices1[w + 8 + state1] * sum1;
                    
                    
                    sum1  = matrices2[w + 12] * partials2[v];
                    sum1 += matrices2[w + 13] * partials2[v + 1];
                    sum1 += matrices2[w + 14] * partials2[v + 2];
                    sum1 += matrices2[w + 15] * partials2[v + 3];
                    
                    *pPartials++ = matrices1[w + 12 + state1] * sum1;
                    
                }
                else if ( state1 < 15) {
                    const double *partials1 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state1];
                    sum1  = matrices1[w]   * partials1[0];
                    sum2  = matrices2[w]   * partials2[v];
                    
                    sum1 += matrices1[w+1] * partials1[1];
                    sum2 += matrices2[w+1] * partials2[v + 1];
                    
                    sum1 += matrices1[w+2] * partials1[2];
                    sum2 += matrices2[w+2] * partials2[v + 2];
                    
                    sum1 += matrices1[w+3] * partials1[3];
                    sum2 += matrices2[w+3] * partials2[v + 3];
                    
                    *pPartials++ = sum1 * sum2;
                    
                    
                    sum1  = matrices1[w+4] * partials1[0];
                    sum2  = matrices2[w+4] * partials2[v];
                    
                    sum1 += matrices1[w+5] * partials1[1];
                    sum2 += matrices2[w+5] * partials2[v + 1];
                    
                    sum1 += matrices1[w+6] * partials1[2];
                    sum2 += matrices2[w+6] * partials2[v + 2];
                    
                    sum1 += matrices1[w+7] * partials1[3];
                    sum2 += matrices2[w+7] * partials2[v + 3];
                    
                    *pPartials++ = sum1 * sum2;
                    
                    
                    sum1  = matrices1[w+8]  * partials1[0];
                    sum2  = matrices2[w+8]  * partials2[v];
                    
                    sum1 += matrices1[w+9]  * partials1[1];
                    sum2 += matrices2[w+9]  * partials2[v + 1];
                    
                    sum1 += matrices1[w+10] * partials1[2];
                    sum2 += matrices2[w+10] * partials2[v + 2];
                    
                    sum1 += matrices1[w+11] * partials1[3];
                    sum2 += matrices2[w+11] * partials2[v + 3];
                    
                    *pPartials++ = sum1 * sum2;
                    
                    
                    sum1  = matrices1[w+12] * partials1[0];
                    sum2  = matrices2[w+12] * partials2[v];
                    
                    sum1 += matrices1[w+13] * partials1[1];
                    sum2 += matrices2[w+13] * partials2[v + 1];
                    
                    sum1 += matrices1[w+14] * partials1[2];
                    sum2 += matrices2[w+14] * partials2[v + 2];
                    
                    sum1 += matrices1[w+15] * partials1[3];
                    sum2 += matrices2[w+15] * partials2[v + 3];
                    
                    *pPartials++ = sum1 * sum2;
                    
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
        
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_states_and_undefined_noexp_integrate_4( const SingleTreeLikelihood *tlk, int idx1, const double * exps1, const double * restrict apartials2, const double * exps2, double *partials3){
#else
void partials_states_and_undefined_noexp_integrate_4( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT exps1, const double * RESTRICT partials2, const double * RESTRICT exps2, double * RESTRICT partials3){
#endif
    double sum1,sum2;
    int v = 0;
    int k;
    int state1;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
#endif
    
    double **evec = tlk->sm->m->eigendcmp->evec;
    double **invevec = tlk->sm->m->eigendcmp->Invevec;
    double *pPartials = partials3;
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            
            if ( state1 < 4) {
                
                for ( int x = 0; x < 4 ; x++ ) {
                    sum1 = 0;
                    sum2 = 0;
                    
                    for ( int i = 0; i < 4 ; i++ ) {
                        sum1 += evec[x][i] * exps1[l*4+i] * invevec[i][state1];
                        sum2 += evec[x][i] * exps2[l*4+i] * partials2[v+i];
                    }
                    *pPartials++ = sum1 * sum2;
                }
        
            }
            else {
                // Child 1 has a gap or unknown state so don't use it
                
                for ( int x = 0; x < 4 ; x++ ) {
                    *pPartials = 0;
                    for ( int i = 0; i < 4 ; i++ ) {
                        *pPartials += evec[x][i] * exps2[l*4+i] * partials2[v+i];
                    }
                    pPartials++;
                }
                
            }
            
            v += 4;
        }
    }
}

    // Auto-vectorization
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _partials_tip_and_internal_4_ancestral( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, int idx2, const double * restrict apartials2, const double * restrict amatrices2, int idx3, double *partials3){
#else
static void _partials_tip_and_internal_4_ancestral( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT partials2, const double * RESTRICT matrices2, int idx3, double * RESTRICT partials3){
#endif
    int v = 0;
    int k;
    int w = 0;
    int state1,state2,state3;
    const int nstate = tlk->sm->nstate;
    
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    
    double *pPartials = partials3;
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            state3 = tlk->sp->patterns[k][idx3];
            
            int p = w + nstate * state3;
            
            if ( state1 < 4) {
                
                *pPartials++ = matrices1[p + state1] * matrices2[p + state2] * partials2[v];
                
            }
            else {
                // Child 1 has a gap or unknown state so don't use it
                
                *pPartials++ = matrices2[p + state2] * partials2[v];
            }
            v++;
        }
        w += 16;
    }
}

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _partials_tip_and_internal_4_ancestral2( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, int idx2, const double * restrict apartials2, const double * restrict amatrices2, double *partials3){
#else
static void _partials_tip_and_internal_4_ancestral2( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT partials2, const double * RESTRICT matrices2, double * RESTRICT partials3){
#endif
    int v = 0;
    int k,i;
    int w = 0;
    int state1,state2;
    const int nstate = tlk->sm->nstate;
    
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    
    double *pPartials = partials3;
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            
            w = l * 16;
            
            if ( state1 < 4) {
                for ( i = 0; i < nstate; i++ ) {
                    *pPartials++ = matrices1[w + state1] * matrices2[w + state2] * partials2[v];
                    w += 4;
                }
            }
            else {
                // Child 1 has a gap or unknown state so don't use it
                for ( i = 0; i < nstate; i++ ) {
                    *pPartials++ = matrices2[w + state2] * partials2[v];
                    w += 4;
                }
            }
            v++;
        }
    }
}

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_undefined_and_undefined_4( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double * restrict partials3 ){
#else
void partials_undefined_and_undefined_4( const SingleTreeLikelihood *tlk, const double * RESTRICT partials1, const double * RESTRICT matrices1, const double * RESTRICT partials2, const double * RESTRICT matrices2, double * RESTRICT partials3){
#endif
    double sum1, sum2;
    
    int v = 0;
    int k;
    int w = 0;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    double *pPartials = partials3;
    int cat_count = tlk->sm->cat_count;
    int sp_count = tlk->sp->count;
    
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

// exps contains exp(lambda_i*bl). dim = cat_count*nstate
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_undefined_and_undefined_noexp_integrate_4( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict exps1, const double * restrict apartials2, const double * restrict exps2, double * restrict partials3 ){
#else
void partials_undefined_and_undefined_noexp_integrate_4( const SingleTreeLikelihood *tlk, const double * RESTRICT partials1, const double * RESTRICT exps1, const double * RESTRICT partials2, const double * RESTRICT exps2, double * RESTRICT partials3){
#endif
    double sum1, sum2;
    
    int v = 0;
    int k;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
#endif
    double *pPartials = partials3;
    double **evec = tlk->sm->m->eigendcmp->evec;
    int cat_count = tlk->sm->cat_count;
    int sp_count = tlk->sp->count;
    
    for ( int l = 0; l < cat_count; l++ ) {
        
        for ( k = 0; k < sp_count; k++ ) {
            
            for ( int x = 0; x < 4 ; x++ ) {
                sum1 = 0;
                sum2 = 0;
                
                for ( int i = 0; i < 4 ; i++ ) {
                    sum1 += evec[x][i] * exps1[l*4+i] * partials1[v+i];
                    sum2 += evec[x][i] * exps2[l*4+i] * partials2[v+i];
                }
                *pPartials++ = sum1 * sum2;
            }
            
            v += 4;
        }
    }
    
}

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _partials_internal_and_internal_4_ancestral( const SingleTreeLikelihood *tlk, int idx1, const double * restrict apartials1, const double * restrict amatrices1, int idx2, const double * restrict apartials2, const double * restrict amatrices2, int idx3, double * restrict partials3 ){
#else
static void _partials_internal_and_internal_4_ancestral( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT partials1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT partials2, const double * RESTRICT matrices2, int idx3, double * RESTRICT partials3){
#endif
        
    int v = 0;
    int k;
    int w = 0;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    double *pPartials = partials3;
    int cat_count = tlk->sm->cat_count;
    int sp_count = tlk->sp->count;
    const int nstate = tlk->sm->nstate;
    
    int state1,state2,state3;
    
    for ( int l = 0; l < cat_count; l++ ) {
        
        
        for ( k = 0; k < sp_count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            state3 = tlk->sp->patterns[k][idx3];
            
            int p = w + nstate * state3;
            
            *pPartials++ = matrices1[p + state1]   * partials1[v] * matrices2[p + state2]   * partials2[v];
            
            v++;
        }
        w += 16;
    }
    
}

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
static void _partials_internal_and_internal_4_ancestral2( const SingleTreeLikelihood *tlk, int idx1, const double * restrict apartials1, const double * restrict amatrices1, int idx2, const double * restrict apartials2, const double * restrict amatrices2, double * restrict partials3 ){
#else
static void _partials_internal_and_internal_4_ancestral2( const SingleTreeLikelihood *tlk, int idx1, const double * RESTRICT partials1, const double * RESTRICT matrices1, int idx2, const double * RESTRICT partials2, const double * RESTRICT matrices2, double * RESTRICT partials3){
#endif
    
    int v = 0;
    int k,i;
    int w = 0;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    double *pPartials = partials3;
    int cat_count = tlk->sm->cat_count;
    int sp_count = tlk->sp->count;
    const int nstate = tlk->sm->nstate;
    
    int state1,state2;
    
    for ( int l = 0; l < cat_count; l++ ) {
        
        for ( k = 0; k < sp_count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            
            w = l * 16;
            
            for ( i = 0; i < nstate; i++ ) {
                *pPartials++ = matrices1[w + state1] * partials1[v] * matrices2[w + state2] * partials2[v];
                w += 4;
            }
            v++;
        }
    }
    
}
    
void update_partials_4_ancestral( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    // internal is known
    if( tlk->mapping[nodeIndex3] != -1 ){
        // leaves are always known
        if( tlk->partials[nodeIndex1] == NULL && tlk->partials[nodeIndex2] == NULL ){
            _partials_tip_and_tip_4_ancestral(tlk, tlk->mapping[nodeIndex1], tlk->matrices[nodeIndex1], tlk->mapping[nodeIndex2], tlk->matrices[nodeIndex2], tlk->mapping[nodeIndex3], tlk->partials[nodeIndex3]);
        }
        else if( tlk->partials[nodeIndex1] == NULL ){
            _partials_tip_and_internal_4_ancestral(tlk, tlk->mapping[nodeIndex1], tlk->matrices[nodeIndex1], tlk->mapping[nodeIndex2], tlk->partials[nodeIndex2], tlk->matrices[nodeIndex2], tlk->mapping[nodeIndex3], tlk->partials[nodeIndex3]);
        }
        else if( tlk->partials[nodeIndex2] == NULL ){
            _partials_tip_and_internal_4_ancestral(tlk, tlk->mapping[nodeIndex2], tlk->matrices[nodeIndex2], tlk->mapping[nodeIndex1], tlk->partials[nodeIndex1], tlk->matrices[nodeIndex1], tlk->mapping[nodeIndex3], tlk->partials[nodeIndex3]);
        }
        else {
            _partials_internal_and_internal_4_ancestral(tlk, tlk->mapping[nodeIndex1], tlk->partials[nodeIndex1], tlk->matrices[nodeIndex1], tlk->mapping[nodeIndex2], tlk->partials[nodeIndex2], tlk->matrices[nodeIndex2], tlk->mapping[nodeIndex3], tlk->partials[nodeIndex3]);
        }
    }
    else {
        if( tlk->partials[nodeIndex1] == NULL ){
            _partials_tip_and_internal_4_ancestral2(tlk, tlk->mapping[nodeIndex1], tlk->matrices[nodeIndex1], tlk->mapping[nodeIndex2], tlk->partials[nodeIndex2], tlk->matrices[nodeIndex2], tlk->partials[nodeIndex3]);
        }
        else if( tlk->partials[nodeIndex2] == NULL ){
            _partials_tip_and_internal_4_ancestral2(tlk, tlk->mapping[nodeIndex2], tlk->matrices[nodeIndex2], tlk->mapping[nodeIndex1], tlk->partials[nodeIndex1], tlk->matrices[nodeIndex1], tlk->partials[nodeIndex3]);
        }
        else {
            _partials_internal_and_internal_4_ancestral2(tlk, tlk->mapping[nodeIndex1], tlk->partials[nodeIndex1], tlk->matrices[nodeIndex1], tlk->mapping[nodeIndex2], tlk->partials[nodeIndex2], tlk->matrices[nodeIndex2], tlk->partials[nodeIndex3]);
        }
    }
    
    if ( tlk->scale ) {
        SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
    }
}

#pragma mark -
#pragma mark OpenMP

#ifdef _OPENMP

void update_partials_4_openmp( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_4_openmp(tlk,
                                               tlk->partials[nodeIndex1],
                                               tlk->matrices[nodeIndex1],
                                               tlk->partials[nodeIndex2],
                                               tlk->matrices[nodeIndex2],
                                               tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_4_openmp(tlk,
                                            tlk->mapping[nodeIndex2],
                                            tlk->matrices[nodeIndex2],
                                            tlk->partials[nodeIndex1],
                                            tlk->matrices[nodeIndex1],
                                            tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_4_openmp(tlk,
                                            tlk->mapping[nodeIndex1],
                                            tlk->matrices[nodeIndex1],
                                            tlk->partials[nodeIndex2],
                                            tlk->matrices[nodeIndex2],
                                            tlk->partials[nodeIndex3]);
            
        }
        else{
            partials_states_and_states_4_openmp(tlk,
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
    
    
void update_partials_uncertainty_4_openmp( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    partials_undefined_and_undefined_4_openmp(tlk,
                                              tlk->partials[nodeIndex1],
                                              tlk->matrices[nodeIndex1],
                                              tlk->partials[nodeIndex2],
                                              tlk->matrices[nodeIndex2],
                                              tlk->partials[nodeIndex3]);
    
    if ( tlk->scale ) {
        SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
    }
}

void partials_states_and_states_4_openmp( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
	
    int nThreads = tlk->nthreads;
    
    #pragma omp parallel for schedule(dynamic,1) num_threads(nThreads)
    for ( int lk = 0; lk < tlk->sp->count*tlk->sm->cat_count; lk++ ) {
        int l = lk / tlk->sp->count;
        int k = lk % tlk->sp->count;
        
        int state1 = tlk->sp->patterns[k][idx1];
        int state2 = tlk->sp->patterns[k][idx2];
        
        int w = l * tlk->matrix_size;
        
        double *pPartials = partials + (l*tlk->sp->count + k)*4;
        
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
    for ( int lk = 0; lk < tlk->sp->count*tlk->sm->cat_count; lk++ ) {
        int l = lk / tlk->sp->count;
        int k = lk % tlk->sp->count;
        
        int w = l * tlk->matrix_size;
        int v = (l*tlk->sp->count + k) * 4;
        
        int state1 = tlk->sp->patterns[k][idx1];
        
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
    for ( int lk = 0; lk < tlk->sp->count*tlk->sm->cat_count; lk++ ) {
        int l = lk / tlk->sp->count;
        int k = lk % tlk->sp->count;
        
        int w = l * tlk->matrix_size;
        int v = (l*tlk->sp->count + k) * 4;
        
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
    if( tlk->sm->cat_count == 1 ){
        memcpy(outPartials, inPartials, tlk->sp->count*4*sizeof(double));
    }
    else {
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            in = _mm_load_pd(pInPartials);
            _mm_store_pd(pPartials, _mm_mul_pd(in, pn));
            pPartials   += 2;
            pInPartials += 2;
            
            in = _mm_load_pd(pInPartials);
            _mm_store_pd(pPartials, _mm_mul_pd(in, pn));
            pPartials   += 2;
            pInPartials += 2;
        }
    
	
	
        for ( int l = 1; l < tlk->sm->cat_count; l++ ) {
            prop[0] = prop[1] = proportions[l];
            pn = _mm_load_pd(prop);
            pPartials = outPartials;
            pSSE = (__m128d*)outPartials;
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
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
	
	for ( int k = 0; k < tlk->sp->count; k++ ) {
		
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

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_states_and_states_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, int idx2, const double * restrict amatrices2, double *partials ){
#else
void partials_states_and_states_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
#endif
    
	int k;
	int state1, state2;
	int u = 0;
	int w;
	
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    
	const double *m1 = matrices1;
	const double *m2 = matrices2;
    double *pPartials = partials;
    
	__m128d m1v0, m1v2,m2v0, m2v2;
	
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			state2 = tlk->sp->patterns[k][idx2];
			
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
    
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_states_and_states_uncertainty_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, int idx2, const double * restrict amatrices2, double *partials ){
#else
void partials_states_and_states_uncertainty_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, int idx2, const double *matrices2, double *partials ){
#endif
    
    int k;
    int state1, state2;
    int w=0;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    
    const double *m1 = matrices1;
    const double *m2 = matrices2;
    double *pPartials = partials;
    
    __m128d p1v0, p1v2, p2v0, p2v2, m1v0, m1v2,m2v0, m2v2;
    double temp[2] __attribute__ ((aligned (16)));
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            state2 = tlk->sp->patterns[k][idx2];
            
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
            else if (state1 < 4 && state2 < 15) {
                const double *partials2 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state2];
                m1 = &matrices1[w + state1*4];
                
                p2v0 = _mm_load_pd(&partials2[0]);
                p2v2 = _mm_load_pd(&partials2[2]);
                
                
                m2v0 = _mm_load_pd(&matrices2[w]);
                m2v2 = _mm_load_pd(&matrices2[w+2]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 * (temp[0]+temp[1]);
                m1++;
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+4]);
                m2v2 = _mm_load_pd(&matrices2[w+6]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 * (temp[0]+temp[1]);
                m1++;
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+8]);
                m2v2 = _mm_load_pd(&matrices2[w+10]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 *  (temp[0]+temp[1]);
                m1++;
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+12]);
                m2v2 = _mm_load_pd(&matrices2[w+14]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 * (temp[0]+temp[1]);
                m1++;
            }
            else if (state2 < 4 && state1 < 15 ) {
                const double *partials1 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state1];
                m2 = &matrices2[w + state2*4];
                
                p1v0 = _mm_load_pd(&partials1[0]);
                p1v2 = _mm_load_pd(&partials1[2]);
                
                
                m1v0 = _mm_load_pd(&matrices1[w]);
                m1v2 = _mm_load_pd(&matrices1[w+2]);
                
                m1v0 = _mm_mul_pd(m1v0, p1v0);
                m1v2 = _mm_mul_pd(m1v2, p1v2);
                
                m1v0 = _mm_add_pd(m1v0, m1v2);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = *m2 * (temp[0]+temp[1]);
                m2++;
                
                
                
                m1v0 = _mm_load_pd(&matrices1[w+4]);
                m1v2 = _mm_load_pd(&matrices1[w+6]);
                
                m1v0 = _mm_mul_pd(m1v0, p1v0);
                m1v2 = _mm_mul_pd(m1v2, p1v2);
                
                m1v0 = _mm_add_pd(m1v0, m1v2);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = *m2 * (temp[0]+temp[1]);
                m2++;
                
                
                
                m1v0 = _mm_load_pd(&matrices1[w+8]);
                m1v2 = _mm_load_pd(&matrices1[w+10]);
                
                m1v0 = _mm_mul_pd(m1v0, p1v0);
                m1v2 = _mm_mul_pd(m1v2, p1v2);
                
                m1v0 = _mm_add_pd(m1v0, m1v2);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = *m2 *  (temp[0]+temp[1]);
                m2++;
                
                
                
                m1v0 = _mm_load_pd(&matrices1[w+12]);
                m1v2 = _mm_load_pd(&matrices1[w+14]);
                
                m1v0 = _mm_mul_pd(m1v0, p1v0);
                m1v2 = _mm_mul_pd(m1v2, p1v2);
                
                m1v0 = _mm_add_pd(m1v0, m1v2);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = *m2 * (temp[0]+temp[1]);
                m2++;
                
            }
            else if (state1 < 15 && state2 < 15) {
                const double *partials1 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state1];
                const double *partials2 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state2];
                
                p1v0 = _mm_load_pd(&partials1[0]);
                p1v2 = _mm_load_pd(&partials1[2]);
                
                p2v0 = _mm_load_pd(&partials2[0]);
                p2v2 = _mm_load_pd(&partials2[2]);
                
                
                
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
                
                m1v0 = _mm_hadd_pd(m1v0, m2v0);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                
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
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                
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
                
                m1v0 = _mm_hadd_pd(m1v0, m2v0);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                
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
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];

            }
            else {
                // both children have a gap or unknown state so set partials to 1
                *pPartials++ = 1.0;
                *pPartials++ = 1.0;
                *pPartials++ = 1.0;
                *pPartials++ = 1.0;
            }
        }
        w += 16;
    }
}
        
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_states_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3 ){
#else
void partials_states_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
#endif
	
	int v = 0;
	int k;
	int w = 0;
	int state1;
	
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    
	double *pPartials = partials3;
	const double *m1 = matrices1;
	
	__m128d p2v0, p2v2, m2v0, m2v2;
	double temp[2] __attribute__ ((aligned (16)));
	
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
			state1 = tlk->sp->patterns[k][idx1];
			
			
			p2v0 = _mm_load_pd(&partials2[v]);
			p2v2 = _mm_load_pd(&partials2[v+2]);
			
			if ( state1 < 4) {
				
				m1 = &matrices1[w + state1*4];
				
				
				m2v0 = _mm_load_pd(&matrices2[w]);
				m2v2 = _mm_load_pd(&matrices2[w+2]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = *m1 * (temp[0]+temp[1]);
				m1++;
				
				
				
				m2v0 = _mm_load_pd(&matrices2[w+4]);
				m2v2 = _mm_load_pd(&matrices2[w+6]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = *m1 * (temp[0]+temp[1]);
				m1++;
				
				
				
				m2v0 = _mm_load_pd(&matrices2[w+8]);
				m2v2 = _mm_load_pd(&matrices2[w+10]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = *m1 *  (temp[0]+temp[1]);
				m1++;
				
				
				
				m2v0 = _mm_load_pd(&matrices2[w+12]);
				m2v2 = _mm_load_pd(&matrices2[w+14]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = *m1 * (temp[0]+temp[1]);
				m1++;
				
			}
			else {
				// Child 1 has a gap or unknown state so don't use it
				
				m2v0 = _mm_load_pd(&matrices2[w]);
				m2v2 = _mm_load_pd(&matrices2[w+2]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = temp[0]+temp[1];
				
				
				
				m2v0 = _mm_load_pd(&matrices2[w+4]);
				m2v2 = _mm_load_pd(&matrices2[w+6]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = temp[0]+temp[1];
				
				
				
				m2v0 = _mm_load_pd(&matrices2[w+8]);
				m2v2 = _mm_load_pd(&matrices2[w+10]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = temp[0]+temp[1];
				
				
				
				m2v0 = _mm_load_pd(&matrices2[w+12]);
				m2v2 = _mm_load_pd(&matrices2[w+14]);
				
				m2v0 = _mm_mul_pd(m2v0, p2v0);
				m2v2 = _mm_mul_pd(m2v2, p2v2);
				
				m2v0 = _mm_add_pd(m2v0, m2v2);
				
				_mm_store_pd(temp, m2v0);
				
				*pPartials++ = temp[0]+temp[1];
				
			}
			v += 4;
		}
		w += 16;
	}
}

#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_states_and_undefined_uncertainty_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3 ){
#else
void partials_states_and_undefined_uncertainty_4_SSE( const SingleTreeLikelihood *tlk, int idx1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
#endif
    
    int v = 0;
    int k;
    int w = 0;
    int state1;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
    
    double *pPartials = partials3;
    const double *m1 = matrices1;
    
    __m128d p1v0, p1v2, m1v0, m1v2, p2v0, p2v2, m2v0, m2v2;
    double temp[2] __attribute__ ((aligned (16)));
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state1 = tlk->sp->patterns[k][idx1];
            
            
            p2v0 = _mm_load_pd(&partials2[v]);
            p2v2 = _mm_load_pd(&partials2[v+2]);
            
            if ( state1 < 4) {
                
                m1 = &matrices1[w + state1*4];
                
                
                m2v0 = _mm_load_pd(&matrices2[w]);
                m2v2 = _mm_load_pd(&matrices2[w+2]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 * (temp[0]+temp[1]);
                m1++;
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+4]);
                m2v2 = _mm_load_pd(&matrices2[w+6]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 * (temp[0]+temp[1]);
                m1++;
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+8]);
                m2v2 = _mm_load_pd(&matrices2[w+10]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 *  (temp[0]+temp[1]);
                m1++;
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+12]);
                m2v2 = _mm_load_pd(&matrices2[w+14]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = *m1 * (temp[0]+temp[1]);
                m1++;
                
            }
            if ( state1 < 15) {
                const double *partials1 = NUCLEOTIDE_AMBIGUITY_STATES_DOUBLE[state1];
                
                p1v0 = _mm_load_pd(&partials1[0]);
                p1v2 = _mm_load_pd(&partials1[2]);
                
                
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
                
                m1v0 = _mm_hadd_pd(m1v0, m2v0);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                
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
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                
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
                
                m1v0 = _mm_hadd_pd(m1v0, m2v0);
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];
                
                
                
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
                
                _mm_store_pd(temp, m1v0);
                
                *pPartials++ = temp[0]*temp[1];
                
            }
            else {
                // Child 1 has a gap or unknown state so don't use it
                
                m2v0 = _mm_load_pd(&matrices2[w]);
                m2v2 = _mm_load_pd(&matrices2[w+2]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = temp[0]+temp[1];
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+4]);
                m2v2 = _mm_load_pd(&matrices2[w+6]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = temp[0]+temp[1];
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+8]);
                m2v2 = _mm_load_pd(&matrices2[w+10]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = temp[0]+temp[1];
                
                
                
                m2v0 = _mm_load_pd(&matrices2[w+12]);
                m2v2 = _mm_load_pd(&matrices2[w+14]);
                
                m2v0 = _mm_mul_pd(m2v0, p2v0);
                m2v2 = _mm_mul_pd(m2v2, p2v2);
                
                m2v0 = _mm_add_pd(m2v0, m2v2);
                
                _mm_store_pd(temp, m2v0);
                
                *pPartials++ = temp[0]+temp[1];
                
            }
            v += 4;
        }
        w += 16;
    }
}

            // Auto-vectorization
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
void partials_undefined_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, const double * restrict apartials1, const double * restrict amatrices1, const double * restrict apartials2, const double * restrict amatrices2, double *partials3 ){
#else
void partials_undefined_and_undefined_4_SSE( const SingleTreeLikelihood *tlk, const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 ){
#endif
	
	int v = 0;
	int k;
	int w = 0;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
    const double *partials1 = __builtin_assume_aligned(apartials1, 16);
    const double *partials2 = __builtin_assume_aligned(apartials2, 16);
    
    const double *matrices1 = __builtin_assume_aligned(amatrices1, 16);
    const double *matrices2 = __builtin_assume_aligned(amatrices2, 16);
#endif
	
	double *pPartials = partials3;
	
	__m128d m1v0, m1v2, m2v0, m2v2, p1v0, p1v2, p2v0, p2v2;
	double temp[2] __attribute__ ((aligned (16)));
	
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
            
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
			
			m1v0 = _mm_hadd_pd(m1v0, m2v0);
			
			_mm_store_pd(temp, m1v0);
			
			*pPartials++ = temp[0]*temp[1];
			
			
			
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
			
			_mm_store_pd(temp, m1v0);
			
			*pPartials++ = temp[0]*temp[1];
			
			
			
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
			
			m1v0 = _mm_hadd_pd(m1v0, m2v0);
			
			_mm_store_pd(temp, m1v0);
			
			*pPartials++ = temp[0]*temp[1];
			
			
			
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
			
			_mm_store_pd(temp, m1v0);
			
			*pPartials++ = temp[0]*temp[1];
			
			
			v += 4;
		}
		w += 16;
	}
}
    
    
    void update_partials_4_SSE( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
        
        if( tlk->mapping[nodeIndex1] == -1 ){
            if(  tlk->mapping[nodeIndex2] == -1 ){
                partials_undefined_and_undefined_4_SSE(tlk,
                                                       tlk->partials[nodeIndex1],
                                                       tlk->matrices[nodeIndex1],
                                                       tlk->partials[nodeIndex2],
                                                       tlk->matrices[nodeIndex2],
                                                       tlk->partials[nodeIndex3]);
            }
            else {
                partials_states_and_undefined_4_SSE(tlk,
                                                    tlk->mapping[nodeIndex2],
                                                    tlk->matrices[nodeIndex2],
                                                    tlk->partials[nodeIndex1],
                                                    tlk->matrices[nodeIndex1],
                                                    tlk->partials[nodeIndex3]);
            }
            
        }
        else{
            if(  tlk->mapping[nodeIndex2] == -1 ){
                partials_states_and_undefined_4_SSE(tlk,
                                                    tlk->mapping[nodeIndex1],
                                                    tlk->matrices[nodeIndex1],
                                                    tlk->partials[nodeIndex2],
                                                    tlk->matrices[nodeIndex2],
                                                    tlk->partials[nodeIndex3]);
                
            }
            else{
                partials_states_and_states_4_SSE(tlk,
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
    
    
void update_partials_uncertainty2_4_SSE( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    partials_undefined_and_undefined_4_SSE(tlk,
                                           tlk->partials[nodeIndex1],
                                           tlk->matrices[nodeIndex1],
                                           tlk->partials[nodeIndex2],
                                           tlk->matrices[nodeIndex2],
                                           tlk->partials[nodeIndex3]);
    
    if ( tlk->scale ) {
        SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
    }
}

    
void update_partials_uncertainty_4_SSE( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_4_SSE(tlk,
                                                   tlk->partials[nodeIndex1],
                                                   tlk->matrices[nodeIndex1],
                                                   tlk->partials[nodeIndex2],
                                                   tlk->matrices[nodeIndex2],
                                                   tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_uncertainty_4_SSE(tlk,
                                                tlk->mapping[nodeIndex2],
                                                tlk->matrices[nodeIndex2],
                                                tlk->partials[nodeIndex1],
                                                tlk->matrices[nodeIndex1],
                                                tlk->partials[nodeIndex3]);
//            partials_undefined_and_undefined_4_SSE(tlk,
//                                                   tlk->partials[nodeIndex1],
//                                                   tlk->matrices[nodeIndex1],
//                                                   tlk->partials[nodeIndex2],
//                                                   tlk->matrices[nodeIndex2],
//                                                   tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_uncertainty_4_SSE(tlk,
                                                tlk->mapping[nodeIndex1],
                                                tlk->matrices[nodeIndex1],
                                                tlk->partials[nodeIndex2],
                                                tlk->matrices[nodeIndex2],
                                                tlk->partials[nodeIndex3]);
//            partials_undefined_and_undefined_4_SSE(tlk,
//                                                   tlk->partials[nodeIndex1],
//                                                   tlk->matrices[nodeIndex1],
//                                                   tlk->partials[nodeIndex2],
//                                                   tlk->matrices[nodeIndex2],
//                                                   tlk->partials[nodeIndex3]);
        }
        else{
            partials_states_and_states_uncertainty_4_SSE(tlk,
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
    
    for ( k = 0; k < tlk->sp->count; k++ ) {
        in = _mm256_load_pd(pInPartials);
        _mm256_store_pd(pPartials, _mm256_mul_pd(in, pn));
        
        pPartials   += 4;
        pInPartials += 4;
    }
    
    
    
    for ( int l = 1; l < tlk->sm->cat_count; l++ ) {
        pn = _mm256_set1_pd(proportions[l]);
        
        pPartials = outPartials;
        //pSSE = (__m256d*)outPartials;
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
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
        
        for ( int k = 0; k < tlk->sp->count; k+=2 ) {
            
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
        
        if( tlk->sp->count & 1 ){
            in0 = _mm256_load_pd(pInPartials);
            
            in0 = _mm256_mul_pd(in0, f0);
            
            _mm256_store_pd(temp, in0);
            
            *pOutPartials = log(temp[0]+temp[1]+temp[2]+temp[3]);
            
            if ( tlk->scale ) {
                *pOutPartials += getLogScalingFactor( tlk, tlk->sp->count-1);
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
        
        for ( int k = 0; k < tlk->sp->count; k++ ) {
            
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
        
        for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
                state1 = tlk->sp->patterns[k][idx1];
                state2 = tlk->sp->patterns[k][idx2];
                
                w = u;
                
                if (state1 < tlk->sm->nstate && state2 < tlk->sm->nstate) {
                    
                    m1v0 = _mm256_load_pd(&matrices1[w+4*state1]);
                    
                    m2v0 = _mm256_load_pd(&matrices2[w+4*state2]);
                    
                    _mm256_store_pd(pPartials, _mm256_mul_pd(m1v0, m2v0));
                    pPartials += 4;
                    
                }
                else if (state1 < tlk->sm->nstate) {
                    // child 1 has a gap or unknown state so treat it as unknown
                    m1 = &matrices1[w+4*state1];
                    
                    *pPartials++ = *m1++;
                    *pPartials++ = *m1++;
                    *pPartials++ = *m1++;
                    *pPartials++ = *m1;
                }
                else if (state2 < tlk->sm->nstate ) {
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
        
        for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
                state1 = tlk->sp->patterns[k][idx1];
                
                
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
        
        for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
                state1 = tlk->sp->patterns[k][idx1];
                
                
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
        for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
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
        
        for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
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
        
        for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
            
            for ( k = 0; k < tlk->sp->count; k++ ) {
                
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


void update_partials_4_AVX( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {

    if( tlk->mapping[nodeIndex1] == -1 ){
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_undefined_and_undefined_4_AVX(tlk,
                                                   tlk->partials[nodeIndex1],
                                                   tlk->matrices[nodeIndex1],
                                                   tlk->partials[nodeIndex2],
                                                   tlk->matrices[nodeIndex2],
                                                   tlk->partials[nodeIndex3]);
        }
        else {
            partials_states_and_undefined_4_AVX(tlk,
                                                tlk->mapping[nodeIndex2],
                                                tlk->matrices[nodeIndex2],
                                                tlk->partials[nodeIndex1],
                                                tlk->matrices[nodeIndex1],
                                                tlk->partials[nodeIndex3]);
        }
        
    }
    else{
        if(  tlk->mapping[nodeIndex2] == -1 ){
            partials_states_and_undefined_4_AVX(tlk,
                                                tlk->mapping[nodeIndex1],
                                                tlk->matrices[nodeIndex1],
                                                tlk->partials[nodeIndex2],
                                                tlk->matrices[nodeIndex2],
                                                tlk->partials[nodeIndex3]);
            
        }
        else{
            partials_states_and_states_4_AVX(tlk,
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

    
void update_partials_uncertainty_4_AVX( SingleTreeLikelihood *tlk, int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
    
    partials_undefined_and_undefined_4_AVX(tlk,
                                           tlk->partials[nodeIndex1],
                                           tlk->matrices[nodeIndex1],
                                           tlk->partials[nodeIndex2],
                                           tlk->matrices[nodeIndex2],
                                           tlk->partials[nodeIndex3]);

    
    if ( tlk->scale ) {
        SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
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
    int w,k;
    //int w2;
    int v = 0;
    int state;
    double *pPartials = partials;
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *matrix_upper   = __builtin_assume_aligned(amatrix_upper, 16);
	const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
    const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
#endif
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
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
    
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
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
    
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( int k = 0; k < tlk->sp->count; k++ ) {

            state1 = tlk->sp->patterns[k][idx1];
            
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
    
    
void update_partials_uncertainty_upper_4( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);
    
    if( Node_isroot(parent) ){
        _update_upper_partials_root_and_undefined(tlk, tlk->partials[ Node_id(sibling) ],  tlk->matrices[ Node_id(sibling) ],  tlk->sm->m->_freqs, tlk->partials_upper[ Node_id(node) ] );
        
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
    int w,k;
    int v = 0;
    double p,sum;
    
    memset(pattern_lk, 0, tlk->sp->count*sizeof(double));
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
	const double *partials_lower = __builtin_assume_aligned(apartials_lower, 16);
	const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
#endif
    
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
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
    
    memset(pattern_lk, 0, tlk->sp->count*sizeof(double));
    
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
	const double *partials_upper = __builtin_assume_aligned(apartials_upper, 16);
	const double *matrix_lower   = __builtin_assume_aligned(amatrix_lower, 16);
#endif
    
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			state = tlk->sp->patterns[k][idx];
            
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
        _partial_lower_upper(tlk, tlk->partials_upper[node_index], tlk->partials[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
    else {
        _partial_lower_upper_leaf(tlk, tlk->partials_upper[node_index], tlk->mapping[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
}

#pragma mark -
#pragma mark Upper Likelihood SSE

#ifdef SSE3_ENABLED

// Called by a node whose parent is NOT the root and the node's sibling is a leaf
// matrix_lower is transposed
static void _update_upper_partials_state_sse2( SingleTreeLikelihood *tlk, const double *matrix_upper, const double *partials_upper, const double *matrix_lower, int sibling_index, double *partials ){
    int w,w2,k;
    int v = 0;
    int state;
    double *pPartials = partials;
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
            state = tlk->sp->patterns[k][ sibling_index ];
            
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
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
    
    for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
        
        for ( k = 0; k < tlk->sp->count; k++ ) {
            
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
    
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( k = 0; k < tlk->sp->count; k++ ) {
			
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
    
	for ( int l = 0; l < tlk->sm->cat_count; l++ ) {
		
		for ( int k = 0; k < tlk->sp->count; k++ ) {
			
            state1 = tlk->sp->patterns[k][idx1];
            
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

    
    
void update_partials_uncertainty_upper_sse_4( SingleTreeLikelihood *tlk, Node *node ){
    Node *parent = Node_parent(node);
    Node *sibling = Node_sibling(node);
    
    if( Node_isroot(parent) ){
        _update_upper_partials_root_and_undefined_sse(tlk, tlk->partials[Node_id(sibling) ],  tlk->matrices[Node_id(sibling) ],  tlk->sm->m->_freqs, tlk->partials_upper[Node_id(node) ] );
    }
    else {
        _update_upper_partials_undefined_sse(tlk, tlk->matrices[ Node_id(parent) ], tlk->partials_upper[ Node_id(parent) ], tlk->matrices[ Node_id(sibling) ], tlk->partials[ Node_id(sibling) ], tlk->partials_upper[ Node_id(node) ]);
    }
}

static void _partial_lower_upper_sse( const SingleTreeLikelihood *tlk, const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ){
    int w,k;
    int v = 0;
    double p;
    int cat_count = tlk->sm->cat_count;
    int sp_count = tlk->sp->count;
    
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

void node_log_likelihoods_upper_sse_4( const SingleTreeLikelihood *tlk, Node *node ){
	int node_index = Node_id(node);
    
    if ( !Node_isleaf(node) ) {
        _partial_lower_upper_sse(tlk, tlk->partials_upper[node_index], tlk->partials[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
    else {
        _partial_lower_upper_leaf_sse(tlk, tlk->partials_upper[node_index], tlk->mapping[node_index], tlk->matrices[node_index], tlk->sm->get_proportions(tlk->sm), tlk->node_pattern_lk );
    }
}
#endif
