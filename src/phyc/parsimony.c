/*
 *  parsimony.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 9/04/2014.
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

#include "parsimony.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


#include "tree.h"
#include "sitepattern.h"
#include "matrix.h"
#include "utils.h"


#include "datatype.h"

#ifdef __SSE2__
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2

#ifdef __SSE4_1__
#include <smmintrin.h> //SSE4.1
#endif


static double _score_fitch_4_sse(Parsimony *parsimony);
static double _score_fitch_4_sse_slow(Parsimony *parsimony);
#endif

static double _score_fitch(Parsimony *parsimony);
static double _score_fitch_slow(Parsimony *parsimony);

static void _update(Parsimony *parsimony);

static void _reconstruct( Parsimony *parsimony );

#ifdef __SSE2__
static inline __m128i muly(const __m128i a, const __m128i b){
#ifdef __SSE4_1__  // modern CPU - use SSE 4.1
    return _mm_mullo_epi32(a, b);
#else               // old CPU - use SSE 2
    __m128i tmp1 = _mm_mul_epu32(a,b); /* mul 2,0*/
    __m128i tmp2 = _mm_mul_epu32( _mm_srli_si128(a,4), _mm_srli_si128(b,4)); /* mul 3,1 */
    return _mm_unpacklo_epi32(_mm_shuffle_epi32(tmp1, _MM_SHUFFLE (0,0,2,0)), _mm_shuffle_epi32(tmp2, _MM_SHUFFLE (0,0,2,0))); /* shuffle results to [63..0] and pack */
#endif
}
#endif

Parsimony * new_Parsimony( SitePattern *sp, Tree *tree ){
    int i,len;
    
    Parsimony *parsimony = (Parsimony*)malloc(sizeof(Parsimony));
    assert(parsimony);
    parsimony->sp = sp;
    parsimony->tree = tree;
    
    len = Tree_node_count(tree)-Tree_tip_count(tree);
    
    parsimony->stateSets = (int8_t**)malloc(Tree_node_count(tree)*sizeof(int8_t*));
    assert(parsimony->stateSets);
    
    parsimony->reconstruct = _reconstruct;
    parsimony->states = NULL;
    
#ifdef INT_WEIGHT
    parsimony->local_scores = (int32_t**)malloc(sizeof(int32_t*)*len);
    assert(parsimony->local_scores);
    parsimony->weights = aligned16_malloc( sizeof(int32_t) * sp->count );
    assert(parsimony->weights);
    
    for ( i = 0; i < sp->count; i++ ) {
        parsimony->weights[i] = (int)sp->weights[i];
    }
#else
    parsimony->local_scores = (double**)malloc(len*sizeof(double*));
    assert(Parsimony->local_scores);
#endif
    
#ifdef __SSE2__
    for ( i = 0; i < Tree_node_count(tree); i++ ) {
        parsimony->stateSets[i] = aligned16_malloc( parsimony->sp->nstate * sp->count * sizeof(int8_t) );
    }
    
#ifdef INT_WEIGHT
    for ( i = 0; i < len; i++ ) {
        parsimony->local_scores[i] = aligned16_malloc( sp->count * sizeof(int32_t) );
        assert(parsimony->local_scores[i]);
    }
    parsimony->scores = aligned16_malloc( sp->count * sizeof(int32_t) );
#else
    for ( i = 0; i < len; i++ ) {
        parsimony->local_scores[i] = aligned16_malloc( sp->count * sizeof(double) );
        assert(parsimony->local_scores[i]);
    }
    parsimony->scores = aligned16_malloc( sp->count * sizeof(double) );
    assert(parsimony->scores);
#endif
    parsimony->calculate = _score_fitch_4_sse;

#else // else __SSE2__
    
    for ( i = 0; i < Tree_node_count(tree); i++ ) {
        parsimony->stateSets[i]    = (int8_t*)calloc( parsimony->sp->nstate*sp->count, sizeof(int8_t));
        assert(parsimony->stateSets[i]);
    }
#ifdef INT_WEIGHT
    for ( i = 0; i < len; i++ ) {
        parsimony->local_scores[i] = (int32_t*)calloc( sp->count,sizeof(int32_t));
        assert(parsimony->local_scores[i]);
    }
    parsimony->scores = (int32_t*)calloc(sp->count, sizeof(int32_t));
    assert(parsimony->scores);
#else
    for ( i = 0; i < len; i++ ) {
        parsimony->local_scores[i] = dvector(sp->count);
    }
    parsimony->scores = dvector(sp->count);
#endif
    parsimony->calculate = _score_fitch;
#endif // endif __SSE2__
    
    parsimony->update_nodes = bvector(Tree_node_count(tree));
    for ( i = 0; i < Tree_node_count(tree); i++ ) {
        parsimony->update_nodes[i] = true;
    }
    
    parsimony->update  = true;
    parsimony->score = 0;
    
    if( parsimony->sp->nstate != 4 ) parsimony->calculate = _score_fitch;
    
    _update(parsimony);
    
    return parsimony;
}

void free_Parsimony( Parsimony *parsimony ){
    int i,len;
    len = Tree_node_count(parsimony->tree)-Tree_tip_count(parsimony->tree);
    free(parsimony->scores);

    for ( i = 0; i < Tree_node_count(parsimony->tree); i++ ) {
        free(parsimony->stateSets[i]);
    }
    for ( i = 0; i < len; i++ ) {
        free(parsimony->local_scores[i]);
    }
    free(parsimony->stateSets);
    free(parsimony->local_scores);
    
#ifdef INT_WEIGHT
    free(parsimony->weights);
#endif
    
    free(parsimony->update_nodes);
    
    if( parsimony->states != NULL ){
        free_ui8matrix(parsimony->states, parsimony->sp->count);
    }
    
    free(parsimony);
}

void Parsimony_update_node( Parsimony *parsimony, Node *node ){
    parsimony->update_nodes[Node_id(node)] = true;
    parsimony->update = true;
}


void Parsimony_update_all_nodes( Parsimony *parsimony ){
	for (int index = 0; index < Tree_node_count(parsimony->tree); index++) {
		parsimony->update_nodes[index] = true;
	}
	parsimony->update = true;
}


void _update( Parsimony *parsimony ){
    Node **nodes = Tree_get_nodes(parsimony->tree, POSTORDER);
    int sp_count = parsimony->sp->count;
    int nstate = parsimony->sp->nstate;
    uint8_t **patterns = parsimony->sp->patterns;
    
    for ( int i = 0; i < Tree_node_count(parsimony->tree); i++ ) {
        
        if( Node_isleaf(nodes[i]) ){
            int idx = get_sequence_index(parsimony->sp, nodes[i]->name);
//            if( parsimony->sp->datatype->type == DATA_TYPE_NUCLEOTIDE ){
//                for ( int j = 0; j < sp_count; j++ ) {
//                    memcpy(&parsimony->stateSets[ nodes[i]->id ][j*nstate], NUCLEOTIDE_AMBIGUITY_STATES[patterns[j][idx]], 4*sizeof(int8_t));
//                }
//            }
//            else{
            
                memset(parsimony->stateSets[ nodes[i]->id ], 0, sizeof(int8_t)*nstate*sp_count);
                for ( int j = 0; j < sp_count; j++ ) {
                    if( patterns[j][ idx ] < nstate ){
                        parsimony->stateSets[ nodes[i]->id ][j*nstate+patterns[j][idx] ] = 1;
                    }
                    else {
                        for ( int k = 0; k < nstate; k++ ) {
                            parsimony->stateSets[ nodes[i]->id ][j*nstate+k] = 1;
                        }
                    }
                }
//            }
        }
    }
    
    
}

bool first_pass( Parsimony *pars, Node *node ){
    bool updated = pars->update_nodes[node->id];
    
    if( !Node_isleaf(node)) {
        
        bool update_child1 = first_pass( pars, Node_left(node) );
		bool update_child2 = first_pass( pars, Node_right(node) );
        
        if( update_child1 || update_child2 ){
            int8_t *states = pars->stateSets[node->id];
            int8_t *l = pars->stateSets[node->left->id];
            int8_t *r = pars->stateSets[node->right->id];
            int nstate = pars->sp->nstate;
            
            for ( int i = 0; i < pars->sp->count; i++ ) {
                pars->local_scores[node->class_id][i] = 0;
                
                int size = 0;
                // intersection
                for ( int j = 0; j < nstate; j++ ) {
                    states[i*nstate+j] = l[i*nstate+j] & r[i*nstate+j];
                    if( states[i*nstate+j] ){
                        size++;
                    }
                }
                
                if( size == 0 ){
                    // Union
                    for ( int j = 0; j < nstate; j++ ) {
                        states[i*nstate+j] = l[i*nstate+j] | r[i*nstate+j];
                    }
                    pars->local_scores[node->class_id][i] = 1;
                }
                
            }
            
            if( !Node_isleaf(node->left)){
                for ( int i = 0; i < pars->sp->count; i++ ) pars->local_scores[node->class_id][i] += pars->local_scores[node->left->class_id][i];
            }
            if( !Node_isleaf(node->right)){
                for ( int i = 0; i < pars->sp->count; i++ )pars->local_scores[node->class_id][i] += pars->local_scores[node->right->class_id][i];
            }
            updated = true;
        }
    }
    return updated;
}


double _score_fitch(Parsimony *parsimony){
    
    if ( !parsimony->update ) {
        return parsimony->score;
    }
    
    first_pass(parsimony, Tree_root(parsimony->tree));
    
    parsimony->score = 0;
    
    for ( int i = 0; i < parsimony->sp->count; i++ ) {
        parsimony->score += parsimony->local_scores[Tree_root(parsimony->tree)->class_id][i] * parsimony->sp->weights[i];
    }
    
    memset(parsimony->update_nodes, 0, Tree_node_count(parsimony->tree)*sizeof(bool));
    parsimony->update = false;
    
    return parsimony->score;
}

static void _reconstruct_aux( Parsimony *parsimony, Node *node ) {
    
    // We do not reconstruct known sequences
    if( !Node_isleaf(node) ){
        int nstate = parsimony->sp->nstate;
        
        for (int i = 0; i < parsimony->sp->count; i++) {
            
            if ( !Node_isroot(node) ){
                uint8_t parentState = parsimony->states[i][Node_id(Node_parent(node))];
                if( parsimony->stateSets[Node_id(node)][i*nstate+parentState]) {
                    parsimony->states[i][Node_id(node)] = parentState;
                    
                    

                }
                else {
                    
                    int j = 0;
                    for ( ; j < nstate; j++ ) {
                        if( parsimony->stateSets[Node_id(node)][i*nstate+j] ) break;
                    }
                    parsimony->states[i][Node_id(node)] = j;
                    
                    if( i!=parsimony->sp->count-1){
                        int count = 0;
                        for ( int jj = 0; jj < nstate; jj++ ) {
                            count += parsimony->stateSets[Node_id(Node_parent(node))][i*nstate+jj];
                        }
                        //if(count >1)fprintf(stderr, "%d %s => %d\n",i, Node_name(node), count);
                    }
//                    if( i!=parsimony->sp->count-1){
//                        int count = 0;
//                        for ( int jj = 0; jj < nstate; jj++ ) {
//                            count += parsimony->stateSets[Node_id(node)][i*nstate+jj];
//                        }
//                        fprintf(stderr, "%d %s %d => %d\n",i, Node_name(node), j, count);
//                    }
                }
            }
            else {
                int j = 0;
                for ( ; j < nstate; j++ ) {
                    if( parsimony->stateSets[Node_id(node)][i*nstate+j] ) break;
                }
                parsimony->states[i][Node_id(node)] = j;
                //if( i!=parsimony->sp->count-1)fprintf(stderr, "%d %s %d\n",i, Node_name(node),j);
            }
        }
    
    
        _reconstruct_aux(parsimony, Node_left(node) );
        _reconstruct_aux(parsimony, Node_right(node) );
    }
}

// indexing of parsimony->states[count] is different from sp->patterns[count]
// parsimony->states[count] is indexed by node->id while sp->patterns[count] indexing is arbitrary (order of the sequnences given in the constructor)
// the order of the first dimension (number of patterns) does not change
static void _reconstruct( Parsimony *parsimony ) {
    _score_fitch(parsimony);
    
    if( parsimony->states == NULL ){
        parsimony->states = ui8matrix(parsimony->sp->count, Tree_node_count(parsimony->tree) );
        Node **nodes = Tree_nodes(parsimony->tree);
        int sp_count = parsimony->sp->count;
         uint8_t **patterns = parsimony->sp->patterns;
        
        for ( int i = 0; i < Tree_node_count(parsimony->tree); i++ ) {
            
            if( Node_isleaf(nodes[i]) ){
                
                int idx = get_sequence_index(parsimony->sp, nodes[i]->name);
                
                for ( int j = 0; j < sp_count; j++ ) {
                    parsimony->states[j][nodes[i]->id] = patterns[j][idx];
                }
            }
        }
    }

    _reconstruct_aux(parsimony, Tree_root(parsimony->tree));
}

void first_pass_slow( Parsimony *pars, Node *node ){
    
    if( !Node_isleaf(node)) {
        
        first_pass_slow( pars, Node_left(node) );
		first_pass_slow( pars, Node_right(node) );
        
        int8_t *states = pars->stateSets[node->id];
        int8_t *l = pars->stateSets[node->left->id];
        int8_t *r = pars->stateSets[node->right->id];
        int nstate = pars->sp->nstate;
        
        for ( int i = 0; i < pars->sp->count; i++ ) {
            pars->local_scores[node->class_id][i] = 0;
            
            int size = 0;
            // intersection
            for ( int j = 0; j < nstate; j++ ) {
                states[i*nstate+j] = l[i*nstate+j] & r[i*nstate+j];
                if( states[i*nstate+j] ){
                    size++;
                }
            }
            
            if( size == 0 ){
                // Union
                for ( int j = 0; j < nstate; j++ ) {
                    states[i*nstate+j] = l[i*nstate+j] | r[i*nstate+j];
                }
                pars->scores[i]++;
            }
            
        }
        
    }
}


double _score_fitch_slow(Parsimony *parsimony){
    
    if ( !parsimony->update ) {
        return parsimony->score;
    }

#ifdef INT_WEIGHT
    memset(parsimony->scores, 0, parsimony->sp->count*sizeof(double));
#else
    memset(parsimony->scores, 0, parsimony->sp->count*sizeof(int32_t));
#endif
    first_pass_slow(parsimony, Tree_root(parsimony->tree));
    
    parsimony->score = 0;
    
    for ( int i = 0; i < parsimony->sp->count; i++ ) {
        parsimony->score += parsimony->scores[i] * parsimony->sp->weights[i];
    }
    
    memset(parsimony->update_nodes, 0, Tree_node_count(parsimony->tree)*sizeof(bool));
    parsimony->update = false;
    
    return parsimony->score;
}


#ifdef __SSE2__

void _first_pass_sse_4_slow( Parsimony *pars, Node *node ){
    
    if( !Node_isleaf(node) ) {
        
        _first_pass_sse_4_slow(pars, Node_left(node));
        _first_pass_sse_4_slow(pars, Node_right(node));
        
        
        int8_t * states = pars->stateSets[node->id];
        
        int8_t *ll = pars->stateSets[node->left->id];
        int8_t *rr = pars->stateSets[node->right->id];
        
        __m128i *l   = (__m128i*)ll;
        __m128i *r   = (__m128i*)rr;
        
        int8_t temp[16] __attribute__ ((aligned (16)));
        int32_t temp2[4];
        
        int n = (pars->sp->count*4/16)*16;
        int i = 0;
        int k = 0;
        
        for ( ; i < n; i+=16, l++, r++, states+=16 ) {
            
            _mm_store_si128((__m128i*)states, _mm_and_si128(*l, *r));
            
            _mm_store_si128((__m128i*)temp, _mm_or_si128(*l, *r));
            
            memcpy(temp2, states, sizeof(int8_t)*16);
            
            if( temp2[0] == 0 ){
                memcpy(states, temp, 4*sizeof(int8_t));
                pars->scores[k]++;
            }
            k++;
            
            if( temp2[1] == 0 ){
                memcpy(states+4, temp+4, 4*sizeof(int8_t));
                pars->scores[k]++;
            }
            k++;
            
            if( temp2[2] == 0 ){
                memcpy(states+8, temp+8, 4*sizeof(int8_t));
                pars->scores[k]++;
            }
            k++;
            
            if( temp2[3] == 0 ){
                memcpy(states+12, temp+12, 4*sizeof(int8_t));
                pars->scores[k]++;
            }
            k++;
        }
        
        ll += i;
        rr += i;
        
        for ( i /= 4; i < pars->sp->count; i++ ) {
            
            int size = 0;
            // intersection
            states[0] = ll[0] & rr[0];
            states[1] = ll[1] & rr[1];
            states[2] = ll[2] & rr[2];
            states[3] = ll[3] & rr[3];
            
            size = states[0] + states[1] + states[2] + states[3];
            
            
            if( size == 0 ){
                // Union
                states[0] = ll[0] | rr[0];
                states[1] = ll[1] | rr[1];
                states[2] = ll[2] | rr[2];
                states[3] = ll[3] | rr[3];
                
                pars->scores[k]++;
            }
            states += 4;
            ll += 4;
            rr += 4;
            k++;
        }
    }
}

double _score_fitch_4_sse_slow(Parsimony *parsimony){
    
    if ( !parsimony->update ) {
        return parsimony->score;
    }
    
#ifdef INT_WEIGHT
    memset(parsimony->scores, 0, parsimony->sp->count*sizeof(int32_t));
#else
    memset(parsimony->scores, 0, parsimony->sp->count*sizeof(double));
#endif
    
    _first_pass_sse_4_slow(parsimony, Tree_root(parsimony->tree));
    
    parsimony->score = 0;
#ifdef INT_WEIGHT

    int32_t temp[4] __attribute__ ((aligned (16)));
    memset(temp, 0, 4*sizeof(int32_t));
    
    __m128i *sum = (__m128i*)temp;
    __m128i *r1 = (__m128i*)parsimony->scores;
    __m128i *r2 = (__m128i*)parsimony->weights;

    int n = (parsimony->sp->count/4)*4;
    int i = 0;
    for ( ; i < n; i += 4 ) {
        *sum = _mm_add_epi32(*sum, muly(*r1, *r2));
        r1++;r2++;
    }
    
    parsimony->score = temp[0] + temp[1] + temp[2] + temp[3];
    
    for (  ; i < parsimony->sp->count; i++ ) {
        parsimony->score = parsimony->scores[i] * parsimony->weights[i];
    }
    //printf("\n%d %d %d %d = %f\n",temp[0] , temp[1] , temp[2] , temp[3], parsimony->score);exit(0);
    
#else
    __m128d r1,r2;
    int n = (parsimony->sp->count/2)*2;
    __m128d sum = _mm_setzero_pd();
    double temp[2] __attribute__ ((aligned (16)));
    
    int i = 0;
    for ( ; i < n; i+=2 ) {
        r1 = _mm_load_pd(parsimony->scores+i);
        r2 = _mm_load_pd(parsimony->sp->weights+i);
        sum = _mm_add_pd(sum, _mm_mul_pd(r1, r2));
    }
    _mm_store_pd(temp, sum);
    parsimony->score = temp[0] + temp[1];
    
    if( parsimony->sp->count & 1 ){
        i--;
        parsimony->score += parsimony->scores[i] * parsimony->sp->weights[i];
    }
#endif
    
    memset(parsimony->update_nodes, 0, Tree_node_count(parsimony->tree)*sizeof(bool));
    parsimony->update = false;
    
    return parsimony->score;
}


bool _first_pass_sse_4( Parsimony *pars, Node *node ){
    bool updated = pars->update_nodes[node->id];
    
    if( !Node_isleaf(node) ) {
        
        bool update_child1 = _first_pass_sse_4(pars, Node_left(node));
        bool update_child2 = _first_pass_sse_4(pars, Node_right(node));
        
        if( update_child1 || update_child2 ){
            
            int8_t * states = pars->stateSets[node->id];
            
            int8_t *ll = pars->stateSets[node->left->id];
            int8_t *rr = pars->stateSets[node->right->id];
            
            __m128i *l   = (__m128i*)ll;
            __m128i *r   = (__m128i*)rr;
            
            int8_t temp[16] __attribute__ ((aligned (16)));
            int32_t temp2[4];
            
            int n = (pars->sp->count*4/16)*16;
            int i = 0;
            int k = 0;
#ifdef INT_WEIGHT
            memset(pars->local_scores[node->class_id], 0, sizeof(int32_t)*pars->sp->count);
#else
            memset(pars->local_scores[node->class_id], 0, sizeof(double)*pars->sp->count);
#endif
            for ( ; i < n; i+=16, l++, r++, states+=16 ) {
                
                _mm_store_si128((__m128i*)states, _mm_and_si128(*l, *r));
                
                _mm_store_si128((__m128i*)temp, _mm_or_si128(*l, *r));
                
                memcpy(temp2, states, sizeof(int8_t)*16);
                
                if( temp2[0] == 0 ){
                    memcpy(states, temp, 4*sizeof(int8_t));
                    pars->local_scores[node->class_id][k] = 1;
                }
                k++;
                
                if( temp2[1] == 0 ){
                    memcpy(states+4, temp+4, 4*sizeof(int8_t));
                    pars->local_scores[node->class_id][k] = 1;
                }
                k++;
                
                if( temp2[2] == 0 ){
                    memcpy(states+8, temp+8, 4*sizeof(int8_t));
                    pars->local_scores[node->class_id][k] = 1;
                }
                k++;
                
                if( temp2[3] == 0 ){
                    memcpy(states+12, temp+12, 4*sizeof(int8_t));
                    pars->local_scores[node->class_id][k] = 1;
                }
                k++;
            }
            
            ll += i;
            rr += i;
            
            for ( i /= 4; i < pars->sp->count; i++ ) {
                int size = 0;
                // intersection
                states[0] = ll[0] & rr[0];
                states[1] = ll[1] & rr[1];
                states[2] = ll[2] & rr[2];
                states[3] = ll[3] & rr[3];
                
                size = states[0] + states[1] + states[2] + states[3];
                
                
                if( size == 0 ){
                    // Union
                    states[0] = ll[0] | rr[0];
                    states[1] = ll[1] | rr[1];
                    states[2] = ll[2] | rr[2];
                    states[3] = ll[3] | rr[3];
                    
                    pars->local_scores[node->class_id][k] = 1;
                }
                states += 4;
                ll += 4;
                rr += 4;
                k++;
            }
            
#ifdef INT_WEIGHT
            __m128i *r1,*r2;
            
            if( !Node_isleaf(node->left)){
                int n = (pars->sp->count/4)*4;
                r1 = (__m128i *)pars->local_scores[node->class_id];
                r2 = (__m128i *)pars->local_scores[node->left->class_id];
                int i = 0;
                for ( ; i < n; i+=4 ) {
                    *r1  = _mm_add_epi32(*r1, *r2);
                    r1++;r2++;
                }
                
                for ( ; i < pars->sp->count; i++ ) {
                    pars->local_scores[node->class_id][i] += pars->local_scores[node->left->class_id][i];
                }
            }
            if( !Node_isleaf(node->right)){
                int n = (pars->sp->count/4)*4;
                r1 = (__m128i *)pars->local_scores[node->class_id];
                r2 = (__m128i *)pars->local_scores[node->right->class_id];
                int i = 0;
                for ( ; i < n; i+=4 ) {
                    *r1  = _mm_add_epi32(*r1, *r2);
                    r1++;r2++;
                }
                
                for ( ; i < pars->sp->count; i++ ) {
                    pars->local_scores[node->class_id][i] += pars->local_scores[node->right->class_id][i];
                }
            }
#else
            __m128d *r1,*r2;
            if( !Node_isleaf(node->left)){
                r1 = (__m128d *)pars->local_scores[node->class_id];
                r2 = (__m128d *)pars->local_scores[node->left->class_id];
                int i = 0;
                for ( ; i < pars->sp->count; i+=2 ) {
                    *r1  = _mm_add_pd(*r1, *r2);
                    r1++;r2++;
                }
                
                if( pars->sp->count & 1 ){
                    --i;
                    pars->local_scores[node->class_id][i] += pars->local_scores[node->left->class_id][i];
                }
            }
            if( !Node_isleaf(node->right)){
                r1 = (__m128d *)pars->local_scores[node->class_id];
                r2 = (__m128d *)pars->local_scores[node->right->class_id];
                int i = 0;
                for ( ; i < pars->sp->count; i+=2 ) {
                    *r1  = _mm_add_pd(*r1, *r2);
                    r1++;r2++;
                }
                
                if( pars->sp->count & 1 ){
                    --i;
                    pars->local_scores[node->class_id][i] += pars->local_scores[node->right->class_id][i];
                }
            }
#endif
            updated = true;
            
        }
        
    }
    return updated;
}



double _score_fitch_4_sse(Parsimony *parsimony){
    
    if ( !parsimony->update ) {
        return parsimony->score;
    }
    
    _first_pass_sse_4(parsimony, Tree_root(parsimony->tree));
    
    
    parsimony->score = 0;
    int root_id = Tree_root(parsimony->tree)->class_id;
    
    
#ifdef INT_WEIGHT
    
    int n = (parsimony->sp->count/4)*4;
    int32_t temp[4] __attribute__ ((aligned (16)));
    memset(temp, 0, 4*sizeof(int32_t));
    
    __m128i *sum = (__m128i*)temp;

    
    __m128i *r1 = (__m128i*)parsimony->local_scores[root_id];
    __m128i *r2 = (__m128i*)parsimony->weights;
    int i = 0;
    for ( ; i < n; i += 4 ) {
        *sum = _mm_add_epi32(*sum, muly(*r1, *r2));
        r1++;r2++;
    }
    
    parsimony->score = temp[0] + temp[1] + temp[2] + temp[3];
    
    
    for (  ; i < parsimony->sp->count; i++ ) {
        parsimony->score += parsimony->local_scores[root_id][i] * parsimony->weights[i];
    }
    
#else
    __m128d r1,r2;
    int n = (parsimony->sp->count/2)*2;
    __m128d sum = _mm_setzero_pd();
    double temp[2] __attribute__ ((aligned (16)));
    
    int i = 0;
    double *scores = parsimony->local_scores[root_id];
    for ( ; i < n; i+=2 ) {
        
        r1 = _mm_load_pd(scores+i);
        r2 = _mm_load_pd(parsimony->sp->weights+i);
        sum = _mm_add_pd(sum, _mm_mul_pd(r1, r2));
    }
    _mm_store_pd(temp, sum);
    parsimony->score = temp[0] + temp[1];
    
    if( parsimony->sp->count & 1 ){
        i--;
        parsimony->score += scores[i] * parsimony->sp->weights[i];
    }
#endif
    
    memset(parsimony->update_nodes, 0, Tree_node_count(parsimony->tree)*sizeof(bool));
    parsimony->update = false;
    
    return parsimony->score;
}
#endif

#pragma mark -
#pragma mark ParsimonyModel

void _parsimony_model_handle_change( Model *self, Model *model, int index ){
	Parsimony *parsimony = (Parsimony*)self->obj;
	if ( strcmp(model->type, "tree") == 0 ) {
		parsimony->update_nodes[index] = true;
		parsimony->update = true;
	}
	else {
		fprintf(stderr, "%s of type %s\n", model->name, model->type);
		error("Unknown change in Parsimony\n");
	}
}

static double _parsimony_model_logP(Model *self){
	self->lp = ((Parsimony*)self->obj)->calculate((Parsimony*)self->obj);
	return self->lp;
}

static void _parsimony_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free parsimony model %s\n", self->name);
		Parsimony* parsimony = (Parsimony*)self->obj;
		Model* mtree = (Model*)self->data;
		
		mtree->free(mtree);
		
		free_SitePattern(parsimony->sp);
		free_Parsimony(parsimony);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _parsimony_model_clone(Model* self, Hashtable* hash){
	Model* mtree = self->data;
	Model *mtreeclone = NULL;

	if(Hashtable_exists(hash, mtree->name)){
		mtreeclone = Hashtable_get(hash, mtree->name);
		mtreeclone->ref_count++;
	}
	else{
		mtreeclone = mtree->clone(mtree, hash);
		Hashtable_add(hash, mtreeclone->name, mtreeclone);
	}
	
	Parsimony* parsimony = (Parsimony*)self->obj;
	Parsimony* clonetlk = new_Parsimony(parsimony->sp, mtreeclone->obj);
	
	Model* clone = new_ParsimonyModel(self->name, clonetlk, mtreeclone);
	Hashtable_add(hash, clone->name, clone);
	mtreeclone->free(mtreeclone);

	return clone;
}

static void _parsimony_model_get_free_parameters(Model* model, Parameters* parameters){
	Model* mtree = model->data;
	mtree->get_free_parameters(mtree, parameters);
}

Model * new_ParsimonyModel(char* name, Parsimony* parsimony, Model* tree){
	Model *model = new_Model("parsimony", name, parsimony);
	tree->listeners->add( tree->listeners, model );
	model->data = tree;
	tree->ref_count++;
	model->update = _parsimony_model_handle_change;
	model->logP = _parsimony_model_logP;
	model->free = _parsimony_model_free;
	model->clone = _parsimony_model_clone;
	model->get_free_parameters = _parsimony_model_get_free_parameters;
	return model;
}

Model * new_ParsimonyModel_from_json(json_node*node, Hashtable*hash){
	json_node* patterns_node = get_json_node(node, "sitepattern");
	json_node* tree_node = get_json_node(node, "tree");
	SitePattern* patterns = NULL;
	Model* mtree = NULL;
	
	if(patterns_node->node_type == MJSON_STRING){
		char* ref = (char*)patterns_node->value;
		// check it starts with a &
		patterns = Hashtable_get(hash, ref+1);
		patterns->ref_count++;
	}
	else{
		char* id = get_json_node_value_string(patterns_node, "id");
		patterns = new_SitePattern_from_json(patterns_node, hash);
		Hashtable_add(hash, id, patterns);
	}
	
	if (tree_node->node_type == MJSON_STRING) {
		char* ref = (char*)tree_node->value;
		// check it starts with a &
		mtree = Hashtable_get(hash, ref+1);
		mtree->ref_count++;
	}
	else{
		char* id = get_json_node_value_string(tree_node, "id");
		mtree = new_TreeModel_from_json(tree_node, hash);
		Hashtable_add(hash, id, mtree);
	}
	
	Parsimony* parsimony = new_Parsimony(patterns, mtree->obj);
	char* id = get_json_node_value_string(node, "id");
	Model* model = new_ParsimonyModel(id, parsimony, mtree);
	mtree->free(mtree);
	return model;
}
