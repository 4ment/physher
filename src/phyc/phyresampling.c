/*
 *  phyclustering.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/01/2015.
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

#include "phyresampling.h"

#include <stdlib.h>
#include <assert.h>

#include "matrix.h"
#include "random.h"



// should only be used with aligned sequences (same length)
Sequences * Sequences_bootstrap( const Sequences *sequences ){
    Sequences *sample = new_Sequences( sequences->size );
    DataType *dataType = sequences->datatype;
    assert(dataType);
    for ( int i = 0; i < sequences->size; i++ ) {
        Sequences_add(sample, new_Sequence2(sequences->seqs[i]->name, sequences->length));
    }
    sample->datatype = clone_DataType(dataType);
    
    int *random_pos = ivector( sequences->length );
    for ( int j = 0; j < sequences->length; j++ ) {
        random_pos[j] = random_int(sequences->length-1);
    }
    
    for ( int i = 0; i < sequences->size; i++ ) {
        for ( int j = 0; j < sequences->length; j++ ) {
            sample->seqs[i]->seq[j] = sequences->seqs[i]->seq[random_pos[j]];
            sample->seqs[i]->length++;
        }
    }
    free(random_pos);
    return sample;
}

// should only be used with aligned sequences (same length)
Sequences * Sequences_jackknife( const Sequences *sequences, int index ){
    Sequences *sample = new_Sequences( sequences->size );
    DataType *dataType = sequences->datatype;
    assert(dataType);
    for ( int i = 0; i < sequences->size; i++ ) {
        Sequence *seq = Sequence_clone(sequences->seqs[i]);
        if( index != sequences->length-1 ){
            memmove(&seq->seq[index], &seq->seq[index+1], (sequences->length-index-1)*sizeof(char));
        }
        seq->length--;
        seq->seq[seq->length] = '\0';
        Sequences_add(sample, seq);
        
    }
    sample->datatype = clone_DataType(dataType);
    return sample;
}

Sequences * Sequences_jackknife_n( const Sequences *sequences, int n ){
    Sequences *sample = new_Sequences( sequences->size );
    DataType *dataType = sequences->datatype;
    assert(dataType);
    
    int *order = ivector(sequences->length);
    for ( int i = 0; i < sequences->length; i++ ) {
        order[i] = i;
    }
    randomize_ivector(order, sequences->length);
    qsort(order, n, sizeof(int), qsort_desc_ivector);
    
    for ( int i = 0; i < sequences->size; i++ ) {
        Sequence *seq = Sequence_clone(sequences->seqs[i]);
        for ( int j = 0; j < n; j++ ) {
            int index = order[j];
            if( index != sequences->length-1 ){
                memmove(&seq->seq[index], &seq->seq[index+1], (sequences->length-index-1)*sizeof(char));
            }
            seq->length--;
            seq->seq[seq->length] = '\0';
        }
        Sequences_add(sample, seq);
        
    }
    sample->datatype = clone_DataType(dataType);
    free(order);
    return sample;
}


#pragma mark *** SitePattern ***

void printit(SitePattern *sp){
    for ( int i = 0; i < sp->count; i++ ) {
        for ( int j = 0; j < sp->size; j++ ) {
            printf("%d,", sp->patterns[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

SitePattern * SitePattern_bootstrap( const SitePattern *original ){
    SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
    assert(sp);
    
    assert(original->datatype);
    sp->datatype = clone_DataType(original->datatype);
    
    assert(original->count > 0);
    assert(original->nsites > 0);
    
    sp->id = original->id+1;
    sp->patterns = (uint8_t**)malloc(original->count * sizeof(uint8_t *) );
    assert(sp->patterns);
    sp->indexes = NULL;
    sp->weights = dvector(original->count);
    sp->count = 0;
    sp->size  = original->size;
    sp->nstate = original->nstate;
    
    double tot= 0;
    for ( int i = 0; i < original->count; i++) {
        tot += original->weights[i];
    }
    
    for ( int i = 0; i < original->nsites; i++ ) {
        int rpos = roulette_wheel2(original->weights, original->count, tot);
        int p = 0;
        for ( ; p < sp->count; p++ ) {
            
            if( memcmp(original->patterns[rpos], sp->patterns[p], original->size*sizeof(uint8_t)) == 0 ){
                sp->weights[p] += 1.;
                //sp->indexes[i] = p;
                break;
            }
        }
        if( p == sp->count ){
            sp->weights[sp->count] = 1.;
            sp->patterns[sp->count] = (uint8_t*)malloc(original->size * sizeof(uint8_t) );
            assert(sp->patterns);
            memcpy( sp->patterns[sp->count], original->patterns[rpos], original->size * sizeof(uint8_t) );
            //sp->indexes[i] = sp->count;
            sp->count++;
        }
        
    }
    
    sp->nsites = original->nsites;
    
    sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
    sp->patterns = realloc(sp->patterns, sp->count * sizeof(uint8_t*) );
    
    sp->names = (char**)malloc( original->size * sizeof(char*) );
    assert(sp->names);
    for ( int i = 0; i < original->size; i++) {
        sp->names[i] = String_clone( original->names[i] );
    }
    return sp;
}


SitePattern * SitePattern_jackknife( const SitePattern *original, int index ){
    SitePattern *sp = clone_SitePattern(original);
    assert(sp);
    sp->datatype = clone_DataType(original->datatype);
    
    if( original->weights[index] > 1 ){
        sp->weights[index]--;
        sp->nsites--;
    }
    else {
        free(sp->patterns[index]);
        sp->patterns[index] = NULL;

        int p = index+1;
        for ( ; p < sp->count; p++ ) {
            sp->patterns[p-1] = sp->patterns[p];
            sp->weights[p-1] = sp->weights[p];
        }
        
        sp->nsites--;
        sp->count--;
        sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
        sp->patterns = realloc(sp->patterns, sp->count * sizeof(uint8_t*) );
//        index--;
//        for ( p = 0; p < original->nsites; p++ ) {
//            if( sp->indexes[p] > index ){
//                sp->indexes[p]--;
//            }
//        }
       
    }
    return sp;
}

SitePattern * SitePattern_jackknife_n( SitePattern *original, int n){
    int *order = ivector(original->nsites);
    int k = 0;
    for ( int i = 0; i < original->count; i++ ) {
        for ( int j = 0; j < original->weights[i]; j++ ) {
            order[k++] = i;
        }
    }
    randomize_ivector(order, original->nsites);
    
    print_ivector(order, original->nsites);
    
    double *weights = clone_dvector(original->weights, original->count);
    for ( int i = 0; i < n; i++ ) {
        weights[order[i]]--;
    }
    
    SitePattern *sp = SitePattern_reweight(original, weights);
    
    free(weights);
    free(order);
    return sp;
}


// if a weight == 0 then it is removed
SitePattern * SitePattern_reweight( const SitePattern *original, const double *weights ){
    SitePattern *sp = clone_SitePattern(original);
    
    memcpy(sp->weights, weights, original->count*sizeof(double));
    int count = sp->count;

    for ( int i = count-1; i >= 0; i-- ) {
        if( sp->weights[i] == 0.0 ){
            
            free(sp->patterns[i]);
            sp->patterns[i] = NULL;
            int p = i+1;
            
            while ( p < sp->count ) {
                sp->patterns[p-1] = sp->patterns[p];
                sp->weights[p-1]  = sp->weights[p];
                p++;
            }
            sp->count--;
            sp->nsites--;
        }
    }
    
    sp->nsites = 0;
    
    for ( int i = 0; i < sp->count; i++ ) {
        sp->nsites += sp->weights[i];
    }
    
    if ( count > sp->count ) {
        sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
        sp->patterns = realloc(sp->patterns, sp->count * sizeof(uint8_t*) );
    }
    
    return sp;
}

