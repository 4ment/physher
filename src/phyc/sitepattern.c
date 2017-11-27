/*
 *  sitepattern.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/15/11.
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
#include "sitepattern.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "matrix.h"
#include "mstring.h"
#include "parser.h"
#include "sequenceio.h"
#include "mstring.h"
#include "geneticcode.h"
#include "random.h"

#include "datatype.h"

static bool compare_patterns( const uint8_t *a, const uint8_t *b, const int n );

static void _make_patterns( const Sequences *a, uint8_t **patterns, int start, int end, int every );
static void _make_patterns_codon( const Sequences *aln, uint8_t **patterns, int start, int end );
static void _make_patterns_generic( const Sequences *aln, uint8_t **patterns );

SitePattern * new_SitePattern( const Sequences *aln ){
	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
	assert(sp);
	sp->id = 0;
	
	sp->datatype = clone_DataType(aln->datatype);
	uint8_t **temp_patterns = NULL;
	sp->nstate = sp->datatype->state_count(sp->datatype);
	
	if ( aln->datatype->type == DATA_TYPE_CODON) {
		sp->nsites = aln->length/3;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns_codon(aln, temp_patterns, 0, aln->length);
	}
	else if ( aln->datatype->type == DATA_TYPE_GENERIC) {
		sp->nsites = 1;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns_generic( aln, temp_patterns );
	}
	else {
		sp->nsites = aln->length;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns( aln, temp_patterns, 0, aln->length, 1 );
	}
	
	assert(sp->nsites);
	sp->patterns = (uint8_t**)malloc(sp->nsites * sizeof(uint8_t *) );
	assert(sp->patterns);
	sp->indexes = ivector(sp->nsites);
	sp->weights = dvector(sp->nsites);
	sp->count = 0;
	sp->size = aln->size;
	
	
	int site, p;
	for( site = 0; site < sp->nsites; site++ ){
		for ( p = 0; p < sp->count; p++ ) {
			if( compare_patterns( temp_patterns[site], sp->patterns[p], aln->size ) ){
				sp->weights[p] += 1.;
				sp->indexes[site] = p;
				break;
			}
		}
		if( p == sp->count ){
			sp->weights[sp->count] = 1.;
			sp->patterns[sp->count] = (uint8_t*)malloc(aln->size * sizeof(uint8_t) );
			assert(sp->patterns[sp->count]);
			memcpy( sp->patterns[sp->count], temp_patterns[site], aln->size * sizeof(uint8_t) );
			sp->indexes[site] = sp->count;
			sp->count++;
		}
	}
	assert(sp->count);
	sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
	sp->patterns = realloc(sp->patterns, sp->count * sizeof(uint8_t*) );
	
	//	//
	//	sp->patterns2 = (uint8_t**)malloc(aln->size * sizeof(uint8_t *) );
	//	for ( int i = 0; i < aln->size; i++ ) {
	//		sp->patterns2[i] = (uint8_t*)calloc(sp->count, sizeof(uint8_t) );
	//	}
	//	for ( int i = 0; i < sp->count; i++ ) {
	//		for ( int j = 0; j < aln->size ; j++ ) {
	//			sp->patterns2[j][i] = sp->patterns[i][j];
	//		}
	//	}
	//	//
	
	sp->names = (char**)malloc( aln->size * sizeof(char*) );
	assert(sp->names);
	
	for ( int i = 0; i < aln->size; i++) {
		sp->names[i] = String_clone( aln->seqs[i]->name );
	}
	
	for ( int i = 0; i < sp->nsites; i++) {
		free(temp_patterns[i]);
	}
	free(temp_patterns);
	sp->ref_count = 1;
	
	return sp;
}


SitePattern * new_SitePattern2( const Sequences *aln, int start, int length, int every ){
	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
	assert(sp);
	sp->id = 0;
	
	sp->datatype = clone_DataType(aln->datatype);
	uint8_t **temp_patterns = NULL;
	sp->nstate = sp->datatype->state_count(sp->datatype);
	
	if ( aln->datatype->type == DATA_TYPE_CODON) {
		sp->nsites = length/3;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns_codon(aln, temp_patterns, start, length);
	}
	else if ( aln->datatype->type == DATA_TYPE_GENERIC) {
		sp->nsites = 1;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns_generic( aln, temp_patterns );
	}
	else {
		sp->nsites = length/every;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns( aln, temp_patterns, start, length, every );
	}
	
	assert(sp->nsites);
	sp->patterns = (uint8_t**)malloc(sp->nsites * sizeof(uint8_t *) );
	assert(sp->patterns);
	sp->indexes = ivector(sp->nsites);
	sp->weights = dvector(sp->nsites);
	sp->count = 0;
	sp->size = aln->size;
	
	
	int site, p;
	for( site = 0; site < sp->nsites; site++ ){
		for ( p = 0; p < sp->count; p++ ) {
			if( compare_patterns( temp_patterns[site], sp->patterns[p], aln->size ) ){
				sp->weights[p] += 1.;
				sp->indexes[site] = p;
				break;
			}
		}
		if( p == sp->count ){
			sp->weights[sp->count] = 1.;
			sp->patterns[sp->count] = (uint8_t*)malloc(aln->size * sizeof(uint8_t) );
			assert(sp->patterns[sp->count]);
			memcpy( sp->patterns[sp->count], temp_patterns[site], aln->size * sizeof(uint8_t) );
			sp->indexes[site] = sp->count;
			sp->count++;
		}
	}
	assert(sp->count);
	sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
	sp->patterns = realloc(sp->patterns, sp->count * sizeof(uint8_t*) );
	
	//	//
	//	sp->patterns2 = (uint8_t**)malloc(aln->size * sizeof(uint8_t *) );
	//	for ( int i = 0; i < aln->size; i++ ) {
	//		sp->patterns2[i] = (uint8_t*)calloc(sp->count, sizeof(uint8_t) );
	//	}
	//	for ( int i = 0; i < sp->count; i++ ) {
	//		for ( int j = 0; j < aln->size ; j++ ) {
	//			sp->patterns2[j][i] = sp->patterns[i][j];
	//		}
	//	}
	//	//
	
	sp->names = (char**)malloc( aln->size * sizeof(char*) );
	assert(sp->names);
	
	for ( int i = 0; i < aln->size; i++) {
		sp->names[i] = String_clone( aln->seqs[i]->name );
	}
	
	for ( int i = 0; i < sp->nsites; i++) {
		free(temp_patterns[i]);
	}
	free(temp_patterns);
	sp->ref_count = 1;
	
	return sp;
}

SitePattern * new_SitePattern3( const Sequences *aln, int start, int length, int a, int b ){
	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
	assert(sp);
	sp->id = 0;
	
	sp->datatype = clone_DataType(aln->datatype);
	uint8_t **temp_patterns = NULL;
	sp->nstate = sp->datatype->state_count(sp->datatype);
	
	if ( aln->datatype->type == DATA_TYPE_CODON) {
		sp->nsites = length/3;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns_codon(aln, temp_patterns, start, length);
	}
	else if ( aln->datatype->type == DATA_TYPE_GENERIC) {
		sp->nsites = 1;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns_generic( aln, temp_patterns );
	}
	else {
		sp->nsites = length/3*2;
		temp_patterns = ui8matrix( sp->nsites, aln->size);
		_make_patterns( aln, temp_patterns, start+a, length, 3 );
		_make_patterns( aln, temp_patterns+length/3, start+b, length, 3 );
	}
	
	assert(sp->nsites);
	sp->patterns = (uint8_t**)malloc(sp->nsites * sizeof(uint8_t *) );
	assert(sp->patterns);
	sp->indexes = ivector(sp->nsites);
	sp->weights = dvector(sp->nsites);
	sp->count = 0;
	sp->size = aln->size;
	
	
	int site, p;
	for( site = 0; site < sp->nsites; site++ ){
		for ( p = 0; p < sp->count; p++ ) {
			if( compare_patterns( temp_patterns[site], sp->patterns[p], aln->size ) ){
				sp->weights[p] += 1.;
				sp->indexes[site] = p;
				break;
			}
		}
		if( p == sp->count ){
			sp->weights[sp->count] = 1.;
			sp->patterns[sp->count] = (uint8_t*)malloc(aln->size * sizeof(uint8_t) );
			assert(sp->patterns[sp->count]);
			memcpy( sp->patterns[sp->count], temp_patterns[site], aln->size * sizeof(uint8_t) );
			sp->indexes[site] = sp->count;
			sp->count++;
		}
	}
	assert(sp->count);
	sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
	sp->patterns = realloc(sp->patterns, sp->count * sizeof(uint8_t*) );
	
	//	//
	//	sp->patterns2 = (uint8_t**)malloc(aln->size * sizeof(uint8_t *) );
	//	for ( int i = 0; i < aln->size; i++ ) {
	//		sp->patterns2[i] = (uint8_t*)calloc(sp->count, sizeof(uint8_t) );
	//	}
	//	for ( int i = 0; i < sp->count; i++ ) {
	//		for ( int j = 0; j < aln->size ; j++ ) {
	//			sp->patterns2[j][i] = sp->patterns[i][j];
	//		}
	//	}
	//	//
	
	sp->names = (char**)malloc( aln->size * sizeof(char*) );
	assert(sp->names);
	
	for ( int i = 0; i < aln->size; i++) {
		sp->names[i] = String_clone( aln->seqs[i]->name );
	}
	
	for ( int i = 0; i < sp->nsites; i++) {
		free(temp_patterns[i]);
	}
	free(temp_patterns);
	sp->ref_count = 1;
	
	return sp;
}

void free_SitePattern( SitePattern *sp ){
	if(sp->ref_count == 1){
		if(sp->patterns != NULL ){
			free_ui8matrix(sp->patterns, sp->count);
			free(sp->weights);
		}
		if(sp->indexes != NULL)free(sp->indexes);
		
		if(sp->names != NULL ){
			for (int i = 0; i < sp->size; i++) {
				free(sp->names[i]);
			}
			free(sp->names);
		}
		free_DataType(sp->datatype);
		free(sp);
	}
	else{
		sp->ref_count--;
	}
}

SitePattern * clone_SitePattern( const SitePattern *sp ){
    SitePattern *newsp = (SitePattern *)malloc( sizeof(SitePattern) );
    assert(newsp);
    
    newsp->id = sp->id;
    newsp->id++;
    newsp->datatype = clone_DataType(sp->datatype);
    
    newsp->patterns = ui8matrix(sp->count, sp->size);
    for ( int p = 0; p < sp->count; p++ ) {
        memcpy(newsp->patterns[p], sp->patterns[p], sp->size * sizeof(uint8_t));
    }
    
    newsp->indexes = clone_ivector( sp->indexes, sp->nsites);
    
    newsp->weights = clone_dvector( sp->weights, sp->count );
    
    newsp->count = sp->count;
    newsp->size  = sp->size;
    
    newsp->nsites  = sp->nsites;
    
    newsp->nstate = sp->nstate;
    
    newsp->names = (char**)malloc( sp->size * sizeof(char*) );
    assert(newsp->names);
    for ( int i = 0; i < sp->size; i++) {
        newsp->names[i] = String_clone( sp->names[i] );
    }
    newsp->ref_count = 1;
    return newsp;
}

static bool _compare_patterns_indexes( const uint8_t *a, const uint8_t *b, int n, const int *indexes ){
	for ( int i = 0; i < n; i++) {
		if( a[ indexes[i] ] != b[i] ) return false;
	}
	return true;
}


//SitePattern * SitePattern_subset( const SitePattern *original, const int *indexes, int size ){
//	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
//	assert(sp);
//	sp->id = original->id+1;
//	sp->datatype = original->datatype;
//	
//	sp->patterns = (uint8_t**)malloc(original->count * sizeof(uint8_t *) );
//	//sp->indexes = ivector(original->nsites);
//	sp->weights = dvector(original->count);
//	sp->count = 0;
//	sp->size  = size;
//    sp->nstate = original->nstate;
//	
//	for ( int i = 0; i < original->nsites; i++ ) {
//		int pos = original->indexes[i];
//		int p = 0;
//		for ( ; p < sp->count; p++ ) {
//			if( _compare_patterns_indexes( original->patterns[pos], sp->patterns[p], sp->size, indexes ) ){
//				sp->weights[p] += 1.;
//				sp->indexes[i] = p;
//				break;
//			}
//		}
//		if( p == sp->count ){
//			sp->weights[sp->count] = 1.;
//			sp->patterns[sp->count] = (uint8_t*)malloc(sp->size * sizeof(uint8_t) );
//			for ( int j = 0; j < sp->size; j++ ) {
//				sp->patterns[sp->count][j] = original->patterns[pos][ indexes[j] ];
//			}
//			//memcpy( sp->patterns[sp->count], original->patterns[rpos], original->size * sizeof(uint8_t) );
//			sp->indexes[i] = sp->count++;
//		}
//		
//	}
//	
//	sp->nsites = original->nsites;
//	
//	sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
//	sp->patterns = realloc(sp->patterns, sp->count * sizeof(char*) );
//	
//	sp->names = (char**)malloc( sp->size * sizeof(char*) );
//	assert(sp->names);
//	
//	for ( int i = 0; i < sp->size; i++) {
//		sp->names[i] = String_clone( original->names[ indexes[i] ] );
//	}
//	return sp;
//}

SitePattern* SitePattern_merge( const SitePattern* sitePattern1, const SitePattern* sitePattern2  ){
	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
	assert(sp);
	sp->id = 0;
	
	sp->datatype = clone_DataType(sitePattern1->datatype);
	sp->nstate = sp->datatype->state_count(sp->datatype);
	
	sp->size = sitePattern1->size;
	sp->nsites = sitePattern1->nsites + sitePattern2->nsites;
	
	assert(sp->nsites);
	sp->patterns = (uint8_t**)malloc(sp->nsites * sizeof(uint8_t *) );
	assert(sp->patterns);
	sp->indexes = ivector(sp->nsites);
	sp->weights = dvector(sp->nsites);
	sp->count = 0;
	
	for( int site = 0; site < sitePattern1->nsites; site++ ){
		sp->weights[site] = sitePattern1->weights[site];
		sp->patterns[site] = (uint8_t*)malloc(sp->size * sizeof(uint8_t) );
		assert(sp->patterns[site]);
		memcpy( sp->patterns[site], sitePattern1->patterns[site], sp->size * sizeof(uint8_t) );
		sp->indexes[site] = site; // does not make sense
		sp->count++;
	}

	for( int site = 0; site < sitePattern2->nsites; site++ ){
		int p = 0;
		for ( ; p < sp->count; p++ ) {
			if( compare_patterns( sitePattern2->patterns[site], sp->patterns[p], sp->size ) ){
				sp->weights[p] += 1.;
				sp->indexes[site] = p;// does not make sense
				break;
			}
		}
		if( p == sp->count ){
			sp->weights[sp->count] = 1.;
			sp->patterns[sp->count] = (uint8_t*)malloc(sp->size * sizeof(uint8_t) );
			assert(sp->patterns[sp->count]);
			memcpy( sp->patterns[sp->count], sitePattern2->patterns[site], sp->size * sizeof(uint8_t) );
			sp->indexes[site] = sp->count; // does not make sense
			sp->count++;
		}
	}
	assert(sp->count);
	sp->weights  = realloc(sp->weights, sp->count * sizeof(double) );
	sp->patterns = realloc(sp->patterns, sp->count * sizeof(uint8_t*) );
	
	sp->names = (char**)malloc( sitePattern1->size * sizeof(char*) );
	assert(sp->names);
	
	for ( int i = 0; i < sp->size; i++) {
		sp->names[i] = String_clone( sitePattern1->names[i] );
	}
	
	sp->ref_count = 1;
	
	return sp;
}

// MARK: SitePattern_split not tested
SitePattern ** SitePattern_split( const SitePattern *sitePattern, const int count  ){
	SitePattern **sps = (SitePattern **)malloc( count*sizeof(SitePattern*) );
    assert(sps);
	
	unsigned chunk_size = sitePattern->count / count;
	unsigned last_chunk_size = sitePattern->count - ( chunk_size * count );
	int p = 0;
	
	for (int i = 0; i < count; i++) {
		sps[i] = (SitePattern *)malloc( sizeof(SitePattern) );
        assert(sps[i]);
		sps[i]->id = i+1;
		sps[i]->size = sitePattern->size;
        
        sps[i]->datatype = clone_DataType(sitePattern->datatype);
        sps[i]->nstate = sitePattern->nstate;
		
		if ( i == count -1 ) {
			sps[i]->count = last_chunk_size;
		}
		else {
			sps[i]->count = count;
		}
		sps[i]->size = sitePattern->size;
		
		sps[i]->patterns =  (uint8_t**)malloc(sps[i]->count * sizeof(uint8_t *) );
        assert(sps[i]->patterns);
		for ( int pp = 0; pp < sps[i]->count; pp++ ) {
			sps[i]->patterns[pp] = (uint8_t*)malloc(sps[i]->size * sizeof(uint8_t) );
            assert(sps[i]->patterns[pp]);
			memcpy( sps[i]->patterns[pp], sitePattern->patterns[p], sitePattern->size * sizeof(uint8_t) );
		}
		
		sps[i]->weights  = dvector(sps[i]->count);
		memcpy(sps[i]->weights, sitePattern->weights+p, sps[i]->count*sizeof(double) );
		
		sps[i]->names = (char**)malloc( sitePattern->size * sizeof(char*) );
		assert(sps[i]->names);
		
		for ( int j = 0; j < sitePattern->size; j++) {
			sps[i]->names[j] = String_clone(sitePattern->names[j]);
		}
		
		sps[i]->nsites = sitePattern->nsites;
		sps[i]->indexes = ivector(sitePattern->nsites);
		p += sps[i]->count;
		sps[i]->ref_count = 1;
	}	
	
	return sps;
}


char * SitePattern_stringify( const SitePattern *sp ){
	StringBuffer *buffer = new_StringBuffer(100);
	
	SitePattern_bufferize(buffer, sp);
	
	char *final = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return final;
	
}

StringBuffer * SitePattern_bufferize( StringBuffer *buffer, const SitePattern *sp ){
	StringBuffer_append_format(buffer, "(SitePattern:\n(id:\"sp%d\")\n",sp->id);
	return buffer;
}

void * SitePattern_SML_to_object( ObjectStore *store, SMLNode node ){
	fprintf(stderr, "SitePattern_SML_to_object\n");
	SitePattern *sp = NULL;
	
	sp = new_SitePattern(readFasta( SML_get_data_of_child( node, "file") ) );
	
	char *id = SML_get_data_of_child( node, "id");
	if ( id != NULL ) {
		while( *id < '0' || *id > '9' ){
			id++;
		}
		sp->id = atoi( id );
	}
	else sp->id = 0;
	
	return sp;
}



void SitePattern_save( const SitePattern *sitepattern, const char *output ){
	Sequences *sequences = SitePattern_to_Sequences(sitepattern);
	Sequences_save_fasta(sequences, output);
	free_Sequences(sequences);
	
}

Sequences * SitePattern_to_Sequences( const SitePattern *sp ){
    Sequences *sequences = new_Sequences(sp->size);
    StringBuffer *buff = new_StringBuffer(sp->nsites);
    
    // that's an inefficient loop
    for ( int j = 0; j < sp->size; j++ ) {
        StringBuffer_empty(buff);
        for ( int s = 0; s < sp->count; s++ ) {
            
            for ( int c = 0; c < sp->weights[s]; c++ ) {
                StringBuffer_append_string( buff, sp->datatype->state_string(sp->datatype,sp->patterns[s][j]) );
            }
        }
        Sequence *seq = new_Sequence(sp->names[j], buff->c);
        Sequences_add(sequences, seq);
    }
    sequences->aligned  = true;
    sequences->datatype = clone_DataType(sp->datatype);
    
    free_StringBuffer(buff);
    return sequences;
}

//Sequences * SitePattern_to_Sequences( const SitePattern *sp ){
//    Sequences *sequences = new_Sequences(sp->size);
//    StringBuffer *buff = new_StringBuffer(sp->nsites);
//    
//    // that's an inefficient loop
//    for ( int j = 0; j < sp->size; j++ ) {
//        StringBuffer_empty(buff);
//        for ( int s = 0; s < sp->nsites; s++ ) {
//            int i = sp->indexes[s];
//            if ( sp->datatype->type == DATA_TYPE_CODON){
//                char codon[3];
//                encoding_to_codon(sp->patterns[i][j],  sp->datatype->genetic_code, codon);
//                StringBuffer_append_char( buff, codon[0] );
//                StringBuffer_append_char( buff, codon[1] );
//                StringBuffer_append_char( buff, codon[2] );
//            }
//            else {
//                StringBuffer_append_char( buff, sp->datatype->state(sp->datatype,sp->patterns[i][j]) );
//            }
//            
//        }
//        Sequence *seq = new_Sequence(sp->names[j], buff->c);
//        Sequences_add(sequences, seq);
//    }
//    sequences->aligned  = true;
//    sequences->datatype = sp->datatype;
//    
//    free_StringBuffer(buff);
//    return sequences;
//}



bool compare_patterns( const uint8_t *a, const uint8_t *b, const int n ){
	int i;
	for ( i = 0; i < n; i++) {
		if( a[i] != b[i] ) return false;
	}
	return true;
}

int get_sequence_index( const SitePattern *sp, const char *name ){
	int i;
	for (i = 0; i < sp->size; i++) {
		if( strcmp(name, sp->names[i] ) == 0 ) return i;
	}
	return -1;
}

int frequencies( const SitePattern *sp, double *freqs, const uint8_t nstate ){
	int j;
	int ambiguity = 0;
	for (int i = 0; i < sp->count; i++) {
		for ( j = 0; j < sp->size; j++) {
			if( sp->patterns[i][j] < nstate ) freqs[sp->patterns[i][j]]++;
			else ambiguity++;
		}
	}
	return ambiguity;
}

void compare_sitepattern( const SitePattern *sp1, const SitePattern *sp2 ){
	assert(sp1->id == sp2->id);
	assert(sp1->size == sp2->size);
	assert(sp1->count == sp2->count);
	
	assert( memcmp(sp1->weights, sp2->weights, sp1->count*sizeof(double) ) == 0 );
	//assert( memcmp(sp1->indexes, sp2->indexes, sp1->count*sizeof(int) ) == 0 );
	
	for (int i = 0; i < sp1->count; i++) {
		assert( memcmp(sp1->patterns[i], sp2->patterns[i], sp1->size*sizeof(char) ) == 0 );
	}
	
	for (int i = 0; i < sp1->size; i++) {
		assert( memcmp(sp1->names[i], sp2->names[i], sp1->size*sizeof(char) ) == 0 );
	}
}

void SitePattern_print( const SitePattern *sp ){
	printf("%d site patterns of length %d\n", sp->count, sp->size );
	for (int i = 0; i < sp->count; i++) {
		for (int j = 0; j < sp->size; j++) {
			//assert( sp->patterns[i][j] <=4 );
			printf("%d%c", sp->patterns[i][j], (j==sp->size-1 ? '\n' : ' ') );
			
		}
	}
}

double unconstrained_lk( const SitePattern *sp, const unsigned int count){
	double ulk = 0.0;
	
	for (int i = 0; i < sp->count; i++) {
		ulk += sp->weights[i] * log(sp->weights[i]);
	}
	
	ulk -= count*log(count);
	return ulk;
}

void _make_patterns( const Sequences *aln, uint8_t **patterns, int start, int length, int every ){
	int j = 0;
	for( int site = start; site < length; site+=every ){
		for( int i = 0; i < aln->size; i++ ){
			patterns[j][i] = aln->datatype->encoding(aln->datatype, aln->seqs[i]->seq[site]);
		}
		j++;
	}
}

void _make_patterns_generic( const Sequences *aln, uint8_t **patterns ){
    char *state = cvector(aln->length+1);
    state[aln->length] = '\0';
    for( int i = 0; i < aln->size; i++ ){
        memcpy(state, aln->seqs[i]->seq, aln->length*sizeof(char));
        patterns[0][i] = aln->datatype->encoding_string(aln->datatype, state);
    }
    free(state);
}


void _make_patterns_codon( const Sequences *aln, uint8_t **patterns, int start, int length ){
	int i;
	assert(length % 3 == 0);
	int pos = 0;
    
    DataType *nuctype = nucleotide_datatype();//SINGLETON_DATATYPE_NUCLEOTIDE;
    
	for( int site = start; site < length; site+=3 ){
		for( i = 0; i < aln->size; i++ ){
            int n1 = nuctype->encoding(nuctype, aln->seqs[i]->seq[site]);// nucleotide_to_encoding( aln->seqs[i]->seq[site] );
			int n2 = nuctype->encoding(nuctype, aln->seqs[i]->seq[site+1]);// nucleotide_to_encoding( aln->seqs[i]->seq[site+1] );
			int n3 = nuctype->encoding(nuctype, aln->seqs[i]->seq[site+2]);// nucleotide_to_encoding( aln->seqs[i]->seq[site+2] );
			
			/*if( GENETIC_CODE_TABLES[aln->datatype->genetic_code][n1*16+n2*4+n3] == '*' ) {
                fprintf(stderr, "%d %d %d %c%c%c\n",n1,n2,n3,aln->seqs[i]->seq[site],aln->seqs[i]->seq[site+1],aln->seqs[i]->seq[site+2]);
				fprintf(stderr, "_make_patterns_codon: There is a codon stop at codon position %d (%d) %c%c%c (%s)\n", pos,site, aln->datatype->state(aln->datatype,n1), aln->datatype->state(aln->datatype,n2), aln->datatype->state(aln->datatype,n3),aln->seqs[i]->name);
				exit(1);
			}*/
			
			if ( n1 > 3 || n2 > 3 || n3 > 3 ) {
				patterns[pos][i] = 65;
			}
			else {
				patterns[pos][i] = n1*16+n2*4+n3;
				
				for ( int j = 0; j < n1*16+n2*4+n3; j++ ) {
					if ( GENETIC_CODE_TABLES[aln->datatype->genetic_code][j] == '*' ) {
						patterns[pos][i]--;
					}
				}
				
			}
		}
		pos++;
	}
	
//	for( int j = 0; j < aln->size; j++ ){
//		for ( int i = 0; i < aln->length/3; i++ ) {
//			fprintf(stdout, "%d ", patterns[i][j]);
//		}
//		fprintf(stdout, "\n");
//	}
}

void SitePattern_sort( SitePattern *sitepattern, int *ordering ){
	uint8_t **patterns = sitepattern->patterns;
	uint8_t * temp = NULL;
	for ( int i = 0; i < sitepattern->count; i++ ) {
		if( i == ordering[i] ) continue;
		temp = patterns [i];
		patterns[i] = patterns[ ordering[i] ];
		patterns[ ordering[i] ] = temp;
	}
}

double ** SitePattern_rates( const SitePattern *sp ){
    double **mat = dmatrix(sp->nstate, sp->nstate);
    for ( int i = 0; i < sp->size; i++ ) {
        for ( int j = i+1; j < sp->size; j++ ) {
    
            for ( int k = 0; k < sp->count; k++ ) {
				int p1 = sp->patterns[k][i];
				int p2 = sp->patterns[k][j];
				if(p1 < sp->nstate && p2 < sp->nstate ){
					if( p1 != p2 ){
						mat[p1][p2] += sp->weights[k];
					}
					mat[p1][p1] += sp->weights[k];
				}
            }
        }
    }
	int sum = 0;
	for (int i = 0; i < sp->nstate; i++) {
		sum += mat[i][i];
	}
	for (int i = 0; i < sp->nstate; i++) {
		mat[i][i] /= sum;
	}
    return mat;
}

unsigned SitePattern_polymorphic_count(SitePattern *sp){
    unsigned polymorphisms = 0;
    for ( int i = 0; i < sp->count; i++ ) {
        uint8_t p = sp->patterns[i][0];
        int j = 1;
        for ( ; j < sp->size; j++ ) {
            if( p != sp->patterns[i][j] ) break;
        }
        if( j != sp->size ){
            polymorphisms += sp->weights[i];
        }
    }
    return polymorphisms;
}


SitePattern* new_SitePattern_from_json(json_node* node, Hashtable* hash){
	json_node* alignment_node = get_json_node(node, "alignment");
	json_node* datatype_node = get_json_node(node, "datatype");
	json_node* start_node = get_json_node(node, "every");
	json_node* length_node = get_json_node(node, "length");
	json_node* every_node = get_json_node(node, "every");
	SitePattern* patterns = NULL;
	
	
	if (alignment_node != NULL) {
		DataType* datatype = new_DataType_from_json(datatype_node, hash);
		Sequences* sequences = new_Sequences_from_json(alignment_node, hash);
		sequences->datatype = datatype;
		if(every_node != NULL || start_node != NULL || length_node != NULL){
			const char* every_string = (char*)every_node->value;
			int start = get_json_node_value_size_t(node, "start", 0);
			int length = get_json_node_value_size_t(node, "length", sequences->length);
			if(strlen(every_string) == 1){
				int every = atoi(every_string);
				patterns = new_SitePattern2(sequences, start, length, every);
			}
			else{
				char temp[2];
				temp[1] = '\0';
				temp[0] = every_string[0];
				int a = atoi(temp);
				temp[0] = every_string[1];
				int b = atoi(temp);
				patterns = new_SitePattern3(sequences, start, length, a, b);
			}
		}
		else{
			patterns = new_SitePattern(sequences);
		}
		free_Sequences(sequences);
	}
	else{
		exit(1);
	}
	return patterns;
}

