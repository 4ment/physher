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
#include <strings.h>
#include <math.h>
#include <ctype.h>

#include "matrix.h"
#include "mstring.h"
#include "sequenceio.h"
#include "mstring.h"
#include "geneticcode.h"
#include "random.h"
#include "hashfunctions.h"
#include "datatype.h"

static void _make_patterns( const Sequences *a, Hashtable* hash, int start, int end, int every );
static void _make_patterns_codon( const Sequences *aln, Hashtable* hash, int start, int end );
static void _make_patterns_generic( const Sequences *aln, Hashtable* hash, int start, int length );

static void _sitepattern_get_partials(const SitePattern* sp, size_t taxonIndex, double* partials);
static void _sitepattern_get_partials_from_partials(const SitePattern* sp, size_t taxonIndex, double* partials);

typedef struct double_vec_t{
	double* values;
	size_t length;
}double_vec_t;


unsigned int hashtable_hash_double_t(const void* data) {
	const double_vec_t* vec = data;
	const double* values = vec->values;
	unsigned int hash = values[0]/71;
	for(int i = 1; i < vec->length; i++){
//		hash += values[i]/71;
		hash += ELFHash((const char *)&values[i], sizeof(values[i]));
	}
	return hash;
}

static bool hashtable_compare_double_t( const void *key1, const void *key2 ){
	const double_vec_t* k1 = (double_vec_t*)key1;
	const double_vec_t* k2 = (double_vec_t*)key2;
	return memcmp(k1->values, k2->values, k1->length*sizeof(double)) == 0;
}

typedef struct uint8vec{
	uint8_t* values;
	size_t length;
}uint8vec;

unsigned int hashtable_hash_uint8_t(const void* data) {
	const uint8vec* vec = data;
	const uint8_t* values = vec->values;
	unsigned int hash = values[0];
	for(int i = 1; i < vec->length; i++){
		hash ^= values[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
	}
	return hash;
}

static bool hashtable_compare_uint8_t( const void *key1, const void *key2 ){
	const uint8vec* k1 = (uint8vec*)key1;
	const uint8vec* k2 = (uint8vec*)key2;
	return memcmp(k1->values, k2->values, k1->length*sizeof(uint8_t)) == 0;
}

SitePattern * new_SitePattern( const Sequences *aln ){
	return new_SitePattern2(aln, 0, aln->length, 1);
}

void _make_partials( double** partials, Hashtable* hash, int size, int length ){
	double_vec_t* pattern = malloc(sizeof(double_vec_t));
	pattern->length = size*4;
	pattern->values = malloc(sizeof(double)*pattern->length);
	
	for( int site = 0; site < length; site++ ){
		for( int i = 0; i < size; i++ ){
			memcpy(pattern->values + i*4, partials[i] + site*4, sizeof(double)*4);
		}
		HashEntry* entry = Hashtable_get_entry(hash, pattern);
		if (entry != NULL) {
			int* counter = HashEntry_value(entry);
			(*counter)++;
		}
		else{
			double_vec_t* vec = malloc(sizeof(double_vec_t));
			vec->length = pattern->length;
			vec->values = clone_dvector(pattern->values, pattern->length);
			Hashtable_add(hash, vec, new_Int(1));
		}
	}
	free(pattern->values);
	free(pattern);
}

// length: number of site in alignment
// size: number of sequences
SitePattern * new_SitePartials( char** names, double **partials, int length, int size, DataType* datatype, bool compress ){
	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
	assert(sp);
	sp->id = 0;
	
	sp->datatype = datatype;
	sp->datatype->ref_count++;
	sp->nstate = sp->datatype->state_count(sp->datatype);
	sp->get_partials = _sitepattern_get_partials_from_partials;
	sp->count = length;
	sp->nsites = length;
	sp->size = size;

	if(compress){
		Hashtable* hash = new_Hashtable( 100, &hashtable_hash_double_t, hashtable_compare_double_t, HASHTABLE_KEY_REFERENCE);
		hashtable_set_key_ownership(hash, true);
		hashtable_set_value_ownership(hash, true);
		_make_partials(partials, hash, size, length);
		
		int compressed_count = Hashtable_length(hash);
		sp->count = compressed_count;
		sp->indexes = ivector(compressed_count);
		sp->weights = dvector(compressed_count);
		
		Hashtable_init_iterator(hash);
		HashEntry *entry = NULL;
		
		double** partials2 = malloc(sizeof(double*)*size);
		for(size_t i = 0; i < size; i++){
			partials2[i] = malloc(sizeof(double)*compressed_count*4);
		}
		
		size_t index = 0;
		while ( (entry = Hashtable_next(hash) ) != NULL ) {
			const double_vec_t *vec = HashEntry_key(entry);
			const double *partial = vec->values;
			const int *counter = HashEntry_value(entry);
			sp->weights[index] = *counter;
			for(size_t i = 0; i < size; i++){
				memcpy(partials2[i] + index*4, partial + i*4, sizeof(double)*4);
			}
			sp->indexes[index] = -1;
			index++;
			free(vec->values);
		}
		for(int i = 0; i < size; i++)
			free(partials[i]);
		free(partials);
		free_Hashtable(hash);
			sp->partials = partials2;
	}
	else{
		sp->indexes = ivector(sp->count);
		sp->weights = dvector(sp->count);
		for (int i = 0; i < sp->count; i++) {
			sp->indexes[i] = i;
			sp->weights[i] = 1.0;
		}
		sp->partials = partials;
	}
	
	sp->patterns = NULL;
	sp->names = names;
	
	sp->ref_count = 1;
	return sp;
}

SitePattern * new_SitePattern2( const Sequences *aln, int start, int length, int every ){
	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
	assert(sp);
	sp->id = 0;
	
	sp->datatype = aln->datatype;
	sp->datatype->ref_count++;
	sp->nstate = sp->datatype->state_count(sp->datatype);
	sp->get_partials = _sitepattern_get_partials;
	sp->partials = NULL;
	
	Hashtable* hash = new_Hashtable( 100, &hashtable_hash_uint8_t, hashtable_compare_uint8_t, HASHTABLE_KEY_REFERENCE);
	hashtable_set_key_ownership(hash, true);
	hashtable_set_value_ownership(hash, true);
	
	if ( aln->datatype->type == DATA_TYPE_CODON) {
		sp->nsites = length/3;
		_make_patterns_codon(aln, hash, start, length);
	}
	else if ( aln->datatype->type == DATA_TYPE_GENERIC) {
		sp->nsites = aln->length/aln->datatype->symbolLength;
		_make_patterns_generic( aln, hash, start, length );
	}
	else {
		sp->nsites = length/every;
		_make_patterns( aln, hash, start, length, every );
	}
	
	sp->count = Hashtable_length(hash);
	sp->size = aln->size;
	sp->indexes = ivector(sp->count);
	sp->weights = dvector(sp->count);
	
	Hashtable_init_iterator(hash);
	HashEntry *entry = NULL;
	
	sp->patterns = malloc(sizeof(uint8_t*)*sp->size);
	for(size_t i = 0; i < aln->size; i++){
		sp->patterns[i] = malloc(sizeof(uint8_t)*sp->count);
	}
	
	size_t index = 0;
	while ( (entry = Hashtable_next(hash) ) != NULL ) {
		const uint8vec *vec = HashEntry_key(entry);
		const uint8_t *pattern = vec->values;
		const int *counter = HashEntry_value(entry);
		sp->weights[index] = *counter;
		for(size_t i = 0; i < aln->size; i++){
			sp->patterns[i][index] = pattern[i];
		}
		sp->indexes[index] = -1;
		index++;
		free(vec->values);
	}
	
	sp->names = (char**)malloc( aln->size * sizeof(char*) );
	assert(sp->names);
	
	for ( int i = 0; i < aln->size; i++) {
		sp->names[i] = String_clone( aln->seqs[i]->name );
	}
	free_Hashtable(hash);
	sp->ref_count = 1;
	
	return sp;
}

SitePattern * new_SitePattern3( const Sequences *aln, int start, int length, int a, int b ){
	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
	assert(sp);
	sp->id = 0;
	
	sp->datatype = aln->datatype;
	sp->datatype->ref_count++;
	sp->nstate = sp->datatype->state_count(sp->datatype);
	sp->get_partials = _sitepattern_get_partials;
	sp->partials = NULL;
	
	Hashtable* hash = new_Hashtable( 100, &hashtable_hash_uint8_t, hashtable_compare_uint8_t, HASHTABLE_KEY_REFERENCE);
	hashtable_set_key_ownership(hash, true);
	hashtable_set_value_ownership(hash, true);
	
	if ( aln->datatype->type == DATA_TYPE_CODON) {
		sp->nsites = length/3;
		_make_patterns_codon(aln, hash, start, length);
	}
	else if ( aln->datatype->type == DATA_TYPE_GENERIC) {
		sp->nsites = aln->length/aln->datatype->symbolLength;
		_make_patterns_generic( aln, hash, start, length );
	}
	else {
		sp->nsites = length/3*2;
		_make_patterns( aln, hash, start+a, length, 3 );
		_make_patterns( aln, hash, start+b, length, 3 );
	}
	
	sp->count = Hashtable_length(hash);
	sp->size = aln->size;
	sp->indexes = ivector(sp->count);
	sp->weights = dvector(sp->count);
	
	Hashtable_init_iterator(hash);
	HashEntry *entry = NULL;
	
	sp->patterns = malloc(sizeof(uint8_t*)*sp->size);
	for(size_t i = 0; i < aln->size; i++){
		sp->patterns[i] = malloc(sizeof(uint8_t)*sp->count);
	}
	
	size_t index = 0;
	while ( (entry = Hashtable_next(hash) ) != NULL ) {
		const uint8vec *vec = HashEntry_key(entry);
		const uint8_t *pattern = vec->values;
		const int *counter = HashEntry_value(entry);
		sp->weights[index] = *counter;
		for(size_t i = 0; i < aln->size; i++){
			sp->patterns[i][index] = pattern[i];
		}
		sp->indexes[index] = -1;
		index++;
		free(vec->values);
	}
	
	sp->names = (char**)malloc( aln->size * sizeof(char*) );
	assert(sp->names);
	
	for ( int i = 0; i < aln->size; i++) {
		sp->names[i] = String_clone( aln->seqs[i]->name );
	}
	free_Hashtable(hash);
	sp->ref_count = 1;
	
	return sp;
}

void free_SitePattern( SitePattern *sp ){
	if(sp->ref_count == 1){
		if(sp->patterns != NULL ){
			free_ui8matrix(sp->patterns, sp->size);
			free(sp->weights);
		}
		if(sp->partials != NULL ){
			free_dmatrix(sp->partials, sp->size);
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
    newsp->datatype = sp->datatype;
	newsp->datatype->ref_count++;
	newsp->get_partials = sp->get_partials;
	newsp->patterns = NULL;
	newsp->partials = NULL;
	
	if(sp->patterns != NULL){
		newsp->patterns = ui8matrix(sp->size, sp->count);
		for ( int p = 0; p < sp->size; p++ ) {
			memcpy(newsp->patterns[p], sp->patterns[p], sp->count * sizeof(uint8_t));
		}
	}
    else if(sp->partials != NULL){
		newsp->partials = dmatrix(sp->size, sp->count * sp->nstate);
		for ( int p = 0; p < sp->size; p++ ) {
			memcpy(newsp->partials[p], sp->partials[p], sp->count * sp->nstate * sizeof(double));
		}
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


//SitePattern * SitePattern_subset( const SitePattern *original, const int *indexes, int size ){
//	SitePattern *sp = (SitePattern *)malloc( sizeof(SitePattern) );
//	assert(sp);
//	sp->id = original->id+1;
//	sp->datatype = original->datatype;
//	sp->get_partials = _sitepattern_get_partials;
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
	sp->get_partials = _sitepattern_get_partials;
	
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
			if(memcmp(sitePattern2->patterns[site], sp->patterns[p], sizeof(uint8_t)*sp->size) == 0){
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
		sps[i]->get_partials = _sitepattern_get_partials;
		
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


int get_sequence_index( const SitePattern *sp, const char *name ){
	int i;
	for (i = 0; i < sp->size; i++) {
		if( strcmp(name, sp->names[i] ) == 0 ) return i;
	}
	return -1;
}

void _sitepattern_get_partials(const SitePattern* sp, size_t taxonIndex, double* partials){
	size_t patternLength = sp->count;
	uint8_t* taxonPatterns = sp->patterns[taxonIndex];
	double* pPartials = partials;
	for (size_t i = 0; i < patternLength; i++) {
		sp->datatype->partial(sp->datatype, taxonPatterns[i], pPartials);
		pPartials += sp->datatype->stateCount;
	}
}

void _sitepattern_get_partials_from_partials(const SitePattern* sp, size_t taxonIndex, double* partials){
	memcpy(partials, sp->partials[taxonIndex], sizeof(double)*sp->count*sp->datatype->stateCount);
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

void _make_patterns( const Sequences *aln, Hashtable* hash, int start, int length, int every ){
	uint8vec* pattern = malloc(sizeof(uint8vec));
	pattern->length = aln->size;
	pattern->values = malloc(sizeof(uint8_t)*pattern->length);
	
	for( int site = start; site < length; site+=every ){
		for( int i = 0; i < aln->size; i++ ){
			pattern->values[i] = aln->datatype->encoding(aln->datatype, aln->seqs[i]->seq[site]);
		}
		HashEntry* entry = Hashtable_get_entry(hash, pattern);
		if (entry != NULL) {
			int* counter = HashEntry_value(entry);
			(*counter)++;
		}
		else{
			uint8vec* vec = malloc(sizeof(uint8vec));
			vec->length = pattern->length;
			vec->values = clone_u8ivector(pattern->values, pattern->length);
			Hashtable_add(hash, vec, new_Int(1));
		}
	}
	free(pattern->values);
	free(pattern);
}

void _make_patterns_generic( const Sequences *aln, Hashtable* hash, int start, int length ){
	size_t symbolLength = aln->datatype->symbolLength;
	uint8vec* pattern = malloc(sizeof(uint8vec));
	pattern->length = aln->size;
	pattern->values = malloc(sizeof(uint8_t)*pattern->length);
    char *state = cvector(symbolLength+1);
    state[symbolLength] = '\0';
	
	for( int site = start; site < length; site+=symbolLength ){
		for( int i = 0; i < aln->size; i++ ){
			memcpy(state, aln->seqs[i]->seq+site, symbolLength*sizeof(char));
			pattern->values[i] = aln->datatype->encoding_string(aln->datatype, state);
		}
		HashEntry* entry = Hashtable_get_entry(hash, pattern);
		if (entry != NULL) {
			int* counter = HashEntry_value(entry);
			(*counter)++;
		}
		else{
			uint8vec* vec = malloc(sizeof(uint8vec));
			vec->length = pattern->length;
			vec->values = clone_u8ivector(pattern->values, pattern->length);
			Hashtable_add(hash, vec, new_Int(1));
		}
	}
	free(state);
	free(pattern->values);
	free(pattern);
}


void _make_patterns_codon( const Sequences *aln, Hashtable* hash, int start, int length ){
	int i;
	assert(length % 3 == 0);
	uint8vec* pattern = malloc(sizeof(uint8vec));
	pattern->length = aln->size;
	pattern->values = malloc(sizeof(uint8_t)*pattern->length);
    
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
				pattern->values[i] = 65;
			}
			else {
				pattern->values[i] = n1*16+n2*4+n3;
				
				for ( int j = 0; j < n1*16+n2*4+n3; j++ ) {
					if ( GENETIC_CODE_TABLES[aln->datatype->genetic_code][j] == '*' ) {
						pattern->values[i]--;
					}
				}
				
			}
		}
		HashEntry* entry = Hashtable_get_entry(hash, pattern);
		if (entry != NULL) {
			int* counter = HashEntry_value(entry);
			(*counter)++;
		}
		else{
			uint8vec* vec = malloc(sizeof(uint8vec));
			vec->length = pattern->length;
			vec->values = clone_u8ivector(pattern->values, pattern->length);
			Hashtable_add(hash, vec, new_Int(1));
		}
	}
	free(pattern->values);
	free(pattern);
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
				int p1 = sp->patterns[i][k];
				int p2 = sp->patterns[j][k];
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
        uint8_t p = sp->patterns[0][i];
		int j = 0;
		while (j < sp->size && sp->patterns[j][i] >= sp->nstate) {
			p = sp->patterns[j][i];
			j++;
		}
        for ( ; j < sp->size; j++ ) {
            if( p != sp->patterns[j][i] && sp->patterns[j][i] < sp->nstate) break;
        }
        if( j != sp->size ){
            polymorphisms += sp->weights[i];
        }
    }
    return polymorphisms;
}

DataType* parse_datatype(json_node* datatype_node, Hashtable* hash){
	DataType* datatype = NULL;
	if(datatype_node->node_type == MJSON_STRING && (strcasecmp((char*)datatype_node->value, "nucleotide") == 0 ||
	   strcasecmp((char*)datatype_node->value, "codon") == 0 || strcasecmp((char*)datatype_node->value, "aa") == 0)){
		datatype = new_DataType_from_json(datatype_node, hash);
	}
	else if (datatype_node->node_type == MJSON_STRING && Hashtable_exists(hash, (char*)datatype_node->value)) {
		datatype = Hashtable_get(hash, (char*)datatype_node->value);
		datatype->ref_count++;
	}
	else{
		datatype = new_DataType_from_json(datatype_node, hash);
		Hashtable_add(hash, datatype->name, datatype);
	}
	return datatype;
}

SitePattern* new_SitePattern_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"alignment",
		"datatype",
		"every",
		"length",
		"partials",
		"start",
		"verbose"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* alignment_node = get_json_node(node, "alignment");
	json_node* partials_node = get_json_node(node, "partials");
	json_node* datatype_node = get_json_node(node, "datatype");
	json_node* start_node = get_json_node(node, "start");
	json_node* length_node = get_json_node(node, "length");
	json_node* every_node = get_json_node(node, "every");
	SitePattern* patterns = NULL;
	
	DataType* datatype = NULL;
	if(datatype_node->node_type == MJSON_STRING && (strcasecmp((char*)datatype_node->value, "nucleotide") == 0 ||
	   strcasecmp((char*)datatype_node->value, "codon") == 0 || strcasecmp((char*)datatype_node->value, "aa") == 0)){
		datatype = new_DataType_from_json(datatype_node, hash);
	}
	else if (datatype_node->node_type == MJSON_STRING && Hashtable_exists(hash, (char*)datatype_node->value)) {
		datatype = Hashtable_get(hash, (char*)datatype_node->value);
		datatype->ref_count++;
	}
	else{
		datatype = new_DataType_from_json(datatype_node, hash);
		Hashtable_add(hash, datatype->name, datatype);
	}
	
	if (alignment_node != NULL) {
		Sequences* sequences = new_Sequences_from_json(alignment_node, hash);
		sequences->datatype = datatype;
		datatype->ref_count++;
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
		
		json_node* seq_node = get_json_node(alignment_node, "sequences");
		if(seq_node != NULL){
			size_t j = 0;
			for ( ; j < alignment_node->child_count; j++) {
				if(alignment_node->children[j] == seq_node) break;
			}
			for ( ; j < alignment_node->child_count; j++) {
				alignment_node->children[j] = alignment_node->children[j+1];
			}
			alignment_node->child_count--;
			json_free_tree(seq_node);
		}
	}
	else if(partials_node != NULL){
		double** partials = malloc(partials_node->child_count*sizeof(double*));
		char** names = malloc(partials_node->child_count*sizeof(char*));
		int length = 0;
		for (int i = 0; i < partials_node->child_count; i++) {
			json_node* child = partials_node->children[i];
			names[i] = String_clone((char*)child->key);
			char* sequence = child->value;
			int l = 0;
			String_trim(sequence);
			partials[i] = String_split_char_double( sequence, ',', &l );
			length = l;
		}
		bool compress = get_json_node_value_bool(node, "compress", true);
		patterns = new_SitePartials(names, partials, length/datatype->stateCount, partials_node->child_count, datatype, compress);
	}
	else{
		fprintf(stderr, "No `alignment' provided in sitepattern object\n");
		exit(1);
	}
	free_DataType(datatype);
	
	bool verbose = get_json_node_value_bool(node, "verbose", true);
	if(verbose){
		double unc_lk = unconstrained_lk(patterns, patterns->nsites);
		fprintf(stdout, "Number of sequences: %d\n", patterns->size );
		fprintf(stdout, "Alignment length: %d\n", patterns->nsites );
		if(patterns->patterns != NULL){
			int polymorphisms = SitePattern_polymorphic_count(patterns);
			fprintf(stdout, "Number of polymorphic sites: %d/%d (%f)\n", polymorphisms, patterns->nsites,((double)polymorphisms/patterns->nsites) );
		}
		fprintf(stdout, "Number of patterns: %d\n", patterns->count );
		fprintf(stdout, "Unconstrained likelihood: %f\n\n", unc_lk );
	}
	return patterns;
}

