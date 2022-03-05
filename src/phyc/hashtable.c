/*
 *  hashtable.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/5/10.
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


#include "hashtable.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "matrix.h"
#include "hashfunctions.h"
#include "mstring.h"

#define DEFAULT_LOAD_FACTOR 0.65

struct _HashEntry{
    void *key;
	void *value;
    unsigned int hash; // used to compare quickly entries in the same bucket
	HashEntry *next;
};

struct _Hashtable {
    HashEntry **table;
	hashtable_key_type type;
    unsigned int size;   // Capacity of the hashtable
    unsigned int length; // Number of entries
    unsigned int loadlimit;
    unsigned int minloadlimit;
    unsigned int primeindex;
	double load_factor;
	hashtable_hashfn hashfn;
    hashtable_cmpfn cmpfn;
	
	bool own_key;
	bool own_value;
	bool override_value;
	
	HashEntry *iterator;
	HashEntry *previous;
	unsigned iter_pos;
};

static const unsigned int primes[] = {
	5,
	53, 97, 193, 389,
	769, 1543, 3079, 6151,
	12289, 24593, 49157, 98317,
	196613, 393241, 786433, 1572869,
	3145739, 6291469, 12582917, 25165843,
	50331653, 100663319, 201326611, 402653189,
	805306457, 1610612741
};
const unsigned int prime_table_length = sizeof(primes)/sizeof(primes[0]);


/*static int *generate_primes( int n ){
	int *p = ivector(n);
	int i,j,count=0;
	for ( i = 2; i <n; j++ ) {
		for ( j = 2; j < i; i++ ) {
			if( i % j == 0) break;;
		}
		if( i == j ) p[count++] = i;
	}
	return p;
}*/

static inline unsigned int hash_to_index(unsigned int hash_size, unsigned int hashvalue) {
    return (hashvalue % hash_size);
}

static bool Hashtable_compare_references(const void *x, const void *y) {
	return x != y;
}

static bool Hashtable_compare_strings( const void *key1, const void *key2 ){
	const char* k1 = (char*)key1;
	const char* k2 = (char*)key2;
	if ( strcmp( k1, k2 ) == 0 ) return true;
	return false;
}

static bool Hashtable_compare_ints( const void *key1, const void *key2 ){
	const int* k1 = (int*)key1;
	const int* k2 = (int*)key2;
	return ( k1 == k2 ? true : false );
}

static bool Hashtable_compare_doubles( const void *key1, const void *key2 ){
	const double* k1 = (double*)key1;
	const double* k2 = (double*)key2;
	return ( k1 == k2 ? true : false );
}

static bool Hashtable_compare_floats( const void *key1, const void *key2 ){
	const float* k1 = (float*)key1;
	const float* k2 = (float*)key2;
	return ( k1 == k2 ? true : false );
}

Hashtable * new_Hashtable( const unsigned int size, hashtable_hashfn hashfn, hashtable_cmpfn cmpfn, hashtable_key_type type ) {
    Hashtable *hash = NULL;
    unsigned int pindex;
    
    if ( size > primes[prime_table_length-1] ) return NULL;
    
    for ( pindex = 0; primes[pindex] < size; pindex++);
	
    hash = (Hashtable *)malloc(sizeof(struct _Hashtable));
    assert(hash);
	
    hash->table = (HashEntry **)calloc( primes[pindex], sizeof(HashEntry*) );
    assert(hash->table);
	
    hash->loadlimit    = (unsigned int) ceil(primes[pindex] * DEFAULT_LOAD_FACTOR);
    hash->minloadlimit = 0;
	hash->size  = primes[pindex];
    hash->primeindex   = pindex;
    hash->length   = 0;
	hash->load_factor = DEFAULT_LOAD_FACTOR;
	hash->type = type;
    hash->hashfn       = hashfn;
    hash->cmpfn        = cmpfn;
	
	hash->own_key   = true;
	hash->own_value = true;
	hash->override_value = true;
	
	hash->iterator = hash->table[0];
	hash->previous = NULL;
	hash->iter_pos = 0;
    return hash;
}


Hashtable * new_Hashtable_string( unsigned int size ){
	return new_Hashtable( size, &JSHash2, Hashtable_compare_strings, HASHTABLE_KEY_STRING);
}

/*! Free a hashtable.
 */
void free_Hashtable( Hashtable *hash ){
	unsigned int i;
	HashEntry *e, *f;
	HashEntry **table = hash->table;
	
	for (i = 0; i < hash->size; i++){
		e = table[i];
		while ( e != NULL ){
			f = e;
			e = e->next;

            if ( hash->own_key){
				free(f->key);
				f->key = NULL; 
			}
            if ( hash->own_value){
				free(f->value);
				f->value = NULL;
			}
			
			free(f); f = NULL;
		}
	}
    
    free(hash->table); hash->table = NULL;
    free(hash); hash = NULL;
}

static unsigned int hashfn( Hashtable *hash, const void *k){
    /* Aim to protect against poor hash functions by adding logic here
     * - logic taken from java 1.4 hashtable source */
	unsigned int i = hash->hashfn(k);
    i += ~(i << 9);
    i ^=  ((i >> 14) | (i << 18)); /* >>> */
    i +=  (i << 4);
    i ^=  ((i >> 10) | (i << 22)); /* >>> */
    return i;
}

static bool Hashtable_expand( Hashtable *hash ){
    HashEntry **newtable;
    HashEntry  *e = NULL;
    HashEntry  **pE = NULL;
    unsigned int newsize, i, index;

    if (hash->primeindex > prime_table_length ) return false;
    newsize = primes[++(hash->primeindex)];
	
    newtable = (HashEntry **)malloc(sizeof(HashEntry*) * newsize);
    if (NULL != newtable){
        memset(newtable, 0, newsize * sizeof(HashEntry*));
        /* This algorithm is not 'stable'. ie. it reverses the list
         * when it transfers entries between the tables */
        for (i = 0; i < hash->size; i++) {
            while (NULL != (e = hash->table[i])) {
                hash->table[i] = e->next;
                index = hash_to_index(newsize,e->hash);
                e->next = newtable[index];
                newtable[index] = e;
            }
        }
        free(hash->table);
        hash->table = newtable;
    }
    else {
        newtable = (HashEntry **)realloc(hash->table, newsize * sizeof(HashEntry *));
        if (NULL == newtable) {
			(hash->primeindex)--; 
			return false;
		}
        hash->table = newtable;
        memset(newtable[hash->size], 0, newsize - hash->size);
        for (i = 0; i < hash->size; i++) {
            for (pE = &(newtable[i]), e = *pE; e != NULL; e = *pE) {
                index = hash_to_index(newsize,e->hash);
                if (index == i) {
                    pE = &(e->next);
                }
                else {
                    *pE = e->next;
                    e->next = newtable[index];
                    newtable[index] = e;
                }
            }
        }
    }
    hash->size = newsize;
    hash->loadlimit   = (unsigned int) ceil(newsize * hash->load_factor);
    return true;
}

bool Hashtable_add_array_unsigned( Hashtable *hash, int *key, const unsigned int n, void *value ){
	unsigned *key_copy = uivector(n);
	memcpy( key_copy, key, n * sizeof(unsigned) );
	return Hashtable_add( hash, key_copy, value);
}

/*! Add a key with its value
 * If the key already exists and override_value and own_value are true then the old value is freed and the new value is used.
 * Otherwise the value is ignored. The user should check if the return function if override_value==false and own_value==true in order to avoid leaks
 * IT is a bit fucked here there is a difference in ownership for keys and values
 */
bool Hashtable_add( Hashtable *hash, const void *key, void *value ){
    unsigned int index;
    HashEntry *entry = NULL;
    if ( hash->length == hash->loadlimit ){
        Hashtable_expand(hash);
    }
	
	HashEntry *p = Hashtable_get_entry( hash, key );
	if( p != NULL){
		if ( hash->override_value && hash->own_value ) {
			free(p->value);
			p->value = value;
		}
		else if( hash->override_value ){
			p->value = value;
		}
		return true;
	}
	
	entry = (HashEntry *)malloc(sizeof(HashEntry));
	assert(entry);
//	if( hash->own_key ) {
//		switch ( hash->type ) {
//			case HASHTABLE_KEY_STRING:{
//				char *key_string = (char*)key;
//				entry->key = String_clone(key_string);
//				break;
//			}
//				
//			default:
//				entry->key = key;
//				break;
//		}
//	}
//	else {
		entry->key = key;
//	}

	entry->value = value;
	entry->hash  = hashfn(hash,key);
	index = hash_to_index(hash->size,entry->hash);
	entry->next  = hash->table[index];
	hash->table[index] = entry;
	hash->length++;
	
    return false;
}


/*! Get entry associated with key from Hashtable
 * return entry
 */
HashEntry * Hashtable_get_entry( Hashtable *hash, const void *key ){
	HashEntry *entry = NULL;
    unsigned int hashvalue = hashfn(hash,key);
    unsigned int index = hash_to_index(hash->size,hashvalue);
    
	for ( entry = hash->table[index]; entry ; entry = entry->next ){
        if ( hashvalue == entry->hash && hash->cmpfn(key, entry->key) ) return entry;
	}
    return NULL;
}

/*! Get value associated with key from Hashtable
 * return reference of value
 */
void * Hashtable_get( Hashtable *hash, const void *key ){
	HashEntry *entry = (HashEntry *)Hashtable_get_entry( hash, key );
    return (entry == NULL ? NULL : entry->value);
}

static void Hashtable_compact( Hashtable *hash ){
	
}

/*! Remove entry from Hashtable
 */
void Hashtable_remove( Hashtable *hash, void *key ){
    HashEntry **pe  = NULL;
	HashEntry *e    = NULL;
	
	if( hash->length < hash->minloadlimit ){
		Hashtable_compact( hash );
	}
	
	unsigned int hashvalue = hashfn(hash,key);
	unsigned int  index = hash_to_index(hash->size, hashfn(hash,key));
    
	for ( pe = &(hash->table[index]); *pe; pe = &(*pe)->next ) {
        if ( hashvalue == (*pe)->hash && hash->cmpfn(key, (*pe)->key) ){
			e = *pe;
            *pe = e->next;
            if ( hash->own_value) free(e->value);
            free(e);
            hash->length--;
            break;
        }
    }
}

void Hashtable_empty( Hashtable *hash ){
	unsigned int i;
	HashEntry *e, *f;
	HashEntry **table = hash->table;
	
	for (i = 0; i < hash->size; i++){
		e = table[i];
		while ( e != NULL ){
			f = e;
			e = e->next;
			
            if ( hash->own_key){
				free(f->key);
				f->key = NULL; 
			}
            if ( hash->own_value){
				free(f->value);
				f->value = NULL;
			}
			
			free(f); f = NULL;
		}
        hash->table[i] = NULL;
	}
	hash->length = 0;
}


bool Hashtable_exists( Hashtable *hash, const void *key ){
	return (Hashtable_get(hash, key) == NULL ? false : true);
}

unsigned int Hashtable_size( Hashtable *hash ){
	return hash->size;
}

unsigned int Hashtable_length( Hashtable *hash ){
	return hash->length;
}

void hashtable_set_key_ownership( Hashtable *hash, bool ownership ){
	hash->own_key = ownership;
}

void hashtable_set_value_ownership( Hashtable *hash, bool ownership ){
	hash->own_value = ownership;
}

void hashtable_set_override( Hashtable *hash, bool override ){
	hash->override_value = override;
}

void Hashtable_init_iterator( Hashtable *hash ){
	if ( hash->length == 0 ) {
		return;
	}
	hash->iter_pos = 0;
	while (hash->table[hash->iter_pos] == NULL) {
		hash->iter_pos++;
	}
	hash->iterator = hash->table[hash->iter_pos];
	hash->previous = NULL;
}

HashEntry * Hashtable_next( Hashtable *hash ){
	if ( hash->length == 0 ) {
		hash->iterator = NULL;
	}
	else if ( hash->previous == NULL ) {
		hash->previous = hash->iterator;
	}
	else {
		hash->previous = hash->iterator;
		if ( hash->iterator->next == NULL) {
			
			while ( ++hash->iter_pos != hash->size && hash->table[hash->iter_pos] == NULL  ) {}
			
			if ( hash->iter_pos == hash->size ){
				hash->iterator = NULL;
				return NULL;
			}
			hash->iterator = hash->table[hash->iter_pos];
		}
		else {
			hash->iterator = hash->iterator->next;
		}
	}
	
	return hash->iterator;
}



void print_Hashtable( FILE *pf, const Hashtable *hash ){
	unsigned int i;
	HashEntry *e, *f;
	HashEntry **table = hash->table;
	
	for ( i = 0; i < hash->size; i++){
		e = table[i];
		fprintf(pf, "\n%d ----------------------------------------------\n", i);
		while ( e != NULL ){
			f = e;
			e = e->next;
		
            //fprintf(pf, "%s %f %d\n", (char*)f->key, *((double*)f->value), f->hash);
			//fprintf(pf, "%s %d %d\n", (char*)f->key, *((int*)f->value), f->hash);
			fprintf(pf, "%s %s %d\n", (char*)f->key, (char*)f->value, f->hash);
		}
	}
}

double * new_Double( const double value){
	double *p = (double *)malloc( sizeof(double) );
	assert(p);
	*p = value;
	return p;
}

int * new_Int( const int value){
	int *p = (int *)malloc( sizeof(int) );
	assert(p);
	*p = value;
	return p;
}

float * new_Float( const float value){
	float *p = (float *)malloc( sizeof(float) );
	assert(p);
	*p = value;
	return p;
}


void * HashEntry_value( const HashEntry *entry ){
	return entry->value;
}

const void * HashEntry_key( const HashEntry *entry ){
	return entry->key;
}
