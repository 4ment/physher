/*
 *  hashtable.h
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

#ifndef _HASHTABLE_H_
#define _HASHTABLE_H_

#include "utils.h"

typedef enum hashtable_key_type{
	HASHTABLE_KEY_STRING,
	HASHTABLE_KEY_INT,
	HASHTABLE_KEY_DOUBLE,
	HASHTABLE_KEY_FLOAT,
	HASHTABLE_KEY_REFERENCE
} hashtable_key_type;

typedef unsigned int (*hashtable_hashfn)( const char *, const unsigned int );

typedef bool (*hashtable_cmpfn) (const void*, const void*);

typedef struct _HashEntry HashEntry;

typedef struct _Hashtable Hashtable;



Hashtable * new_Hashtable( const unsigned int size, hashtable_hashfn hashfn, hashtable_cmpfn cmpfn, hashtable_key_type type );

Hashtable * new_Hashtable_string( unsigned int size );

Hashtable * new_Hashtable_int( unsigned int size );

void free_Hashtable( Hashtable *hash );



bool Hashtable_add( Hashtable *hash, void *key, void *value );

void * Hashtable_get( Hashtable *hash, const void *key );

HashEntry * Hashtable_get_entry( Hashtable *hash, const void *key );

void Hashtable_remove( Hashtable *hash, void *key );

void Hashtable_empty( Hashtable *hash );

bool Hashtable_exists( Hashtable *hash, const void *key );

unsigned int Hashtable_size( Hashtable *hash );

unsigned int Hashtable_length( Hashtable *hash );

void hashtable_set_key_ownership( Hashtable *hash, bool ownership );

void hashtable_set_value_ownership( Hashtable *hash, bool ownership );

void hashtable_set_override( Hashtable *hash, bool override );

double * new_Double( const double value);

int * new_Int( const int value);

float * new_Float( const float value);


void print_Hashtable( FILE *pf, const Hashtable *hash );

HashEntry * Hashtable_next( Hashtable *hash );

void Hashtable_init_iterator( Hashtable *hash );

const void * HashEntry_value( const HashEntry *entry );

const void * HashEntry_key( const HashEntry *entry );

#endif
