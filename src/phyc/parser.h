/*
 *  parser.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/8/11.
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

#ifndef _PARSER_H_
#define _PARSER_H_

#include "utils.h"

#include "mstring.h"



// positive values >1 are reserved for the number of cells in an array
typedef enum sml_token {
	SML_OPEN      = -1,
	SML_CLOSE     = -2, 
	SML_TAG       = -3,
	
	SML_STRING    = -4,
	SML_REFERENCE = -5,
	SML_NODE      = -6,
	SML_NUMBER    =  1,
	
	SML_EOS       = -9
} sml_token;

typedef struct Tokenizer{
	char *string;
	const char *pstr;
	StringBuffer *buffer;
}Tokenizer;


struct _SMLNode;
typedef struct _SMLNode *SMLNode;

struct _ObjectStore;
typedef struct _ObjectStore ObjectStore;


typedef void * (*pParser)(ObjectStore *, SMLNode);

typedef struct ParserStore{
	char    **name;
	pParser *parser;
	unsigned int count;
	unsigned int capacity;
}ParserStore;

#pragma mark -
#pragma mark ParserStore


ParserStore * new_ParserStore();

void free_ParserStore( ParserStore *ref );

ParserStore * ParserStore_add( ParserStore *ref, const pParser func, const char *name );

pParser ParserStore_get( const ParserStore *ref,  const char *name );

#pragma mark -
#pragma mark ObjectStore

ObjectStore * new_ObjectStore( const SMLNode root );

void free_ObjectStore( ObjectStore *store );

void ObjectStore_add( ObjectStore *store, char *id, void *obj );

void * ObjectStore_get_object( const ObjectStore *store, const char *id );

void * ObjectStore_get( ObjectStore *store, const char *id );

bool ObjectStore_exists( const ObjectStore *store, const char *id );

#pragma mark -
#pragma mark Tokenizer

Tokenizer * new_Tokenizer( const char *string );

void free_Tokenizer( Tokenizer *tokenizer );

sml_token sml_next_token( Tokenizer *tokenizer);

#pragma mark -
#pragma mark SMLNode

SMLNode new_SMLNode( const SMLNode parent );

void free_SMLNode( SMLNode node );

#pragma mark -
#pragma mark SMLTree

SMLNode new_SMLTree( Tokenizer *tok );

void free_SMLTree( SMLNode root);

char *SML_get_data( const SMLNode node );

char *SML_get_tag( const SMLNode node );

char *SML_get_data_of_child( SMLNode node, const char *tag );

SMLNode SML_get_element( SMLNode node, const char *tag );

SMLNode SML_get_element_by_index( SMLNode node, const int index );

int SML_get_child_count( const SMLNode node );

SMLNode SML_find_from_id( SMLNode node, const char *id );

void SML_insert_buffer( SMLNode node, const StringBuffer* buffer );

void SML_insert_string( SMLNode node, const char *string );

void print_SMLTree( FILE *pfile, const SMLNode root );

StringBuffer * SMLTree_to_buffer( SMLNode root, bool pretty );

double * sml_string_to_double_array( char *string, const unsigned int n );

#endif
