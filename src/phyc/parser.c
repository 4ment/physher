/*
 *  parser.c
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

#include "parser.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mstring.h"
#include "matrix.h"
#include "hashtable.h"

// for parsers
#include "tree.h"
#include "sitepattern.h"

struct _SMLNode{
	char *tag;
	int type;
	struct _SMLNode **child;
	struct _SMLNode *parent;
	int child_count;
	int child_capacity;
	char *data;
};

// could be more efficient than an arry but the size should be small
struct _ObjectStore{
	Hashtable *hash;
	SMLNode root;
	ParserStore *p;
};

#pragma mark ParserStore


ParserStore * new_ParserStore(){
	ParserStore *store = (ParserStore*)malloc(sizeof(ParserStore) ) ;
	assert(store);
	store->capacity = 5;
	store->count = 0;
	
	store->parser = (pParser*)malloc( store->capacity * sizeof(pParser) );
	assert(store->parser);
	
	store->name = (char**)malloc( store->capacity * sizeof(char*) );
	assert(store->name);
	
	ParserStore_add( store, Tree_SML_to_object, "Tree" );
	ParserStore_add( store, SitePattern_SML_to_object, "SitePattern");
	
	return store;
}

void free_ParserStore( ParserStore *ref ){
	for (int i = 0; i < ref->count; i++) {	
		free(ref->name[i]);
	}
	free(ref->name);
	free(ref);
}

ParserStore * ParserStore_add( ParserStore *ref, pParser const func, const char *name ){
	if( ref->count == ref->capacity ){
		++ref->capacity;
		ref->parser = realloc( ref->parser, ref->capacity * sizeof(pParser) );
		assert(ref->parser);
		ref->name = realloc( ref->name, ref->capacity * sizeof(char*) );
		
	}
	
	ref->parser[ref->count]   = func;
	
	ref->name[ref->count] = (char*)malloc( (strlen(name)+1) * sizeof(char) );
	assert(ref->name[ref->count]);
	strcpy( ref->name[ref->count], name);
	
	ref->count++;
	return ref;
}

pParser ParserStore_get( const ParserStore *ref,  const char *name ){
	for (int i = 0; i < ref->count; i++) {
		if ( strcmp( ref->name[i], name ) == 0 ){
			return ref->parser[i];
		}
	}
	return NULL;
}



#pragma mark -
#pragma mark ObjectStore

static void _checkObjectStore( const SMLNode node, int *count ){
	if ( node != NULL ) {
		if ( node->type == SML_REFERENCE ) 
			(*count)++;
		
		for (int i = 0; i < node->child_count; i++) {
			_checkObjectStore( node->child[i], count );
		}
	}
}

ObjectStore * new_ObjectStore( const SMLNode root ){
	int capacity = 0;
	_checkObjectStore( root, &capacity);
	ObjectStore *store = (ObjectStore*)malloc( sizeof(ObjectStore) );
	assert(store);
	
	
	store->hash = new_Hashtable_string(capacity);
	hashtable_set_value_ownership( store->hash, false );
	
	store->root = root;
	
	store->p = new_ParserStore();
	
	return store;
}

void free_ObjectStore( ObjectStore *store ){
	free_Hashtable( store->hash );
	
	if( store->p != NULL ) free_ParserStore( store->p );
	
	free( store );
	store = NULL;
}

// should check if the hash work, i have changed the Hashtable_add function
void ObjectStore_add( ObjectStore *store, char *id, void *obj ){
	bool duplicate = Hashtable_add( store->hash, id, obj);
	if ( duplicate ) error("Duplicates in ObjectStore\n");
}


void * ObjectStore_get_object( const ObjectStore *store, const char *id ){
	return Hashtable_get( store->hash, id );
}

void * ObjectStore_get( ObjectStore *store, const char *id ){
	void * obj = Hashtable_get( store->hash, id );
	
	if ( obj == NULL ) {
		SMLNode node = SML_find_from_id( store->root, id );
		assert(node);
		pParser parser = ParserStore_get( store->p, SML_get_tag(node) );
		assert(parser);
		obj = (*parser)( store, node);
		assert(obj);
	}
	return obj;
}

bool ObjectStore_exists( const ObjectStore *store, const char *id ){
	void * obj = Hashtable_get( store->hash, id );
	return ( obj == NULL ? true : false );
}

#pragma mark -
#pragma mark Tokenizer

Tokenizer * new_Tokenizer( const char *string ){
	Tokenizer *tok = (Tokenizer*)malloc( sizeof(Tokenizer) );
	assert(tok);
	tok->string = (char*)malloc( (strlen(string)+1) * sizeof(char) );
	assert(tok->string);
	strcpy( tok->string, string );
	tok->pstr = &tok->string[0];
	tok->buffer = new_StringBuffer( 1000 );
	return tok;
}

void free_Tokenizer( Tokenizer *tokenizer ){
	free_StringBuffer( tokenizer->buffer );
	free(tokenizer->string);
	tokenizer->string = NULL;
	tokenizer->pstr = NULL;
	free(tokenizer);
	tokenizer = NULL;
}

sml_token sml_next_token( Tokenizer *tok ){
	StringBuffer_empty(tok->buffer);
	
	while ( tok->pstr[0] != '\0' ) {
		
		if( tok->pstr[0] == '\"'){
			tok->pstr++;
			while ( tok->pstr[0] != '\"' ){
				tok->buffer = StringBuffer_append_char(tok->buffer, tok->pstr[0]);
				tok->pstr++;
			}
			tok->pstr++;
			return SML_STRING;
		}
		else if( tok->pstr[0] == '(' ){
			tok->buffer = StringBuffer_append_char(tok->buffer, tok->pstr[0]);
			tok->pstr++;
			return SML_OPEN;
		}
		// TAG
		else if( (tok->pstr[0] >= 'a' && tok->pstr[0] <= 'z') || (tok->pstr[0] >= 'A' && tok->pstr[0] <= 'Z') ){
			while ( tok->pstr[0] != ':'  ) {
				tok->buffer = StringBuffer_append_char(tok->buffer, tok->pstr[0]);
				tok->pstr++;
			}
			tok->pstr++;
			return SML_TAG;
		}
		// Number or array of numbers
		else if ( (tok->pstr[0] >= 48 && tok->pstr[0] <= 57) || tok->pstr[0] == '.' || tok->pstr[0] == '-' ){
			int count = 0;
			while ( tok->pstr[0] != ')' ){
				if(count != 0) tok->buffer = StringBuffer_append_char(tok->buffer, ' ');
				while( (tok->pstr[0] >= 48 && tok->pstr[0] <= 57) || tok->pstr[0] == '.' || tok->pstr[0] == '-'
					  || tok->pstr[0] == 'e' || tok->pstr[0] == 'E' ){
					tok->buffer = StringBuffer_append_char(tok->buffer, tok->pstr[0]);
					tok->pstr++;
				}
				count++;
				while ( tok->pstr[0] == 9 || tok->pstr[0] == 32 ){tok->pstr++;}
			}
			return count;
		}
		else if( tok->pstr[0] == ')' ){
			tok->buffer = StringBuffer_append_char(tok->buffer, tok->pstr[0]);
			tok->pstr++;
			return SML_CLOSE;
		}
		// Reference
		else if( tok->pstr[0] == '&' ){
			tok->pstr++;
			while ( tok->pstr[0] != ')' ){
				tok->buffer = StringBuffer_append_char(tok->buffer, tok->pstr[0]);
				tok->pstr++;
			}
			return SML_REFERENCE;
		}
		// space, tab, \n, \r
		else {
			tok->pstr++;
		}

	}
	return SML_EOS;
}


#pragma mark -
#pragma mark SMLNode

SMLNode sml_set_tag_node( SMLNode node, const char *tag );

SMLNode new_SMLNode( const SMLNode parent ){
	SMLNode node = (struct _SMLNode*)malloc( sizeof(struct _SMLNode) );
	assert(node);
	
	node->tag = NULL;
	
	node->child = (struct _SMLNode**)malloc( 2 * sizeof(struct _SMLNode*) );
	assert(node->child);
	for (int i = 0; i < 2; i++) {
		node->child[i] = NULL;
	}
	node->child_capacity = 2;
	node->child_count = 0;
	
	node->data = NULL;
	
	node->parent = parent;
	
	node->type = 0;
	
	return node;
}

void free_SMLNode( SMLNode node ){
	free(node->tag);
	node->tag = NULL;
	
	free( node->data);
	node->data = NULL;
	
	free(node->child);
	node->child = NULL;
	
	free(node);
	node = NULL;
}

#pragma mark -
#pragma mark SMLTree

SMLNode sml_add_node( SMLNode dst, SMLNode src ){
	if( dst->child_capacity == dst->child_count ){
		dst->child = realloc( dst->child, (++dst->child_capacity) *sizeof(SMLNode*));
		dst->child_capacity++;
	}
	dst->child[dst->child_count] = src;
	dst->child_count++;
	return dst;
}

SMLNode sml_set_tag_node( SMLNode node, const char *tag ){
	if( tag != NULL){
		node->tag = (char*)malloc( (strlen(tag)+1) * sizeof(char) );
		assert(node->tag);
		strcpy(node->tag, tag);
	}
	return node;
}

SMLNode sml_set_data_node( SMLNode node, const char *data, const int type ){
	node->data = (char*)malloc( (strlen(data)+1) * sizeof(char) );
	assert(node->data);
	strcpy(node->data, data);
	node->type = type;
	return node;
}

SMLNode new_SMLTree( Tokenizer *tok ){
	int type = 0;
	SMLNode current = NULL;
	SMLNode root    = NULL;
	int open = 0;
	int close = 0;
	while( (type = sml_next_token(tok)) != SML_EOS ){
		switch ( type ) {
			case SML_OPEN:{
				open++;
				SMLNode n = new_SMLNode( current );
				if (current != NULL ) current = sml_add_node( current, n);
				else root = n;
				current = n;
				break;
			}
			case SML_CLOSE:{
				close++;
				current = current->parent;
				break;
			}
			case SML_TAG:{
				current = sml_set_tag_node(current, tok->buffer->c);
				break;
			}
			default:{
				current = sml_set_data_node(current, tok->buffer->c, type );
				break;
			}
		}
	}
	assert(open == close);
	return root;
}

static void free_SMLTree_aux( SMLNode node){
	if (node == NULL) return;
	for (int i = 0; i < node->child_count; i++) {
		free_SMLTree_aux( node->child[i]);
	}
	free_SMLNode(node);
	node = NULL;
}

void free_SMLTree( SMLNode root){
	free_SMLTree_aux(root);
}


char * SML_get_data( const SMLNode node ){
	return node->data;
}

char * SML_get_tag( const SMLNode node ){
	return node->tag;
}


char *SML_get_data_of_child( const SMLNode node, const char *tag ){
	for (int i = 0; i < node->child_count; i++) {
		if( strcmp( tag, node->child[i]->tag) == 0 ){
			return node->child[i]->data;
		}
	}
	return NULL;
}

SMLNode SML_get_element( SMLNode node, const char *tag ){
	for (int i = 0; i < node->child_count; i++) {
		if( strcmp( tag, node->child[i]->tag) == 0 ){
			return node->child[i];
		}
	}
	return NULL;
}

SMLNode SML_get_element_by_index( const SMLNode node, const int index ){
	if( index < node->child_count && index >= 0) {
		return node->child[index];
	}
	return NULL;
}

int SML_get_child_count( const SMLNode node ){
	return node->child_count;
}

static void SML_find_from_id_aux( SMLNode node, const char *id, SMLNode *thenode ){
	if (node != NULL ) {
		SMLNode test = NULL;
		if ( (test = SML_get_element(node, "id") ) != NULL) {
			char *test_id = SML_get_data(test);
			if( test != NULL && strcmp( test_id, id) == 0 ){
				*thenode = node;
			}
		}
		
		for (int i = 0; i < node->child_count; i++) {
			SML_find_from_id_aux( node->child[i], id, thenode );
		}
	}
	
}

SMLNode SML_find_from_id( SMLNode node, const char *id ){
	SMLNode thenode = NULL;
	SML_find_from_id_aux( node, id, &thenode);
	return thenode;
}


static void padit( FILE *pfile, const int n ){
	for (int i = 0; i < n; i++) {
		fprintf(pfile, "  ");
	}
}

static void print_SMLTree_aux( FILE *pfile, const SMLNode node, int indent ){
	if( node == NULL ) return;
	padit(pfile, indent);
	if( node->child_count != 0 ){
		fprintf(pfile, "(%s: \n", node->tag);
	}
	else{
		fprintf(pfile, "(%s: %s)\n", node->tag, node->data);
	}
	indent++;
	for (int i = 0; i < node->child_count; i++) {
		print_SMLTree_aux( pfile, node->child[i], indent );
	}
	indent--;
	if( node->child_count != 0 ){
		padit(pfile, indent);
		fprintf(pfile, ")\n");
	}
}

void SML_insert_buffer( SMLNode node, const StringBuffer* buffer ){
	SML_insert_string( node, buffer->c );
}

void SML_insert_string( SMLNode node, const char *string ){
	Tokenizer *tok = new_Tokenizer( string);
	SMLNode new = new_SMLTree(tok);
	sml_add_node( node, new );
	new->parent = node;
	free_Tokenizer(tok);
}

void print_SMLTree( FILE *pfile, const SMLNode root ){
	int indent = 0;
	print_SMLTree_aux( pfile, root, indent );
	
}

static void pad_buffer( StringBuffer *buffer, const int n ){
	for (int i = 0; i < n; i++) {
		StringBuffer_append_string(buffer, "  " );
	}
}
static void SMLTree_to_buffer_aux( SMLNode node, StringBuffer *buffer, int indent, bool pretty ){
	if ( node == NULL ) return;
	if( pretty ) pad_buffer(buffer, indent);
	if( node->child_count != 0 ){
		StringBuffer_append_format(buffer, "(%s:", node->tag );
		if( pretty ) StringBuffer_append_format(buffer, " \n" );
	}
	else{
		StringBuffer_append_format(buffer, "(%s: %s)", node->tag, node->data );
		if( pretty ) StringBuffer_append_char(buffer, '\n' );
	}
	indent++;
	for (int i = 0; i < node->child_count; i++) {
		SMLTree_to_buffer_aux( node->child[i], buffer, indent, pretty );
	}
	indent--;
	if( node->child_count != 0 ){
		if( pretty ) pad_buffer(buffer, indent);
		buffer = StringBuffer_append_char(buffer, ')' );
		if( pretty ) StringBuffer_append_char(buffer, '\n' );
	}
}

StringBuffer * SMLTree_to_buffer( SMLNode root, bool pretty ){
	StringBuffer *buffer = new_StringBuffer(10);
	int indent = 0;
	SMLTree_to_buffer_aux( root, buffer, indent, pretty );
	
	return buffer;
}

double * sml_string_to_double_array( char *string, const unsigned int n ){
	double *array = NULL;
	char *pstr = &string[0];
	StringBuffer *buffer = new_StringBuffer(10);
	int count = n;
	
	if( count == 0){
		while ( *pstr != '\0' ) {
			if( (*pstr > '0' && *pstr< '9') || *pstr == '.' || *pstr == '-' ){
				while ( (*pstr > '0' && *pstr< '9') || *pstr == '.' ) {
					pstr++;
				}
				count++;
			}
		}
	}
	
	array = dvector( count );
	
	pstr = &string[0];
	count = 0;
	while ( *pstr != '\0' ) {
		if( isspace(*pstr) && buffer->length != 0 ){
			array[count++] = atof( buffer->c );
			StringBuffer_empty( buffer );
		}
		else {
			StringBuffer_append_char(buffer, *pstr);
		}
		
		pstr++;
	}
	
	if ( buffer->length != 0) {
		array[count] = atof( buffer->c );
	}
	free_StringBuffer(buffer);
	return array;
}
