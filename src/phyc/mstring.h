/*
 *  mstring.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/2/10.
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

#ifndef _MSTRING_H_
#define _MSTRING_H_

#include "utils.h"

#include <stdbool.h>

typedef struct StringBuffer{
	char *c;
	size_t length;
	size_t capacity;
} StringBuffer;


typedef unsigned char u_char;

#pragma mark -
#pragma mark String

char * String_new( const size_t capacity );

char * String_append_char( char *dst, const char src, size_t *capacity );

char * String_append_string( char *dst, const char *src );

char * String_compact( char *dst );

char * String_clone( const char *src );

char * String_ltrim( char *s );
char * String_rtrim( char *s );
char * String_trim( char *s );

bool String_is_empty( const char *s );

double * String_to_double_array( const char *string, const char separator, unsigned *length );

int * String_to_int_array( const char *string, const char separator, unsigned *length );

unsigned * String_to_uint_array( const char *string, const char separator, unsigned *length );

char ** String_to_string_array( const char *string, const char separator, unsigned *length );

bool String_to_boolean( const char *str );

bool String_start_with( const char *str, const char *substr, bool ignore_space );
bool String_i_start_with( const char *str, const char *substr, bool ignore_space );

bool String_end_with( const char *str, const char *substr, bool ignore_space );

bool String_equals( const char *str, const char *substr, bool ignore_space );

char ** String_split( const char *str, const char *pattern, int *count );

char ** String_split_char( const char *str, const char pattern, int *count );

double * String_split_char_double( const char *str, const char pattern, int *count );

bool String_contains( const char *str, char car );

bool String_contains_str( const char *str, const char *chars );

int String_contains_str_count( const char *str, const char *chars );

int String_index_of_str( const char *str, const char *chars );


int String_strcasecmp( const char *s1, const char *s2 ) ;

int String_strncasecmp( const char *s1, const char *s2, size_t n );

bool String_iequals( const char *s1, const char *s2 );

bool String_iequals_n( const char *s1, const char *s2, size_t n );

void String_print_join( FILE *pf, const char *sep, const char **array, int n);

#pragma mark -
#pragma mark StringBuffer

StringBuffer * new_StringBuffer( const size_t capacity );

void free_StringBuffer( StringBuffer *buffer );

StringBuffer * StringBuffer_clone( StringBuffer *buffer );

void StringBuffer_prepend_string( StringBuffer *buffer, const char *src );

char * StringBuffer_assign( const StringBuffer *buffer, char *dst );

StringBuffer * StringBuffer_reserve( StringBuffer *buffer, const size_t size );

StringBuffer * StringBuffer_set_string( StringBuffer *buffer, const char *src );

StringBuffer * StringBuffer_append_char( StringBuffer *buffer, const char src );

StringBuffer * StringBuffer_append_string( StringBuffer *buffer, const char *src );

StringBuffer * StringBuffer_append_nstring( StringBuffer *buffer, const char *src, const size_t n );

StringBuffer * StringBuffer_append_strings( StringBuffer *buffer, const int n, ... );

StringBuffer * StringBuffer_append_substring( StringBuffer *buffer, const char *src, const int n );

StringBuffer * StringBuffer_append_range( StringBuffer *buffer, const char *src, const int start, const int n );

StringBuffer * StringBuffer_append_format( StringBuffer *buffer, const char *format, ... );

StringBuffer * StringBuffer_append_StringBuffer( StringBuffer *buffer, StringBuffer *src );

StringBuffer * StringBuffer_compact( StringBuffer *buffer );

void StringBuffer_empty( StringBuffer *buffer );

char * StringBuffer_tochar( const StringBuffer *buffer );

void StringBuffer_trim( StringBuffer *buffer );

void StringBuffer_chop( StringBuffer *buffer );

void StringBuffer_cut( StringBuffer *buffer );

void StringBuffer_assign_lowercase( const StringBuffer *src, StringBuffer *dst );

void StringBuffer_assign_uppercase( const StringBuffer *src, StringBuffer *dst );

void StringBuffer_lowercase( StringBuffer *src );

void StringBuffer_uppercase( StringBuffer *src );


#endif
