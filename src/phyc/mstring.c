/*
 *  mstring.c
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


#include "mstring.h"

#include <string.h>
#include <strings.h>
#include <assert.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "matrix.h"


#pragma mark String

char * String_new( const size_t capacity ){
	assert(capacity != 0);
	char *string = (char *)calloc( capacity, sizeof(char));
	assert(string);
	string[0] = '\0';
	return string;
}

char * String_clone( const char *src ){
	char *string = (char *)calloc( (strlen(src)+1), sizeof(char));
	assert(string);
	strcpy(string, src);
	return string;
}

char * String_append_char( char *dst, const char src, size_t *capacity ){
	size_t length_dst = strlen(dst);
	
	if( length_dst == *capacity-1 ){
		*capacity *= 2;
		dst = realloc( dst,  *capacity *sizeof(char) );
		assert(dst);
	}
	dst[length_dst]   = src;
	dst[length_dst+1] = '\0';
	return dst;
}

char * String_append_string( char *dst, const char *src ){
	size_t length_dst = strlen(dst);
	size_t length_src = strlen(src);
	
	dst = realloc( dst, (length_dst + length_src + 1) * sizeof(char) );
	assert(dst);
	
	strcpy(&dst[length_dst], src);
	dst[length_src+length_dst] = '\0';
	return dst;
}

char * String_compact( char *dst ){
	dst = realloc( dst, (strlen(dst)+1) * sizeof(char));
	assert(dst);
	return dst;
}


char * String_ltrim( char *s ) {
	char *ptr = s;
    while(isspace(*ptr)) ptr++;
	if( ptr != s ) memmove(s, ptr, strlen(ptr)+1);
    return s;
}

// This is one dangerous because we cannot free it since the first cell of the array can be different
char * String_ltrimX( char *s ) {
	while(isspace(*s)) s++;
	return s;
}

char * String_rtrim( char *s ) {
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char * String_trim( char *s ){
	if (!s)
        return NULL;   // handle NULL string
    if (!*s)
        return s;      // handle empty string
	
	s = String_rtrim(s);
    return String_ltrim(s);
}

// A string is empty it only contains spaces or tabs
bool String_is_empty( const char *s ){
	const char *ptr = s;
	while ( *ptr != '\0' ) {
		if ( *ptr != ' ' && *ptr != '\t' ) {
			break;
		}
		ptr++;
	}
	return 	*ptr == '\0' ? true : false;
}

double* String_to_double_array( const char *string, const char separator, unsigned *length ){
    double *array = NULL;
    const char *pstr = &string[0];
    StringBuffer *buffer = new_StringBuffer(10);
    int count = 0;
    
    while ( *pstr != '\0' ) {
        if( *pstr == separator ){
            count++;
        }
        pstr++;
    }
    count++;
    array = dvector( count );
    
    pstr = &string[0];
    count = 0;
    while ( *pstr != '\0' ) {
        if( *pstr == separator && buffer->length != 0 ){
            StringBuffer_trim(buffer);
            array[count++] = atof( buffer->c );
            StringBuffer_empty( buffer );
        }
        else {
            StringBuffer_append_char(buffer, *pstr);
        }
        pstr++;
    }
    
    if ( buffer->length != 0) {
        array[count++] = atof( buffer->c );
    }
    free_StringBuffer(buffer);
    *length = count;
    
    return array;
}

int* String_to_int_array( const char *string, const char separator, unsigned *length ){
    int *array = NULL;
    const char *pstr = &string[0];
    StringBuffer *buffer = new_StringBuffer(10);
    int count = 0;
    
    while ( *pstr != '\0' ) {
        if( *pstr == separator ){
            count++;
        }
        pstr++;
    }
    
    
    array = ivector( count+1 );
    
    pstr = &string[0];
    count = 0;
    while ( *pstr != '\0' ) {
        if( *pstr == separator && buffer->length != 0 ){
            StringBuffer_trim(buffer);
            array[count++] = atoi( buffer->c );
            StringBuffer_empty( buffer );
        }
        else {
            StringBuffer_append_char(buffer, *pstr);
        }
        
        pstr++;
    }
    
    if ( buffer->length != 0) {
        array[count++] = atoi( buffer->c );
    }
    free_StringBuffer(buffer);
    *length = count;
    
    return array;
}

unsigned * String_to_uint_array( const char *string, const char separator, unsigned *length ){
    unsigned *array = NULL;
    const char *pstr = &string[0];
    StringBuffer *buffer = new_StringBuffer(10);
    int count = 0;
    
    while ( *pstr != '\0' ) {
        if( *pstr == separator ){
            count++;
        }
        pstr++;
    }
    
    
    array = uivector( count+1 );
    
    pstr = &string[0];
    count = 0;
    while ( *pstr != '\0' ) {
        if( *pstr == separator && buffer->length != 0 ){
            StringBuffer_trim(buffer);
            array[count++] = atoi( buffer->c );
            StringBuffer_empty( buffer );
        }
        else {
            StringBuffer_append_char(buffer, *pstr);
        }
        
        pstr++;
    }
    
    if ( buffer->length != 0) {
        array[count++] = atoi( buffer->c );
    }
    free_StringBuffer(buffer);
    *length = count;
    
    return array;
}

char ** String_to_string_array( const char *string, const char separator, unsigned *length ){
    char **array = NULL;
    const char *pstr = &string[0];
    StringBuffer *buffer = new_StringBuffer(10);
    int count = 0;
    
    while ( *pstr != '\0' ) {
        if( *pstr == separator ){
            count++;
        }
        pstr++;
    }
    
    
    array = (char**)malloc( sizeof(char*)*(count+1) );
    assert(array);
    
    pstr = &string[0];
    count = 0;
    while ( *pstr != '\0' ) {
        if( *pstr == separator && buffer->length != 0 ){
            StringBuffer_trim(buffer);
            array[count++] = String_clone(buffer->c);
            StringBuffer_empty( buffer );
        }
        else {
            StringBuffer_append_char(buffer, *pstr);
        }
        
        pstr++;
    }
    
    if ( buffer->length != 0) {
        array[count++] = String_clone(buffer->c);
    }
    free_StringBuffer(buffer);
    *length = count;
    
    return array;
}

bool String_to_boolean( const char *str ){
	if ( strcmp(str, "true") == 0 || strcmp(str, "TRUE") == 0 ) {
		return true;
	}
	return false;
}

bool String_start_with( const char *str, const char *substr, bool ignore_space ){
	assert ( str != NULL && substr != NULL );
	
	
	if ( strlen(str) < strlen(substr) ) {
		return false;
	}
	
	const char *ptr = str;
	if ( ignore_space ) {
		while (  *ptr != '\0' && isspace(*ptr) ){
			++ptr;
		}
	}
	
	if ( strlen(ptr) < strlen(substr) ) {
		return false;
	}
	
	return strncmp(ptr, substr, strlen(substr)) == 0;
}

bool String_i_start_with( const char *str, const char *substr, bool ignore_space ){
	assert ( str != NULL && substr != NULL );
	
	
	if ( strlen(str) < strlen(substr) ) {
		return false;
	}
	
	const char *ptr = str;
	if ( ignore_space ) {
		while (  *ptr != '\0' && isspace(*ptr) ){
			++ptr;
		}
	}
	
	if ( strlen(ptr) < strlen(substr) ) {
		return false;
	}
	
	const char *ptr2 = substr;
	
	while (  *ptr2 != '\0' ){
		if ( toupper(*ptr) != toupper(*ptr2) ) {
			return false;
		}
		ptr++;
		ptr2++;
	}
	return true;
}

bool String_end_with( const char *str, const char *substr, bool ignore_space ){
	assert ( str != NULL && substr != NULL );
	
	
	if ( strlen(str) < strlen(substr) ) {
		return false;
	}
	
	
	const char *ptr = &str[strlen(str)-1];
	int nBlanks = 0;
	if ( ignore_space ) {
		while (  *ptr && isspace(*ptr) ){
			--ptr;
			++nBlanks;
		}
	}
	
	if ( strlen(str)-nBlanks < strlen(substr) ) {
		return false;
	}
	
	return strncmp(str + (strlen(str)-nBlanks - strlen(substr)), substr, strlen(substr)) == 0;
}

bool String_equals( const char *str, const char *substr, bool ignore_space ){
	assert ( str != NULL && substr != NULL );
	
	
	if ( strlen(str) < strlen(substr) ) {
		return false;
	}
	
	if ( !ignore_space ) {
		return strncmp(str, substr, strlen(substr)) == 0;
	}

	const char *ptr = str;
	while (  *ptr != '\0' && isspace(*ptr) ){
		++ptr;
	}
	
	if ( strlen(ptr) < strlen(substr) ) {
		return false;
	}
	return strncmp(ptr, substr, strlen(substr)) == 0;
}

int String_strcasecmp( const char *s1, const char *s2 ) {
	const u_char *us1 = (const u_char *)s1, *us2 = (const u_char *)s2;
	
	while (tolower(*us1) == tolower(*us2)) {
		if (*us1++ == '\0')
			return 0;
		us2++;
	}
	return (tolower(*us1) - tolower(*us2));
}

int String_strncasecmp( const char *s1, const char *s2, size_t n ) {
	
	if (n != 0) {
		const u_char *us1 = (const u_char *)s1;
		const u_char *us2 = (const u_char *)s2;
		
		do {
			if ( tolower(*us1) != tolower(*us2) ){
				return (tolower(*us1) - tolower(*us2));
			}
			if (*us1++ == '\0') {
				break;
			}
			us2++;
		} while ( --n != 0 );
	}
	return 0;
}

bool String_iequals( const char *s1, const char *s2 ) {
	const u_char *us1 = (const u_char *)s1, *us2 = (const u_char *)s2;
	
	while (tolower(*us1) == tolower(*us2)) {
		if (*us1++ == '\0')
			return true;
		us2++;
	}
	return false;
}

bool String_iequals_n( const char *s1, const char *s2, size_t n ) {
	
	if (n != 0) {
		const u_char *us1 = (const u_char *)s1;
		const u_char *us2 = (const u_char *)s2;
		
		do {
			if ( tolower(*us1) != tolower(*us2) ){
				return false;
			}
			if (*us1++ == '\0') {
				break;
			}
			us2++;
		} while ( --n != 0 );
	}
	return true;
}


bool String_contains( const char *str, char car ){
	for ( int i = 0; i < strlen(str); i++ ) {
		if ( str[i] == car ) {
			return true;
		}
	}
	return false;
}

bool String_contains_str( const char *str, const char *chars ){
	if( str == NULL ){
		return false;
	}
	const char *pch = str;
	size_t len = strlen(str);
	size_t len2 = strlen(chars);
	
	while ( len >= len2 ) {
		if ( strncmp( pch, chars, len2 ) == 0 ){
			return true;
		}
		pch++;
		assert(len > 0);
		len--;
	}
	return false;
}

int String_contains_str_count( const char *str, const char *chars ){
	const char *pch = str;
	size_t len = strlen(str);
	size_t len2 = strlen(chars);
	int count = 0;
	while ( len >= len2 ) {
		if ( strncmp( pch, chars, len2 ) == 0 ){
			count++;
		}
		pch++;
		len--;
	}
	return count;
}

int String_index_of_str( const char *str, const char *chars ){
	const char *pch = str;
	size_t len = strlen(str);
	size_t len2 = strlen(chars);
	int pos = 0;
	
	while ( len >= len2 ) {
		if ( strncmp( pch, chars, len2 ) == 0 ){
			return pos;
		}
		pos++;
		pch++;
		len--;
	}
	return -1;
}


// strok modify str
char ** String_split( const char *str, const char *pattern, int *count ){
	*count = 0;
	int capacity = 10;
	char **substrs = NULL;
    
    size_t len = strlen(str);
    size_t len_pattern = strlen(pattern);
    
	if ( len_pattern > len || !String_contains_str(str, pattern) ) {
		return NULL;
	}
		
	substrs = (char**)calloc( capacity, sizeof(char*) );
    assert(substrs);
	StringBuffer *buffer = new_StringBuffer(100);
	
	const char *pch = str;
	
	while ( len >= len_pattern ) {
		if ( strncmp( pch, pattern, len_pattern ) == 0 ){
			if( *count == capacity ){
				capacity *= 2;
				substrs = realloc(substrs, capacity * sizeof(char*));
			}
			substrs[*count] = String_clone(buffer->c);
			(*count)++;
			StringBuffer_empty(buffer);
			pch += len_pattern;
			//len -= len_pattern;
		}
		else {
			StringBuffer_append_char(buffer, *pch);
			pch++;
		}
		len--;
	}
	
	if( buffer->length != 0 ){
		substrs[*count] = String_clone(buffer->c);
		(*count)++;
	}
	if( *count > 0 && *count != capacity ){
		substrs = realloc(substrs, *count * sizeof(char*) );
	}
	free_StringBuffer(buffer);
	return substrs;
}

char ** String_split_char( const char *str, const char pattern, int *count ){
	*count = 0;
	int c = 0;
	const char *pch = str;
    
	while ( *pch != '\0') {
		if ( *pch == pattern ) c++;
		pch++;
	}
	
	if ( c == 0 ) {
		char **substrs = malloc( sizeof(char*) );
		substrs[0] = String_clone(str);
		*count = 1;
		return substrs;
	}
	
	// c+1 because we count the number separators
	// if str == "a=b" then c==1 but we need space  for c+1
	char **substrs = (char**)calloc( c+1, sizeof(char*) );
    assert(substrs);
	
	StringBuffer *buffer = new_StringBuffer(strlen(str));
	
	pch = str;
	while ( *pch != '\0' ) {
		if ( *pch != pattern ) {
			StringBuffer_append_char(buffer, *pch);
		}
		else {
			substrs[*count] = String_clone(buffer->c);
			StringBuffer_empty(buffer);
			(*count)++;
		}
		pch++;
	}
	
	if (buffer->length != 0 ) {
		substrs[*count] = String_clone(buffer->c);
		(*count)++;
	}
    
	free_StringBuffer(buffer);
	
	return substrs;
}

double * String_split_char_double( const char *str, const char pattern, int *count ){
	*count = 0;
	int c = 0;
	const char *pch = str;
    
	while ( *pch != '\0') {
		if ( *pch == pattern ) c++;
		pch++;
	}
    
// it could be an array of length 1
//	if ( c == 0 ) {
//		return NULL;
//	}
	
	// c+1 because we count the number separators
	// if str == "a=b" then c==1 but we need space  for c+1
	double *substrs = dvector(c+1);

	
	StringBuffer *buffer = new_StringBuffer(strlen(str));
	
	pch = str;
	while ( *pch != '\0' ) {
		if ( *pch != pattern ) {
			StringBuffer_append_char(buffer, *pch);
		}
		else {
			substrs[*count] = atof(buffer->c);
			StringBuffer_empty(buffer);
            (*count)++;
		}
		pch++;
	}
	
	if (buffer->length != 0 ) {
		substrs[*count] = atof(buffer->c);
		(*count)++;
	}
    
	free_StringBuffer(buffer);
	
	return substrs;
}

void String_print_join( FILE *pf, const char *sep, const char **array, int n){
    for ( int i = 0; i < n; i++ ) {
        fprintf(pf, "%s%s",array[i],(i==n-1?"":sep));
    }
}

#pragma mark -
#pragma mark StringBuffer

StringBuffer * new_StringBuffer( const size_t capacity ){
	assert(capacity != 0);
	StringBuffer *buffer = (StringBuffer *)malloc(sizeof(StringBuffer));
	assert(buffer);
	buffer->c = String_new(capacity);
	//buffer->start = NULL;
	buffer->capacity = capacity;
	buffer->length = 0;
	return buffer;
}

void free_StringBuffer( StringBuffer *buffer ){
	free(buffer->c );
	buffer->c = NULL;
	free(buffer);
	buffer = NULL;
}

StringBuffer * StringBuffer_clone( StringBuffer *buffer ){
	StringBuffer *aBuffer = (StringBuffer *)malloc(sizeof(StringBuffer));
	assert(aBuffer);
	
	aBuffer->c = StringBuffer_tochar(buffer);
	aBuffer->length = buffer->length;
	aBuffer->capacity = buffer->length + 1;
	return aBuffer;
}

// only expand the buffer, do nothing if the buffer is bigger than size + 1
StringBuffer * StringBuffer_reserve( StringBuffer *buffer, const size_t size ){
	if ( buffer->capacity < size+1 ) {
		buffer->c = realloc(buffer->c, (size+1)*sizeof(char) );
		buffer->capacity = size+1;
	}
	return buffer;
}


#pragma mark -

void StringBuffer_prepend_string( StringBuffer *buffer, const char *src ){
	size_t src_len = strlen(src);
	if ( src_len + buffer->length >= buffer->capacity ) {
		buffer->c = realloc( buffer->c, (src_len + buffer->length + 1) * sizeof(char) );
		assert(buffer->c);
		buffer->capacity += src_len;
	}
	memmove(buffer->c + src_len, buffer->c, (buffer->length+1) *sizeof(char) ); // move \0 with it
	strncpy(buffer->c, src, src_len * sizeof(char)); // ignore \0
	buffer->length += src_len; 
}

// does not change the size of the string if it is bigger than the buffer
char * StringBuffer_assign( const StringBuffer *buffer, char *dst ){
	if ( dst == NULL ) {
		dst = String_clone(buffer->c);
	}
	else if( buffer->length > strlen(dst) ){
        dst = realloc( dst, (buffer->length+1) *sizeof(char) );
        strcpy(dst, buffer->c);
	}
	return dst;
}



StringBuffer * StringBuffer_set_string( StringBuffer *buffer, const char *src ){
	StringBuffer_empty(buffer);
	
	size_t length_src = strlen(src);
	
	if( length_src == 0 ){
        buffer->c[0] = '\0';
		return buffer;
	}
	
	if( length_src >= buffer->capacity ){
		buffer->c = realloc( buffer->c, (length_src + 1) * sizeof(char) );
		assert(buffer->c);
		buffer->capacity = length_src + 1;
	}
	strcpy(buffer->c, src);
	buffer->length = length_src;
	buffer->c[buffer->length] = '\0';
	return buffer;
}

StringBuffer * StringBuffer_append_char( StringBuffer *buffer, const char src ){
	if( buffer->length == buffer->capacity-1 ){
		buffer->capacity *= 2;
		buffer->c = realloc( buffer->c, buffer->capacity *sizeof(char) );
		assert(buffer->c);
	}
	buffer->c[buffer->length++]   = src;
	buffer->c[buffer->length] = '\0';
	return buffer;
}

StringBuffer * StringBuffer_append_StringBuffer( StringBuffer *buffer, StringBuffer *src ){
	return StringBuffer_append_string(buffer, src->c);
}

StringBuffer * StringBuffer_append_string( StringBuffer *buffer, const char *src ){
	size_t length_src = strlen(src);
	
	if( length_src == 0 ){
		return buffer;
	}
	
	if( buffer->length +length_src >= buffer->capacity ){
		buffer->c = realloc( buffer->c, (buffer->length + length_src + 1) * sizeof(char) );
		assert(buffer->c);
		buffer->capacity = buffer->length + length_src + 1;
	}
	strcpy(&buffer->c[buffer->length], src);
	buffer->length += length_src;
	buffer->c[buffer->length] = '\0';
	return buffer;
}

StringBuffer * StringBuffer_append_nstring( StringBuffer *buffer, const char *src, const size_t n ){
	size_t length_src = stmin(strlen(src), n);
	
	if( buffer->length +length_src >= buffer->capacity ){
		buffer->c = realloc( buffer->c, (buffer->length + length_src + 1) * sizeof(char) );
		assert(buffer->c);
		buffer->capacity = buffer->length + length_src + 1;
	}
	strncpy(&buffer->c[buffer->length], src, n);
	buffer->length += length_src;
	buffer->c[buffer->length] = '\0';
	return buffer;
}

StringBuffer * StringBuffer_append_strings( StringBuffer *buffer, const int n, ... ){
	va_list ap;
	va_start(ap,n);
	const char *str = NULL;
	for (int i = 0; i < n; i++) {
		str = va_arg(ap, const char*);
		StringBuffer_append_string( buffer, str );
	}
	va_end(ap);
	
	return buffer;
}

StringBuffer * StringBuffer_append_substring( StringBuffer *buffer, const char *src, const int n ){
	size_t length_src = strlen(src);
	
	if( length_src < n ) error("StringBuffer.c: StringBuffer_append_substring: error");
	
	if( buffer->length +n >= buffer->capacity ){
		buffer->c = realloc( buffer->c, (buffer->length + n + 1) * sizeof(char) );
		assert(buffer->c);
		buffer->capacity = buffer->length + n + 1;
	}
	
	strncpy( &buffer->c[buffer->length], src, n );
	buffer->length += n;
	buffer->c[ buffer->length ] = '\0';
	return buffer;
}

StringBuffer * StringBuffer_append_range( StringBuffer *buffer, const char *src, const int start, const int n ){
	size_t length_src = strlen(src);
	int length_substr = n - start;
	
	if( length_src < length_substr ) error("StringBuffer.c: StringBuffer_append_range: error");
	
	if( buffer->length +length_substr >= buffer->capacity ){
		buffer->c = realloc( buffer->c, (buffer->length + length_substr + 1) * sizeof(char) );
		assert(buffer->c);
		buffer->capacity = buffer->length + length_substr + 1;
	}
	
	strncpy( &buffer->c[buffer->length], &src[start], n );
	buffer->length += n;
	buffer->c[ buffer->length ] = '\0';
	return buffer;
}

StringBuffer * StringBuffer_append_format( StringBuffer *buffer, const char *format, ... ){
	va_list ap;
	va_start(ap,format);
	const char *cp = format;
	char buf[150];
	StringBuffer *buffer2 = new_StringBuffer(10);
	
	while( *cp != '\0'){
		if( *cp == '%' ){
			// escaped %
			if ( *(cp+1) == '%' ) {
				StringBuffer_append_char(buffer, '%');
				cp++;
			}
			else {
				StringBuffer_empty(buffer2);
                StringBuffer_append_char(buffer2, '%');
                
                if( *(cp+1) == '.' ){
                    StringBuffer_append_char(buffer2, *(++cp));
                    StringBuffer_append_char(buffer2, *(++cp));
                }
                
				char specifier = *(++cp);
                StringBuffer_append_char(buffer2, specifier);
                
				if( specifier == 'd' || specifier == 'i' ){
					int var = va_arg(ap,  int);
					sprintf(buf, buffer2->c, var);
					buffer = StringBuffer_append_string( buffer, buf );
				}
				else if( specifier == 's' ){
					char *str = va_arg(ap, char*);
					buffer = StringBuffer_append_string( buffer, str );
					
				}
				else if( specifier == 'c' ){
					char c = (char)va_arg(ap, int);
					buffer = StringBuffer_append_char( buffer, c );
					
				}
				else if( specifier == 'f' || specifier == 'e' || specifier == 'E' ){
					double var = va_arg(ap, double);
					sprintf(buf, buffer2->c, var);
					buffer = StringBuffer_append_string( buffer, buf );
				}
				else if( specifier == 'u' ){
					unsigned int var = va_arg(ap, unsigned int);
					sprintf(buf, buffer2->c, var);
					buffer = StringBuffer_append_string( buffer, buf );
				}
				else if( specifier == 'l' ){
					specifier = *(++cp);
					StringBuffer_append_char(buffer2, *cp);
					
					if( *cp == 'u' ){
						unsigned long var = va_arg(ap, unsigned long);
						sprintf(buf, buffer2->c, var);
						buffer = StringBuffer_append_string( buffer, buf );
					}
					else if( cp[1] == 'd' || specifier == 'i' ){
						long int var = va_arg(ap, long int);
						sprintf(buf, buffer2->c, var);
						buffer = StringBuffer_append_string( buffer, buf );
					}
					else {
						printf("specifier l %c %c\n", cp[1], cp[2]);
						exit(1);
					}
				}
				else if( specifier == 'z' ){
					specifier = *(++cp);
					StringBuffer_append_char(buffer2, *cp);
					
					if( *cp == 'u' ){
						size_t var = va_arg(ap, size_t);
						sprintf(buf, buffer2->c, var);
						buffer = StringBuffer_append_string( buffer, buf );
					}
					else {
						printf("specifier z %c %c\n", cp[1], cp[2]);
						exit(1);
					}
				}
			}
		}
		else{
			StringBuffer_append_char(buffer, *cp);
		}
		cp++;
	}
	va_end(ap);
	free_StringBuffer(buffer2);
	return buffer;
}

StringBuffer * StringBuffer_compact( StringBuffer *buffer ){
	buffer->c = realloc( buffer->c, (buffer->length + 1) * sizeof(char) );
	assert(buffer->c);
	return buffer;
}

char * StringBuffer_tochar( const StringBuffer *buffer ){
	char *c = String_clone(buffer->c);
	return c;
}

void StringBuffer_empty( StringBuffer *buffer ){
	buffer->c[0] = '\0';
	buffer->length = 0;
}

void StringBuffer_trim( StringBuffer *buffer ){
	buffer->c = String_trim( buffer->c );
	buffer->length = strlen(buffer->c);
}

void StringBuffer_chop( StringBuffer *buffer ){
    if( buffer->length != 0 ){
        buffer->c[buffer->length-1] = '\0';
        buffer->length--;
    }
}

void StringBuffer_cut( StringBuffer *buffer ){
    if( buffer->length != 0 ){
        memmove(buffer->c, &buffer->c[1], (buffer->length-1)* sizeof(char));
        buffer->c[buffer->length-1] = '\0';
        buffer->length--;
    }
}

void StringBuffer_assign_lowercase( const StringBuffer *src, StringBuffer *dst ){
	char *ptr = src->c;
	if( src->length +1 > dst->capacity ){
		dst->c = realloc( dst->c, (src->length + 1) * sizeof(char) );
		dst->capacity = src->length + 1;
	}
	char *ptr2 = dst->c;
	while ( *ptr != '\0' ) {
		*ptr2 = tolower(*ptr);
		ptr++;
		ptr2++;
	}
	dst->length = src->length;
}

void StringBuffer_assign_uppercase( const StringBuffer *src, StringBuffer *dst ){
	char *ptr = src->c;
	if( src->length +1 > dst->capacity ){
		dst->c = realloc( dst->c, (src->length + 1) * sizeof(char) );
		dst->capacity = src->length + 1;
	}
	char *ptr2 = dst->c;
	while ( *ptr != '\0' ) {
		*ptr2 = toupper(*ptr);
		ptr++;
		ptr2++;
	}
	dst->length = src->length;
}

void StringBuffer_lowercase( StringBuffer *src ){
	char *ptr = src->c;
	while ( *ptr != '\0' ) {
		*ptr = tolower(*ptr);
		ptr++;
	}
}

void StringBuffer_uppercase( StringBuffer *src ){
	char *ptr = src->c;
	while ( *ptr != '\0' ) {
		*ptr = toupper(*ptr);
		ptr++;
	}
}

