/*
 *  filereader.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/14/11.
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

#include "filereader.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"

static bool _FileReader_readline_fget( FileReader *reader );
static bool _FileReader_readline_fread( FileReader *reader );

FileReader * new_FileReader( const char *filename, const unsigned buffer_size ){
	FileReader *reader = malloc( sizeof(FileReader) );
	assert(reader);
	reader->line = NULL;
	reader->file = fopen(filename, "r");
	if( reader->file == NULL ){
		fprintf(stderr, "File does not exist: %s\n", filename);
		exit(1);
	}
	reader->buffer = new_StringBuffer(buffer_size);
	reader->line = reader->buffer->c;
    reader->b[0] = '\0';
    reader->read_line = _FileReader_readline_fread;
    reader->index = 0;
	return reader;
}

FileReader * new_FileReader_with_mode( const char *filename, const unsigned buffer_size, const char *mode ){
	FileReader *reader = malloc( sizeof(FileReader) );
	assert(reader);
	reader->line = NULL;
	reader->file = fopen(filename, mode);
	if( reader->file == NULL ){
		fprintf(stderr, "File does not exist: %s\n", filename);
		exit(1);
	}
	reader->buffer = new_StringBuffer(buffer_size);
    reader->line = reader->buffer->c;
    reader->index = 0;
    reader->b[0] = '\0';
    reader->read_line = _FileReader_readline_fread;
	return reader;
}

void free_FileReader( FileReader *reader ){
	free_StringBuffer(reader->buffer);
	reader->line = NULL;
	fclose(reader->file);
	free(reader);
	reader = NULL;
}

bool _FileReader_readline_fget( FileReader *reader ){
    char c;
    
    if ( feof(reader->file)  ) {
        return false;
    }
    
    StringBuffer_empty(reader->buffer);

    while ( 1 ) {
        c = fgetc(reader->file);
        
        if ( c == '\n' || c == '\r' ) {
            // check if it is windows
            if ( c == '\r' ) {
                fpos_t pos;
                fgetpos(reader->file, &pos);
                char cc = fgetc(reader->file);
                if( cc != '\n' ){
                    fsetpos(reader->file, &pos);
                }
            }
            break;
        }
        else {
            if( c == EOF ){
                break;
            }
            else {
                StringBuffer_append_char(reader->buffer, c);
            }
        }
    }
    reader->line = reader->buffer->c;
    return true;
}


bool _FileReader_readline_fread2( FileReader *reader ){
    
    int ret,i,len;
    char c;
    char *b;
    
    if ( feof(reader->file) && reader->b[0] == '\e'  ) {
        return false;
    }
    
    StringBuffer_empty(reader->buffer);
    len = strlen(&reader->b[reader->index]) + reader->index;
    
    for ( i = reader->index; i < len; i++ ) {
        b = &reader->b[i];
        
        if( *b == '\n' || *b == '\r' ){
            *b = '\0';
            StringBuffer_append_string(reader->buffer, &reader->b[reader->index]);
            
            if( i < len-1){
                if( b[0] == '\r' && b[1] == '\n'){
                    i++;
                }
            }
            reader->index = i+1;
            reader->line = reader->buffer->c;
            return true;
        }
    }
    
    // last line of the file without \n and \r at the end
    if( i == strlen(reader->b) && feof(reader->file) ){
        StringBuffer_append_string(reader->buffer, reader->b);
        reader->b[0] = '\e';
        return true;
    }
    
    StringBuffer_set_string(reader->buffer, reader->b+reader->index);
    reader->b[0] = '\0';
    reader->index = 0;
    
    while ( 1 ) {
        ret = fread(reader->b, sizeof(char), 512, reader->file);
        if( ret == 0 ) return false;
        
        for ( i = 0; i < ret; i++ ) {
            c = reader->b[i];
            if( c == '\n' || c == '\r' ){
                reader->b[i] = '\0';
                StringBuffer_append_string(reader->buffer, reader->b);
                
                if( i < ret-1){
                    if( c == '\r' && reader->b[i+1] == '\n'){
                        i++;
                    }
                }
                reader->index = i+1;
                reader->line = reader->buffer->c;
                return true;
            }
        }
        reader->b[i] = '\0';
        StringBuffer_append_string(reader->buffer, reader->b);
        
    }
    
    reader->line = reader->buffer->c;
    return true;
}

bool _FileReader_readline_fread( FileReader *reader ){
    
    int ret,i;
    char c;
    unsigned len;
    
    if ( feof(reader->file) && reader->b[0] == '\e'  ) {
        return false;
    }
    
    StringBuffer_empty(reader->buffer);
    len = strlen(reader->b);
    
    for ( i = 0; i < len; i++ ) {
        c = reader->b[i];
        if( c == '\n' || c == '\r' ){
            reader->b[i] = '\0';
            StringBuffer_append_string(reader->buffer, reader->b);
            
            if( i < len-1){
                if( c == '\r' && reader->b[i+1] == '\n'){
                    i++;
                }
            }
            memmove(&reader->b[0], &reader->b[i+1], (len-i-1)*sizeof(char));
            reader->b[len-i-1] = '\0';
            
            reader->line = reader->buffer->c;
            return true;
        }
    }
    
    // last line of the file without \n and \r at the end
    if( i == strlen(reader->b) && feof(reader->file) ){
        StringBuffer_append_string(reader->buffer, reader->b);
        reader->b[0] = '\e';
        return true;
    }
    
    StringBuffer_set_string(reader->buffer, reader->b);
    reader->b[0] = '\0';
    
    while ( 1 ) {
        ret = fread(reader->b, sizeof(char), 512, reader->file);
        if( ret == 0 ){
            reader->b[0] = '\e';
            break;
        }
        
        for ( i = 0; i < ret; i++ ) {
            c = reader->b[i];
            if( c == '\n' || c == '\r' ){
                reader->b[i] = '\0';
                StringBuffer_append_string(reader->buffer, reader->b);
                
                if( i < ret-1){
                    if( c == '\r' && reader->b[i+1] == '\n'){
                        i++;
                    }
                }
                memmove(&reader->b[0], &reader->b[i+1], (ret-i-1)*sizeof(char));
                reader->b[ret-i-1] = '\0';
                
                reader->line = reader->buffer->c;
                return true;
            }
        }
        reader->b[i] = '\0';
        StringBuffer_append_string(reader->buffer, reader->b);
        
    }
    
    reader->line = reader->buffer->c;
    return true;
}

double ** FileReader_csv_double( const char *filename, int nrow, int ncol ){
    int count = 0;
    char *ptr = NULL;
    double *temp = NULL;
    int l,i;
    double **matrix = dmatrix(nrow, ncol);
    
    FileReader *reader = new_FileReader(filename, 1000);
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);

		if ( reader->buffer->length == 0){
            continue;
        }
        
        if( count == ncol ){
            fprintf(stderr, "FileReader_csv_double: ncol is different from the number of rows in file %s\n", filename);
            exit(1);
        }
        
        ptr = reader->line;
        l = 0;
        temp = String_split_char_double( ptr, ',', &l );
        for ( i = 0; i < l; i++ ) {
            matrix[i][count] = temp[i];
        }
        free(temp);
		count++;
		if(l != nrow ){
            fprintf(stderr, "FileReader_csv_double: nrow is different from the number of columns in file %s (nrow: %d l: %d)\n", filename, nrow, l);
            exit(1);
        }
	}
	free_FileReader(reader);
    return  matrix;
}


char* load_file(const char *filename){
	FileReader* reader = new_FileReader(filename, 1000);
	StringBuffer* buffer = new_StringBuffer(100);
	while ( reader->read_line(reader) ) {
		StringBuffer_append_string(buffer, reader->buffer->c);
	}
	free_FileReader(reader);
	char* content = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return content;
}
