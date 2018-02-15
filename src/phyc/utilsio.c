//
//  utilsio.c
//  physher
//
//  Created by Mathieu Fourment on 16/02/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "utilsio.h"

#include <string.h>
#include "filereader.h"
#include "matrix.h"
#include "mstring.h"

Vector* read_log_column_with_id( const char *filename, size_t burnin, const char* id ){
	int count = 0;
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	Vector* vec = new_Vector(1000);
	
	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);// discard header
	ptr = reader->line;
	l = 0;
	char** header = String_split_char( ptr, '\t', &l );
	int i = 0;
	for (; i < l; i++) {
		if(strcmp(header[i], id) == 0){
			break;
		}
	}
	if(i == l){
		fprintf(stderr, "Could not find ID `%s`\n", id);
		exit(1);
	}
	for (int j = 0; j < l; j++) {
		free(header[j]);
	}
	free(header);
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		
		if ( reader->buffer->length == 0){
			continue;
		}
		if ( count >= burnin){
			ptr = reader->line;
			l = 0;
			temp = String_split_char_double( ptr, '\t', &l );
			Vector_push(vec, temp[i]);
			free(temp);
		}
		count++;
	}
	free_FileReader(reader);
	Vector_pack(vec);
	return  vec;
}

Vector** read_log_column_with_ids( const char *filename, size_t burnin, const char** tags, size_t tag_count ){
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	size_t sample = 0;
	
	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);
	ptr = reader->line;
	l = 0;
	char** header = String_split_char(ptr, '\t', &l);
	bool* in = bvector(l);
	
	Vector** vecs = malloc(2*sizeof(Vector*));
	
	for (int j = 0; j < tag_count; j++) {
		vecs[j] = new_Vector(1000);
		for (int i = 0; i < l; i++) {
			if(strcmp(tags[j], header[i]) == 0){
				in[i] = true;
				break;
			}
		}
	}
	
	for (int i = 0; i < l; i++) {
		free(header[i]);
	}
	free(header);
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		
		if ( reader->buffer->length == 0){
			continue;
		}
		if ( sample >= burnin){
			ptr = reader->line;
			l = 0;
			temp = String_split_char_double( ptr, '\t', &l );
			size_t index = 0;
			for (int i = 0; i < l; i++) {
				if(in[i]){
					Vector_push(vecs[index], temp[i]);
					index++;
				}
			}
			free(temp);
		}
		sample++;
	}
	free_FileReader(reader);
	for (int j = 0; j < tag_count; j++) {
		Vector_pack(vecs[j]);
	}
	free(in);
	return  vecs;
}
