//
//  parametersio.c
//  physher
//
//  Created by Mathieu Fourment on 6/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "parametersio.h"

#include <string.h>
#include <strings.h>

#include "filereader.h"

Vector** read_log_for_parameters( const char *filename, size_t burnin, size_t* count, Parameters* params ){
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	size_t capacity = 1000;
	*count = 0;
	size_t paramCount = Parameters_count(params);
	Vector** vecs = malloc(capacity*sizeof(Vector*));
	size_t sample = 0;
	
	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);// discard header
	ptr = reader->line;
	l = 0;
	char** header = String_split_char(ptr, '\t', &l);
	bool* in = bvector(l);
	
	for (int j = 0; j < paramCount; j++) {
		for (int i = 0; i < l; i++) {
			if(strcmp(Parameters_name(params, j), header[i]) == 0){
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
			if(*count == capacity){
				capacity *= 2;
				vecs = realloc(vecs, capacity*sizeof(Vector*));
			}
			ptr = reader->line;
			l = 0;
			temp = String_split_char_double( ptr, '\t', &l );
			vecs[*count] = new_Vector(paramCount);
			for (int i = 0; i < l; i++) {
				if(in[i]){
					Vector_push(vecs[*count], temp[i]);
				}
			}
			free(temp);
			(*count)++;
		}
		sample++;
	}
	free_FileReader(reader);
	vecs = realloc(vecs, *count*sizeof(Vector*));
	free(in);
	return  vecs;
}

Vector** read_log_for_parameters_t( const char *filename, size_t burnin, Parameters* params ){
	size_t paramCount = Parameters_count(params);
	
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	Vector** vec = malloc(sizeof(Vector*)*paramCount);
	for(int i = 0; i < paramCount; i++){
		vec[i] = new_Vector(1000);
	}
	size_t sample = 0;
	
	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);
	ptr = reader->line;
	l = 0;
	char** header = String_split_char(ptr, '\t', &l);
	bool* in = bvector(l);
	
	for (int j = 0; j < paramCount; j++) {
		for (int i = 0; i < l; i++) {
			if(strcmp(Parameters_name(params, j), header[i]) == 0){
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
			int index = 0;
			for (int i = 0; i < l; i++) {
				if(in[i]){
					Vector_push(vec[index], temp[i]);
					index++;
				}
			}
			free(temp);
		}
		sample++;
	}
	free_FileReader(reader);
	free(in);
	return vec;
}

Vector** read_log_for_names_t( const char *filename, size_t burnin, char** params, size_t paramCount ){
	
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	Vector** vec = malloc(sizeof(Vector*)*paramCount);
	for(int i = 0; i < paramCount; i++){
		vec[i] = new_Vector(1000);
	}
	size_t sample = 0;
	
	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);
	ptr = reader->line;
	l = 0;
	char** header = String_split_char(ptr, '\t', &l);
	bool* in = bvector(l);
	
	for (int j = 0; j < paramCount; j++) {
		for (int i = 0; i < l; i++) {
			if(strcmp(params[j], header[i]) == 0){
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
			int index = 0;
			for (int i = 0; i < l; i++) {
				if(in[i]){
					Vector_push(vec[index], temp[i]);
					index++;
				}
			}
			free(temp);
		}
		sample++;
	}
	free_FileReader(reader);
	free(in);
	return vec;
}

