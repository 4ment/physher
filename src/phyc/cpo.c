//
//  cpo.c
//  physher
//
//  Created by Mathieu Fourment on 22/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "cpo.h"

#include <float.h>

#include "mstring.h"
#include "matrix.h"
#include "filereader.h"

void _cpo_calculate(struct CPO* cpo){
	char *ptr = NULL;
	double *temp = NULL;
	int l = 0;
	size_t capacity = 1000;
	size_t count = 0;
	Vector** vecs = malloc(capacity*sizeof(Vector*));
	size_t sample = 0;
	double* weights = NULL;
	FileReader *reader = new_FileReader(cpo->filename, 1000);
	reader->read_line(reader);
	weights = String_split_char_double( reader->line+1, '\t', &l );
	reader->read_line(reader);// discard header
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		
		if ( reader->buffer->length == 0){
			continue;
		}
		if ( sample >= cpo->burnin){
			if(count == capacity){
				capacity *= 2;
				vecs = realloc(vecs, capacity*sizeof(Vector*));
			}
			ptr = reader->line;
			l = 0;
			temp = String_split_char_double( ptr, '\t', &l );
			vecs[count] = new_Vector(l-1);
			for (int i = 1; i < l; i++) {
				Vector_push(vecs[count], temp[i]);
			}
			free(temp);
			count++;
		}
		sample++;
	}
	free_FileReader(reader);
	
	size_t n = Vector_length(vecs[0]);
	double logCPO = 0;
	for (int i = 0; i < n; i++) {
		double sum = -DBL_MAX;
		double min = Vector_at(vecs[0], i);
		for (int j = 1; j < count; j++) {
			min = dmin(Vector_at(vecs[j], i), min);
		}
		for (int j = 0; j < count; j++) {
			sum = logaddexp(sum, min-Vector_at(vecs[j], i));
		}
		logCPO += (log(count) + min - sum)*weights[i];
	}
	printf("logCPO: %f\n", logCPO);
	for (int i = 0; i < count; i++) {
		free_Vector(vecs[i]);
	}
	free(vecs);
	free(weights);
}

void _free_cpo(CPO* cpo){
	free(cpo->filename);
	free(cpo);
}


CPO* new_CPO_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"burnin",
		"filename"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* filename = get_json_node_value_string(node, "filename");
	CPO* cpo = malloc(sizeof(CPO));
	cpo->filename = String_clone(filename);
	cpo->burnin = get_json_node_value_size_t(node, "burnin", 0);
	cpo->calculate = _cpo_calculate;
	cpo->free = _free_cpo;
	return cpo;
}
