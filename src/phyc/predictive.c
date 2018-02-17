//
//  predictive.c
//  physher
//
//  Created by Mathieu Fourment on 16/02/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "predictive.h"

#include <float.h>

#include "mstring.h"
#include "matrix.h"
#include "filereader.h"

// computed log pointwise predictive density
void _predictive_calculate(Predictive* predictive){
	char *ptr = NULL;
	double *temp = NULL;
	int l = 0;
	size_t capacity = 1000;
	size_t count = 0;
	Vector** vecs = malloc(capacity*sizeof(Vector*));
	size_t sample = 0;
	double* weights = NULL;
	FileReader *reader = new_FileReader(predictive->filename, 1000);
	reader->read_line(reader);
	weights = String_split_char_double( reader->line+1, '\t', &l );
	reader->read_line(reader);// discard header
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		
		if ( reader->buffer->length == 0){
			continue;
		}
		if ( sample >= predictive->burnin){
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
	double llpd = 0;
	for (int i = 0; i < n; i++) {
		double sum = Vector_at(vecs[0], i);
		for (int j = 1; j < count; j++) {
			sum = logaddexp(sum, Vector_at(vecs[j], i));
		}
		llpd += (sum - log(count))*weights[i];
	}
	printf("lpd: %f\n", llpd);
	for (int i = 0; i < count; i++) {
		free_Vector(vecs[i]);
	}
	free(vecs);
	free(weights);
}

void _free_predictive(Predictive* predictive){
	free(predictive->filename);
	free(predictive);
}


Predictive* new_Predictive_from_json(json_node* node, Hashtable* hash){
	char* filename = get_json_node_value_string(node, "filename");
	Predictive* predictive = malloc(sizeof(Predictive));
	predictive->filename = String_clone(filename);
	predictive->burnin = get_json_node_value_size_t(node, "burnin", 0);
	predictive->calculate = _predictive_calculate;
	predictive->free = _free_predictive;
	return predictive;
}