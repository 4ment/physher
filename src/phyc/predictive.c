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
#include "statistics.h"

// computed log pointwise predictive density
void _predictive_calculate(Predictive* predictive){
	char *ptr = NULL;
	double *temp = NULL;
	int l = 0;
	size_t count = 0;
	size_t sample = 0;
	double* weights = NULL;
	FileReader *reader = new_FileReader(predictive->filename, 1000);
	reader->read_line(reader);
	weights = String_split_char_double( reader->line+1, '\t', &l );
	size_t siteCount = l;
	Vector** vecs = malloc(siteCount*sizeof(Vector*));
	for (int i = 0; i < siteCount; i++) {
		vecs[i] = new_Vector(100);
	}
	reader->read_line(reader);// discard header
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		
		if ( reader->buffer->length == 0){
			continue;
		}
		if ( sample >= predictive->burnin){
			ptr = reader->line;
			l = 0;
			//printf("=%s=\n", ptr);
			temp = String_split_char_double( ptr, '\t', &l );
			for (int i = 0; i < siteCount; i++) {
				Vector_push(vecs[i], temp[i+1]);
			}
			free(temp);
			count++;
		}
		sample++;
	}
	free_FileReader(reader);
	
	double llpd = 0;
	for (int i = 0; i < siteCount; i++) {
		double sum = Vector_at(vecs[i], 0);
		for (int j = 1; j < count; j++) {
			sum = logaddexp(sum, Vector_at(vecs[i], j));
		}
		llpd += (sum - log(count))*weights[i];
	}
	printf("lppd: %f\n", llpd);
	
	double pwaic = 0;
	for (int i = 0; i < siteCount; i++) {
		double m = mean(Vector_data(vecs[i]), count);
		pwaic += weights[i]*variance(Vector_data(vecs[i]), count, m);
	}
	printf("pwaic: %f\n", pwaic);
	
	for (int i = 0; i < siteCount; i++) {
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
	char* allowed[] = {
		"burnin",
		"filename"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* filename = get_json_node_value_string(node, "filename");
	Predictive* predictive = malloc(sizeof(Predictive));
	predictive->filename = String_clone(filename);
	predictive->burnin = get_json_node_value_size_t(node, "burnin", 0);
	predictive->calculate = _predictive_calculate;
	predictive->free = _free_predictive;
	return predictive;
}