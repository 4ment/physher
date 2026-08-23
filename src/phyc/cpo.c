// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "cpo.h"

#include <float.h>
#include <string.h>

#include "mstring.h"
#include "matrix.h"
#include "filereader.h"
#include "treelikelihood.h"

// Locate the per-pattern log-likelihood columns in a trace file's header.
//
// The logger that writes them is an ordinary column logger, so the file may
// carry anything else beside them (the posterior, a tree length, another
// likelihood). Columns are therefore matched by *name* rather than by position:
// the loggable interface names them "<model>.patternLogLikelihood.<i>"
// (_TreeLikelihoodModel_log_name), so every pattern of the requested model is
// found wherever it sits and everything else is skipped.
//
// Returns a malloc'd array of `pattern_count` field indices, or dies if any
// pattern is missing.
static int* _cpo_map_columns(const char* header, const char* model_name,
                             int pattern_count){
	int field_count = 0;
	char** fields = String_split_char(header, '\t', &field_count);
	int* map = ivector(pattern_count);
	for (int i = 0; i < pattern_count; i++) map[i] = -1;

	StringBuffer* buffer = new_StringBuffer(64);
	for (int i = 0; i < pattern_count; i++) {
		StringBuffer_empty(buffer);
		StringBuffer_append_format(buffer, "%s.patternLogLikelihood.%d",
		                           model_name, i);
		for (int j = 0; j < field_count; j++) {
			if (strcmp(fields[j], buffer->c) == 0) {
				map[i] = j;
				break;
			}
		}
		if (map[i] == -1) {
			fprintf(stderr, "cpo: no column named '%s' in the trace file. The "
			                "logger must carry {\"ref\": \"@%s\", \"quantity\": "
			                "\"patternLogLikelihood\"}\n", buffer->c, model_name);
			exit(1);
		}
	}
	free_StringBuffer(buffer);
	for (int j = 0; j < field_count; j++) free(fields[j]);
	free(fields);
	return map;
}

void _cpo_calculate(struct CPO* cpo){
	SingleTreeLikelihood* tlk = cpo->model->obj;
	const int pattern_count = tlk->sp->count;
	const double* weights = tlk->sp->weights;

	size_t capacity = 1000;
	size_t count = 0;
	Vector** vecs = malloc(capacity*sizeof(Vector*));
	size_t sample = 0;
	FileReader *reader = new_FileReader(cpo->filename, 1000);

	// The weights used to be a comment line at the top of the file; they are now
	// read off the site pattern the log was written from, which is the same
	// object the logger used and cannot disagree with it.
	if (!reader->read_line(reader)) {
		fprintf(stderr, "cpo: '%s' is empty\n", cpo->filename);
		exit(1);
	}
	StringBuffer_trim(reader->buffer);
	int* map = _cpo_map_columns(reader->line, cpo->model->name, pattern_count);

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
			int l = 0;
			double* temp = String_split_char_double( reader->line, '\t', &l );
			vecs[count] = new_Vector(pattern_count);
			for (int i = 0; i < pattern_count; i++) {
				if (map[i] >= l) {
					fprintf(stderr, "cpo: row %zu of '%s' has %d fields, too few "
					                "for column %d\n", sample, cpo->filename, l,
					        map[i]);
					exit(1);
				}
				Vector_push(vecs[count], temp[map[i]]);
			}
			free(temp);
			count++;
		}
		sample++;
	}
	free_FileReader(reader);
	free(map);

	if (count == 0) {
		fprintf(stderr, "cpo: no samples left in '%s' after a burnin of %zu\n",
		        cpo->filename, cpo->burnin);
		exit(1);
	}

	double logCPO = 0;
	for (int i = 0; i < pattern_count; i++) {
		double sum = -DBL_MAX;
		double min = Vector_at(vecs[0], i);
		for (int j = 1; j < count; j++) {
			min = fmin(Vector_at(vecs[j], i), min);
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
}

void _free_cpo(CPO* cpo){
	free(cpo->filename);
	free(cpo);
}


CPO* new_CPO_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"burnin", JSON_OPTIONAL, JSON_NUMBER},
	    {"filename", JSON_REQUIRED, JSON_STRING},
	    {"model", JSON_REQUIRED, JSON_STRING},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	
	char* filename = get_json_node_value_string(node, "filename");
	char* ref = get_json_node_value_string(node, "model");
	// The likelihood the trace was logged from: it supplies the pattern weights
	// and the column names to look for, so the file needs to carry neither.
	Model* model = Hashtable_get(hash, ref + 1);
	if (model == NULL || model->type != MODEL_TREELIKELIHOOD) {
		json_die(node, "\"model\" must reference a tree likelihood (e.g. \"@treelikelihood\"): '%s'", ref);
	}

	CPO* cpo = malloc(sizeof(CPO));
	cpo->filename = String_clone(filename);
	cpo->model = model;
	cpo->burnin = get_json_node_value_size_t(node, "burnin", 0);
	cpo->calculate = _cpo_calculate;
	cpo->free = _free_cpo;
	return cpo;
}
