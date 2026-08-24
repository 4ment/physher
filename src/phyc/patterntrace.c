// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "patterntrace.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "filereader.h"
#include "matrix.h"
#include "mstring.h"
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
static int* _map_columns(const char* header, const char* model_name,
                         int pattern_count, const char* action){
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
			fprintf(stderr, "%s: no column named '%s' in the trace file. The "
			                "logger must carry {\"ref\": \"@%s\", \"quantity\": "
			                "\"patternLogLikelihood\"}\n", action, buffer->c,
			        model_name);
			exit(1);
		}
	}
	free_StringBuffer(buffer);
	for (int j = 0; j < field_count; j++) free(fields[j]);
	free(fields);
	return map;
}

double** read_pattern_log_likelihoods(Model* model, const char* filename,
                                      size_t burnin, const char* action,
                                      size_t* sample_count){
	SingleTreeLikelihood* tlk = model->obj;
	const int pattern_count = tlk->sp->count;

	Vector** vecs = malloc(pattern_count*sizeof(Vector*));
	for (int i = 0; i < pattern_count; i++) {
		vecs[i] = new_Vector(1000);
	}

	size_t count = 0;
	size_t sample = 0;
	FileReader *reader = new_FileReader(filename, 1000);

	// The weights used to be a comment line at the top of the file; they are now
	// read off the site pattern the log was written from, which is the same
	// object the logger used and cannot disagree with it.
	if (!reader->read_line(reader)) {
		fprintf(stderr, "%s: '%s' is empty\n", action, filename);
		exit(1);
	}
	StringBuffer_trim(reader->buffer);
	int* map = _map_columns(reader->line, model->name, pattern_count, action);

	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);

		if ( reader->buffer->length == 0){
			continue;
		}
		if ( sample >= burnin){
			int l = 0;
			double* temp = String_split_char_double( reader->line, '\t', &l );
			for (int i = 0; i < pattern_count; i++) {
				if (map[i] >= l) {
					fprintf(stderr, "%s: row %zu of '%s' has %d fields, too few "
					                "for column %d\n", action, sample, filename,
					        l, map[i]);
					exit(1);
				}
				Vector_push(vecs[i], temp[map[i]]);
			}
			free(temp);
			count++;
		}
		sample++;
	}
	free_FileReader(reader);
	free(map);

	if (count == 0) {
		fprintf(stderr, "%s: no samples left in '%s' after a burnin of %zu\n",
		        action, filename, burnin);
		exit(1);
	}

	double** trace = malloc(pattern_count*sizeof(double*));
	for (int i = 0; i < pattern_count; i++) {
		trace[i] = malloc(count*sizeof(double));
		memcpy(trace[i], Vector_data(vecs[i]), count*sizeof(double));
		free_Vector(vecs[i]);
	}
	free(vecs);

	*sample_count = count;
	return trace;
}
