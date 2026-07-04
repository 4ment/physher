// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef logger_h
#define logger_h

#include <stdio.h>

#include "parameters.h"

// One ordered output column. Exactly one of the two pointers is non-NULL:
//   model     != NULL -> log this Model's logP        (JSON ref "@id")
//   parameter != NULL -> log this Parameter's value(s) (JSON ref "&id" or "%id")
typedef struct LogColumn{
	Model* model;
	Parameter* parameter;
}LogColumn;

// Parse the "columns" array (each item "@id", "&id" or "%id") into a freshly
// allocated ordered array; returns the number of columns (0 and *columns==NULL
// when the node is absent). "%id"/vector refs expand to one column per element.
size_t get_columns_from_json(json_node* node, Hashtable* hash, LogColumn** columns);

typedef struct Logger{
	Parameters* parameters;
	Model** models;
	size_t model_count;
	LogColumn* columns;
	size_t column_count;
	void (*log)(struct Logger*);
	FILE* file;
	char* filename;
	char* format;
	bool tree;
	bool internal; // show internal node name
    char sep; //separator
}Logger;

struct Dumper{
    Parameters** parameters;
    size_t parameter_count;
    Model** models;
    size_t model_count;
    void (*dump)(struct Dumper*);
    void (*free)(struct Dumper*);
    FILE* file;
    char* filename;
};


struct Logger* new_logger_from_json(json_node* node, Hashtable* hash);

void free_Logger(struct Logger* logger);


struct Dumper* new_Dumper_from_json(json_node* node, Hashtable* hash);

void free_Dumper(struct Dumper* logger);

#endif /* logger_h */
