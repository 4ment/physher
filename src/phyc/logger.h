// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef logger_h
#define logger_h

#include <stdio.h>

#include "parameters.h"

typedef struct Logger{
	Parameters* parameters;
	Model** models;
	size_t model_count;
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
