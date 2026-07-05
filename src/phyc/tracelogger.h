// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef tracelogger_h
#define tracelogger_h

#include <stdio.h>
#include <sys/time.h>

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

typedef struct Trace{
	Model** models;
	size_t model_count;
	LogColumn* columns;
	size_t column_count;
	// tree annotation (tree loggers only)
	Model** trait_models;   // BranchModel sources for per-branch traits
	char** trait_tags;      // output key written per trait, e.g. "rate"
	char** trait_formats;   // printf format per trait, e.g. "%e" (default) or "%.10f"
	size_t trait_count;
	Model** scalar_models;  // whole-tree likelihoods, logged as [&name=logP] (nexus only)
	char** scalar_formats;  // printf format per annotation, e.g. "%e" (default)
	size_t scalar_count;
	FILE* file;
	char* filename;
	size_t every;
	bool append;
	void(*initialize)(struct Trace* logger);
	void(*finalize)(struct Trace* logger);
	void(*write)(struct Trace* logger, size_t);
	void(*free)(struct Trace*);
	bool cpo;
	bool tree;
	char* format;
	bool force;// force calculation of model (i.e. do not use stored lnl)
	struct timeval start;
	struct timeval end;
}Trace;
	
 Trace* new_Trace_from_json(json_node* node, Hashtable* hash);

#endif /* tracelogger_h */
