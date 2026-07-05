// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef tracelogger_h
#define tracelogger_h

#include <stdio.h>
#include <sys/time.h>

#include "parameters.h"
#include "logger.h"

typedef struct Trace{
	Parameters* x;
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
	void(*write_with)(struct Trace* logger, size_t, const char*);
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
