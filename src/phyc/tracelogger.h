// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef tracelogger_h
#define tracelogger_h

#include <stdio.h>
#include <sys/time.h>

#include "parameters.h"

// One ordered output column. Exactly one of the two pointers is non-NULL:
//   model     != NULL -> log this Model's logP, unless `quantity` is set, in
//                        which case log the model's named loggable quantity via
//                        its log_* interface  (JSON ref "@id")
//   parameter != NULL -> log this Parameter's value(s) (JSON ref "&id" or "%id")
// format is an optional per-column printf conversion; NULL falls back to the
// logger-level format. name is an optional header/label override; NULL falls
// back to the model/parameter's own name. quantity (model columns only) names a
// derived value to log instead of logP; loggable_count caches the model's
// log_count for that quantity, resolved at parse time.
typedef struct LogColumn{
	Model* model;
	Parameter* parameter;
	char* format;
	char* name;
	char* quantity;
	size_t loggable_count;
}LogColumn;

// Parse the "columns" array into a freshly allocated ordered array; returns the
// number of columns (0 and *columns==NULL when the node is absent). Each item is
// either a bare ref string ("@id", "&id" or "%id") or an object
// {"ref": "@id", "name": "lnL", "format": "%f"} with an optional per-column
// format and header/label name. "%id"/vector refs expand to one column per
// element (sharing the item's format); a "name" override requires a single ref.
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
	void(*report)(struct Trace* logger);
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
