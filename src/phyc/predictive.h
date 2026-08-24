// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef predictive_h
#define predictive_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "model.h"

#define JSON_PREDICTIVE "predictive"

struct _Predictive;
typedef struct _Predictive Predictive;

struct _Predictive{
	char* filename;
	Model* model;  // the tree likelihood the trace was logged from
	size_t burnin;
	void(*calculate)(Predictive*);
	void(*free)(Predictive*);
};

Predictive* new_Predictive_from_json(json_node* node, Hashtable* hash);

#endif /* predictive_h */
