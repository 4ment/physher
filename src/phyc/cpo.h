// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef cpo_h
#define cpo_h

#include "mjson.h"
#include "hashtable.h"
#include "model.h"

typedef struct CPO{
	char* filename;
	Model* model;  // the tree likelihood the trace was logged from
	size_t burnin;
	void(*calculate)(struct CPO*);
	void(*free)(struct CPO*);
	// something print
}CPO;

CPO* new_CPO_from_json(json_node* node, Hashtable* hash);

//void print_hessian(Hessian* hessian);

#endif /* cpo_h */
