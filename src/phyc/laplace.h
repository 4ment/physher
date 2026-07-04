// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef laplace_h
#define laplace_h

#include <stdio.h>

#include "parameters.h"

typedef struct Laplace{
	Model* model;
	Parameters* parameters;
	Model* refdist;
    Model* empirical;
	double(*calculate)(struct Laplace*);
	void(*free)(struct Laplace*);
	// something print
}Laplace;

Laplace* new_Laplace_from_json(json_node* node, Hashtable* hash);

#endif /* laplace_h */
