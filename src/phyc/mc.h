// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef mc_h
#define mc_h

#include <stdio.h>

#include "mjson.h"
#include "parameters.h"

typedef struct MC{
	Model* likelihood;
	Model* prior;
	Parameters* parameters;
	size_t samples;
	double(*calculate)(struct MC*);
	void(*free)(struct MC*);
}MC;

MC* new_MonteCarlo_from_json(json_node* node, Hashtable* hash);

#endif /* mc_h */
