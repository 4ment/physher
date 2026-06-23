//
//  operator.h
//  physher
//
//  Created by Mathieu Fourment on 4/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef operator_h
#define operator_h

#include <stdio.h>

#include "parameters.h"

#include <gsl/gsl_rng.h>

typedef struct Operator{
	char* name;
	Parameters* x;
	Model** models;
	size_t model_count;
	double* parameters;
	//int* indexes;
	double weight;
	bool all;
	size_t rejected_count;
	size_t accepted_count;
	size_t failure_count;
	size_t tuning_delay;
	size_t accepted_at_delay; // accepted_count snapshot when tuning began
	size_t count_at_delay;    // accepted+rejected snapshot when tuning began
	bool tuning_started;
	double target;
	bool (*propose)(struct Operator*, double*);
	void (*optimize)(struct Operator*, double);
	gsl_rng* rng;
}Operator;

Operator* new_Operator_from_json(json_node* node, Hashtable* hash);

// Post-delay tuning statistics. Returns false while still inside the tuning
// delay window. On the first call past the delay it snapshots the baseline
// counts, so *count and *accepted measure only proposals made since tuning
// started (an unbiased acceptance ratio is *accepted / *count).
bool operator_tuning_stats(Operator* op, long* count, long* accepted);

void free_Operator(Operator* op);

#endif /* operator_h */
