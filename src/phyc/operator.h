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
	bool (*propose)(struct Operator*, double*);
	void (*optimize)(struct Operator*, double);
	gsl_rng* rng;
}Operator;

Operator* new_Operator_from_json(json_node* node, Hashtable* hash);

void free_Operator(Operator* op);

#endif /* operator_h */
