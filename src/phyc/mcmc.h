//
//  mcmc.h
//  physher
//
//  Created by Mathieu Fourment on 2/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef mcmc_h
#define mcmc_h

#include <stdio.h>

#include "parameters.h"

typedef struct Log{
	Parameters* x;
	Model** models;
	size_t model_count;
	Model** simplexes;
	size_t simplex_count;
	FILE* file;
	char* filename;
	size_t every;
	bool append;
	void(*write)(struct Log* logger, size_t);
}Log;

typedef struct Operator{
	char* name;
	Parameters* x;
	Model* simplex;
	int index;
	double* storage;
	double* parameters;
	double weight;
	size_t rejected_count;
	size_t accepted_count;
	bool (*propose)(struct Operator*, double*);
	void (*store)(struct Operator*);
	void (*restore)(struct Operator*);
	void (*optimize)(struct Operator*, double);
}Operator;

typedef struct MCMC{
	Model* model;
	void (*run)(struct MCMC*);
	Operator** operators;
	size_t operator_count;
	size_t chain_length;
	Log** logs;
	size_t log_count;
	int verbose;
} MCMC;

MCMC* new_MCMC_from_json(json_node* node, Hashtable* hash);

void free_MCMC(MCMC* mcmc);

#endif /* mcmc_h */
