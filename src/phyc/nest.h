//
//  nest.h
//  physher
//
//  Created by Mathieu Fourment on 22/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef nest_h
#define nest_h

#include <stdio.h>

#include "mcmc.h"

typedef struct NEST{
	Model* prior;
	Model* likelihood;
	Operator** operators;
	size_t operator_count;
	long chain_length;
	double precision;
	size_t steps;
	size_t burnin;
	size_t N;
	void (*run)(struct NEST*);
//	char* log_file;
	Parameters* x;
	gsl_rng* rng;
} NEST;

NEST* new_NEST_from_json(json_node* node, Hashtable* hash);
void free_NEST(NEST* nest);
#endif /* nest_h */
