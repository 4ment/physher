//
//  mcmc.h
//  physher
//
//  Created by Mathieu Fourment on 2/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef mcmc_h
#define mcmc_h

#include "parameters.h"
#include "operator.h"
#include "logmcmc.h"

typedef struct MCMC{
	Model* model;
	void (*run)(struct MCMC*);
	void (*free)(struct MCMC*);
	Operator** operators;
	size_t operator_count;
	size_t chain_length;
	double chain_temperature;
	Log** logs;
	size_t log_count;
	size_t tuning_frequency;
	int verbose;
	bool gss;
} MCMC;

MCMC* new_MCMC_from_json(json_node* node, Hashtable* hash);

#endif /* mcmc_h */
