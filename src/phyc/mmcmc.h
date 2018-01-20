//
//  mmcmc.h
//  physher
//
//  Created by Mathieu Fourment on 17/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef mmcmc_h
#define mmcmc_h

#include "mcmc.h"

#include <gsl/gsl_rng.h>

typedef struct MMCMC{
	MCMC* mcmc;
	double* temperatures;
	size_t temperature_count;
	size_t burnin;
	void (*run)(struct MMCMC*);
	// for marginal likelihood calculation
	char* log_file;
	Parameters* x; // for bridge sampling
	gsl_rng* rng;
} MMCMC;


MMCMC* new_MMCMC_from_json(json_node* node, Hashtable* hash);

void free_MMCMC(MMCMC* mmcmc);

#endif /* mmcmc_h */
