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

typedef struct MMCMC{
	MCMC* mcmc;
	double* temperatures;
	size_t temperature_count;
	void (*run)(struct MMCMC*);
	void (*free)(struct MMCMC*);
	bool gss;
	bool bf;
	size_t prior_samples;
} MMCMC;


MMCMC* new_MMCMC_from_json(json_node* node, Hashtable* hash);

#endif /* mmcmc_h */
