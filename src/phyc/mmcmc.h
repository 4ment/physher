// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
	int start;
	size_t prior_samples;
} MMCMC;


MMCMC* new_MMCMC_from_json(json_node* node, Hashtable* hash);

#endif /* mmcmc_h */
