//
//  marginal.h
//  physher
//
//  Created by Mathieu Fourment on 18/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef marginal_h
#define marginal_h

#include <stdio.h>

#include "matrix.h"
#include "mjson.h"
#include "hashtable.h"

typedef struct MarginaLikelihood{
	char* file;
	double* temperatures;
	size_t temperature_count;
	size_t burnin;
	char* likelihood_tag;
	char* refdist_tag;
	void (*run)(struct MarginaLikelihood*);
	void (*free)(struct MarginaLikelihood*);
	
}MarginaLikelihood;

double log_arithmetic_mean(const Vector* vecvalues);

double log_harmonic_mean(const Vector* values);

double log_stablilized_harmonic_mean(double guess, const Vector* values, double delta);

double log_marginal_stepping_stone(const Vector** values, size_t temp_count, const double* temperatures, double* lrssk);

double log_marginal_path_sampling(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk);

double log_marginal_path_sampling_modified(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk);

MarginaLikelihood* new_MarginaLikelihood_from_json(json_node* node, Hashtable* hash);

#endif /* marginal_h */
