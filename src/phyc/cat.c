//
//  cat.c
//  physher
//
//  Created by Mathieu Fourment on 4/6/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "cat.h"

#include "parameters.h"
#include "treelikelihood.h"
#include "utils.h"
#include "matrix.h"


void fasttree_cat(SingleTreeLikelihood* tlk, bool verbosity){
	SiteModel* sm = tlk->sm;
	int cat_count = sm->cat_count;
	sm->cat_count = 1;
	int* counts = ivector(cat_count);
	double* rates = log_spaced_spaced_vector(1.0/cat_count, cat_count, cat_count);
	print_dvector(rates, cat_count);
	double* likelihoods = malloc(sizeof(double)*tlk->sp->count*cat_count);
	tlk->calculate(tlk);// make sure everything is up-to-date
	sm->need_update = false;
	for (int i = 0; i < cat_count; i++) {
		sm->cat_rates[0] = rates[i];
		SingleTreeLikelihood_update_all_nodes(tlk);
		tlk->calculate(tlk);
		
		memcpy(likelihoods+(i*tlk->sp->count), tlk->pattern_lk, sizeof(double)*tlk->sp->count);
	}
	
	for (int i = 0; i < tlk->sp->count; i++) {
		int best = 0;
		double bestLnl = -DBL_MAX;
		for (int j = 0; j < cat_count; j++) {
			double logP = likelihoods[j*tlk->sp->count+i];// + 2.0*log(rates[j]) - 3.0*rates[j]; // gamma prior as in fasttree
			if(logP > bestLnl){
				bestLnl = logP;
				best = j;
			}
		}
		sm->site_category[i] = best;
		if(verbosity > 0) printf("CAT pattern %d rate index %d rate %f [%f\n", i, best, rates[best], bestLnl);
		counts[best]++;
	}
	Parameters_set_value(sm->rates, 0, rates[0]);
	for (int i = 1; i < cat_count; i++) {
		Parameters_set_value_quietly(sm->rates, i, rates[i]);
	}
	if(verbosity > 0){
		for (int i = 0; i < cat_count; i++) {
			if(counts[i] == 0) fprintf(stdout, "CAT rate %f not used\n", rates[i]);
		}
	}
	SingleTreeLikelihood_use_rescaling(tlk, false);
	SingleTreeLikelihood_update_all_nodes(tlk);
	
	sm->cat_count = cat_count;
	free(likelihoods);
	free(rates);
	free(counts);
}

void cat_estimator_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"model",
		"rates",
		"verbosity"
	};
	//json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* ref = get_json_node_value_string(node, "model");
	Model* mtlk = Hashtable_get(hash, ref+1);
	int verbosity = get_json_node_value_int(node, "verbosity", 0);
	SingleTreeLikelihood *tlk = mtlk->obj;
	fasttree_cat(tlk, verbosity);
}
