//
//  mc.c
//  physher
//
//  Created by Mathieu Fourment on 3/5/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "mc.h"

double _naive_monte_carlo_calculate(MC* mc){
	size_t samples = mc->samples;
	Model* prior = mc->prior;
	Model* likelihood = mc->likelihood;
	double sum = -DBL_MAX;
	for (size_t i = 0; i < samples; i++) {
		prior->sample(prior);
		sum = logaddexp(sum, likelihood->logP(likelihood));
	}
	double logP = sum - log(samples);
	printf("Monte Carlo: %f\n", logP);
	return logP;
}

void _free_MC(MC* mc){
	free_Parameters(mc->parameters);
	mc->likelihood->free(mc->likelihood);
	mc->prior->free(mc->prior);
	free(mc);
}

MC* new_MonteCarlo_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"likelihood",
		"parameters",
		"prior",
		"samples"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	MC* mc = malloc(sizeof(MC));
	char* ref_lk = get_json_node_value_string(node, "likelihood");
	char* ref_prior = get_json_node_value_string(node, "prior");
	mc->samples = get_json_node_value_size_t(node, "samples", 1000);
	mc->likelihood = Hashtable_get(hash, ref_lk+1);
	mc->likelihood->ref_count++;
	mc->prior = Hashtable_get(hash, ref_prior+1);
	mc->prior->ref_count++;
	mc->calculate = _naive_monte_carlo_calculate;
	mc->free = _free_MC;
	mc->parameters = new_Parameters(1);
	get_parameters_references(node, hash, mc->parameters);
	return mc;
}
