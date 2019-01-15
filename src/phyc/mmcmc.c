//
//  mmcmc.c
//  physher
//
//  Created by Mathieu Fourment on 17/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "mmcmc.h"

#include <string.h>
#include <strings.h>

#include "matrix.h"
#include "beta.h"
#include "compoundmodel.h"

void mmcmc_run(MMCMC* mmcmc){
	StringBuffer* buffer = new_StringBuffer(10);
	MCMC* mcmc = mmcmc->mcmc;
	char** filenames = malloc(sizeof(char*)*mcmc->log_count);
	
	for (int j = 0; j < mcmc->log_count; j++) {
		filenames[j] = NULL;
		if(mcmc->logs[j]->filename != NULL){
			filenames[j] = String_clone(mcmc->logs[j]->filename);
		}
	}
	int i = 0;
	if (mmcmc->gss) {
		mcmc->gss = true;
		i = 1;
	}
	// Check that we can sample directly from prior withtou MCMC
	bool sampleable = mmcmc->prior_samples > 0;
	CompoundModel* cm = (CompoundModel*)(mcmc->model->obj);
	for (int i = 1; i < cm->count; i++) {
		sampleable &= cm->models[i]->samplable;
	}
	
	// temperatures should be in decreasing order
	for (; i < mmcmc->temperature_count; i++) {
		for (int j = 0; j < mcmc->log_count; j++) {
			if(filenames[j] != NULL){
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "%d%s",i, filenames[j]);
				free(mcmc->logs[j]->filename);
				mcmc->logs[j]->filename = StringBuffer_tochar(buffer);
			}
		}
		printf("Temperature: %e - %s\n", mmcmc->temperatures[i],buffer->c);
		mcmc->chain_temperature = mmcmc->temperatures[i];
		
		if(sampleable && mcmc->chain_temperature == 0){
			// Initialize loggers
			for (int j = 0; j < mcmc->log_count; j++) {
				mcmc->logs[j]->initialize(mcmc->logs[j]);
			}
			
			for (size_t s = 0; s < mmcmc->prior_samples; s++) {
				for (int j = 1; j < cm->count; j++) {
					cm->models[j]->sample_evaluate(cm->models[j]);
				}
				mcmc->model->logP(mcmc->model);
				for (int j = 0; j < mcmc->log_count; j++) {
					// log every sample to file but not to stdout
					if(filenames[j] != NULL || (s % mcmc->logs[j]->every == 0 && filenames[j] == NULL)){
						mcmc->logs[j]->write(mcmc->logs[j], s);
					}
				}
			}
			for (int j = 0; j < mcmc->log_count; j++) {
				mcmc->logs[j]->finalize(mcmc->logs[j]);
			}
		}
		else{
			mcmc->run(mcmc);
		}
	}
	
	// Leave it as a standard mcmc with original loggers
	mcmc->chain_temperature = -1;
	for (int j = 0; j < mcmc->log_count; j++) {
		if(filenames[j] != NULL){
			free(mcmc->logs[j]->filename);
			mcmc->logs[j]->filename = String_clone(filenames[j]);
			free(filenames[j]);
		}
	}
	mcmc->gss = false;
	
	free(filenames);
	free_StringBuffer(buffer);
}


static void _free_MMCMC(MMCMC* mmcmc){
	mmcmc->mcmc->free(mmcmc->mcmc);
	free(mmcmc->temperatures);
	free(mmcmc);
}

MMCMC* new_MMCMC_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"distribution",
		"gss",
		"mcmc",
		"samples",
		"steps",
		"temperatures"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	MMCMC* mmcmc = malloc(sizeof(MMCMC));
	json_node* mcmc_node = get_json_node(node, "mcmc");
	json_node* temp_node = get_json_node(node, "temperatures");
	json_node* steps_node = get_json_node(node, "steps");
	mmcmc->gss = get_json_node_value_bool(node, "gss", false);
	mmcmc->prior_samples = get_json_node_value_size_t(node, "samples", 0);
	
	if (temp_node != NULL) {
		if (temp_node->node_type != MJSON_ARRAY) {
			fprintf(stderr, "attribute `temperatures` should be an array\n\n");
			exit(1);
		}
		mmcmc->temperature_count = temp_node->child_count;
		mmcmc->temperatures = dvector(mmcmc->temperature_count);
		for (int i = 0; i < temp_node->child_count; i++) {
			mmcmc->temperatures[i] = atof((char*)temp_node->children[i]->value);
		}
	}
	else if (steps_node != NULL){
		char* dist_string = get_json_node_value_string(node, "distribution");
		mmcmc->temperature_count = get_json_node_value_size_t(node, "steps", 100);
		mmcmc->temperatures = dvector(mmcmc->temperature_count);
		mmcmc->temperatures[0] = 1;
		mmcmc->temperatures[mmcmc->temperature_count-1] = 0;
		
		// temperatures are in descreasing order
		if(strcasecmp(dist_string, "beta") == 0){
			double alpha = get_json_node_value_double(node, "alpha", 0.3);
			double beta = get_json_node_value_double(node, "beta", 1.0);
			double value = 1;
			double incr = 1.0/(mmcmc->temperature_count-1);
			for (size_t i = 1; i < mmcmc->temperature_count-1; i++) {
				value -= incr;
				mmcmc->temperatures[i] = invbetai(value, alpha, beta);
			}
		}
		else if(strcasecmp(dist_string, "uniform") == 0){
			double incr = 1.0/(mmcmc->temperature_count-1);
			for (size_t i = 1; i < mmcmc->temperature_count-1; i++) {
				mmcmc->temperatures[i] = mmcmc->temperatures[i+1] - incr;
			}
		}
		else{
			fprintf(stderr, "Attribute `distribution` should be specified (`uniform` or `beta`)\n\n");
			exit(1);
		}
	}
	else{
		fprintf(stderr, "Attribute `temperatures` or `steps` should be specified\n\n");
		exit(1);
	}
	mmcmc->mcmc = new_MCMC_from_json(mcmc_node, hash);
	mmcmc->run = mmcmc_run;
	mmcmc->free = _free_MMCMC;
	
	if (mmcmc->gss) {
		CompoundModel* cm = mmcmc->mcmc->model->obj;
		if (cm->count != 2) {
			fprintf(stderr, "Exactly 2 models (posterior and reference distributions) should be specified for GSS\n");
			exit(2);
		}
	}
	
	return mmcmc;
}
