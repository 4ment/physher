//
//  mmcmc.c
//  physher
//
//  Created by Mathieu Fourment on 17/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "mmcmc.h"

#include "matrix.h"

void mmcmc_run(MMCMC* mmcmc){
	StringBuffer* buffer = new_StringBuffer(10);
			MCMC* mcmc = mmcmc->mcmc;
	for (int i = 0; i < mmcmc->temperature_count; i++) {
		for (int j = 0; j < mcmc->log_count; j++) {
			if(mcmc->logs[j]->filename != NULL){
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "%f-", mmcmc->tempratures[i]);
				StringBuffer_append_string(buffer, mcmc->logs[j]->filename);
				fclose(mcmc->logs[j]->file);
				mcmc->logs[j]->file = fopen(buffer->c, "w");
			}
		}
		mcmc->chain_temperature = mmcmc->tempratures[i];
		mcmc->run(mcmc);
	}
	// Leave it as a standard mcmc
	mcmc->chain_temperature = -1;
	free_StringBuffer(buffer);
}

MMCMC* new_MMCMC_from_json(json_node* node, Hashtable* hash){
	MMCMC* mmcmc = malloc(sizeof(MMCMC));
	json_node* mcmc_node = get_json_node(node, "mcmc");
	json_node* temp_node = get_json_node(node, "temperatures");
	
	if (temp_node->node_type == MJSON_ARRAY) {
		mmcmc->temperature_count = temp_node->child_count;
		mmcmc->tempratures = dvector(mmcmc->temperature_count);
		for (int i = 0; i < temp_node->child_count; i++) {
			mmcmc->tempratures[i] = atof((char*)temp_node->children[i]->value);
		}
	}
	mmcmc->mcmc = new_MCMC_from_json(mcmc_node, hash);
	mmcmc->run = mmcmc_run;
	return mmcmc;
}

void free_MMCMC(MMCMC* mmcmc){
	free_MCMC(mmcmc->mcmc);
	free(mmcmc->tempratures);
	free(mmcmc);
}