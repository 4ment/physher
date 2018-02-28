//
//  mcmc.c
//  physher
//
//  Created by Mathieu Fourment on 2/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "mcmc.h"

#include <strings.h>
#include <stdio.h>

#include "random.h"
#include "matrix.h"
#include "compoundmodel.h"
#include "treelikelihood.h"


void run(MCMC* mcmc){
	Model* model = mcmc->model;
	
	size_t iter = 0;
	
	int tuneFrequency = 10;
	
	double* weights = dvector(mcmc->operator_count);
	double sum = 0;
	for (int i = 0; i < mcmc->operator_count; i++) {
		sum += mcmc->operators[i]->weight;
	}
	for (int i = 0; i < mcmc->operator_count; i++) {
		weights[i] = mcmc->operators[i]->weight/sum;
	}
	
	StringBuffer *buffer = new_StringBuffer(20);
	double logP;
	
	// Initialize loggers
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->initialize(mcmc->logs[i]);
	}
	
	// loggers log headers
	if(mcmc->chain_temperature < 0){
		logP = model->logP(model);
	}
	else if(mcmc->gss){
		model->logP(model);
		CompoundModel* cm = (CompoundModel*)model->obj;
		// Un-normalized posterior
		double posteriorLogP = cm->models[0]->logP(cm->models[0]) * mcmc->chain_temperature;
		
		// Reference distribution
		double referenceLogP = cm->models[1]->logP(cm->models[1]) * (1.0 - mcmc->chain_temperature);
		
		logP = posteriorLogP + referenceLogP;
	}
	else{
		model->logP(model);
		CompoundModel* cm = (CompoundModel*)model->obj;
		logP = cm->models[0]->logP(cm->models[0]) * mcmc->chain_temperature;
		for (int i = 1; i < cm->count; i++) {
			logP += cm->models[i]->logP(cm->models[i]);
		}
	}
	
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->write(mcmc->logs[i], iter);
	}
	
	while ( iter < mcmc->chain_length) {
		
		// Choose operator
		int op_index = roulette_wheel(weights, mcmc->operator_count);
		Operator* op = mcmc->operators[op_index];
		
		double logHR = 0;
		
		model->store(model);
		
		op->store(op);
		bool success = op->propose(op, &logHR);
		
		if (success) {
			double proposed_logP;
			if(mcmc->chain_temperature < 0){
				proposed_logP = model->logP(model);
			}
			else if(mcmc->gss){
				model->logP(model);
				CompoundModel* cm = (CompoundModel*)model->obj;
				double posteriorLogP = cm->models[0]->logP(cm->models[0]) * mcmc->chain_temperature;
				
				double referenceLogP = cm->models[1]->logP(cm->models[1]) * (1.0 - mcmc->chain_temperature);
				
				proposed_logP = posteriorLogP + referenceLogP;
			}
			else{
				model->logP(model);
				CompoundModel* cm = (CompoundModel*)model->obj;
				proposed_logP = cm->models[0]->logP(cm->models[0]) * mcmc->chain_temperature;
				for (int i = 1; i < cm->count; i++) {
					proposed_logP += cm->models[i]->logP(cm->models[i]);
				}
			}
			
			double alpha = proposed_logP - logP + logHR;
			
			// accept
			if ( alpha >=  0 || alpha > log(random_double()) ) {
	//			printf("%zu %f %f\n", iter, logP, proposed_logP);
				logP = proposed_logP;
				op->accepted_count++;
			}
			// reject
			else {
	//			printf("%zu %f %f *\n", iter, logP, proposed_logP);
				model->restore(model);
				op->restore(op);
				op->rejected_count++;
			}
			
			if(op->optimize != NULL){// && iter % tuneFrequency == 0){
				op->optimize(op, alpha);
			}
		}
		iter++;
		
		for (int i = 0; i < mcmc->log_count; i++) {
			if(iter % mcmc->logs[i]->every == 0){
				mcmc->logs[i]->write(mcmc->logs[i], iter);
			}
		}
	}
	
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->finalize(mcmc->logs[i]);
	}
	
	if( mcmc->verbose > 0){
		for (int i = 0; i < mcmc->operator_count; i++) {
			Operator* op = mcmc->operators[i];
			printf("Acceptance ratio %s: %f %f (failures: %zu)\n", op->name, ((double)op->accepted_count/(op->accepted_count+op->rejected_count)), op->parameters[0], op->failure_count);
		}
	}
	free(weights);
	free_StringBuffer(buffer);
}

void _free_MCMC(MCMC* mcmc){
	mcmc->model->free(mcmc->model);
	for (int i = 0; i < mcmc->operator_count; i++) {
		free_Operator(mcmc->operators[i]);
	}
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->free(mcmc->logs[i]);
	}
	free(mcmc->logs);
	free(mcmc->operators);
	free(mcmc);
}

MCMC* new_MCMC_from_json(json_node* node, Hashtable* hash){
	MCMC* mcmc = malloc(sizeof(MCMC));
	json_node* model_node = get_json_node(node, "model");
	mcmc->operator_count = 0;
	json_node* ops = get_json_node(node, "operators");
	json_node* logs = get_json_node(node, "log");
	mcmc->chain_length = get_json_node_value_size_t(node, "length", 100000);
	mcmc->verbose = get_json_node_value_int(node, "verbose", 1);
	mcmc->run = run;
	
	if (ops->child_count == 0) {
		fprintf(stderr, "Please specify at least one operator\n");
		exit(1);
	}
	
	if(model_node->node_type == MJSON_OBJECT){
		char* type = get_json_node_value_string(model_node, "type");
		char* id = get_json_node_value_string(model_node, "id");
		
		if (strcasecmp(type, "compound") == 0) {
			mcmc->model = new_CompoundModel_from_json(model_node, hash);
		}
		else if (strcasecmp(type, "treelikelihood") == 0) {
			mcmc->model = new_TreeLikelihoodModel_from_json(model_node, hash);
		}
		Hashtable_add(hash, id, mcmc->model);
		// could be a parametric distribution
	}
	// ref
	else if(model_node->node_type == MJSON_STRING){
		char* model_string = model_node->value;
		mcmc->model = Hashtable_get(hash, model_string+1);
		mcmc->model->ref_count++;
	}
	
	for (int i = 0; i < ops->child_count; i++) {
		json_node* child = ops->children[i];
		Operator* op = new_Operator_from_json(child, hash);
		char* id_string = get_json_node_value_string(child, "id");
		Hashtable_add(hash, id_string, op);
		if(mcmc->operator_count == 0){
			mcmc->operators = malloc(sizeof(Operator*));
		}
		else{
			mcmc->operators = realloc(mcmc->operators, sizeof(Operator*)*(mcmc->operator_count+1));
		}
		mcmc->operators[mcmc->operator_count] = op;
		mcmc->operator_count++;
	}
	
	mcmc->log_count = 0;
    if(logs != NULL){
        for (int i = 0; i < logs->child_count; i++) {
            json_node* child = logs->children[i];
            if(mcmc->log_count == 0){
                mcmc->logs = malloc(sizeof(Log*));
            }
            else{
                mcmc->logs = realloc(mcmc->logs, sizeof(Log*)*(mcmc->log_count+1));
            }
            mcmc->logs[i] = new_Log_from_json(child, hash);
            mcmc->log_count++;
        }
    }
	mcmc->chain_temperature = -1;
	mcmc->gss = get_json_node_value_bool(node, "gss", false);
	mcmc->free = _free_MCMC;
	return mcmc;
}

