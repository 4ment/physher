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


void log_log(Log* logger, size_t iter){
	fprintf(logger->file, "%zu", iter);
	for (int i = 0; i < logger->model_count; i++) {
		fprintf(logger->file, ",%e", logger->models[i]->logP(logger->models[i]));
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		fprintf(logger->file, ",%e", Parameters_value(logger->x, i));
	}
	for (int i = 0; i < logger->simplex_count; i++) {
		Simplex* simplex = logger->simplexes[i]->obj;
		for (int j = 0; j < simplex->K; j++) {
			fprintf(logger->file, ",%e", simplex->get_value(simplex, j));
		}
	}
	
	fprintf(logger->file, "\n");
}

Log* new_Log_from_json(json_node* node, Hashtable* hash){
	Log* logger = malloc(sizeof(Log));
	logger->x = new_Parameters(1);
	json_node* x_node = get_json_node(node, "x");
	if (x_node != NULL) {
		get_parameters_references2(node, hash, logger->x, "x");
	}

	logger->every = get_json_node_value_size_t(node, "every", 1000);
	json_node* filename_node = get_json_node(node, "file");
	logger->write = log_log;
	logger->file = stdout;
	logger->filename = NULL;
	if (filename_node != NULL) {
		char* filename = (char*)filename_node->value;
		if (strcasecmp(filename, "stderr") == 0) {
			logger->file = stderr;
		}
		else if (strcasecmp(filename, "stdout") != 0) {
			logger->filename = String_clone(filename);
			char a[2] = "w";
			bool append = get_json_node_value_bool(node, "append", false);
			if(append) a[0] = 'a';
			logger->file = fopen(logger->filename, a);
		}
	}
	
	json_node* models_node = get_json_node(node, "models");
	logger->model_count = 0;
	logger->models = NULL;
	if (models_node != NULL) {
		if (models_node->node_type == MJSON_ARRAY) {
			for (int i = 0; i < models_node->child_count; i++) {
				json_node* child = models_node->children[i];
				char* child_string = child->value;
				Model* m = Hashtable_get(hash, child_string+1);
				if(logger->model_count == 0){
					logger->models = malloc(sizeof(Model*));
				}
				else{
					logger->models = realloc(logger->models, sizeof(Model*)*(logger->model_count+1));
				}
				logger->models[i] = m;
				logger->model_count++;
			}
		}
		else if (models_node->node_type == MJSON_STRING) {
			char* ref = models_node->value;
			logger->models = malloc(sizeof(Model*));
			logger->models[0] = Hashtable_get(hash, ref+1);;
			logger->model_count++;
		}
	}
	
	json_node* simplexes_node = get_json_node(node, "simplexes");
	logger->simplex_count = 0;
	logger->simplexes = NULL;
	if (simplexes_node != NULL) {
		if (simplexes_node->node_type == MJSON_ARRAY) {
			for (int i = 0; i < simplexes_node->child_count; i++) {
				json_node* child = simplexes_node->children[i];
				char* child_string = child->value;
				Model* m = Hashtable_get(hash, child_string+1);
				if(logger->simplex_count == 0){
					logger->simplexes = malloc(sizeof(Model*));
				}
				else{
					logger->simplexes = realloc(logger->simplexes, sizeof(Model*)*(logger->simplex_count+1));
				}
				logger->simplexes[i] = m;
				logger->simplex_count++;
			}
		}
		else if (simplexes_node->node_type == MJSON_STRING) {
			char* ref = simplexes_node->value;
			logger->simplexes = malloc(sizeof(Model*));
			logger->simplexes[0] = Hashtable_get(hash, ref+1);
			logger->simplex_count++;
		}
	}
	
	fprintf(logger->file, "iter");
	for (int i = 0; i < logger->model_count; i++) {
		fprintf(logger->file, ",%s", logger->models[i]->name);
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		fprintf(logger->file, ",%s", Parameters_name(logger->x, i));
	}
	StringBuffer* buffer = new_StringBuffer(10);
	for (int i = 0; i < logger->simplex_count; i++) {
		Simplex* simplex = logger->simplexes[i]->obj;
		for (int j = 0; j < simplex->K; j++) {
			StringBuffer_set_string(buffer, logger->simplexes[i]->name);
			StringBuffer_append_format(buffer, "%d", (j+1));
			fprintf(logger->file, ",%s", buffer->c);
		}
	}
	free_StringBuffer(buffer);
	fprintf(logger->file, "\n");
	return logger;
}

void free_Log(Log* logger){
	free_Parameters(logger->x);
	if (logger->filename != NULL) {
		free(logger->filename);
		fclose(logger->file);
	}
	if(logger->model_count > 0){
		free(logger->models);
	}
	free(logger);
}

void run(MCMC* mcmc){
	Model* model = mcmc->model;
	
//	double logTau = logTopologyPrior(Tree_tip_count(tlk->tree));

	double logP = model->logP(model);
	
	size_t iter = 0;
	
	int operator = 1;
	
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
	
	while ( iter < mcmc->chain_length) {
		
		// Choose operator
		int op_index = roulette_wheel(weights, mcmc->operator_count);
		Operator* op = mcmc->operators[op_index];
		
		double logHR = 0;
		
		op->store(op);
		bool success = op->propose(op, &logHR);
		
		// operator did not propose a valid value
		if (success == false) {
			iter++;
			continue;
		}
		double proposed_logP = model->logP(model);
		
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
			op->restore(op);
			op->rejected_count++;
		}
		
		iter++;
		
		if(op->optimize != NULL){
			op->optimize(op, alpha);
		}
		
		for (int i = 0; i < mcmc->log_count; i++) {
			if(iter % mcmc->logs[i]->every == 0){
				mcmc->logs[i]->write(mcmc->logs[i], iter);
			}
		}
	}
	if( mcmc->verbose > 0){
		for (int i = 0; i < mcmc->operator_count; i++) {
			Operator* op = mcmc->operators[i];
			printf("Acceptance ratio %s: %f %f\n", op->name, ((double)op->accepted_count/(op->accepted_count+op->rejected_count)), op->parameters[0]);
		}
	}
	free(weights);
	free_StringBuffer(buffer);
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
	return mcmc;
}

void free_MCMC(MCMC* mcmc){
	free_Model(mcmc->model);
	for (int i = 0; i < mcmc->operator_count; i++) {
		free_Operator(mcmc->operators[i]);
	}
	for (int i = 0; i < mcmc->log_count; i++) {
		free_Log(mcmc->logs[i]);
	}
	free(mcmc->logs);
	free(mcmc->operators);
	free(mcmc);
}
