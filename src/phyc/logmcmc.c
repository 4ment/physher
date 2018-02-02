//
//  logmcmc.c
//  physher
//
//  Created by Mathieu Fourment on 30/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "logmcmc.h"

#include "treelikelihood.h"

void _log_write_header(Log* logger, bool save_cold_ll){
	StringBuffer* buffer = new_StringBuffer(10);
	if(logger->cpo){
		Model* treelikelihood = logger->models[0];
		SingleTreeLikelihood* tlk = treelikelihood->obj;
		
		for (int i = 0; i < tlk->sp->count; i++) {
			StringBuffer_empty(buffer);
			StringBuffer_append_format(buffer, "%f", tlk->sp->weights[i]);
			fprintf(logger->file, "%c%s", (i == 0 ? '#': '\t'), buffer->c);
		}
		fprintf(logger->file, "\niter");
		for (int i = 0; i < tlk->sp->count; i++) {
			StringBuffer_set_string(buffer, "p");
			StringBuffer_append_format(buffer, "%d", i);
			fprintf(logger->file, "\t%s", buffer->c);
		}
	}
	else{
		fprintf(logger->file, "iter");
		for (int i = 0; i < logger->model_count; i++) {
			fprintf(logger->file, "\t%s", logger->models[i]->name);
		}
		for (int i = 0; i < Parameters_count(logger->x); i++) {
			fprintf(logger->file, "\t%s", Parameters_name(logger->x, i));
		}
		for (int i = 0; i < logger->simplex_count; i++) {
			Simplex* simplex = logger->simplexes[i]->obj;
			for (int j = 0; j < simplex->K; j++) {
				StringBuffer_set_string(buffer, logger->simplexes[i]->name);
				StringBuffer_append_format(buffer, "%d", (j+1));
				fprintf(logger->file, "\t%s", buffer->c);
			}
		}
		if(save_cold_ll){
			fprintf(logger->file, "\t%s", "coldll");
		}
	}
	free_StringBuffer(buffer);
	fprintf(logger->file, "\n");
}

void log_log(Log* logger, size_t iter){
	fprintf(logger->file, "%zu", iter);
	for (int i = 0; i < logger->model_count; i++) {
		fprintf(logger->file, "\t%e", logger->models[i]->logP(logger->models[i]));
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		fprintf(logger->file, "\t%e", Parameters_value(logger->x, i));
	}
	for (int i = 0; i < logger->simplex_count; i++) {
		Simplex* simplex = logger->simplexes[i]->obj;
		for (int j = 0; j < simplex->K; j++) {
			fprintf(logger->file, "\t%e", simplex->get_value(simplex, j));
		}
	}
	
	fprintf(logger->file, "\n");
}

void log_log_cpo(Log* logger, size_t iter){
	fprintf(logger->file, "%zu", iter);
	Model* treelikelihood = logger->models[0];
	SingleTreeLikelihood* tlk = treelikelihood->obj;
	tlk->calculate(tlk);// update partials
	for (int i = 0; i < tlk->sp->count; i++) {
		fprintf(logger->file, "\t%e", tlk->pattern_lk[i]);
	}
	fprintf(logger->file, "\n");
}

void log_log_with(Log* logger, size_t iter, const char* more){
	fprintf(logger->file, "%zu", iter);
	for (int i = 0; i < logger->model_count; i++) {
		fprintf(logger->file, "\t%e", logger->models[i]->logP(logger->models[i]));
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		fprintf(logger->file, "\t%e", Parameters_value(logger->x, i));
	}
	for (int i = 0; i < logger->simplex_count; i++) {
		Simplex* simplex = logger->simplexes[i]->obj;
		for (int j = 0; j < simplex->K; j++) {
			fprintf(logger->file, "\t%e", simplex->get_value(simplex, j));
		}
	}
	
	fprintf(logger->file, "\t%s\n", more);
}

void log_initialize(Log* logger){
	if (logger->filename != NULL) {
		char a[2] = "w";
		if(logger->append) a[0] = 'a';
		logger->file = fopen(logger->filename, a);
	}
}

void log_finalize(Log* logger){
	if (logger->filename != NULL) {
		fclose(logger->file);
	}
}

void _free_Log(Log* logger){
	//printf("free logger");
	free_Parameters(logger->x);
	if (logger->filename != NULL) {
		//		printf("close logger");
		free(logger->filename);
	}
	if(logger->file != NULL){
		fclose(logger->file);
	}
	if(logger->model_count > 0){
		free(logger->models);
	}
	free(logger);
}

Log* new_Log_from_json(json_node* node, Hashtable* hash){
	Log* logger = malloc(sizeof(Log));
	logger->x = new_Parameters(1);
	json_node* header_node = get_json_node(node, "header");
	json_node* x_node = get_json_node(node, "x");
	if (x_node != NULL) {
		get_parameters_references2(node, hash, logger->x, "x");
	}
	
	logger->every = get_json_node_value_size_t(node, "every", 1000);
	json_node* filename_node = get_json_node(node, "file");
	logger->write = log_log;
	logger->write_with = log_log_with;
	logger->write_header = _log_write_header;
	logger->file = stdout;
	logger->filename = NULL;
	if (filename_node != NULL) {
		char* filename = (char*)filename_node->value;
		if (strcasecmp(filename, "stderr") == 0) {
			logger->file = stderr;
		}
		else if (strcasecmp(filename, "stdout") != 0) {
			logger->filename = String_clone(filename);
			logger->append = get_json_node_value_bool(node, "append", false);
			logger->file = NULL;
		}
	}
	
	json_node* models_node = get_json_node(node, "models");
	logger->cpo = get_json_node_value_bool(node, "cpo", false);
	
	if (logger->cpo) {
		logger->write = log_log_cpo;
	}
	
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
	
	logger->initialize = log_initialize;
	logger->finalize = log_finalize;
	logger->free = _free_Log;
	
	return logger;
}
