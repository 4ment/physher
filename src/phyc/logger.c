//
//  logger.c
//  physher
//
//  Created by Mathieu Fourment on 26/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "logger.h"

#include <string.h>

#include "treeio.h"

void _log(struct Logger* logger){
	for (int i = 0; i < logger->model_count; i++) {
		Model* model = logger->models[i];
		fprintf(logger->file, "%s: %f\n", model->name, model->logP(model));
	}
	
	for (int i = 0; i < Parameters_count(logger->parameters); i++) {
		fprintf(logger->file, "%s: %f\n", Parameters_name(logger->parameters, i), Parameters_value(logger->parameters, i));
	}
	for (int i = 0; i < logger->simplexCount; i++) {
		Model* msimplex = logger->simplexes[i];
		Simplex* simplex = msimplex->obj;
		fprintf(logger->file, "%s", msimplex->name);
		for (int j = 0; j < simplex->K; j++) {
			fprintf(logger->file, " %f", simplex->get_value(simplex, j));
		}
		fprintf(logger->file, "\n");
	}
	if (logger->tree) {
		Tree_print_newick(logger->file, logger->tree->obj, true);
	}
	fprintf(logger->file, "\n");
}

void get_references(json_node* node, Hashtable* hash, struct Logger* logger){
	json_node* x_node = get_json_node(node, "parameters");
	json_node* tree_node = get_json_node(node, "tree");
	json_node* models_node = get_json_node(node, "models");
	logger->simplexCount = 0;
	logger->model_count = 0;
	
	if (x_node != NULL) {
		if(x_node->node_type == MJSON_ARRAY){
			for (int i = 0; i < x_node->child_count; i++) {
				json_node* child = x_node->children[i];
				char* ref = (char*)child->value;
				// it's a ref
				if (child->node_type == MJSON_STRING) {
					if (ref[0] == '&') {
						Parameter* p = Hashtable_get(hash, ref+1);
						Parameters_add(logger->parameters, p);
						
					}
					// tree
					else if (ref[0] == '%') {
						Parameters* ps = Hashtable_get(hash, ref+1);
						Parameters_add_parameters(logger->parameters, ps);
					}
					// simplex
					else if (ref[0] == '$') {
						if(logger->simplexes == NULL){
							logger->simplexes = malloc(sizeof(Model*));
						}
						else{
							logger->simplexes = realloc(logger->simplexes, sizeof(Model*)*(logger->simplexCount+1));
						}
						Model* msimplex = Hashtable_get(hash, ref+1);
						logger->simplexes[logger->simplexCount] = msimplex;
						logger->simplexCount++;
					}
				}
				else{
					exit(1);
				}
			}
		}
		// it's a ref
		else if(x_node->node_type == MJSON_STRING){
			char* ref = (char*)x_node->value;
			if (ref[0] == '&') {
				Parameter* p = Hashtable_get(hash, ref+1);
				Parameters_add(logger->parameters, p);
			}
			else if (ref[0] == '%') {
				Parameters* ps = Hashtable_get(hash, ref+1);
				Parameters_add_parameters(logger->parameters, ps);
			}
			// simplex
			else if (ref[0] == '$') {
				if(logger->simplexes == NULL){
					logger->simplexes = malloc(sizeof(Model*));
				}
				else{
					logger->simplexes = realloc(logger->simplexes, sizeof(Model*)*(logger->simplexCount+1));
				}
				Model* msimplex = Hashtable_get(hash, ref+1);
				logger->simplexes[logger->simplexCount] = msimplex;
				logger->simplexCount++;
			}
		}
		else{
			exit(1);
		}
	}
	//TODO: could be an array with one elements
	else if(tree_node != NULL){
		char* ref = (char*)tree_node->value;
		logger->tree = Hashtable_get(hash, ref+1);
		
	}
	
	if (models_node != NULL) {
		if(models_node->node_type == MJSON_ARRAY){
			for (int i = 0; i < models_node->child_count; i++) {
				json_node* child = models_node->children[i];
				char* ref = (char*)child->value;
				// it's a ref
				if (child->node_type == MJSON_STRING && ref[0] == '&') {
					if (logger->model_count == 0) {
						logger->models = malloc(sizeof(Model*));
					}
					else{
						logger->models = realloc(logger->models, sizeof(Model*)*(logger->model_count+1));
					}
					logger->models[logger->model_count] = Hashtable_get(hash, ref+1);
					logger->model_count++;
				}
				else{
					exit(1);
				}
			}
		}
		// it's a ref
		else if(models_node->node_type == MJSON_STRING){
			char* ref = (char*)models_node->value;
			if (ref[0] == '&') {
				logger->models = malloc(sizeof(Model*));
				logger->models[logger->model_count] = Hashtable_get(hash, ref+1);
				logger->model_count++;
			}
		}
		else{
			exit(1);
		}
	}
}

struct Logger* new_logger_from_json(json_node* node, Hashtable* hash){
	struct Logger* logger = malloc(sizeof(struct Logger));
	logger->parameters = new_Parameters(1);
	logger->simplexCount = 0;
	logger->model_count = 0;
	logger->models = NULL;
	logger->simplexes = NULL;
	logger->tree = NULL;
	get_references(node, hash, logger);
	
	json_node* file_node = get_json_node(node, "file");
	logger->file = stdout;
	logger->filename = NULL;
	if(file_node != NULL){
		char* filename = file_node->key;
		if (strcmp(filename, "stderr") != 0 && strcmp(filename, "stdout") != 0) {
			logger->filename = String_clone(file_node->key);
			logger->file = fopen(logger->filename, "w");
		}
		else if(strcmp(filename, "stderr") == 0){
			logger->file = stderr;
		}
		else{
			logger->file = stdout;
		}
	}
	logger->log = _log;
	json_node* format_node = get_json_node(node, "format");
	logger->format = NULL;
	if (format_node != NULL) {
		logger->format = String_clone(format_node->key);
	}
	return logger;
}

void free_Logger(struct Logger* logger){
	free_Parameters(logger->parameters);
	if(logger->simplexes != NULL) free(logger->simplexes);
	if(logger->filename != NULL){
		free(logger->filename);
		fclose(logger->file);
	}
	if(logger->model_count > 0) free(logger->models);
	if (logger->format != NULL) free(logger->format);
	free(logger);
}
