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
#include "simplex.h"
#include "discreteparameter.h"


void _log(struct Logger* logger){
	for (int i = 0; i < logger->model_count; i++) {
		Model* model = logger->models[i];
		fprintf(logger->file, "%s:", model->name);
		if(model->type == MODEL_DISCRETE_PARAMETER){
			DiscreteParameter* dp = model->obj;
			for (int j = 0; j < dp->length; j++) {
				fprintf(logger->file, " %d", dp->values[j]);
			}
		}
		else if (model->type == MODEL_SIMPLEX) {
			Simplex* simplex = model->obj;
			for (int j = 0; j < simplex->K; j++) {
				fprintf(logger->file, " %f", simplex->get_value(simplex, j));
			}
			fprintf(logger->file, "\n");
		}
		else{
			fprintf(logger->file, " %f\n", model->logP(model));
		}
	}
	
	for (int i = 0; i < Parameters_count(logger->parameters); i++) {
		fprintf(logger->file, "%s: %f\n", Parameters_name(logger->parameters, i), Parameters_value(logger->parameters, i));
	}
	fprintf(logger->file, "\n");
}

void _log_tree(struct Logger* logger){
	Tree* tree = logger->models[0]->obj;
	if(strcasecmp(logger->format, "newick") == 0){
		Tree_print_newick(logger->file, tree, false);
	}
	else if(strcasecmp(logger->format, "nexus") == 0){
		fprintf(logger->file, "#NEXUS\n\n");
		Tree_print_nexus_taxa_block(logger->file, tree);
		Tree_print_nexus_header_figtree_BeginTrees(logger->file, tree);

		char root_tag = 'U';
		if(Tree_rooted(tree)){
			root_tag = 'R';
		}
		fprintf(logger->file, "tree = [&%c] ", root_tag);
		Tree_print_nexus(logger->file, tree);
		fprintf(logger->file, "\nend;\n");
	}
	fprintf(logger->file, "\n");
}

void get_references(json_node* node, Hashtable* hash, struct Logger* logger){
	json_node* x_node = get_json_node(node, "parameters");
	json_node* models_node = get_json_node(node, "models");
	logger->model_count = 0;
	
	if (x_node != NULL) {
		get_parameters_references(node, hash, logger->parameters);
	}
	else if (models_node != NULL) {
		if(models_node->node_type == MJSON_ARRAY){
			for (int i = 0; i < models_node->child_count; i++) {
				json_node* child = models_node->children[i];
				char* ref = (char*)child->value;
				// it's a ref
				if (child->node_type == MJSON_STRING && (ref[0] == '&' || ref[0] == '$')) {
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
			fprintf(stderr, "Logger %s cannot read value of key %s\n", get_json_node_value_string(node, "id"), models_node->key);
			exit(2);
		}
	}
}

struct Logger* new_logger_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"file",
		"format",
		"models",
		"parameters",
		"tree"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	struct Logger* logger = malloc(sizeof(struct Logger));
	logger->parameters = new_Parameters(1);
	logger->model_count = 0;
	logger->models = NULL;
	logger->tree = get_json_node_value_bool(node, "tree", false);
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
	if (logger->tree) {
		logger->log = _log_tree;
	}
	
	char* format = get_json_node_value_string(node, "format");
	if(format != NULL){
		logger->format = String_clone(format);
	}
	else{
		logger->format = String_clone("newick");
	}
	return logger;
}

void free_Logger(struct Logger* logger){
	free_Parameters(logger->parameters);
	if(logger->filename != NULL){
		free(logger->filename);
		fclose(logger->file);
	}
	if(logger->model_count > 0) free(logger->models);
	free(logger->format);
	free(logger);
}
