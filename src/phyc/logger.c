//
//  logger.c
//  physher
//
//  Created by Mathieu Fourment on 26/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "logger.h"

void _log(struct Logger* logger){
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
}

void get_references(json_node* node, Hashtable* hash, struct Logger* logger){
	json_node* x_node = get_json_node(node, "parameters");
	logger->simplexCount = 0;
	
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

struct Logger* new_logger_from_json(json_node* node, Hashtable* hash){
	struct Logger* logger = malloc(sizeof(struct Logger));
	logger->parameters = new_Parameters(1);
	logger->simplexCount = 0;
	logger->simplexes = NULL;
	get_references(node, hash, logger);
	logger->log = _log;
	logger->file = stdout;
	logger->filename = NULL;
	logger->format = NULL;
	return logger;
}

void free_Logger(struct Logger* logger){
	free_Parameters(logger->parameters);
	if(logger->simplexes != NULL)free(logger->simplexes);
	free(logger);
}