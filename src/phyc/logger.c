// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "logger.h"

#include <string.h>
#include <strings.h>

#include "treeio.h"
#include "discreteparameter.h"


void _log(struct Logger* logger){
	for (int i = 0; i < logger->model_count; i++) {
		Model* model = logger->models[i];
		//only print name of model with stderr or stdout
		if(logger->filename == NULL) {
			fprintf(logger->file, "%s:", model->name);
		}
		if(model->type == MODEL_DISCRETE_PARAMETER){
			DiscreteParameter* dp = model->obj;
			for (int j = 0; j < dp->length; j++) {
				fprintf(logger->file, " %d", dp->values[j]);
			}
		}
		else if (model->print != NULL) {
			//only print name of model with stderr or stdout
			if(logger->filename == NULL) {
				fprintf(logger->file, "\n");
			}
			model->print(model, logger->file);
		}
		else{
			fprintf(logger->file, " %.10f\n", model->logP(model));
		}
	}
	
	for (int i = 0; i < Parameters_count(logger->parameters); i++) {
		Parameter* parameter = Parameters_at(logger->parameters, i);
		fprintf(logger->file, "%s:", Parameter_name(parameter));
		for (size_t j = 0; j < Parameter_size(parameter); j++) {
			fprintf(logger->file, " %f", Parameter_value_at(parameter, j));
		}
		fprintf(logger->file, "\n");
	}
	fprintf(logger->file, "\n");
}

void _log_columns(struct Logger* logger){
	for (size_t i = 0; i < logger->column_count; i++) {
		if (logger->columns[i].model != NULL) {
			Model* model = logger->columns[i].model;
			//only print name of model with stderr or stdout
			if(logger->filename == NULL) {
				fprintf(logger->file, "%s:", model->name);
			}
			if(model->type == MODEL_DISCRETE_PARAMETER){
				DiscreteParameter* dp = model->obj;
				for (int j = 0; j < dp->length; j++) {
					fprintf(logger->file, " %d", dp->values[j]);
				}
			}
			else if (model->print != NULL) {
				//only print name of model with stderr or stdout
				if(logger->filename == NULL) {
					fprintf(logger->file, "\n");
				}
				model->print(model, logger->file);
			}
			else{
				fprintf(logger->file, " %.10f\n", model->logP(model));
			}
		}
		else{
			Parameter* parameter = logger->columns[i].parameter;
			fprintf(logger->file, "%s:", Parameter_name(parameter));
			for (size_t j = 0; j < Parameter_size(parameter); j++) {
				fprintf(logger->file, " %f", Parameter_value_at(parameter, j));
			}
			fprintf(logger->file, "\n");
		}
	}
	fprintf(logger->file, "\n");
}

void _log_tree(struct Logger* logger){
	Tree* tree = logger->models[0]->obj;
	if(strcasecmp(logger->format, "newick") == 0){
		Tree_print_newick(logger->file, tree, logger->internal, 12);
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

size_t get_columns_from_json(json_node* node, Hashtable* hash, LogColumn** columns){
	json_node* columns_node = get_json_node(node, "columns");
	*columns = NULL;
	if (columns_node == NULL) return 0;

	// normalize to an array of ref strings
	json_node** items;
	size_t item_count;
	if (columns_node->node_type == MJSON_ARRAY) {
		items = columns_node->children;
		item_count = columns_node->child_count;
	}
	else if (columns_node->node_type == MJSON_STRING) {
		items = &columns_node;
		item_count = 1;
	}
	else{
		fprintf(stderr, "\"columns\" must be a string or an array of strings\n");
		exit(1);
	}

	LogColumn* cols = NULL;
	size_t count = 0;
	for (size_t i = 0; i < item_count; i++) {
		char* ref = (char*)items[i]->value;
		if (ref[0] == '@') {
			cols = realloc(cols, sizeof(LogColumn) * (count + 1));
			cols[count].model = Hashtable_get(hash, ref + 1);
			cols[count].parameter = NULL;
			count++;
		}
		else if (ref[0] == '&' || ref[0] == '%') {
			// resolve the (possibly multi-element) reference, then splay it into
			// one column per Parameter so the header/value order stays aligned
			Parameters* tmp = new_Parameters(1);
			get_parameter_reference(ref, hash, tmp);
			for (size_t j = 0; j < Parameters_count(tmp); j++) {
				cols = realloc(cols, sizeof(LogColumn) * (count + 1));
				cols[count].model = NULL;
				cols[count].parameter = Parameters_at(tmp, j);
				count++;
			}
			free_Parameters_weak(tmp);
		}
		else{
			fprintf(stderr, "column reference '%s' must start with '@' (Model), '&' (Parameter) or '%%' (Parameters)\n", ref);
			exit(1);
		}
	}
	*columns = cols;
	return count;
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
			fprintf(stderr, "Logger %s cannot read value of key %s\n", get_json_node_value_string(node, "id"), models_node->key);
			exit(2);
		}
	}
}

struct Logger* new_logger_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"columns", JSON_OPTIONAL, JSON_ANY},
	    {"file", JSON_OPTIONAL, JSON_STRING},
	    {"format", JSON_OPTIONAL, JSON_STRING},
	    {"internal", JSON_OPTIONAL, JSON_BOOL},
	    {"models", JSON_OPTIONAL, JSON_ANY},
	    {"parameters", JSON_OPTIONAL, JSON_ANY},
	    {"tree", JSON_OPTIONAL, JSON_BOOL},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

	struct Logger* logger = malloc(sizeof(struct Logger));
	logger->parameters = new_Parameters(1);
	logger->model_count = 0;
	logger->models = NULL;
	logger->tree = get_json_node_value_bool(node, "tree", false);
	// "columns" is the ordered superset of "parameters"+"models"; when present it
	// wins and the legacy keys are ignored.
	logger->columns = NULL;
	logger->column_count = get_columns_from_json(node, hash, &logger->columns);
	if (logger->column_count == 0) {
		get_references(node, hash, logger);
	}
	
	logger->internal = get_json_node_value_bool(node, "internal", false);
	
	json_node* file_node = get_json_node(node, "file");
	logger->file = stdout;
	logger->filename = NULL;
	if(file_node != NULL){
		char* filename = file_node->value;
		if (strcmp(filename, "stderr") != 0 && strcmp(filename, "stdout") != 0) {
			logger->filename = String_clone(filename);
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
	if (logger->column_count > 0) {
		logger->log = _log_columns;
	}
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
    char* sep = get_json_node_value_string(node, "sep");
    logger->sep = ',';
    if (sep != NULL) {
        logger->sep = sep[0];
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
	if(logger->column_count > 0) free(logger->columns);
	free(logger->format);
	free(logger);
}
