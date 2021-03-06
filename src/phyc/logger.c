//
//  logger.c
//  physher
//
//  Created by Mathieu Fourment on 26/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "logger.h"

#include <string.h>
#include <strings.h>

#include "treeio.h"
#include "simplex.h"
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
		else if (model->type == MODEL_SIMPLEX) {
			Simplex* simplex = model->obj;
			for (int j = 0; j < simplex->K; j++) {
				fprintf(logger->file, " %f", simplex->get_value(simplex, j));
			}
			fprintf(logger->file, "\n");
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
		fprintf(logger->file, "%s: %f\n", Parameters_name(logger->parameters, i), Parameters_value(logger->parameters, i));
	}
	fprintf(logger->file, "\n");
}

void _log_tree(struct Logger* logger){
	Tree* tree = logger->models[0]->obj;
	if(strcasecmp(logger->format, "newick") == 0){
		Tree_print_newick(logger->file, tree, logger->internal);
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
		"internal",
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
	free(logger->format);
	free(logger);
}


void dumper_get_references(json_node* node, Hashtable* hash, struct Dumper* dumper){
    json_node* models_node = get_json_node(node, "models");
    dumper->model_count = 0;
    if (models_node != NULL) {
        if(models_node->node_type == MJSON_ARRAY){
            for (int i = 0; i < models_node->child_count; i++) {
                json_node* child = models_node->children[i];
                char* ref = (char*)child->value;
                // it's a ref
                if (child->node_type == MJSON_STRING && ref[0] == '&') {
                    if (dumper->model_count == 0) {
                        dumper->models = malloc(sizeof(Model*));
                    }
                    else{
                        dumper->models = realloc(dumper->models, sizeof(Model*)*(dumper->model_count+1));
                    }
                    dumper->models[dumper->model_count] = Hashtable_get(hash, ref+1);
                    dumper->model_count++;
                }
                else{
                    fprintf(stderr, "Dumper only accepts models: `%s'\n", ref);
                    exit(12);
                }
            }
        }
        // it's a ref
        else if(models_node->node_type == MJSON_STRING){
            char* ref = (char*)models_node->value;
            if (ref[0] == '&') {
                dumper->models = malloc(sizeof(Model*));
                dumper->models[dumper->model_count] = Hashtable_get(hash, ref+1);
                dumper->model_count++;
            }
        }
        else{
            fprintf(stderr, "Dmuper %s cannot read value of key %s\n", get_json_node_value_string(node, "id"), models_node->key);
            exit(2);
        }
    }
}

static void _dump(struct Dumper* dumper){
    json_node* jroot = create_json_node(NULL);
    jroot->node_type = MJSON_OBJECT;
    for (int i = 0; i < dumper->model_count; i++) {
        Model* model = dumper->models[i];
        model->jsonize(model, jroot);
    }
    json_tree_fprint(jroot, dumper->file);
    json_free_tree(jroot);
}

static void _free_Logger(struct Dumper* dumper){
    if(dumper->filename != NULL){
        free(dumper->filename);
        fclose(dumper->file);
    }
    if(dumper->model_count > 0) free(dumper->models);
    free(dumper);
}

struct Dumper* new_Dumper_from_json(json_node* node, Hashtable* hash){
    char* allowed[] = {
        "file",
        "models",
        "parameters",
    };
    json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
    
    struct Dumper* dumper = malloc(sizeof(struct Dumper));
    dumper->parameters = NULL;
    dumper->parameter_count = 0;
    dumper->model_count = 0;
    dumper->models = NULL;
    
    dumper_get_references(node, hash, dumper);
    
    json_node* file_node = get_json_node(node, "file");
    dumper->file = stdout;
    dumper->filename = NULL;
    if(file_node != NULL){
        char* filename = file_node->value;
        if (strcmp(filename, "stderr") != 0 && strcmp(filename, "stdout") != 0) {
            dumper->filename = String_clone(filename);
            dumper->file = fopen(dumper->filename, "w");
        }
        else if(strcmp(filename, "stderr") == 0){
            dumper->file = stderr;
        }
        else{
            dumper->file = stdout;
        }
    }
    dumper->dump = _dump;
    dumper->free = _free_Logger;
    
    return dumper;
}
