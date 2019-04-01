//
//  logmcmc.c
//  physher
//
//  Created by Mathieu Fourment on 30/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "logmcmc.h"

#include <strings.h>

#include "treelikelihood.h"
#include "treeio.h"

static void _log_write_header(Log* logger){
	StringBuffer* buffer = new_StringBuffer(10);
	if(logger->cpo){
		fprintf(logger->file, "#");
		for(int j = 0; j < logger->model_count; j++){
			Model* treelikelihood = logger->models[j];
			SingleTreeLikelihood* tlk = treelikelihood->obj;
			
			for (int i = 0; i < tlk->sp->count; i++) {
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "%f", tlk->sp->weights[i]);
				fprintf(logger->file, "%s%s", (i == 0 && j == 0 ? "": "\t"), buffer->c);
			}
		}
		fprintf(logger->file, "\niter");
		for(int j = 0; j < logger->model_count; j++){
			Model* treelikelihood = logger->models[j];
			SingleTreeLikelihood* tlk = treelikelihood->obj;
			for (int i = 0; i < tlk->sp->count; i++) {
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "%s%s%d", treelikelihood->name, ".p", i);
				fprintf(logger->file, "\t%s", buffer->c);
			}
		}
		fflush(logger->file);
	}
	else{
		fprintf(logger->file, "iter");
		for (int i = 0; i < logger->model_count; i++) {
			Model* model = logger->models[i];
			if(model->type == MODEL_DISCRETE_PARAMETER){
				DiscreteParameter* dp = model->obj;
				for (int j = 0; j < dp->length; j++) {
					fprintf(logger->file, "\t%s.%d", model->name, j+1);
				}
			}
			else if(model->type == MODEL_SIMPLEX){
				Simplex* simplex = model->obj;
				for (int j = 0; j < simplex->K; j++) {
					StringBuffer_set_string(buffer, model->name);
					StringBuffer_append_format(buffer, ".%d", (j+1));
					fprintf(logger->file, "\t%s", buffer->c);
				}
			}
			else{
				fprintf(logger->file, "\t%s", model->name);
			}
		}
		for (int i = 0; i < Parameters_count(logger->x); i++) {
			fprintf(logger->file, "\t%s", Parameters_name(logger->x, i));
		}
	}
	free_StringBuffer(buffer);
	fprintf(logger->file, "\n");
}

void log_tree(Log* logger, size_t iter){
	Tree* tree = logger->models[0]->obj;
	if(strcasecmp(logger->format, "newick") == 0){
		Tree_print_newick(logger->file, tree, false);
	}
	else if(strcasecmp(logger->format, "nexus") == 0){
		char root_tag = 'U';
		if(Tree_rooted(tree)){
			root_tag = 'R';
		}
		fprintf(logger->file, "tree STATE_%lu = [&%c] ", iter, root_tag);
		Tree_print_nexus(logger->file, tree);
	}
	fprintf(logger->file, "\n");
}

void log_log(Log* logger, size_t iter){
	fprintf(logger->file, "%zu", iter);
	for (int i = 0; i < logger->model_count; i++) {
		Model* model = logger->models[i];
		if(model->type == MODEL_DISCRETE_PARAMETER){
			DiscreteParameter* dp = model->obj;
			for (int j = 0; j < dp->length; j++) {
				fprintf(logger->file, "\t%d", dp->values[j]);
			}
		}
		else if(model->type == MODEL_SIMPLEX){
			Simplex* simplex = model->obj;
			for (int j = 0; j < simplex->K; j++) {
				fprintf(logger->file, "\t%e", simplex->get_value(simplex, j));
			}
		}
		else{
			if (logger->force) {
				fprintf(logger->file, "\t%e", model->logP(model));
			}
			else fprintf(logger->file, "\t%e", model->lp);
		}
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		fprintf(logger->file, "\t%e", Parameters_value(logger->x, i));
	}
	
	if (logger->filename == NULL) {
		if(iter > 0){
			gettimeofday(&logger->end, NULL);
			double diff_time = (double)(logger->end.tv_usec - logger->start.tv_usec) / 1000000 + (double)(logger->end.tv_sec - logger->start.tv_sec);
			double speed = diff_time/logger->every*1e6;
			if (speed < 1) {
				fprintf(logger->file, "  %.2f sec/million", speed);
			}
			else{
				fprintf(logger->file, "  %.2f min/million", speed/60);
			}
			
			logger->start = logger->end;
		}
	}

	fprintf(logger->file, "\n");
	fflush(logger->file);
}

void log_log_cpo(Log* logger, size_t iter){
	fprintf(logger->file, "%zu", iter);
	for(int j = 0; j < logger->model_count; j++){
		Model* treelikelihood = logger->models[j];
		SingleTreeLikelihood* tlk = treelikelihood->obj;
		tlk->calculate(tlk);// update partials
		for (int i = 0; i < tlk->sp->count; i++) {
			fprintf(logger->file, "\t%e", tlk->pattern_lk[i]);
		}
	}
	fprintf(logger->file, "\n");
}

void log_log_with(Log* logger, size_t iter, const char* more){
	fprintf(logger->file, "%zu", iter);
	for (int i = 0; i < logger->model_count; i++) {
		Model* model = logger->models[i];
		if(model->type == MODEL_SIMPLEX){
			Simplex* simplex = model->obj;
			for (int j = 0; j < simplex->K; j++) {
				fprintf(logger->file, "\t%e", simplex->get_value(simplex, j));
			}
		}
		else{
			fprintf(logger->file, "\t%e", model->lp);
		}
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		fprintf(logger->file, "\t%e", Parameters_value(logger->x, i));
	}
	
	fprintf(logger->file, "\t%s\n", more);
}

void log_initialize(Log* logger){
	if (logger->filename != NULL) {
		char a[2] = "w";
		if(logger->append) a[0] = 'a';
		logger->file = fopen(logger->filename, a);
	}
	if(logger->tree && strcasecmp(logger->format, "nexus") == 0){
		Tree* tree = logger->models[0]->obj;
		fprintf(logger->file, "#NEXUS\n\n");
		Tree_print_nexus_taxa_block(logger->file, tree);
		Tree_print_nexus_header_figtree_BeginTrees(logger->file, tree);
	}
	else if(!logger->tree || (logger->tree && strcasecmp(logger->format, "newick") != 0)){
		_log_write_header(logger);
	}
}

void log_finalize(Log* logger){
	if(logger->tree && strcasecmp(logger->format, "nexus") == 0){
		fprintf(logger->file, "end;\n");
	}
	if (logger->filename != NULL) {
		fclose(logger->file);
		logger->file = NULL;
	}
}

void _free_Log(Log* logger){
	//printf("free logger");
	free_Parameters(logger->x);
	if (logger->filename != NULL) {
		//		printf("close logger");
		free(logger->filename);
	}
	
	if(logger->model_count > 0){
		free(logger->models);
	}
	if (logger->format != NULL) {
		free(logger->format);
	}
	free(logger);
}

Log* new_Log_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"append",
		"cpo",
		"every",
		"file",
		"force",
		"format",
		"header",
		"models",
		"x"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
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
	logger->file = stdout;
	logger->filename = NULL;
	gettimeofday(&logger->start, NULL);
	logger->tree = false;
	logger->format = NULL;
	char* format = get_json_node_value_string(node, "format");
	logger->force = get_json_node_value_bool(node, "force", false);
	if(format != NULL){
		logger->format = String_clone(format);
	}
	
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
				if (m->type == MODEL_TREE) {
					logger->tree = true;
				}
			}
		}
		else if (models_node->node_type == MJSON_STRING) {
			char* ref = models_node->value;
			logger->models = malloc(sizeof(Model*));
			logger->models[0] = Hashtable_get(hash, ref+1);;
			logger->model_count++;
			if (logger->models[0]->type == MODEL_TREE) {
				logger->tree = true;
			}
		}
	}
	
	if (logger->tree) {
		logger->write = log_tree;
		if(format == NULL)logger->format = String_clone("newick");
	}
	
	logger->initialize = log_initialize;
	logger->finalize = log_finalize;
	logger->free = _free_Log;
	
	return logger;
}

