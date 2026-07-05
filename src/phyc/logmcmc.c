// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "logmcmc.h"

#include <ctype.h>
#include <strings.h>

#include "treelikelihood.h"
#include "treeio.h"
#include "branchmodel.h"
#include "node.h"

// Validate that `fmt` is a single floating-point printf conversion (e.g. "%e" or
// "%.10f"). The trait/annotation values are doubles, so only float conversions
// are allowed; anything else would read the wrong argument type at runtime.
static void _validate_log_format(json_node* node, const char* fmt){
	const char* p = fmt;
	if(*p != '%'){
		json_die(node, "\"format\" must be a printf conversion starting with '%%' (e.g. \"%%.10f\"): '%s'", fmt);
	}
	p++;
	while(*p != '\0' && (isdigit((unsigned char)*p) || *p == '.' || *p == '+' ||
	                     *p == '-' || *p == ' ' || *p == '#' || *p == '0')){
		p++;
	}
	char conv = *p;
	if((conv != 'e' && conv != 'E' && conv != 'f' && conv != 'F' &&
	    conv != 'g' && conv != 'G') || *(p + 1) != '\0'){
		json_die(node, "\"format\" must be a single floating-point printf conversion (e.g. \"%%e\" or \"%%.10f\"): '%s'", fmt);
	}
}

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
			else{
				fprintf(logger->file, "\t%s", model->name);
			}
		}
		for (int i = 0; i < Parameters_count(logger->x); i++) {
			Parameter* parameter = Parameters_at(logger->x, i);
			if(Parameter_size(parameter) == 1){
				fprintf(logger->file, "\t%s", Parameters_name(logger->x, i));
			}
			else{
				for(size_t j = 0; j < Parameter_size(parameter); j++){
					StringBuffer_empty(buffer);
					StringBuffer_append_format(buffer, "%s.%zu", Parameter_name(parameter), j);
					fprintf(logger->file, "\t%s", buffer->c);
				}
			}
		}
	}
	free_StringBuffer(buffer);
	fprintf(logger->file, "\n");
}

// Refresh each node's annotation table from the branch-model traits. Rates
// change every sample, so the previous annotations are cleared first.
static void _log_tree_annotate(Log* logger, Tree* tree){
	if(logger->trait_count == 0) return;
	Node** nodes = Tree_get_nodes(tree, POSTORDER);
	int node_count = Tree_node_count(tree);
	for(int i = 0; i < node_count; i++){
		Node_empty_annotation(nodes[i]);
	}
	char buffer[64];
	for(size_t t = 0; t < logger->trait_count; t++){
		Model* model = logger->trait_models[t];
		const char* tag = logger->trait_tags[t];
		const char* fmt = logger->trait_formats[t];
		// Dispatch on the trait source. Only branch models are supported today;
		// other per-node providers (e.g. discrete traits, node parameters) add
		// their own branch here.
		if(model->type == MODEL_BRANCHMODEL){
			BranchModel* bm = model->obj;
			for(int i = 0; i < node_count; i++){
				Node* n = nodes[i];
				if(Node_isroot(n)) continue; // the root has no branch
				// snprintf, not StringBuffer_append_format: the latter only parses
				// a single-digit precision (so "%.10f" would break).
				snprintf(buffer, sizeof(buffer), fmt, bm->get(bm, n));
				Node_set_annotation(n, tag, buffer);
			}
		}
		else{
			fprintf(stderr, "logger: cannot annotate tree with trait '%s': unsupported model type %d\n", tag, model->type);
			exit(1);
		}
	}
}

void log_tree(Log* logger, size_t iter){
	Tree* tree = logger->models[0]->obj;
	_log_tree_annotate(logger, tree);

	if(strcasecmp(logger->format, "newick") == 0){
		if(logger->trait_count > 0){
			Tree_print_newick_with_annotation(logger->file, tree, false, 12);
		}
		else{
			Tree_print_newick(logger->file, tree, false, 12);
		}
	}
	else if(strcasecmp(logger->format, "nexus") == 0){
		char root_tag = 'U';
		if(Tree_rooted(tree)){
			root_tag = 'R';
		}
		fprintf(logger->file, "tree STATE_%lu ", iter);
		// whole-tree scalar annotations (likelihoods), nexus only
		for(size_t s = 0; s < logger->scalar_count; s++){
			Model* m = logger->scalar_models[s];
			fprintf(logger->file, "%s%s=", (s == 0 ? "[&" : ","), m->name);
			fprintf(logger->file, logger->scalar_formats[s], m->logP(m));
		}
		if(logger->scalar_count > 0) fprintf(logger->file, "] ");
		fprintf(logger->file, "= [&%c] ", root_tag);
		if(logger->trait_count > 0){
			Tree_print_nexus_with_annotation2(logger->file, tree, Tree_is_time_mode(tree));
		}
		else{
			Tree_print_nexus(logger->file, tree);
		}
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
		else{
			// if (logger->force) {
				fprintf(logger->file, "\t%e", model->logP(model));
			// }
			// else fprintf(logger->file, "\t%e", model->lp);
		}
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		Parameter* parameter = Parameters_at(logger->x, i);
		const double* values = Parameter_values(parameter);
		for(size_t j = 0; j < Parameter_size(parameter); j++){
			fprintf(logger->file, "\t%e", values[j]);
		}
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
		fprintf(logger->file, "\t%e", model->lp);
	}
	for (int i = 0; i < Parameters_count(logger->x); i++) {
		Parameter* parameter = Parameters_at(logger->x, i);
		const double* values = Parameter_values(parameter);
		for(size_t j = 0; j < Parameter_size(parameter); j++){
			fprintf(logger->file, "\t%e", values[j]);
		}
	}
	
	fprintf(logger->file, "\t%s\n", more);
}

static void _log_columns_header(Log* logger){
	StringBuffer* buffer = new_StringBuffer(10);
	fprintf(logger->file, "iter");
	for (size_t i = 0; i < logger->column_count; i++) {
		if (logger->columns[i].model != NULL) {
			Model* model = logger->columns[i].model;
			if(model->type == MODEL_DISCRETE_PARAMETER){
				DiscreteParameter* dp = model->obj;
				for (int j = 0; j < dp->length; j++) {
					fprintf(logger->file, "\t%s.%d", model->name, j+1);
				}
			}
			else{
				fprintf(logger->file, "\t%s", model->name);
			}
		}
		else{
			Parameter* parameter = logger->columns[i].parameter;
			if(Parameter_size(parameter) == 1){
				fprintf(logger->file, "\t%s", Parameter_name(parameter));
			}
			else{
				for(size_t j = 0; j < Parameter_size(parameter); j++){
					StringBuffer_empty(buffer);
					StringBuffer_append_format(buffer, "%s.%zu", Parameter_name(parameter), j);
					fprintf(logger->file, "\t%s", buffer->c);
				}
			}
		}
	}
	free_StringBuffer(buffer);
	fprintf(logger->file, "\n");
}

void log_columns(Log* logger, size_t iter){
	fprintf(logger->file, "%zu", iter);
	for (size_t i = 0; i < logger->column_count; i++) {
		if (logger->columns[i].model != NULL) {
			Model* model = logger->columns[i].model;
			if(model->type == MODEL_DISCRETE_PARAMETER){
				DiscreteParameter* dp = model->obj;
				for (int j = 0; j < dp->length; j++) {
					fprintf(logger->file, "\t%d", dp->values[j]);
				}
			}
			else{
				fprintf(logger->file, "\t%e", model->logP(model));
			}
		}
		else{
			Parameter* parameter = logger->columns[i].parameter;
			const double* values = Parameter_values(parameter);
			for(size_t j = 0; j < Parameter_size(parameter); j++){
				fprintf(logger->file, "\t%e", values[j]);
			}
		}
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
		if(logger->column_count > 0){
			_log_columns_header(logger);
		}
		else{
			_log_write_header(logger);
		}
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
	if(logger->column_count > 0){
		free(logger->columns);
	}
	if(logger->trait_count > 0){
		for(size_t i = 0; i < logger->trait_count; i++){
			free(logger->trait_tags[i]);
			free(logger->trait_formats[i]);
		}
		free(logger->trait_tags);
		free(logger->trait_formats);
		free(logger->trait_models);
	}
	if(logger->scalar_count > 0){
		for(size_t i = 0; i < logger->scalar_count; i++){
			free(logger->scalar_formats[i]);
		}
		free(logger->scalar_formats);
		free(logger->scalar_models);
	}
	if (logger->format != NULL) {
		free(logger->format);
	}
	free(logger);
}

Log* new_Log_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"annotate", JSON_OPTIONAL, JSON_ARRAY},
	    {"append", JSON_OPTIONAL, JSON_ANY},
	    {"columns", JSON_OPTIONAL, JSON_ANY},
	    {"cpo", JSON_OPTIONAL, JSON_ANY},
	    {"every", JSON_OPTIONAL, JSON_ANY},
	    {"file", JSON_OPTIONAL, JSON_ANY},
	    {"force", JSON_OPTIONAL, JSON_ANY},
	    {"format", JSON_OPTIONAL, JSON_ANY},
	    {"header", JSON_OPTIONAL, JSON_ANY},
	    {"models", JSON_FORBIDDEN, JSON_ANY, "replaced by 'columns'"},
	    {"traits", JSON_OPTIONAL, JSON_ARRAY},
	    {"x", JSON_FORBIDDEN, JSON_ANY, "replaced by 'columns'"},
	    {"tree", JSON_OPTIONAL, JSON_STRING},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	json_validate_xor(node, "tree", "columns", NULL);

	Log* logger = malloc(sizeof(Log));
	logger->x = new_Parameters(1);
	// A logger is either a tree logger ("tree": "@tree") or a column logger
	// ("columns": [...]); json_validate_xor above guarantees exactly one.
	logger->columns = NULL;
	logger->column_count = get_columns_from_json(node, hash, &logger->columns);

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
	
	logger->cpo = get_json_node_value_bool(node, "cpo", false);

	logger->model_count = 0;
	logger->models = NULL;
	logger->trait_models = NULL;
	logger->trait_tags = NULL;
	logger->trait_formats = NULL;
	logger->trait_count = 0;
	logger->scalar_models = NULL;
	logger->scalar_formats = NULL;
	logger->scalar_count = 0;

	// "tree": "@tree" turns this into a tree logger; the referenced model is
	// stashed in models[0] the same way the tree-printing helpers expect it.
	json_node* tree_node = get_json_node(node, "tree");
	if (tree_node != NULL) {
		char* ref = (char*)tree_node->value;
		Model* m = Hashtable_get(hash, ref + 1);
		if (m == NULL || m->type != MODEL_TREE) {
			json_die(node, "logger \"tree\" must reference a tree model (e.g. \"@tree\"): '%s'", ref);
		}
		logger->models = malloc(sizeof(Model*));
		logger->models[0] = m;
		logger->model_count = 1;
		logger->tree = true;
	}

	if (logger->cpo) {
		logger->write = log_log_cpo;
	}
	if (logger->column_count > 0) {
		logger->write = log_columns;
	}
	if (logger->tree) {
		logger->write = log_tree;
		if(format == NULL)logger->format = String_clone("newick");
	}

	// Tree annotation: per-branch "traits" (branch-model rates, both formats) and
	// whole-tree "annotate" scalars (likelihoods, nexus only). Both require a tree
	// logger. Branch models are optional, so "traits" may legitimately be absent
	// (e.g. unrooted trees have no branch model).
	json_node* traits_node = get_json_node(node, "traits");
	json_node* annotate_node = get_json_node(node, "annotate");
	if (!logger->tree && (traits_node != NULL || annotate_node != NULL)) {
		json_die(node, "\"traits\" and \"annotate\" are only valid on a tree logger");
	}
	if (traits_node != NULL) {
		logger->trait_count = traits_node->child_count;
		logger->trait_models = malloc(sizeof(Model*) * logger->trait_count);
		logger->trait_tags = malloc(sizeof(char*) * logger->trait_count);
		logger->trait_formats = malloc(sizeof(char*) * logger->trait_count);
		for (size_t i = 0; i < logger->trait_count; i++) {
			json_node* child = traits_node->children[i];
			char* tag = get_json_node_value_string(child, "tag");
			char* model_ref = get_json_node_value_string(child, "model");
			if (tag == NULL || model_ref == NULL) {
				json_die(node, "each \"traits\" entry needs a \"tag\" and a \"model\"");
			}
			Model* m = Hashtable_get(hash, model_ref + 1);
			if (m == NULL || m->type != MODEL_BRANCHMODEL) {
				json_die(node, "trait \"model\" must reference a branch model (e.g. \"@clock\"): '%s'", model_ref);
			}
			char* fmt = get_json_node_value_string(child, "format");
			if (fmt != NULL) _validate_log_format(node, fmt);
			logger->trait_models[i] = m;
			logger->trait_tags[i] = String_clone(tag);
			logger->trait_formats[i] = String_clone(fmt != NULL ? fmt : "%e");
		}
	}
	if (annotate_node != NULL) {
		if (logger->format == NULL || strcasecmp(logger->format, "nexus") != 0) {
			json_die(node, "\"annotate\" (likelihood annotations) requires \"format\": \"nexus\"");
		}
		logger->scalar_count = annotate_node->child_count;
		logger->scalar_models = malloc(sizeof(Model*) * logger->scalar_count);
		logger->scalar_formats = malloc(sizeof(char*) * logger->scalar_count);
		for (size_t i = 0; i < logger->scalar_count; i++) {
			json_node* child = annotate_node->children[i];
			char* ref;
			char* fmt = NULL;
			if (child->node_type == MJSON_OBJECT) {
				// { "model": "@posterior", "format": "%.10f" }
				ref = get_json_node_value_string(child, "model");
				if (ref == NULL) {
					json_die(node, "each object in \"annotate\" needs a \"model\"");
				}
				fmt = get_json_node_value_string(child, "format");
				if (fmt != NULL) _validate_log_format(node, fmt);
			}
			else {
				// plain string reference "@posterior"
				ref = (char*)child->value;
			}
			Model* m = Hashtable_get(hash, ref + 1);
			if (m == NULL || (m->type != MODEL_TREELIKELIHOOD && m->type != MODEL_COMPOUND && m->type != MODEL_COALESCENT)) {
				json_die(node, "\"annotate\" must reference a likelihood (treelikelihood, compound or coalescent): '%s'", ref);
			}
			logger->scalar_models[i] = m;
			logger->scalar_formats[i] = String_clone(fmt != NULL ? fmt : "%f");
		}
	}

	logger->initialize = log_initialize;
	logger->finalize = log_finalize;
	logger->free = _free_Log;
	
	return logger;
}

