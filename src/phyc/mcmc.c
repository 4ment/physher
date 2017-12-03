//
//  mcmc.c
//  physher
//
//  Created by Mathieu Fourment on 2/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "mcmc.h"

#include <strings.h>

#include "random.h"
#include "matrix.h"
#include "compoundmodel.h"
#include "treelikelihood.h"
#include "gamma.h"
#include "dirichlet.h"
#include "simplex.h"

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
			logger->simplexes[0] = Hashtable_get(hash, ref+1);;
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

double tune(int accepted, int count, double target, double tuner, double min, double max, bool inverse){
	double delta = 1./sqrt(count);
	//delta = fmin(0.01, delta);
	double logTuner = log(tuner);
	double accRatio = (double)accepted/count;
	if (inverse) {
		delta = -delta;
	}
	if(accRatio > target){
		logTuner += delta;
	}
	else {
		logTuner -= delta;
	}
	double t = exp(logTuner);
	
	if (t <= min || t >= max)
		t = tuner;
	return t;
}

double getDeltaP(int count, double logAlpha, double target){
	return ((1.0 / count) * (exp(fmin(logAlpha, 0.)) - target));
}

double optimizeScaleFactor(double scaleFactor, int count, double logAlpha, double target ){
	double delta = getDeltaP(count, logAlpha, target);
	delta += log(1./scaleFactor-1.0);
	return 1./(exp(delta)+1);
}

void operator_store(Operator* op){
	if(op->x != NULL){
		Parameters_store_value(op->x, op->storage);
	}
	else if(op->simplex != NULL){
		memcpy(op->storage, ((Simplex*)op->simplex->obj)->values, sizeof(double)*((Simplex*)op->simplex->obj)->K);
	}
}

void operator_restore(Operator* op){
	if(op->x != NULL){
		if(op->index >= 0){
			Parameters_set_value(op->x, op->index, op->storage[op->index]);
		}
		else{
			Parameters_restore_value(op->x, op->storage);
		}
	}
	else if(op->simplex != NULL){
		Simplex* simplex = op->simplex->obj;
		simplex->set_values(simplex, op->storage);
	}
}

bool operator_scaler(Operator* op, double* logHR){
	op->index = random_int(Parameters_count(op->x)-1);
	Parameter* p = Parameters_at(op->x, op->index);
	double scaler_scaleFactor = op->parameters[0];
	double v = Parameter_value(p);
	double vv = v;
	double s = (scaler_scaleFactor + (random_double() * ((1.0 / scaler_scaleFactor) - scaler_scaleFactor)));
	vv *= s;
	*logHR = -log(s);
	if ( vv > Parameter_upper(p) || vv < Parameter_lower(p ) ) {
		op->rejected_count++;
		return false;
	}
	Parameter_set_value(p, vv);
	return true;
}

bool operator_slider(Operator* op, double* logHR){
	op->index = random_int(Parameters_count(op->x)-1);
	Parameter* p = Parameters_at(op->x, op->index);
	double v = Parameter_value(p);
	double vv = v;
	double w = (random_double() - 0.5)*op->parameters[0];//slider_delta;
	vv += w;
//	printf("%s %f %f ", Parameter_name(p), w, Parameter_value(p));
	if ( vv > Parameter_upper(p) ) {
		vv = 2.0*Parameter_upper(p) - vv;
	}
	else if ( vv < Parameter_lower(p ) ) {
		vv = 2.0*Parameter_lower(p) - vv;
	}
	*logHR = 0;
	if ( vv > Parameter_upper(p) || vv < Parameter_lower(p) ) {
		op->rejected_count++;
		return false;
	}
	Parameter_set_value(p, vv);
//	printf("%f\n", Parameter_value(p));
	return true;
}

bool operator_dirichlet(Operator* op, double* logHR){
	Simplex* simplex = op->simplex->obj;
	double alpha = op->parameters[0];
	double sum = 0;
	double* scaledOld = dvector(simplex->K);
	double* newScaled = dvector(simplex->K);
	double* newValues = dvector(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		scaledOld[i] = alpha*simplex->values[i];
		newValues[i] = rgamma(scaledOld[i]);
		sum += newValues[i];
	}
	for (int i = 0; i < simplex->K; i++) {
		newValues[i] /= sum;
		newScaled[i] = alpha*newValues[i];
	}
	
	for (int i = 0; i < simplex->K; i++) {
		if (isnan(newValues[i])) {
			free(scaledOld);
			free(newValues);
			free(newScaled);
			return false;
		}
	}
	simplex->set_values(simplex, newValues);
	
	double f = ddirchletln(newValues, simplex->K, scaledOld);
	double b = ddirchletln(simplex->values, simplex->K, newScaled);
	*logHR = b-f;
	free(scaledOld);
	free(newValues);
	free(newScaled);
	return true;
}


void operator_scaler_optimize(Operator* op, double alpha){
	op->parameters[0] = optimizeScaleFactor(op->parameters[0], op->accepted_count+op->rejected_count+1, alpha, 0.24);
}

void operator_slider_optimize(Operator* op, double alpha){
	op->parameters[0] = tune(op->accepted_count, op->accepted_count+op->rejected_count, 0.24, op->parameters[0], 0.0001, 20, false);
}

void operator_dirichlet_optimize(Operator* op, double alpha){
	op->parameters[0] = tune(op->accepted_count, op->accepted_count+op->rejected_count, 0.24, op->parameters[0], 0.1, 10000, true);
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
//			printf("%d %f %f\n", iter, logP, proposed_logP);
			logP = proposed_logP;
			op->accepted_count++;
		}
		// reject
		else {
//			printf("%d %f %f *\n", iter, logP, proposed_logP);
			op->restore(op);
			op->rejected_count++;
		}
		
		iter++;
		
		op->optimize(op, alpha);
		
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

Operator* new_Operator_from_json(json_node* node, Hashtable* hash){
	Operator* op = malloc(sizeof(Operator));
	const char* id_string = get_json_node_value_string(node, "id");
	const char* algorithm_string = get_json_node_value_string(node, "algorithm");
	const char* x_string = get_json_node_value_string(node, "x");
	op->x = NULL;
	op->weight = get_json_node_value_double(node, "weight", 1);
	op->name = String_clone(id_string);
	op->store = operator_store;
	op->restore = operator_restore;
	op->parameters = NULL;
	op->simplex = NULL;

	if (strcasecmp(algorithm_string, "scaler") == 0) {
		op->x = new_Parameters(1);
		get_parameters_references2(node, hash, op->x, "x");
		op->storage = dvector(Parameters_count(op->x));
		op->propose = operator_scaler;
		op->optimize = operator_scaler_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 0.9;
		
	}
	else if (strcasecmp(algorithm_string, "slider") == 0) {
		op->x = new_Parameters(1);
		get_parameters_references2(node, hash, op->x, "x");
		op->storage = dvector(Parameters_count(op->x));
		op->propose = operator_slider;
		op->optimize = operator_slider_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 0.001;
	}
	else if (strcasecmp(algorithm_string, "dirichlet") == 0) {
		char* ref = get_json_node_value_string(node, "x");
		op->simplex = Hashtable_get(hash, ref+1);
		op->simplex->ref_count++;
		op->storage = dvector(((Simplex*)op->simplex->obj)->K);
		op->propose = operator_dirichlet;
		op->optimize = operator_dirichlet_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 1000;
	}
	op->rejected_count = 0;
	op->accepted_count = 0;
	return op;
}

void free_Operator(Operator* op){
	if(op->parameters != NULL) free(op->parameters);
	if(op->simplex != NULL) op->simplex->free(op->simplex);
	free(op->storage);
	free_Parameters(op->x);
	free(op->name);
	free(op);
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
