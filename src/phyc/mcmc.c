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
#include <signal.h>

#include "utilsgsl.h"
#include "matrix.h"
#include "compoundmodel.h"
#include "treelikelihood.h"

bool mcmc_interrupted = false;

void mcmc_signal_callback_handler( int signum ) {
	printf("Caught signal %d\n",signum);
	if ( mcmc_interrupted == true ) {
		exit(SIGINT);
	}
	mcmc_interrupted = true;
}

static double  _calculate_logP(MCMC* mcmc){
	Model* model = mcmc->model;
	double logP = 0;
	if(mcmc->chain_temperature < 0){
		logP = model->logP(model);
	}
	// bayes factor or generalized marginal with geometric tempering
	else if(mcmc->generalized || mcmc->bf){
		model->logP(model);
		CompoundModel* cm = (CompoundModel*)model->obj;
		logP = cm->models[0]->logP(cm->models[0]) * mcmc->chain_temperature; // Un-normalized posterior in marginal GSS
		logP += cm->models[1]->logP(cm->models[1]) * (1.0 - mcmc->chain_temperature); // Reference distribution in marginal GSS
		
		// allows some terms not to be heated
		for (int i = 2; i < cm->count; i++) {
			logP += cm->models[i]->logP(cm->models[i]);
		}
	}
	// standard marginal with geometric tempering
	else{
		model->logP(model);
		CompoundModel* cm = (CompoundModel*)model->obj;
		logP = cm->models[0]->logP(cm->models[0]) * mcmc->chain_temperature;
		for (int i = 1; i < cm->count; i++) {
			logP += cm->models[i]->logP(cm->models[i]);
		}
	}
	return logP;
}

void run(MCMC* mcmc){
	if (mcmc->interruptible) {
		signal(SIGINT, mcmc_signal_callback_handler);
		mcmc_interrupted = false;
	}
	
	Model* model = mcmc->model;
	
	size_t iter = 0;
	
	double* weights = dvector(mcmc->operator_count);
	double sum = 0;
	for (int i = 0; i < mcmc->operator_count; i++) {
		sum += mcmc->operators[i]->weight;
	}
	for (int i = 0; i < mcmc->operator_count; i++) {
		weights[i] = mcmc->operators[i]->weight/sum;
	}
	
	StringBuffer *buffer = new_StringBuffer(20);
	double logP;
	
	// Initialize loggers
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->initialize(mcmc->logs[i]);
	}
	
	for (int i = 0; i < mcmc->operator_count; i++) {
		mcmc->operators[i]->rejected_count = 0;
		mcmc->operators[i]->accepted_count = 0;
		mcmc->operators[i]->failure_count = 0;
	}
	
	logP = _calculate_logP(mcmc);
	
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->write(mcmc->logs[i], iter);
	}

#ifdef DEBUG_OPERATORS
	FILE* ofile = fopen("operators.csv", "w");
	fprintf(ofile, "iter,operator,accepted,rejected,prob,parameter\n");
	for (int i = 0; i < mcmc->operator_count; i++) {
		Operator* op = mcmc->operators[i];
		double prob = (double)op->accepted_count/(op->accepted_count+op->rejected_count);
		fprintf(ofile, "%d,%s,%zu,%zu,%f,%f\n", 0, op->name, op->accepted_count, op->rejected_count, prob,(op->parameters != NULL ? op->parameters[0]: -1));
	}
	fflush(ofile);
#endif

	while ( iter < mcmc->chain_length) {
		
		// Choose operator
		size_t op_index = roulette_wheel_gsl(mcmc->rng, weights, mcmc->operator_count);
		Operator* op = mcmc->operators[op_index];
		
		double logHR = 0;
		
		model->store(model);
		bool success = op->propose(op, &logHR);
		
		if (success) {
			double proposed_logP = _calculate_logP(mcmc);
			
			double alpha = proposed_logP - logP + logHR;
			
			// accept
			if ( alpha >=  0 || alpha > log(gsl_rng_uniform(mcmc->rng)) ) {
	//			printf("%zu %f %f\n", iter, logP, proposed_logP);
				logP = proposed_logP;
				op->accepted_count++;
			}
			// reject
			else {
	//			printf("%zu %f %f *\n", iter, logP, proposed_logP);
				model->restore(model);
				op->rejected_count++;
			}
			
			if(op->optimize != NULL /*&& iter % mcmc->tuning_frequency == 0 && iter >	1000*/){
				op->optimize(op, alpha);
			}
		}
		iter++;
		
		for (int i = 0; i < mcmc->log_count; i++) {
			if(iter % mcmc->logs[i]->every == 0){
				mcmc->logs[i]->write(mcmc->logs[i], iter);
			}
		}
#ifdef DEBUG_OPERATORS
		if(iter % DEBUG_OPERATORS == 0)
		for (int i = 0; i < mcmc->operator_count; i++) {
			Operator* op = mcmc->operators[i];
			double prob = (double)op->accepted_count/(op->accepted_count+op->rejected_count);
			fprintf(ofile, "%zu,%s,%zu,%zu,%f,%f\n", iter, op->name, op->accepted_count, op->rejected_count, prob,(op->parameters != NULL ? op->parameters[0]: -1));
		}
		fflush(ofile);
#endif
		if(mcmc->interruptible && mcmc_interrupted){
			break;
		}
	}
#ifdef DEBUG_OPERATORS
	fclose(ofile);
#endif
	
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->finalize(mcmc->logs[i]);
	}
	
	if( mcmc->verbose > 0){
		for (int i = 0; i < mcmc->operator_count; i++) {
			Operator* op = mcmc->operators[i];
			printf("Acceptance ratio %s: %f", op->name, ((double)op->accepted_count/(op->accepted_count+op->rejected_count)));
			if(op->parameters != NULL) printf(" %f", op->parameters[0]);
			printf(" (failures: %zu; attempts: %zu)\n", op->failure_count, op->accepted_count+op->rejected_count);
		}
	}
	free(weights);
	free_StringBuffer(buffer);
	if (mcmc->interruptible) {
		signal(SIGINT, SIG_DFL);// restore the default handler
	}
}

void _free_MCMC(MCMC* mcmc){
	mcmc->model->free(mcmc->model);
	for (int i = 0; i < mcmc->operator_count; i++) {
		free_Operator(mcmc->operators[i]);
	}
	for (int i = 0; i < mcmc->log_count; i++) {
		mcmc->logs[i]->free(mcmc->logs[i]);
	}
	free(mcmc->logs);
	free(mcmc->operators);
	free(mcmc);
}

MCMC* new_MCMC_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"bf",
		"interruptible",
		"generalized",
		"length",
		"log",
		"model",
		"operators",
		"temperature",
		"tuningfrequency",
		"verbose"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	MCMC* mcmc = malloc(sizeof(MCMC));
	json_node* model_node = get_json_node(node, "model");
	mcmc->operator_count = 0;
	json_node* ops = get_json_node(node, "operators");
	json_node* logs = get_json_node(node, "log");
	mcmc->chain_length = get_json_node_value_size_t(node, "length", 100000);
	mcmc->tuning_frequency = get_json_node_value_size_t(node, "tuningfrequency", 1);
	mcmc->verbose = get_json_node_value_int(node, "verbose", 1);
	mcmc->interruptible = get_json_node_value_bool(node, "interruptible", true);
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
    if(logs != NULL){
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
    }
	mcmc->chain_temperature = get_json_node_value_double(node, "temperature", -1);
	mcmc->generalized = get_json_node_value_bool(node, "generalized", false);
	mcmc->bf = get_json_node_value_bool(node, "bf", false);
	mcmc->free = _free_MCMC;
	mcmc->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return mcmc;
}

