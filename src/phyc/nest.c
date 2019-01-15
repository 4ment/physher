//
//  nest.c
//  physher
//
//  Created by Mathieu Fourment on 22/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "nest.h"

#include "matrix.h"
#include "random.h"

double nest_mcmc(NEST* mcmc, double minLnL){
	Model* prior = mcmc->prior;
	
	long iter = -mcmc->burnin;
	
//	int tuneFrequency = 10;
	
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
//	for (int i = 0; i < mcmc->log_count; i++) {
//		mcmc->logs[i]->initialize(mcmc->logs[i]);
//	}
	
	logP = prior->logP(prior);
//	for (int i = 0; i < mcmc->log_count; i++) {
//		_mcmc_write_header(mcmc->logs[i], false);
//		mcmc->logs[i]->write(mcmc->logs[i], iter);
//	}
	double logLikelihood = INFINITY;
	while ( iter < mcmc->chain_length) {
		
		// Choose operator
		int op_index = roulette_wheel(weights, mcmc->operator_count);
		Operator* op = mcmc->operators[op_index];
		
		double logHR = 0;
		
		prior->store(prior);
		op->store(op);
		bool success = op->propose(op, &logHR);
		
		if (success) {
			double proposed_logPrior;
			double proposed_logLikelihood = INFINITY;
			
			if(!isinf(minLnL)){
				proposed_logLikelihood = mcmc->likelihood->logP(mcmc->likelihood);
			}

			proposed_logPrior = prior->logP(prior);
			
			double alpha = proposed_logPrior - logP + logHR;
			
			// accept
			if ( (alpha >=  0 || alpha > log(random_double())) && proposed_logLikelihood > minLnL ) {
				//			printf("%zu %f %f\n", iter, logP, proposed_logP);
				logP = proposed_logPrior;
				logLikelihood = proposed_logLikelihood;
				op->accepted_count++;
			}
			// reject
			else {
				//			printf("%zu %f %f *\n", iter, logP, proposed_logP);
				prior->restore(prior);
				op->restore(op);
				op->rejected_count++;
			}
			
			if(op->optimize != NULL){
				op->optimize(op, alpha);
			}
			
			iter++;
		}
	
//		for (int i = 0; i < mcmc->log_count; i++) {
//			if(iter % mcmc->logs[i]->every == 0){
//				mcmc->logs[i]->write(mcmc->logs[i], iter);
//			}
//		}
		
	}
	
//	for (int i = 0; i < mcmc->log_count; i++) {
//		mcmc->logs[i]->finalize(mcmc->logs[i]);
//	}
	
//	if( mcmc->verbose > 0){
//		for (int i = 0; i < mcmc->operator_count; i++) {
//			Operator* op = mcmc->operators[i];
//			printf("Acceptance ratio %s: %f %f (failures: %zu)\n", op->name, ((double)op->accepted_count/(op->accepted_count+op->rejected_count)), op->parameters[0], op->failure_count);
//		}
//	}
	free(weights);
	free_StringBuffer(buffer);
	return logLikelihood;
}

void nest_run(NEST* nest){
	double logZ = -DBL_MAX;
	double H = -DBL_MAX;
	size_t N = nest->N;
	double NN = N;
	size_t paramCount = Parameters_count(nest->x);
	double* logL = dvector(nest->N);
	double** params = dmatrix(N, paramCount);
	
	for (int i = 0; i < N; i++) {
		nest_mcmc(nest, -INFINITY);
		Parameters_store_value(nest->x, params[i]);
		logL[i] = nest->likelihood->logP(nest->likelihood);
	}
	
	double logw = log(1.0 - exp(-1.0/N));
	
	for (int i = 1; i <= nest->steps; i++) {
		double logZp = logZ;
		double ii = i;
		size_t worst = which_dmin(logL, N);
		double lw = logL[worst] + logw - (ii-1.0)/NN;
		logZ = logaddexp(logZ, lw);
		H = exp(lw - logZ) * logL[worst] - logZ + exp(logZp - logZ)*(H + logZp);

		if(i % 100 == 0) printf("%d %f %f", i, logZ, logL[worst]);
		
		size_t index = 0;
		if(N > 1){
			index = random_int(N-1);
			while (index == worst) {
				index = random_int(N-1);
			}
		}
		Parameters_restore_value(nest->x, params[index]);
		double logLikelihood = nest_mcmc(nest, logL[worst]);
		Parameters_store_value(nest->x, params[worst]);
		logL[worst] = logLikelihood;
		
		if(i % 100 == 0) printf(" %f\n", logLikelihood);

		if (fabs(logZ-logZp) < nest->precision) {
			break;
		}
	}

	double sum = -DBL_MAX;
	for (int i = 0; i < N; i++) {
		sum = logaddexp(sum, logL[i]);
	}

	sum += -(double)nest->steps/N - log(NN);
	logZ = logaddexp(sum, logZ);
	
	printf("logZ: %f\n", logZ);
	printf("H: %f\n", H);
	printf("stdev: %f\n", sqrt(H/N));
	free(logL);
	for (int i = 0; i < N; i++) {
		free(params[i]);
	}
	free(params);
}

NEST* new_NEST_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"burnin",
		"length",
		"likelihood",
		"N",
		"operators",
		"precision",
		"prior",
		"steps",
		"x"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	NEST* nest = malloc(sizeof(NEST));
	
	json_node* likelihood_node = get_json_node(node, "likelihood");
	json_node* prior_node = get_json_node(node, "prior");
	
	json_node* ops = get_json_node(node, "operators");
	
	nest->chain_length = get_json_node_value_size_t(node, "length", 100000);
	nest->steps = get_json_node_value_size_t(node, "steps", 1000);
	nest->burnin = get_json_node_value_size_t(node, "burnin", 0);
	nest->N = get_json_node_value_size_t(node, "N", 100);
	nest->precision = get_json_node_value_double(node, "precision", 1.e-6);
	
	nest->run = nest_run;
	
	nest->x = new_Parameters(1);
	get_parameters_references2(node, hash, nest->x, "x");
	
	// Prior as a compound
	char* model_string = prior_node->value;
	nest->prior = Hashtable_get(hash, model_string+1);
	nest->prior->ref_count++;
	
	// Treelikelihood (compound or model)
	model_string = likelihood_node->value;
	nest->likelihood = Hashtable_get(hash, model_string+1);
	nest->likelihood->ref_count++;
	
	nest->operator_count = 0;
	for (int i = 0; i < ops->child_count; i++) {
		json_node* child = ops->children[i];
		Operator* op = new_Operator_from_json(child, hash);
		char* id_string = get_json_node_value_string(child, "id");
		Hashtable_add(hash, id_string, op);
		if(nest->operator_count == 0){
			nest->operators = malloc(sizeof(Operator*));
		}
		else{
			nest->operators = realloc(nest->operators, sizeof(Operator*)*(nest->operator_count+1));
		}
		nest->operators[nest->operator_count] = op;
		nest->operator_count++;
	}
	
//	mmcmc->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return nest;
}

void free_NEST(NEST* nest){
	nest->likelihood->free(nest->likelihood);
	nest->prior->free(nest->prior);
	free_Parameters(nest->x);
	free(nest);
}
