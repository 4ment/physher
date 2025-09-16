//
//  is.c
//  physher
//
//  Created by Mathieu Fourment on 19/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "is.h"

#include "matrix.h"
#include "compoundmodel.h"
#include "utilsgsl.h"

double _calculate_importance_sampler( ImportanceSampler* mvb ){
	Model* mvar = mvb->distribution[0];
	Model* posterior = mvb->model;
	size_t n = mvb->samples;
	double lmarg;
	size_t dim = Parameters_count(mvb->parameters);
	double* backup = dvector(dim);
	Parameters_store_value(mvb->parameters, backup);
	if(mvb->normalize == false){
		double sum = -DBL_MAX;
		
		for(size_t i = 0; i < n; i++){
			double logQ;
			if(mvb->distribution_count > 1){
				size_t idx = roulette_wheel_gsl(mvb->rng, mvb->weights, mvb->distribution_count);
				mvar = mvb->distribution[idx];
			}
			// mvar->sample(mvar, samples, &logQ);
			// Parameters_restore_value(mvb->parameters, samples);
			mvar->sample(mvar);
			logQ = mvar->logP(mvar);
			double logP = posterior->logP(posterior);
			sum = logaddexp(sum, logP-logQ);
			//            printf(",%f",logP);
		}
		lmarg = sum-log(n);
	}
	else{
		double sum = -DBL_MAX;
		double denom = -DBL_MAX;
		CompoundModel* cm = posterior->obj;
		Model* lk = cm->models[0];
		Model* prior = cm->models[1];
		for(size_t i = 0; i < n; i++){
			double logQ;
			if(mvb->distribution_count > 1){
				size_t idx = roulette_wheel_gsl(mvb->rng, mvb->weights, mvb->distribution_count);
				mvar = mvb->distribution[idx];
			}
			// mvar->sample(mvar, samples, &logQ);
			// Parameters_restore_value(mvb->parameters, samples);
			mvar->sample(mvar);
			logQ = mvar->logP(mvar);
			posterior->logP(posterior);
			double logP = lk->lp;
			double logPrior = prior->lp;
			double w = logPrior - logQ;
			sum = logaddexp(sum, logP+w);
			denom = logaddexp(denom, w);
			//            printf(",%f",logP);
		}
		lmarg = sum - denom;
	}
	Parameters_restore_value(mvb->parameters, backup);
	
	return lmarg;
	
}
void _free_ImportanceSampler(ImportanceSampler* mvb){
	mvb->model->free(mvb->model);
	for (int i = 0; i < mvb->distribution_count; i++) {
		free(mvb->distribution[i]);
	}
	free(mvb->distribution);
	free(mvb->weights);
	free_Parameters(mvb->parameters);
	free(mvb);
}

ImportanceSampler* new_ImportanceSampler_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"distribution",
		"model",
		"normalize",
		"parameters",
		"samples",
		"weights",
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* ref = get_json_node_value_string(node, "model");
	json_node* dist_node = get_json_node(node, "distribution");
	ImportanceSampler* mvb = malloc(sizeof(ImportanceSampler));
	mvb->model = Hashtable_get(hash, ref+1);
	mvb->model->ref_count++;
	mvb->distribution_count = 0;
	
	if (dist_node->node_type == MJSON_ARRAY) {
		json_node* weights_node = get_json_node(node, "weights");
		mvb->weights = dvector(weights_node->child_count);
		for (int i = 0; i < weights_node->child_count; i++) {
			json_node* child = weights_node->children[i];
			if(child->node_type != MJSON_PRIMITIVE){
				fprintf(stderr, "sadf\n");
				exit(1);
			}
			mvb->weights[i] = atof((char*)child->value);
		}
		
		mvb->distribution = malloc(sizeof(Model*)*dist_node->child_count);
		for (int i = 0; i < dist_node->child_count; i++) {
			char* ref_dist = get_json_node_value_string(dist_node->children[i], "distribution");
			mvb->distribution[i] = Hashtable_get(hash, ref_dist+1);
			mvb->distribution[i]->ref_count++;
			mvb->distribution_count++;
		}
	}
	else{
		char* ref_dist = get_json_node_value_string(node, "distribution");
		mvb->weights = dvector(1);
		mvb->weights[0] = 1;
		mvb->distribution = malloc(sizeof(Model*));
		mvb->distribution[0] = Hashtable_get(hash, ref_dist+1);
		mvb->distribution[0]->ref_count++;
		mvb->distribution_count = 1;
	}
	mvb->samples = get_json_node_value_size_t(node, "samples", 100000);
	mvb->calculate = _calculate_importance_sampler;
	mvb->free = _free_ImportanceSampler;
	mvb->parameters = new_Parameters(1);
	get_parameters_references(node, hash, mvb->parameters);
	mvb->normalize = get_json_node_value_bool(node, "normalize", true);
	mvb->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return mvb;
}


