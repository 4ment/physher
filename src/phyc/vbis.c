//
//  vbis.c
//  physher
//
//  Created by Mathieu Fourment on 19/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "vbis.h"

#include "matrix.h"
#include "compoundmodel.h"

double _calculate_marginal_vb( marginal_vb_t* mvb ){
	Model* mvar = mvb->distribution;
	Model* posterior = mvb->model;
	size_t n = mvb->samples;
	double lmarg;
	size_t dim = Parameters_count(mvb->parameters);
	if(mvb->normalize == false){
		double* samples = dvector(dim);
		double sum = -DBL_MAX;
		
		for(size_t i = 0; i < n; i++){
			double logQ;
			mvar->sample(mvar, samples, &logQ);
			Parameters_restore_value(mvb->parameters, samples);
			double logP = posterior->logP(posterior);
			sum = logaddexp(sum, logP-logQ);
			//            printf(",%f",logP);
		}
		lmarg = sum-log(n);
		free(samples);
	}
	else{
		double* samples = dvector(dim);
		double sum = -DBL_MAX;
		double denom = -DBL_MAX;
		CompoundModel* cm = posterior->obj;
		Model* lk = cm->models[0];
		Model* prior = cm->models[1];
		for(size_t i = 0; i < n; i++){
			double logQ;
			mvar->sample(mvar, samples, &logQ);
			Parameters_restore_value(mvb->parameters, samples);
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
	
	return lmarg;
	
}
void _free_Marginal_VB(marginal_vb_t* mvb){
	mvb->model->free(mvb->model);
	mvb->distribution->free(mvb->distribution);
	free_Parameters(mvb->parameters);
	free(mvb);
}

marginal_vb_t* new_Marginal_VB_from_json(json_node* node, Hashtable* hash){
	char* ref = get_json_node_value_string(node, "model");
	char* ref_dist = get_json_node_value_string(node, "distribution");
	marginal_vb_t* mvb = malloc(sizeof(marginal_vb_t));
	mvb->model = Hashtable_get(hash, ref+1);
	mvb->model->ref_count++;
	mvb->distribution = Hashtable_get(hash, ref_dist+1);
	mvb->distribution->ref_count++;
	mvb->samples = get_json_node_value_size_t(node, "samples", 100000);
	mvb->calculate = _calculate_marginal_vb;
	mvb->free = _free_Marginal_VB;
	mvb->parameters = new_Parameters(1);
	get_parameters_references(node, hash, mvb->parameters);
	mvb->normalize = get_json_node_value_bool(node, "normalize", true);
	return mvb;
}


