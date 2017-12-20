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
	Model* mvar = mvb->model;
	variational_t* var = mvar->obj;
	size_t n = mvb->samples;
	
	Model* posterior = var->posterior;
	double lmarg;
	if(false){
		size_t dim = Parameters_count(var->parameters);
		double* samples = dvector(dim);
		double sum = -DBL_MAX;
		
		for(size_t i = 0; i < n; i++){
			var->sample(var, samples);
			Parameters_restore_value(var->parameters, samples);
			double logP = posterior->logP(posterior);
			double logQ = var->logP(var, samples);
			sum = logaddexp(sum, logP-logQ);
			//            printf(",%f",logP);
		}
		lmarg = sum-log(n);
	}
	else{
		CompoundModel* cm = posterior->obj;
		//TODO: Likleihoods and joint prior should be encapuslated in compound and given as input in json
		Model* lk = cm->models[0];
		size_t dim = Parameters_count(var->parameters);
		double* samples = dvector(dim);
		double sum = -DBL_MAX;
		double denom = -DBL_MAX;
		
		for(size_t i = 0; i < n; i++){
			var->sample(var, samples);
			Parameters_restore_value(var->parameters, samples);
			double logP = lk->logP(lk);
			double logPrior = 0;
			for(int j = 1; j < cm->count; j++){
				logPrior += cm->models[j]->logP(cm->models[j]);
			}
			double logQ = var->logP(var, samples);
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
	free(mvb);
}

marginal_vb_t* new_Marginal_VB_from_json(json_node* node, Hashtable* hash){
	char* ref = get_json_node_value_string(node, "model");
	marginal_vb_t* mvb = malloc(sizeof(marginal_vb_t));
	mvb->model = Hashtable_get(hash, ref+1);
	mvb->model->ref_count++;
	mvb->samples = get_json_node_value_size_t(node, "samples", 10000);
	mvb->calculate = _calculate_marginal_vb;
	mvb->free = _free_Marginal_VB;
	return mvb;
}


