// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "klpq.h"

#include <string.h>

#include "matrix.h"
#include "gaussian.h"
#include "transforms.h"
#include "klqp.h"


#include "bound.h"
#include "compoundmodel.h"
#include "distmodelfactory.h"
#include "treelikelihood.h"
#include "matrix.h"


double _klpq_snis_logP(Bound* klpq){
	size_t infCount = 0;
    Model** distributions = &klpq->variational;
	size_t distCount = 1;
	if(klpq->variational->type == MODEL_COMPOUND){
		CompoundModel* cm = klpq->variational->obj;
		distributions = cm->models;
		distCount = cm->count;
	}

	double* logw = dvector(klpq->samples);
	double sum = -DBL_MAX;
	double estimate = 0;

    for(size_t i = 0; i < klpq->samples; i++){
		klpq->variational->sample(klpq->variational);
		double logQ = klpq->variational->logP(klpq->variational);
		double logP = klpq->joint->logP(klpq->joint);
		if(!isinf(logP) && !isnan(logP)){
			logw[i] = logP - logQ;
		}
		else{
			infCount++;
			i--;
		}
		if (infCount == 10) {
			free(logw);
			return NAN;
		}
		sum = logaddexp(logw[i], sum);
    }
	for (size_t i = 0; i < klpq->samples; i++) {
		estimate += exp(logw[i] - sum) * logw[i];
	}
	free(logw);
    return estimate;
}

double _klpq_snis_gradient(Bound* klpq, Parameters* parameters){
	size_t infCount = 0;
    Model** distributions = &klpq->variational;
	size_t distCount = 1;
	if(klpq->variational->type == MODEL_COMPOUND){
		CompoundModel* cm = klpq->variational->obj;
		distributions = cm->models;
		distCount = cm->count;
	}
	size_t totalSize = 0;
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		totalSize += Parameter_size(Parameters_at(parameters, i));
	}
	double* grads = dvector(totalSize*klpq->samples); 

	double* logw = dvector(klpq->samples);
	double sum = -DBL_MAX;
	double estimate = 0;
	size_t index = 0;

    for(size_t i = 0; i < klpq->samples; i++){
		klpq->variational->sample(klpq->variational);
		double logQ = klpq->variational->logP(klpq->variational);
		double logP = klpq->joint->logP(klpq->joint);

		Parameters_zero_grad(parameters);
		klpq->variational->gradient(klpq->variational, parameters);

		for(size_t k = 0; k < Parameters_count(parameters); k++){
			Parameter* p = Parameters_at(parameters, k);
			for(size_t j = 0; j < Parameter_size(p); j++){
				grads[index++] = p->grad[j];
			}
		}

		if(!isinf(logP) && !isnan(logP)){
			logw[i] = logP - logQ;
		}
		else{
			infCount++;
			i--;
		}
		if (infCount == 10) {
			free(logw);
			return NAN;
		}
		sum = logaddexp(logw[i], sum);
    }

	Parameters_zero_grad(parameters);
	index = 0;
	for (size_t i = 0; i < klpq->samples; i++) {
		double w = exp(logw[i] - sum);
		for(size_t k = 0; k < Parameters_count(parameters); k++){
			Parameter* p = Parameters_at(parameters, k);
			for(size_t j = 0; j < Parameter_size(p); j++){
				p->grad[j] -= grads[index++] * w;
			}
		}
		estimate += w * logw[i];
	}
	free(logw);
	free(grads);
    return -estimate;
}

Model* new_KLpqBound_from_json(json_node* node, Hashtable* hash) {
    Model* model = new_AbstractBoundModel_from_json(node, hash);
    Bound* bound = model->obj;
    bound->logP = _klpq_snis_logP;
    bound->gradient = _klpq_snis_gradient;
    return model;
}