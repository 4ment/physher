// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "compoundmodel.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>

#include "treelikelihood.h"
#include "distmodelfactory.h"
#include "demographicmodels.h"
#include "jacobian.h"
#include "matrix.h"


double _compoundModel_logP(CompoundModel* cm){
	double logP = 0;
	if (cm->weights != NULL) {
		logP = -DBL_MAX;
		const double* weights = Parameter_values(cm->weights);
		for(int i = 0; i < cm->count; i++){
			logP = logaddexp(logP, log(weights[i]) + cm->models[i]->logP(cm->models[i]));
		}
	}
	else{
		for(int i = 0; i < cm->count; i++){
			logP += cm->models[i]->logP(cm->models[i]);
		}
	}
	return logP;
}

double _compoundModel_full_logP(CompoundModel* cm){
	double logP = 0;
	if (cm->weights != NULL) {
		logP = -DBL_MAX;
		const double* weights = Parameter_values(cm->weights);
		for(int i = 0; i < cm->count; i++){
			logP = logaddexp(logP, log(weights[i]) + cm->models[i]->full_logP(cm->models[i]));
		}
	}
	else{
		for(int i = 0; i < cm->count; i++){
			logP += cm->models[i]->full_logP(cm->models[i]);
		}
	}
	return logP;
}

static void _compoundModel_add(CompoundModel* cm, Model*model){
	cm->models = realloc(cm->models, sizeof(Model*)*(cm->count+1));
	cm->models[cm->count] = model;
	model->ref_count++;
	cm->count++;
}

static void _compoundModel_move(CompoundModel* cm, Model*model){
	cm->models = realloc(cm->models, sizeof(Model*)*(cm->count+1));
	cm->models[cm->count] = model;
	cm->count++;
}

static void _compoundModel_remove( CompoundModel* cm, Model*model ){
	int i = 0;
	for ( ; i < cm->count; i++ ) {
		if ( cm->models[i] == model ) {
			break;
		}
	}
	if ( i == cm->count) {
		return;
	}
	i++;
	for ( ; i < cm->count; i++ ) {
		cm->models[i-1] = cm->models[i];
	}
	cm->models[cm->count-1] = NULL;
	cm->count--;
	model->ref_count--;
}

static void _compoundModel_remove_all( CompoundModel* cm ){
	for ( int i = 0; i < cm->count; i++ ) {
		cm->models[i]->ref_count--;
		cm->models[i] = NULL;
	}
	cm->count = 0;
}

static void _free_compound_model(CompoundModel* cm){
	for (int i = 0; i < cm->count; i++) {
		cm->models[i]->free(cm->models[i]);
	}
	if(cm->weights!=NULL) free_Parameter(cm->weights);
	free(cm->models);
	free(cm);
}

CompoundModel* clone_compound_model(CompoundModel* cm){
	CompoundModel* clone = new_CompoundModel();
	clone->add = cm->add;
	clone->move = cm->move;
	clone->remove = cm->remove;
	clone->removeAll = cm->removeAll;
	clone->logP = cm->logP;
	clone->free = cm->free;
    clone->weights = NULL;
	return clone;
}

#pragma mark-
#pragma mark Model

static Model* _compound_model_clone( Model *self, Hashtable* hash ){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	CompoundModel* cm = self->obj;
	CompoundModel* cmclone = clone_compound_model(cm);
	for (int i = 0; i < cm->count; i++) {
		Model* m = cm->models[i];
		Model* mclone = NULL;
		if (Hashtable_exists(hash, m->name)) {
			mclone = Hashtable_get(hash, m->name);
			mclone->ref_count++;
		}
		else{
			mclone = m->clone(m, hash);
			Hashtable_add(hash, mclone->name, mclone);
		}
		cmclone->add(cmclone, mclone);
		mclone->free(mclone);
	}
	if(cm->weights != NULL){
		Parameter* weights_clone = NULL;
		if (Hashtable_exists(hash, Parameter_name(cm->weights))) {
			weights_clone = Hashtable_get(hash, Parameter_name(cm->weights));
			weights_clone->refCount++;
		}
		else{
			weights_clone = clone_Parameter(cm->weights);
			Hashtable_add(hash, Parameter_name(weights_clone), weights_clone);
		}
		cmclone->weights = weights_clone;
	}
	Model* clone = new_CompoundModel2(self->name, cmclone);
	
	Hashtable_add(hash, clone->name, clone);
	clone->store = self->store;
	clone->restore = self->restore;
	clone->storedLogP = self->storedLogP;
	clone->lp = self->lp;
	clone->samplable = self->samplable;
	clone->sample = self->sample;
	clone->full_logP = self->full_logP;
	return clone;
}

CompoundModel* new_CompoundModel(){
	CompoundModel* cm = (CompoundModel*)malloc(sizeof(CompoundModel));
	assert(cm);
	cm->models = (Model**)malloc(sizeof(Model*)*2);
	cm->models[0] = NULL;
	cm->models[1] = NULL;
	assert(cm->models);
	cm->count = 0;
	cm->add = _compoundModel_add;
	cm->move = _compoundModel_move;
	cm->remove = _compoundModel_remove;
	cm->removeAll = _compoundModel_remove_all;
	cm->logP = _compoundModel_logP;
	cm->full_logP = _compoundModel_full_logP;
	cm->free = _free_compound_model;
	return cm;
}

static void _compound_model_free( Model *self ){
	assert(self->ref_count >= 1);
	if(self->ref_count == 1){
		//printf("Free compound model %s\n", self->name);
		CompoundModel* cm = (CompoundModel*)self->obj;
		cm->free(cm);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static void _compoundModel_store(Model* self){
	if(!self->stored){
		self->storedLogP = self->lp;
		CompoundModel* cm = (CompoundModel*)self->obj;
		for (int i = 0; i < cm->count; i++) {
			cm->models[i]->store(cm->models[i]);
		}
		if (cm->weights != NULL) {
			Parameter_store(cm->weights);
		}
		self->stored = true;
	}
}

static void _compoundModel_restore(Model* self){
	if(self->stored){
		self->lp = self->storedLogP;
		CompoundModel* cm = (CompoundModel*)self->obj;
		for (int i = 0; i < cm->count; i++) {
			cm->models[i]->restore(cm->models[i]);
		}
		if (cm->weights != NULL) {
			Parameter_restore(cm->weights);
		}
		self->stored = false;
	}
}

static void _compoundModel_accept(Model* self){
	if(self->stored){
		CompoundModel* cm = (CompoundModel*)self->obj;
		for (size_t i = 0; i < cm->count; i++) {
			cm->models[i]->accept(cm->models[i]);
		}
		if (cm->weights != NULL) {
			Parameter_accept(cm->weights);
		}
		self->stored = false;
	}
}

double _compoundModel_logP2(Model *self){
	CompoundModel* cm = (CompoundModel*)self->obj;
	self->lp = cm->logP(cm);
	return self->lp;
}

double _compoundModel_full_logP2(Model *self){
	CompoundModel* cm = (CompoundModel*)self->obj;
	self->lp = cm->full_logP(cm);
	return self->lp;
}

void _compoundModel_gradient2(Model* self, Parameters* ps){
	CompoundModel* cm = (CompoundModel*)self->obj;

	if (cm->weights == NULL) {
		double logP = 0;
		for (size_t i = 0; i < cm->count; i++) {
			cm->models[i]->gradient(cm->models[i], ps);
		}
		return;
	}

	// Mixture: logP = logsumexp_i(log w_i + logP_i). The gradient is the
	// responsibility-weighted sum  d logP = sum_i r_i * d logP_i  with
	// r_i = exp(log w_i + logP_i - logP). Because submodels accumulate into
	// the shared p->grad, each component must be computed in isolation,
	// scaled by r_i, and summed into a separate accumulator.
	size_t parameterCount = Parameters_count(ps);
	size_t size = Parameters_size(ps);

	const double* weights = Parameter_values(cm->weights);

	// Pass 1: per-component logP_i and the total logP.
	double* logPi = dvector(cm->count);
	double logP = -DBL_MAX;
	for (size_t i = 0; i < cm->count; i++) {
		logPi[i] = cm->models[i]->logP(cm->models[i]);
		logP = logaddexp(logP, log(weights[i]) + logPi[i]);
	}

	// Pass 2: ps must not be zeroed - it may hold the weights, parameters
	// shared with the mixture components, or parameters belonging to sibling
	// models that already accumulated into these buffers. The leaf gradients
	// accumulate (+=), so after each component call the increment over the
	// running total is exactly that component's own gradient g_i. The calls
	// leave the unweighted sum (sum_i g_i) in p->grad; we want the weighted sum
	// (sum_i r_i * g_i), so we collect the correction sum_i (r_i - 1) * g_i and
	// add it at the end, never overwriting the pre-existing gradient.
	double* running = dvector(size);
	size_t index = 0;
	for (size_t j = 0; j < parameterCount; j++) {
		Parameter* p = Parameters_at(ps, j);
		memcpy(running + index, p->grad, sizeof(double) * Parameter_size(p));
		index += Parameter_size(p);
	}

	double* correction = dvector(size);
	for (size_t i = 0; i < cm->count; i++) {
		double r = exp(log(weights[i]) + logPi[i] - logP);
		cm->models[i]->gradient(cm->models[i], ps);
		index = 0;
		for (size_t j = 0; j < parameterCount; j++) {
			Parameter* p = Parameters_at(ps, j);
			for (size_t k = 0; k < Parameter_size(p); k++) {
				double gi = p->grad[k] - running[index];
				correction[index] += (r - 1.0) * gi;
				running[index] = p->grad[k];
				index++;
			}
		}
	}

	// Accumulate the correction to turn the unweighted sum left in p->grad into
	// the responsibility-weighted mixture gradient.
	index = 0;
	for (size_t j = 0; j < parameterCount; j++) {
		Parameter* p = Parameters_at(ps, j);
		for (size_t k = 0; k < Parameter_size(p); k++) {
			p->grad[k] += correction[index++];
		}
	}

	// Gradient wrt the (constrained) weights: d logP / d w_i = exp(logP_i - logP).
	// Only computed when ps contains the weights or their unconstrained version.
	Parameter* weightsx = Parameters_depends(ps, cm->weights);
	if (weightsx != NULL) {
		double* dweights = dvector(cm->count);
		for (size_t i = 0; i < cm->count; i++) {
			dweights[i] = exp(logPi[i] - logP);
		}
		if (weightsx == cm->weights) {
			for (size_t i = 0; i < cm->count; i++) {
				cm->weights->grad[i] += dweights[i];
			}
		}
		else {
			// chain rule into the unconstrained parameter's grad (accumulates).
			cm->weights->transform->backward(cm->weights->transform, dweights);
		}
		free(dweights);
	}

	free(logPi);
	free(running);
	free(correction);
}

void _compound_model_sample(Model *self){
	CompoundModel* cm = (CompoundModel*)self->obj;
	for (size_t i = 0; i < cm->count; i++) {
		cm->models[i]->sample(cm->models[i]);
	}
}

void _compound_model_rsample(Model *self){
	CompoundModel* cm = (CompoundModel*)self->obj;
	for (size_t i = 0; i < cm->count; i++) {
		cm->models[i]->rsample(cm->models[i]);
	}
}

// double _compound_model_sample_evaluate(Model *self){
// 	self->lp = 0;
// 	CompoundModel* cm = (CompoundModel*)self->obj;
// 	for (int i = 0; i < cm->count; i++) {
// 		self->lp += cm->models[i]->sample_evaluate(cm->models[i]);
// 	}
// 	return self->lp;
// }

// For a summed compound (logP = sum_i logP_i) the Hessian is the sum of the
// submodel Hessians, so each submodel's own (analytic or FD) hessian is used
// and accumulated. The mixture case (weighted) has extra rank-1 terms and is
// left to finite differences.
void _compoundModel_hessian(Model* self, const Parameters* parameters, hessian_mode_t mode, double* out){
	CompoundModel* cm = (CompoundModel*)self->obj;
	if(cm->weights != NULL){
		Model_hessian_fd(self, parameters, mode, out);
		return;
	}
	size_t dim = Parameters_size(parameters);
	size_t n = (mode == HESSIAN_FULL) ? dim*dim : dim;
	memset(out, 0, sizeof(double)*n);
	double* tmp = dvector(n);
	for(size_t i = 0; i < cm->count; i++){
		cm->models[i]->hessian(cm->models[i], parameters, mode, tmp);
		for(size_t k = 0; k < n; k++) out[k] += tmp[k];
	}
	free(tmp);
}

Model* new_CompoundModel2(const char* name, CompoundModel* cm){
	Model *model = new_Model(MODEL_COMPOUND, name, cm);
	model->logP = _compoundModel_logP2;
	model->full_logP = _compoundModel_full_logP2;
	model->gradient = _compoundModel_gradient2;
	model->hessian = _compoundModel_hessian;
	model->free = _compound_model_free;
	model->clone = _compound_model_clone;
	model->store = _compoundModel_store;
	model->restore = _compoundModel_restore;
	model->accept = _compoundModel_accept;
	model->samplable = true;
	for (int i = 0; i < cm->count; i++) {
		if (!cm->models[i]->samplable) {
			model->samplable = false;
			break;
		}
	}
//	for(int i = 0; i < cm->count; i++){
//		cm->models[i]->listeners->add(cm->models[i]->listeners, model),
//	}
	if (cm->weights != NULL) {
		cm->weights->listeners->add( cm->weights->listeners, model );
		Parameters_add_recursively(model->parameters, cm->weights);
	}
	model->sample = _compound_model_sample;
	model->rsample = _compound_model_rsample;
	// model->sample_evaluate = _compound_model_sample_evaluate;
	return model;
}

Model* new_CompoundModel_from_json(json_node*node, Hashtable*hash){
	static const json_field schema[] = {
		{"distributions", JSON_REQUIRED, JSON_ARRAY},
		{"weights", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	
	CompoundModel* cm = new_CompoundModel();
	json_node* distributions_node = get_json_node(node, "distributions");
	json_node* simplex_node = get_json_node(node, "weights");
	char* id = get_json_node_value_string(node, "id");
	assert(distributions_node);

	for (int i = 0; i < distributions_node->child_count; i++) {
		json_node* child = distributions_node->children[i];
		
		if (child->node_type == MJSON_STRING) {
			char* ref = (char*)child->value;
			Model* model = safe_get_reference_model(ref, hash, id);
			model->ref_count++;
			cm->add(cm, model);
			model->free(model);
			continue;
		}
		
		char* type = get_json_node_value_string(child, "type");
		model_t model_type = check_model(type);
		if(model_type == MODEL_TREELIKELIHOOD){
			Model* likelihood = NULL;
			if (child->node_type == MJSON_OBJECT) {
				likelihood = new_TreeLikelihoodModel_from_json(child, hash);
				char* id = get_json_node_value_string(child, "id");
				Hashtable_add(hash, id, likelihood);
			}
			else if(child->node_type == MJSON_STRING){
				char* ref = (char*)child->value;
				likelihood = safe_get_reference_model(ref, hash, id);
				likelihood->ref_count++;
			}
			else{
				exit(10);
			}
			cm->add(cm, likelihood);
			likelihood->free(likelihood);
		}
		else if (model_type == MODEL_DISTRIBUTION){
			Model* compound = NULL;
			if (child->node_type == MJSON_OBJECT) {
				compound = new_DistributionModel_from_json(child, hash);
				char* id = get_json_node_value_string(child, "id");
				Hashtable_add(hash, id, compound);
			}
			else if(child->node_type == MJSON_STRING){
				char* ref = (char*)child->value;
				compound = safe_get_reference_model(ref, hash, id);
				compound->ref_count++;
			}
			else{
				exit(10);
			}
			cm->add(cm, compound);
			compound->free(compound);
		}else if (model_type == MODEL_COMPOUND){
			Model* compound = NULL;
			if (child->node_type == MJSON_OBJECT) {
				compound = new_CompoundModel_from_json(child, hash);
				char* id = get_json_node_value_string(child, "id");
				Hashtable_add(hash, id, compound);
			}
			else if(child->node_type == MJSON_STRING){
				char* ref = (char*)child->value;
				compound = safe_get_reference_model(ref, hash, id);
				compound->ref_count++;
			}
			else{
				exit(10);
			}
			cm->add(cm, compound);
			compound->free(compound);
		}
		else if(model_type == MODEL_COALESCENT){
			Model* coalescent = NULL;
			if (child->node_type == MJSON_OBJECT) {
				coalescent = new_CoalescentModel_from_json(child, hash);
				char* id = get_json_node_value_string(child, "id");
				Hashtable_add(hash, id, coalescent);
			}
			else if(child->node_type == MJSON_STRING){
				char* ref = (char*)child->value;
				coalescent = safe_get_reference_model(ref, hash, id);
				coalescent->ref_count++;
			}
			else{
				exit(10);
			}
			cm->add(cm, coalescent);
			coalescent->free(coalescent);
		}
		else if(model_type == MODEL_JACOBIAN){
			Model* jac = NULL;
			if (child->node_type == MJSON_OBJECT) {
				jac = new_JacobianModel_from_json(child, hash);
				char* id = get_json_node_value_string(child, "id");
				Hashtable_add(hash, id, jac);
			}
			else if(child->node_type == MJSON_STRING){
				char* ref = (char*)child->value;
				jac = safe_get_reference_model(ref, hash, id);
				jac->ref_count++;
			}
			else{
				exit(10);
			}
			cm->add(cm, jac);
			jac->free(jac);
		}
		else{
			printf("json CompoundModel unknown: (%s)\n", type);
			exit(1);
		}
	}
	
	cm->weights = NULL;
	
	// it's a mixture
	if (simplex_node != NULL) {
		if (simplex_node->node_type == MJSON_OBJECT) {
			cm->weights = new_Parameter_from_json(simplex_node, hash);
			Hashtable_add(hash, Parameter_name(cm->weights), cm->weights);
		}
		else if(simplex_node->node_type == MJSON_STRING){
			char* ref = (char*)simplex_node->value;
			cm->weights = safe_get_reference_model(ref, hash, id);
			cm->weights->refCount++;
		}
	}
	Model* model = new_CompoundModel2(id, cm);
	return model;
}
