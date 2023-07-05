//
//  compoundmodel.c
//  physher
//
//  Created by Mathieu Fourment on 1/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "compoundmodel.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>

#include "treelikelihood.h"
#include "distmodelfactory.h"
#include "demographicmodels.h"


double _compoundModel_logP(CompoundModel* cm){
	double logP = 0;
	if (cm->weights != NULL) {
		Simplex* s = cm->weights->obj;
		logP = -DBL_MAX;
		for(int i = 0; i < cm->count; i++){
			logP = logaddexp(logP, log(s->get_value(s, i)) + cm->models[i]->logP(cm->models[i]));
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
		Simplex* s = cm->weights->obj;
		logP = -DBL_MAX;
		for(int i = 0; i < cm->count; i++){
			logP = logaddexp(logP, log(s->get_value(s, i)) + cm->models[i]->full_logP(cm->models[i]));
		}
	}
	else{
		for(int i = 0; i < cm->count; i++){
			logP += cm->models[i]->full_logP(cm->models[i]);
		}
	}
	return logP;
}

double _compoundModel_dlogP(CompoundModel* cm, const Parameter* p){
	double dlogP = 0;
	if (cm->weights != NULL) {
		Simplex* s = cm->weights->obj;
		for(int i = 0; i < cm->count; i++){
			dlogP += log(s->get_value(s, i)) + cm->models[i]->dlogP(cm->models[i], p);
		}
	}
	else{
		for(int i = 0; i < cm->count; i++){
			dlogP += cm->models[i]->dlogP(cm->models[i], p);
		}
	}
	return dlogP;
}

double _compoundModel_d2logP(CompoundModel* cm, const Parameter* p){
	double d2logP = 0;
	if (cm->weights != NULL) {
		Simplex* s = cm->weights->obj;
		for(int i = 0; i < cm->count; i++){
			d2logP += log(s->get_value(s, i)) + cm->models[i]->d2logP(cm->models[i], p);
		}
	}
	else{
		for(int i = 0; i < cm->count; i++){
			d2logP += cm->models[i]->d2logP(cm->models[i], p);
		}
	}
	return d2logP;
}

double _compoundModel_ddlogP(CompoundModel* cm, const Parameter* p1, const Parameter* p2){
	double ddlogP = 0;
	if(cm->weights!=NULL){
		fprintf(stderr, "Mixture _compoundModel_ddlogP not yet implemented\n");
		exit(11);
	}
	for(int i = 0; i < cm->count; i++){
		ddlogP += cm->models[i]->ddlogP(cm->models[i], p1, p2);
	}
	return ddlogP;
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
	if(cm->weights!=NULL) cm->weights->free(cm->weights);
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
	clone->dlogP = cm->dlogP;
	clone->d2logP = cm->d2logP;
	clone->ddlogP = cm->ddlogP;
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
		Model* msimplex_rates_clone = NULL;
		if (Hashtable_exists(hash, cm->weights->name)) {
			msimplex_rates_clone = Hashtable_get(hash, cm->weights->name);
			msimplex_rates_clone->ref_count++;
		}
		else{
			msimplex_rates_clone = cm->weights->clone(cm->weights, hash);
			Hashtable_add(hash, msimplex_rates_clone->name, msimplex_rates_clone);
		}
		cmclone->weights = msimplex_rates_clone;
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
	cm->dlogP = _compoundModel_dlogP;
	cm->d2logP = _compoundModel_d2logP;
	cm->ddlogP = _compoundModel_ddlogP;
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
	self->storedLogP = self->lp;
	CompoundModel* cm = (CompoundModel*)self->obj;
	for (int i = 0; i < cm->count; i++) {
		cm->models[i]->store(cm->models[i]);
	}
	if (cm->weights != NULL) {
		cm->weights->store(cm->weights);
	}
}

static void _compoundModel_restore(Model* self){
	self->lp = self->storedLogP;
	CompoundModel* cm = (CompoundModel*)self->obj;
	for (int i = 0; i < cm->count; i++) {
		cm->models[i]->restore(cm->models[i]);
	}
	if (cm->weights != NULL) {
		cm->weights->restore(cm->weights);
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

void _compoundModel_prepare_gradient(Model* self, const Parameters* ps){
	CompoundModel* cm = (CompoundModel*)self->obj;
	for (size_t i = 0; i < cm->count; i++) {
		cm->models[i]->prepare_gradient(cm->models[i], ps);
	}
}

double _compoundModel_dlogP2(Model *self, const Parameter* p){
	CompoundModel* cm = (CompoundModel*)self->obj;
	return cm->dlogP(cm, p);
}

double _compoundModel_d2logP2(Model *self, const Parameter* p){
	CompoundModel* cm = (CompoundModel*)self->obj;
	return cm->d2logP(cm, p);
}

double _compoundModel_ddlogP2(Model *self, const Parameter* p1, const Parameter* p2){
	CompoundModel* cm = (CompoundModel*)self->obj;
	return cm->ddlogP(cm, p1, p2);
}

void _compound_model_sample(Model *self, double* samples, double* logP){
	if (samples != NULL) {
		fprintf(stderr, "_compound_model_sample not implemented\n");
		exit(1);
	}
	CompoundModel* cm = (CompoundModel*)self->obj;
	if(logP != NULL){
		*logP = 0;
	}
	for (size_t i = 0; i < cm->count; i++) {
		if(logP != NULL){
			double model_logP = 0;
			cm->models[i]->sample(cm->models[i], NULL, &model_logP);
			*logP += model_logP;
		}
		else{
			cm->models[i]->sample(cm->models[i], NULL, NULL);
		}
	}
}

double _compound_model_sample_evaluate(Model *self){
	self->lp = 0;
	CompoundModel* cm = (CompoundModel*)self->obj;
	for (int i = 0; i < cm->count; i++) {
		self->lp += cm->models[i]->sample_evaluate(cm->models[i]);
	}
	return self->lp;
}

Model* new_CompoundModel2(const char* name, CompoundModel* cm){
	Model *model = new_Model(MODEL_COMPOUND, name, cm);
	model->logP = _compoundModel_logP2;
	model->full_logP = _compoundModel_full_logP2;
	model->dlogP = _compoundModel_dlogP2;
	model->d2logP = _compoundModel_d2logP2;
	model->ddlogP = _compoundModel_ddlogP2;
	model->free = _compound_model_free;
	model->clone = _compound_model_clone;
	model->store = _compoundModel_store;
	model->restore = _compoundModel_restore;
	model->prepare_gradient = _compoundModel_prepare_gradient;
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
//	if (cm->weights != NULL) {
//		cm->weights->listeners->add( cm->weights->listeners, model );
//	}
	model->sample = _compound_model_sample;
	model->sample_evaluate = _compound_model_sample_evaluate;
	return model;
}

Model* new_CompoundModel_from_json(json_node*node, Hashtable*hash){
	char* allowed[] = {
		"distributions",
		"weights"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	CompoundModel* cm = new_CompoundModel();
	json_node* distributions_node = get_json_node(node, "distributions");
	json_node* simplex_node = get_json_node(node, "weights");
	char* id = get_json_node_value_string(node, "id");
	assert(distributions_node);

	for (int i = 0; i < distributions_node->child_count; i++) {
		json_node* child = distributions_node->children[i];
		
		if (child->node_type == MJSON_STRING) {
			char* ref = (char*)child->value;
			Model* model = Hashtable_get(hash, ref+1);
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
				likelihood = Hashtable_get(hash, ref+1);
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
				compound = Hashtable_get(hash, ref+1);
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
				compound = Hashtable_get(hash, ref+1);
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
				coalescent = Hashtable_get(hash, ref+1);
				coalescent->ref_count++;
			}
			else{
				exit(10);
			}
			cm->add(cm, coalescent);
			coalescent->free(coalescent);
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
			cm->weights = new_SimplexModel_from_json(simplex_node, hash);
			char* id = get_json_node_value_string(simplex_node, "id");
			Hashtable_add(hash, id, simplex_node);
		}
		else if(simplex_node->node_type == MJSON_STRING){
			char* ref = (char*)simplex_node->value;
			// check it starts with a &
			cm->weights = Hashtable_get(hash, ref+1);
			cm->weights->ref_count++;
		}
	}
	Model* model = new_CompoundModel2(id, cm);
	return model;
}
