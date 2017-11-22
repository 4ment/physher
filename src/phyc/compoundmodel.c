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

#include "treelikelihood.h"
#include "distmodel.h"


double _compoundModel_logP(CompoundModel* cm){
	double logP = 0;
	for(int i = 0; i < cm->count; i++){
		logP += cm->models[i]->logP(cm->models[i]);
	}
	return logP;
}

double _compoundModel_dlogP(CompoundModel* cm, const Parameter* p){
	double dlogP = 0;
	for(int i = 0; i < cm->count; i++){
		dlogP += cm->models[i]->dlogP(cm->models[i], p);
	}
	return dlogP;
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
	clone->free = cm->free;
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
	Model* clone = new_CompoundModel2(self->name, cmclone);
	
	Hashtable_add(hash, clone->name, clone);
	return clone;
}

static void _compound_model_get_free_parameters(Model* model, Parameters* parameters){
	CompoundModel* cm = model->obj;
	for(int i = 0; i < cm->count; i++){
		Model* mm = (Model*)cm->models[i];
		mm->get_free_parameters(mm, parameters);
	}
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
	cm->dlogP = _compoundModel_dlogP;
	cm->free = _free_compound_model;
	return cm;
}

static void _compound_model_free( Model *self ){
	assert(self->ref_count >= 1);
	if(self->ref_count == 1){
		printf("Free compound model %s\n", self->name);
		CompoundModel* cm = (CompoundModel*)self->obj;
		cm->free(cm);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}


double _compoundModel_logP2(Model *self){
	CompoundModel* cm = (CompoundModel*)self->obj;
	return cm->logP(cm);
}

double _compoundModel_dlogP2(Model *self, const Parameter* p){
	CompoundModel* cm = (CompoundModel*)self->obj;
	return cm->dlogP(cm, p);
}

Model* new_CompoundModel2(const char* name, CompoundModel* cm){
	Model *model = new_Model(name, cm);
	model->logP = _compoundModel_logP2;
	model->dlogP = _compoundModel_dlogP2;
	model->free = _compound_model_free;
	model->clone = _compound_model_clone;
	model->get_free_parameters = _compound_model_get_free_parameters;
	return model;
}

Model* new_CompoundModel_from_json(json_node*node, Hashtable*hash){
	CompoundModel* cm = new_CompoundModel();
	json_node* distributions_node = get_json_node(node, "distributions");
	char* id = get_json_node_value_string(node, "id");
	assert(distributions_node);

	for (int i = 0; i < distributions_node->child_count; i++) {
		json_node* child = distributions_node->children[i];
		char* type = get_json_node_value_string(child, "type");
		printf("type %s\n", type);
		if(strcasecmp(type, "treelikelihood") == 0){
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
				exit(1);
			}
			cm->add(cm, likelihood);
			likelihood->free(likelihood);
		}
		else if (strcasecmp(type, "distribution") == 0){
			Model* model = new_DistributionModel_from_json(child, hash);
			cm->add(cm, model);
			model->free(model);
			char* id = get_json_node_value_string(child, "id");
			Hashtable_add(hash, id, model);
		}
		else{
			printf("json CompoundModel unknown: (%s)\n", child->type);
			exit(1);
		}
	}
	Model* model = new_CompoundModel2(id, cm);
	return model;
}