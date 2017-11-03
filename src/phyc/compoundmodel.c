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

CompoundModel* new_CompoundModel(){
	CompoundModel* cm = (CompoundModel*)malloc(sizeof(CompoundModel));
	assert(cm);
	cm->models = (Model**)malloc(sizeof(Model*)*2);
	cm->models[0] = NULL;
	cm->models[1] = NULL;
	assert(cm->models);
	cm->count = 0;
	cm->add = _compoundModel_add;
	cm->remove = _compoundModel_remove;
	cm->removeAll = _compoundModel_remove_all;
	cm->logP = _compoundModel_logP;
	cm->dlogP = _compoundModel_dlogP;
	cm->free = _free_compound_model;
	return cm;
}

static void _compound_model_free( Model *self ){
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
	return model;
}