//
//  discreteparameter.c
//  physher
//
//  Created by Mathieu Fourment on 1/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "discreteparameter.h"

#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "matrix.h"

static void _free_DiscreteParameter( DiscreteParameter *p ){
	p->listeners->free(p->listeners);
	free(p->values);
	free(p->stored_values);
	free(p);
}

static DiscreteParameter * _clone_DiscreteParameter( DiscreteParameter *p ){
	DiscreteParameter *pnew = new_DiscreteParameter(p->length);
	memcpy(pnew->values, p->values, sizeof(unsigned)*p->length);
	memcpy(pnew->stored_values, p->stored_values, sizeof(unsigned)*p->length);
	return pnew;
}

static void _set_value_discrete(DiscreteParameter* p, int index, unsigned value) {
    assert(index < p->length);
    p->values[index] = value;
    p->listeners->fire(p->listeners, NULL, NULL, index);
}

static void _set_values_discrete(DiscreteParameter* p, const unsigned* values) {
    memcpy(p->values, values, p->length * sizeof(unsigned));
    p->listeners->fire(p->listeners, NULL, NULL, -1);
}

DiscreteParameter * new_DiscreteParameter_with_postfix_values( const char *postfix, const unsigned* values, size_t dim ){
	DiscreteParameter *p = (DiscreteParameter *)malloc( sizeof(DiscreteParameter) );
	assert(p);
	
	if(values != NULL){
		p->values = clone_uivector(values, dim);
		p->stored_values = clone_uivector(values, dim);
	}
	else{
		p->values = uivector(dim);
		p->stored_values = uivector(dim);
	}
	p->length = dim;
	p->set_value = _set_value_discrete;
	p->set_values = _set_values_discrete;
	p->free = _free_DiscreteParameter;
	p->clone = _clone_DiscreteParameter;
	p->listeners = new_ListenerList(1);
	return p;
}


DiscreteParameter * new_DiscreteParameter_with_values( const unsigned* values, size_t dim ){
	return new_DiscreteParameter_with_postfix_values("", values, dim);
}

DiscreteParameter * new_DiscreteParameter_with_postfix( const char *postfix, size_t dim ){
	return new_DiscreteParameter_with_postfix_values(postfix, NULL, dim);
}

DiscreteParameter * new_DiscreteParameter( size_t dim ){
	return new_DiscreteParameter_with_postfix("", dim);
}

#pragma mark -
#pragma mark Model implementation

static void _discrete_parameter_model_handle_change(Model* self, Model* model,
                                                    Parameter* parameter, int index) {
    self->listeners->fire(self->listeners, self, parameter, index);
}

static void _discrete_parameter_model_handle_restore( Model *self, Model *model, int index ){
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _discrete_parameter_model_store(Model* self){
	DiscreteParameter* dp = (DiscreteParameter*)self->obj;
	memcpy(dp->stored_values, dp->values, dp->length*sizeof(unsigned));
}

static void _discrete_parameter_model_restore(Model* self){
	DiscreteParameter* dp = (DiscreteParameter*)self->obj;
	if (memcmp(dp->values, dp->stored_values, dp->length*sizeof(unsigned)) != 0) {
		memcpy(dp->values, dp->stored_values, dp->length*sizeof(unsigned));
		self->listeners->fire_restore(self->listeners, self, 0);
	}
}

static void _discrete_parameter_model_free( Model *self ){
	if(self->ref_count == 1){
		DiscreteParameter* dp = (DiscreteParameter*)self->obj;
		dp->free(dp);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _discrete_parameter_model_clone(Model* self, Hashtable *hash){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	DiscreteParameter* dp = (DiscreteParameter*)self->obj;
	DiscreteParameter* sclone = dp->clone(dp);
	Model *clone = new_DiscreteParameterModel(self->name, sclone);
	Hashtable_add(hash, clone->name, clone);
	return clone;
}

Model * new_DiscreteParameterModel( const char* name, DiscreteParameter *dp ){
	Model *model = new_Model(MODEL_DISCRETE_PARAMETER, name, dp);
	
	model->update = _discrete_parameter_model_handle_change;
	model->handle_restore = _discrete_parameter_model_handle_restore;
	model->store = _discrete_parameter_model_store;
	model->restore = _discrete_parameter_model_restore;
	model->free = _discrete_parameter_model_free;
	model->clone = _discrete_parameter_model_clone;
	
	dp->listeners->add(dp->listeners, model);
	
	return model;
}

Model* new_DiscreteParameterModel_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"dimension",
		"values"
	};
	
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* values = get_json_node(node, "values");
	int dim = get_json_node_value_int(node, "dimension", 0);
	size_t K = values->child_count;
	if (dim == 0) {
		dim = K;
	}
	if(dim < K){
		fprintf(stderr, "dimension attribute (%d) cannot be smaller than the number of values (%zu)\n", dim, K);
		exit(2);
	}
	unsigned* x = uivector(dim);
	int i = 0;
	while(i != dim) {
		for(int j = 0; j < values->child_count; j++){
			if(i == dim) break;
			json_node* child = values->children[j];
			x[j] = atoi((char*)child->value);
			i++;
		}
	}
	
	DiscreteParameter* dp = new_DiscreteParameter_with_values(x, dim);
	free(x);
	
	json_node* id_node = get_json_node(node, "id");
	Model* model = new_DiscreteParameterModel((char*)id_node->value, dp);
	
	return model;
}
