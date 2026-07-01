/*
 *  model.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/15/10.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "model.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>

#include "mstring.h"
#include "parameters.h"

static void _dummy_update(Model* self, Model* model, Parameter* parameter, int index) {}
static double _logP(Model *model){return 0;}
static double _fulllogP(Model *model){return model->logP(model);}
static void _dummy_gradient(Model *model, Parameters* ps){
    fprintf(stderr, "gradient function not implemented for model %s\n", model->name);
    exit(2);
}
static void _dummy_reset(Model* m){}
static void _dummy_restore(Model* m){}
static void _dummy_accept(Model* m){}
static void _dummy_store(Model* m){}
static void _dummy_sample(Model* m) {
    fprintf(stderr, "Cannot sample from model %s\n", m->name);
    exit(2);
}
static void _dummy_rsample(Model* m) {
    fprintf(stderr, "Cannot rsample from model %s\n", m->name);
    exit(2);
}

static void _dummy_jsonize(Model* m, json_node* node){
    fprintf(stderr, "jsonize function not implemented for model %s\n", m->name);
    exit(2);
}

#pragma mark -


Model * new_Model( model_t type, const char *name, void *obj ){
	Model *model = (Model*)malloc(sizeof(Model));
	assert(model);
	model->name = String_clone(name);
	model->type = type;
	model->obj = obj;
	model->logP = _logP;
	model->full_logP = _fulllogP;
	model->gradient = _dummy_gradient;
	model->hessian = Model_hessian_fd;
	model->update = _dummy_update;
	model->free = free_Model;
	model->data = NULL;
	model->listeners = new_ListenerList(1);
	model->clone = NULL;
	model->ref_count = 1;
	model->reset = _dummy_reset;
	model->restore = _dummy_restore;
	model->store = _dummy_store;
	model->accept = _dummy_accept;
	model->stored = false;
	model->logP = 0;
	model->lp = 0;
	model->sample = _dummy_sample;
	model->rsample = _dummy_rsample;
	model->samplable = false;
	model->print = NULL;
    model->jsonize = _dummy_jsonize;
	model->epsilon = 0.0;
	model->get = NULL;
	model->set = NULL;
	model->parameters = new_Parameters(1);
	return model;
}

void free_Model( Model *model ){
	assert(model->ref_count >= 1);
	if(model->ref_count == 1){
		free(model->name);
		free_Parameters(model->parameters);
		model->listeners->free(model->listeners);
		free(model);
	}
	else{
		model->ref_count--;
	}

}

model_t check_model(const char* type){
	int len = sizeof(model_type_strings)/sizeof(model_type_strings[0]);
	for (int i = 0; i < len; i++) {
		if (strcasecmp(model_type_strings[i], type) == 0) {
			return i;
		}
	}
	return -1;
}

#pragma mark -

static void _free_ListenerList( ListenerList *listeners ){
	free(listeners->models);
	// listeners use weak references of Parameter objects. do not use free_Parameters or free_Parameter
	free_Parameters_weak(listeners->parameters);
	free(listeners);
}

static void _ListenerList_fire(ListenerList* listeners, Model* model,
                               Parameter* parameter, int index) {
    if (listeners->enabled == false) return;
    for (int i = 0; i < listeners->count; i++) {
        listeners->models[i]->update(listeners->models[i], model, parameter, index);
    }

    for (size_t i = 0; i < Parameters_count(listeners->parameters); i++) {
        Parameter* p = Parameters_at(listeners->parameters, i);
        p->update(p, model, parameter, index);
    }
}

static void _ListenerList_remove( ListenerList *listeners, Model* model ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		if ( listeners->models[i] == model ) {
			break;
		}
	}
	if ( i == listeners->count) {
		return;
	}
	i++;
	for ( ; i < listeners->count; i++ ) {
		listeners->models[i-1] = listeners->models[i];
	}
	listeners->models[listeners->count-1] = NULL;
	listeners->count--;
}

static void _ListenerList_remove_all( ListenerList *listeners ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		listeners->models[i] = NULL;
	}
	listeners->count = 0;
}

static void _ListenerList_add( ListenerList *listeners, Model *model ){
	if ( listeners->count == listeners->capacity) {
		listeners->capacity++;
		listeners->models = realloc(listeners->models, listeners->capacity*sizeof(Model*));
	}
	listeners->models[listeners->count] = model;
	listeners->count++;
}

static void _ListenerList_add_parameter(ListenerList* listeners, Parameter* p){
	// we use weak references so no ref increment
	Parameters_move(listeners->parameters, p);
}

ListenerList * new_ListenerList( const unsigned capacity ){
	ListenerList *listeners = (ListenerList*)malloc( sizeof(ListenerList));
	assert(listeners);
	listeners->capacity = (capacity == 0 ? 1 : capacity);
	listeners->count = 0;
	listeners->models = (Model**)malloc( listeners->capacity * sizeof(Model*));
	assert(listeners->models);
    listeners->parameters = new_Parameters(1);
    listeners->free = _free_ListenerList;
	listeners->add = _ListenerList_add;
	listeners->remove = _ListenerList_remove;
	listeners->removeAll = _ListenerList_remove_all;
	listeners->add_parameter = _ListenerList_add_parameter;
	listeners->fire = _ListenerList_fire;
	listeners->enabled = true;
	return listeners;
}

double Model_mixed_derivative( Model *model, Parameter* p1, Parameter* p2 ) {
	double eps = 0.001;

	double v1 = Parameter_value(p1);
	double e1 = eps * v1;
	double v2 = Parameter_value(p2);
	double e2 = eps * v2;

	// + +
	Parameter_set_value(p1, v1 + e1);
	Parameter_set_value(p2, v2 + e2);
	double pp = model->logP(model);

	// - -
	Parameter_set_value(p1, v1 - e1);
	Parameter_set_value(p2, v2 - e2);
	double mm = model->logP(model);

	// + -
	Parameter_set_value(p1, v1 + e1);
	Parameter_set_value(p2, v2 - e2);
	double pm = model->logP(model);

	// - +
	Parameter_set_value(p1, v1 - e1);
	Parameter_set_value(p2, v2 + e2);
	double mp = model->logP(model);

	Parameter_set_value(p1, v1);
	Parameter_set_value(p2, v2);

	return (pp + mm - pm - mp) / (4.0*e1*e2);
}

double Model_second_derivative( Model *model, Parameter* parameter, double* first, double eps ) {
	//double eps = 0.00001;

	double lnl = model->logP(model);

	double v = Parameter_value(parameter);
	double e = eps * v;

	// +
	Parameter_set_value(parameter, v + e);
	double p = model->logP(model);

	// -
	Parameter_set_value(parameter, v - e);
	double m = model->logP(model);

	Parameter_set_value(parameter, v);

	if (first != NULL) {
		*first = (p - m)/(2.0*e);
	}

	return (p + m -2*lnl)/(e*e);
}

// First derivative using central differences
// eps = 1.0e-8;
double Model_first_derivative( Model *model, Parameter* parameter, double eps ) {
	double v = Parameter_value(parameter);
	double e = eps * v;

	Parameter_set_value(parameter, v + e);
	double pp = model->logP(model);

	Parameter_set_value(parameter, v - e);
	double mm = model->logP(model);

	Parameter_set_value(parameter, v);

	return (pp - mm)/(2.0*e);
}

void Model_first_derivatives( Model *model, Parameter* parameter, double eps, double* grad ) {
	const double* values = Parameter_values(parameter);
	for(size_t i = 0; i < Parameter_size(parameter); i++){
		double v = values[i];
		double e = eps * v;

		Parameter_set_value_at(parameter, v + e, i);
		double pp = model->logP(model);

		Parameter_set_value_at(parameter, v - e, i);
		double mm = model->logP(model);

		Parameter_set_value_at(parameter, v, i);

		grad[i] = (pp - mm)/(2.0*e);
	}
}

// Map a flattened element index to its owning Parameter and local index.
static Parameter* _flat_parameter( const Parameters *parameters, size_t k, size_t *local ){
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		Parameter* p = Parameters_at((Parameters*)parameters, i);
		size_t sz = Parameter_size(p);
		if(k < sz){ *local = k; return p; }
		k -= sz;
	}
	*local = 0;
	return NULL;
}

void Model_hessian_fd( Model *model, const Parameters *parameters,
                       hessian_mode_t mode, double *out ){
	const double eps = model->epsilon > 0.0 ? model->epsilon : 1.0e-4;
	size_t dim = Parameters_size(parameters);
	double lnl = model->logP(model);

	for(size_t k = 0; k < dim; k++){
		size_t lk;
		Parameter* pk = _flat_parameter(parameters, k, &lk);
		double vk = Parameter_value_at(pk, lk);
		double ek = eps * (vk != 0.0 ? fabs(vk) : 1.0);

		// diagonal: central second difference
		Parameter_set_value_at(pk, vk + ek, lk);
		double pp = model->logP(model);
		Parameter_set_value_at(pk, vk - ek, lk);
		double mm = model->logP(model);
		Parameter_set_value_at(pk, vk, lk);

		double dkk = (pp + mm - 2.0*lnl)/(ek*ek);
		if(mode == HESSIAN_DIAGONAL){
			out[k] = dkk;
			continue;
		}
		out[k*dim + k] = dkk;

		// off-diagonal: central mixed difference, mirror to lower triangle
		for(size_t l = k+1; l < dim; l++){
			size_t ll;
			Parameter* pl = _flat_parameter(parameters, l, &ll);
			double vl = Parameter_value_at(pl, ll);
			double el = eps * (vl != 0.0 ? fabs(vl) : 1.0);

			Parameter_set_value_at(pk, vk + ek, lk);
			Parameter_set_value_at(pl, vl + el, ll);
			double ppe = model->logP(model);

			Parameter_set_value_at(pl, vl - el, ll);
			double pme = model->logP(model);

			Parameter_set_value_at(pk, vk - ek, lk);
			double mme = model->logP(model);

			Parameter_set_value_at(pl, vl + el, ll);
			double mpe = model->logP(model);

			Parameter_set_value_at(pk, vk, lk);
			Parameter_set_value_at(pl, vl, ll);

			double dkl = (ppe + mme - pme - mpe)/(4.0*ek*el);
			out[k*dim + l] = dkl;
			out[l*dim + k] = dkl;
		}
	}
}

#pragma mark -
#pragma mark CatParameterModel

static void _cat_parameter_model_get(Model* self, double* values){
	Parameters* parameters = (Parameters*)self->obj;
	size_t counter = 0;
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		Parameter* p = Parameters_at(parameters, i);
		memcpy(values+counter, Parameter_values(p), sizeof(double)*Parameter_size(p));
		counter += Parameter_size(p);
	}
}

static void _cat_parameter_model_set(Model* self, const double* values){
	Parameters* parameters = (Parameters*)self->obj;
	size_t counter = 0;
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		Parameter* p = Parameters_at(parameters, i);
		Parameter_set_values(p, values+counter);
		counter += Parameter_size(p);
	}
}

static void _cat_parameter_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free cat parameter model %s\n", self->name);
		Parameters* parameters = (Parameters*)self->obj;
		free_Parameters(parameters);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

Model * new_CatParameterModel( const char* name, Parameters *parameters ){
	Model *model = new_Model(MODEL_PARAMETERS, name, parameters);

	model->free = _cat_parameter_model_free;
	model->get = _cat_parameter_model_get;
	model->set = _cat_parameter_model_set;
	return model;
}

Parameter* new_CatParameter_from_json(json_node*node, Hashtable*hash){
	static const json_field schema[] = {
	    {"parameters", JSON_OPTIONAL, JSON_ANY},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

    char* id = get_json_node_value_string(node, "id");
	json_node* parameters_node = get_json_node(node, "parameters");
	Parameters* parameters = new_Parameters(10);
	grab_parameters(parameters_node, hash, parameters);

	Model* model = new_CatParameterModel(id, parameters);
	Parameter* p = new_ParameterModel(id, NULL, Parameters_count(parameters), NULL, model);
	return p;
}
