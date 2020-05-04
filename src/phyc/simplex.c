//
//  simplex.c
//  physher
//
//  Created by Mathieu Fourment on 7/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "simplex.h"

#include <string.h>

#include "matrix.h"

// free to constrained
void _inverse_transform(const Parameters* parameters, double* values){
	size_t N = Parameters_count(parameters); // n=K-1
		values[N] = 1.0;
		for (int i = 0; i < N; i++) {
			values[N] += Parameters_value(parameters, i);
		}
		 values[N] = 1.0/values[N];
		for (int i = 0; i < N; i++) {
			values[i] = Parameters_value(parameters, i)*values[N];
		}
}

// constrained to free
void _transform(const double* values, Parameters* parameters){
	size_t N = Parameters_count(parameters); // N=K-1
		for (int i = 0; i < N; i++) {
			Parameters_set_value(parameters, i, values[i]/values[N]);
		}
}

// stick: free to constrained
void _inverse_transform_stick(const Parameters* parameters, double* values){
	size_t N = Parameters_count(parameters); // n=K-1
	double stick = 1.0;
	int k = 0;
	for(; k < N; k++){
		values[k] = stick * Parameters_value(parameters, k);
		stick -= values[k];
	}
	values[k] = stick;
}

// stick: constrained to free
void _transform_stick(const double* values, Parameters* parameters){
	size_t N = Parameters_count(parameters); // N=K-1
	double stick = 1.0;
	for (int i = 0; i < N; i++) {
		Parameters_set_value(parameters, i, values[i]/stick);
		stick -= values[i];
	}
}


void free_Simplex(Simplex* simplex){
	free_Parameters(simplex->parameters);
	free(simplex->values);
	free(simplex);
}

Simplex* clone_Simplex(const Simplex* simplex){
	Simplex* clone = (Simplex*)malloc(sizeof(Simplex));
	clone->K = simplex->K;
	if(simplex->parameters != NULL){
		clone->parameters = clone_Parameters(simplex->parameters );
	}
	clone->values = clone_dvector(simplex->values, simplex->K);
	clone->need_update = simplex->need_update;
	clone->get_value = simplex->get_value;
	clone->get_values = simplex->get_values;
	clone->set_parameter_value = simplex->set_parameter_value;
	clone->set_values = simplex->set_values;
	return clone;
}

void set_values(Simplex* simplex, const double* values){
	memcpy(simplex->values, values, sizeof(double)*simplex->K);
	_transform_stick(simplex->values, simplex->parameters);
}

void set_parameter_value(Simplex* simplex, int index, double value){
	Parameters_set_value(simplex->parameters, index, value);
	simplex->need_update  = true;
}

const double* get_values(Simplex* simplex){
	if (simplex->need_update) {
		_inverse_transform_stick(simplex->parameters, simplex->values);
		simplex->need_update  = false;
	}
	return simplex->values;
}

double get_value(Simplex* simplex, int i){
	if (simplex->need_update) {
		_inverse_transform_stick(simplex->parameters, simplex->values);
		simplex->need_update  = false;
	}
	return simplex->values[i];
}

// x is a simplex of dimension K
//Simplex* new_Simplex_with_parameters(const Parameters *x){
//	Simplex* simplex = (Simplex*)malloc(sizeof(Simplex));
//	simplex->K = Parameters_count(x);
//	simplex->cparameters = new_Parameters(simplex->K);
//	for(int i = 0; i < simplex->K; i++){
//		Parameters_move(simplex->cparameters, Parameters_at(x, i));
//	}
//	simplex->values = dvector(simplex->K);
//	simplex->get_values = get_values;
//	simplex->need_update = true;
//	return simplex;
//}

// x is of dimension K
Simplex* new_Simplex_with_values(const char* name, const double *x, size_t K){
	Simplex* simplex = (Simplex*)malloc(sizeof(Simplex));
	simplex->K = K;
	simplex->parameters = new_Parameters_with_name(name, K-1);
	simplex->values = clone_dvector(x, K);
    StringBuffer* buffer = new_StringBuffer(10);
	size_t N = K-1;
	double stick = 1;
    double sum = 0;
    for(int i = 0; i < K; i++){
        sum += x[i];
    }
    if(fabs(sum - 1.0) > 0.0001) {
        fprintf(stderr, "simplex %s does not sum to 1.0 (+- 0.0001)\n", name);
        exit(1);
    }

	for(int i = 0; i < N; i++){
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "%s.phi%d", name, i);
        Constraint* constraint = new_Constraint(0, 1);
        Constraint_set_fupper(constraint, 0.99);
		double phi = simplex->values[i]/stick;
		Parameters_move(simplex->parameters, new_Parameter(buffer->c, phi, constraint));
		Parameters_at(simplex->parameters, i)->id = i;
		stick -= simplex->values[i];
	}
    free_StringBuffer(buffer);
	simplex->get_values = get_values;
	simplex->get_value = get_value;
	simplex->set_values = set_values;
	simplex->set_parameter_value = set_parameter_value;
	simplex->need_update = true;
	return simplex;
}

Simplex* new_Simplex(const char* name, size_t K){
	Simplex* simplex = (Simplex*)malloc(sizeof(Simplex));
	simplex->K = K;
	simplex->parameters = new_Parameters(K-1);
    Parameters_set_name2(simplex->parameters, name);
	simplex->values = dvector(K);
	size_t N = K-1;
    StringBuffer* buffer = new_StringBuffer(10);
	double p = 1.0/K;
	for(int i = 0; i < K; i++){
		simplex->values[i] = p;
	}
	double stick = 1;
	for(int i = 0; i < N; i++){
        Constraint* constraint = new_Constraint(0, 1);
        Constraint_set_fupper(constraint, 0.999);
        Constraint_set_flower(constraint, 0.001);
		double phi = simplex->values[i]/stick;
        Parameters_move(simplex->parameters, new_Parameter(buffer->c, phi, constraint));
		Parameters_at(simplex->parameters, i)->id = i;
		stick -= simplex->values[i];
	}
    free_StringBuffer(buffer);
	simplex->get_values = get_values;
	simplex->get_value = get_value;
	simplex->set_values = set_values;
	simplex->set_parameter_value = set_parameter_value;
	simplex->need_update = true;
	return simplex;
}

static void _simplex_model_handle_change( Model *self, Model *model, int index ){
//	printf("update %d\n", index);
	Simplex* simplex = (Simplex*)self->obj;
	simplex->need_update = true;
	self->listeners->fire( self->listeners, self, index );
}

static void _simplex_model_handle_restore( Model *self, Model *model, int index ){
	Simplex* simplex = (Simplex*)self->obj;
	simplex->need_update = true;
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _simplex_model_store(Model* self){
	Simplex* simplex = self->obj;
	if ( simplex->parameters != NULL ) {
		for (int i = 0; i < Parameters_count(simplex->parameters); i++) {
			Parameter_store(Parameters_at(simplex->parameters, i));
		}
	}
}

static void _simplex_model_restore(Model* self){
	Simplex* simplex = self->obj;
	if ( simplex->parameters != NULL ) {
		bool changed = false;
		Parameter*p = NULL;
		for (int i = 0; i < Parameters_count(simplex->parameters); i++) {
			p = Parameters_at(simplex->parameters, i);
			if (Parameter_changed(p)) {
				changed = true;
				Parameter_restore_quietly(p);
			}
		}
		if (changed) {
			p->listeners->fire_restore(p->listeners, NULL, p->id);
		}
	}
}

static void _simplex_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free simplex model %s\n", self->name);
		Simplex* simplex = (Simplex*)self->obj;
		free_Simplex(simplex);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _simplex_model_clone(Model* self, Hashtable *hash){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	Simplex* simplex = (Simplex*)self->obj;
	Simplex* sclone = clone_Simplex(simplex);
	if ( sclone->parameters != NULL ) {
		for ( int i = 0; i < Parameters_count(sclone->parameters); i++ ) {
			Hashtable_add(hash, Parameters_name(sclone->parameters, i), Parameters_at(sclone->parameters, i));
		}
	}
	Model *clone = new_SimplexModel(self->name, sclone);
	Hashtable_add(hash, clone->name, clone);
	return clone;
}

// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SimplexModel( const char* name, Simplex *simplex ){
	Model *model = new_Model(MODEL_SIMPLEX, name, simplex);
	int i = 0;
	
	if ( simplex->parameters != NULL ) {
		for ( i = 0; i < Parameters_count(simplex->parameters); i++ ) {
			Parameters_at(simplex->parameters, i)->listeners->add( Parameters_at(simplex->parameters, i)->listeners, model );
		}
	}
	
	model->update = _simplex_model_handle_change;
	model->handle_restore = _simplex_model_handle_restore;
	model->store = _simplex_model_store;
	model->restore = _simplex_model_restore;
	model->free = _simplex_model_free;
	model->clone = _simplex_model_clone;
	return model;
}

Model* new_SimplexModel_from_json(json_node*node, Hashtable*hash){
	char* allowed[] = {
		"dimension",
		"values"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
    char* id = get_json_node_value_string(node, "id");
	json_node* values = get_json_node(node, "values");
	Simplex* simplex = NULL;
	
	if(values != NULL){
		double* x = dvector(values->child_count);
        double sum = 0;
		for (int i = 0; i < values->child_count; i++) {
			json_node* child = values->children[i];
			if(child->node_type != MJSON_PRIMITIVE){
				fprintf(stderr, "sadf\n");
				exit(1);
			}
			x[i] = atof((char*)child->value);
            sum += x[i];
            if(x[i] <= 0 || x[i] >= 1){
                fprintf(stderr, "Each value in a simplex should be between 0 and 1 (id: %s)\n", id);
                exit(1);
            }
		}
        if(fabs(sum - 1.0) < 0.0001) {
            for(int i = 0; i < values->child_count; i++){
                 x[i] /= sum;
            }
        }
        else{
            fprintf(stderr, "simplex %s does not sum to 1.0 (+- 0.0001)\n", id);
            exit(1);
        }

		simplex = new_Simplex_with_values(id, x, values->child_count);
		free(x);
	}
	else{
		json_node* dimension = get_json_node(node, "dimension");
		simplex = new_Simplex(id, atoi((char*)dimension->value));
	}
	Model* model = new_SimplexModel(id, simplex);
	return model;
}

