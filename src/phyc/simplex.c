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
#include "transforms.h"

// free to constrained
void _inverse_transform(const Parameter* parameter, double* values){
	size_t N = Parameter_size(parameter); // n=K-1
	const double* p = Parameter_values(parameter);
	values[N] = 1.0;
	for (size_t i = 0; i < N; i++) {
		values[N] += p[i];
	}
	values[N] = 1.0/values[N];
	for (size_t i = 0; i < N; i++) {
		values[i] = p[i]*values[N];
	}
}

// constrained to free
void _transform(const double* values, Parameter* parameter){
	size_t N = Parameter_size(parameter); // N=K-1
	for (size_t i = 0; i < N; i++) {
		Parameter_set_value_at(parameter, values[i]/values[N], i);
	}
}

// stick: free to constrained
void _inverse_transform_stick(const Parameter* parameter, double* values){
	size_t N = Parameter_size(parameter); // n=K-1
	const double* p = Parameter_values(parameter);
	double stick = 1.0;
	size_t k = 0;
	for(; k < N; k++){
		values[k] = stick * p[k];
		stick -= values[k];
	}
	values[k] = stick;
}

// stick: constrained to free
void _transform_stick(const double* values, Parameter* parameter){
	size_t N = Parameter_size(parameter); // N=K-1
	double stick = 1.0;
	for (int i = 0; i < N; i++) {
		Parameter_set_value_at(parameter, values[i]/stick, i);
		stick -= values[i];
	}
}

void _inverse_transform_stick_stan(const Parameter* parameter, double* values){
	size_t N = Parameter_size(parameter); // N=K-1
	const double* p = Parameter_values(parameter);
	double stick = 1.0;
	for(size_t k = 0; k < N; k++){
		double z = inverse_logit(p[k] - log(N-k));
		values[k] = stick * z;
		stick -= values[k];
	}
	values[N] = stick;
}

void _transform_stick_stan(const double* values, Parameter* parameter){
	size_t N = Parameter_size(parameter); // N=K-1
	double sum = 0;
	for (size_t i = 0; i < N; i++) {
		double z = values[i] / (1.0 - sum);
		Parameter_set_value_at(parameter, logit(z) + log(N-i), i);
		sum += values[i];
	}
}


void free_Simplex(Simplex* simplex){
	free_Parameter(simplex->parameter);
	free(simplex->values);
	free(simplex->stored_values);
	free(simplex);
}

Simplex* clone_Simplex(const Simplex* simplex){
	Simplex* clone = (Simplex*)malloc(sizeof(Simplex));
	clone->K = simplex->K;
	clone->parameter = clone_Parameter(simplex->parameter);
	clone->values = clone_dvector(simplex->values, simplex->K);
	clone->stored_values = clone_dvector(simplex->stored_values, simplex->K);
	clone->need_update = simplex->need_update;
	clone->get_value = simplex->get_value;
	clone->get_values = simplex->get_values;
	clone->set_parameter_value = simplex->set_parameter_value;
	clone->set_values = simplex->set_values;
	return clone;
}

void set_values(Simplex* simplex, const double* values){
	memcpy(simplex->values, values, sizeof(double)*simplex->K);
	_transform_stick(simplex->values, simplex->parameter);
}

void set_parameter_value(Simplex* simplex, int index, double value){
	Parameter_set_value_at(simplex->parameter, index, value);
	simplex->need_update  = true;
}

const double* get_values(Simplex* simplex){
	if (simplex->need_update) {
		_inverse_transform_stick(simplex->parameter, simplex->values);
		simplex->need_update  = false;
	}
	return simplex->values;
}

double get_value(Simplex* simplex, int i){
	if (simplex->need_update) {
		_inverse_transform_stick(simplex->parameter, simplex->values);
		simplex->need_update  = false;
	}
	return simplex->values[i];
}

void set_values_stan(Simplex* simplex, const double* values){
	memcpy(simplex->values, values, sizeof(double)*simplex->K);
	_transform_stick_stan(simplex->values, simplex->parameter);
}

const double* get_values_stan(Simplex* simplex){
	if (simplex->need_update) {
		_inverse_transform_stick_stan(simplex->parameter, simplex->values);
		simplex->need_update  = false;
	}
	return simplex->values;
}

double get_value_stan(Simplex* simplex, int i){
	if (simplex->need_update) {
		_inverse_transform_stick_stan(simplex->parameter, simplex->values);
		simplex->need_update  = false;
	}
	return simplex->values[i];
}

// dp/dphi_index
void _simplex_gradient(Simplex* simplex, size_t index, double* gradient){
	const double* freqs = simplex->get_values(simplex);
	memset(gradient, 0, sizeof(double)*simplex->K);
	
	double p = Parameter_value_at(simplex->parameter, index);
	gradient[index] = freqs[index]/p;
	double cum = gradient[index];
	for (size_t i = index+1; i <simplex->K-1; i++) {
		gradient[i] = -cum*Parameter_value_at(simplex->parameter, i);
		cum += gradient[i];
	}
	gradient[simplex->K-1] = -cum;
}

// dp/dphi_index
void _simplex_gradient_stan(Simplex* simplex, size_t index, double* gradient){
	const double* freqs = simplex->get_values(simplex);
	memset(gradient, 0, sizeof(double)*simplex->K);
	
	double p = Parameter_value_at(simplex->parameter, index);
	double stickRemaining = 1.0;
	for(size_t i = 0; i < index; i++){
		stickRemaining -= freqs[i];
	}
	gradient[index] = stickRemaining * grad_inverse_logit(p - log(simplex->K-1-index));
	double cum = gradient[index];
	for (size_t i = index+1; i <simplex->K-1; i++) {
		gradient[i] = -cum*inverse_logit(Parameter_value_at(simplex->parameter, i) - log(simplex->K-1-i));
		cum += gradient[i];
	}
	gradient[simplex->K-1] = -cum;
}

void Simplex_use_stan_transform(Simplex* simplex, bool use_stan){
	if(use_stan){
		simplex->get_values = get_values_stan;
		simplex->get_value = get_value_stan;
		simplex->set_values = set_values_stan;
		simplex->gradient = _simplex_gradient_stan;
	}
	else{
		simplex->get_values = get_values;
		simplex->get_value = get_value;
		simplex->set_values = set_values;
		simplex->gradient = _simplex_gradient;
	}
	simplex->need_update = true;
}

// Simplex uses the reparameterization of Stan
// If unconstrained parameters are all equal to zero then constrained values are all equal
Simplex* new_Simplex_with_values(const char* name, const double *x, size_t K){
	size_t N = K-1;
    double sum = 0;
    for(int i = 0; i < K; i++){
        sum += x[i];
    }
    if(fabs(sum - 1.0) > 0.0001) {
        fprintf(stderr, "simplex %s does not sum to 1.0 (+- 0.0001) %f\n", name, sum);
        exit(1);
    }

	Parameter* parameter = new_Parameter2(name, x, K-1, new_Constraint(-INFINITY, INFINITY));
	StringBuffer* buffer = new_StringBuffer(10);	
    StringBuffer_append_format(buffer, "%s.phi", name);
	Parameter_set_name(parameter, buffer->c);
    free_StringBuffer(buffer);

	Simplex* simplex = new_Simplex_with_parameter(name, parameter);
	simplex->set_values(simplex, x);
	return simplex;
}

Simplex* new_Simplex_with_parameter(const char* name, Parameter* parameter){
	Simplex* simplex = (Simplex*)malloc(sizeof(Simplex));
	simplex->K = Parameter_size(parameter) + 1;
	simplex->parameter = parameter;
	simplex->values = dvector(simplex->K);
	simplex->stored_values = dvector(simplex->K);
	simplex->get_values = get_values_stan;
	simplex->get_value = get_value_stan;
	simplex->set_values = set_values_stan;
	simplex->set_parameter_value = set_parameter_value;
	simplex->gradient = _simplex_gradient_stan;
	simplex->need_update = true;
	return simplex;
}

Simplex* new_Simplex(const char* name, size_t K){
	double* values = dvector(K);
	for(size_t i = 0; i < K; i++){
		values[i] = 1.0/K;
	}
	Simplex* simplex = new_Simplex_with_values(name, values, K);
	free(values);
	return simplex;
}

static void _simplex_model_handle_change( Model *self, Model *model, Parameter* parameter, int index ){
	Simplex* simplex = (Simplex*)self->obj;
	simplex->need_update = true;
	self->listeners->fire( self->listeners, self, parameter, index );
	// printf("_simplex_model_handle_change %s  %d\n", self->name, self->listeners->count);
}

static void _simplex_model_handle_restore( Model *self, Model *model, int index ){
	Simplex* simplex = (Simplex*)self->obj;
	simplex->need_update = true;
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _simplex_model_store(Model* self){
	Simplex* simplex = self->obj;
	Parameter_store(simplex->parameter);
	memcpy(simplex->stored_values, simplex->values, sizeof(double)*simplex->K);
}

static void _simplex_model_restore(Model* self){
	Simplex* simplex = self->obj;
	Parameter_restore(simplex->parameter);
	memcpy(simplex->values, simplex->stored_values, sizeof(double)*simplex->K);
	simplex->need_update = false;
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
	if ( sclone->parameter != NULL ) {
		Hashtable_add(hash, Parameter_name(sclone->parameter), sclone->parameter);
	}
	Model *clone = new_SimplexModel(self->name, sclone);
	Hashtable_add(hash, clone->name, clone);
	return clone;
}

static void _simplex_model_get(Model* self, double* values){
	Simplex* simplex = (Simplex*)self->obj;
	const double* constrained = simplex->get_values(simplex);
	memcpy(values, constrained, sizeof(double)*simplex->K);
}

static void _simplex_model_set(Model* self, const double* values){
	Simplex* simplex = (Simplex*)self->obj;
	simplex->set_values(simplex, values);
}

// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SimplexModel( const char* name, Simplex *simplex ){
	Model *model = new_Model(MODEL_SIMPLEX, name, simplex);

	if ( simplex->parameter != NULL ) {
		simplex->parameter->listeners->add(simplex->parameter->listeners, model );
	}
	
	model->update = _simplex_model_handle_change;
	model->handle_restore = _simplex_model_handle_restore;
	model->store = _simplex_model_store;
	model->restore = _simplex_model_restore;
	model->free = _simplex_model_free;
	model->clone = _simplex_model_clone;
	model->get = _simplex_model_get;
	model->set = _simplex_model_set;
	Parameters_add(model->parameters, simplex->parameter);
	return model;
}

Parameter* new_SimplexParameter_from_json(json_node*node, Hashtable*hash){
	Model* msimplex = new_SimplexModel_from_json(node, hash);
	Simplex* simplex = (Simplex*)msimplex->obj;
	Constraint* cnstr = new_Constraint(0, 1);
	char* id = get_json_node_value_string(node, "id");
	Parameter* parameter =  new_ParameterModel(id, NULL, simplex->K, cnstr, msimplex);
	return parameter;
}

Model* new_SimplexModel_from_json(json_node*node, Hashtable*hash){
	char* allowed[] = {
		"centered",
		"dimension",
		"parameter",
		"values",
		"x"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
    char* id = get_json_node_value_string(node, "id");
	json_node* values = get_json_node(node, "x");
	if(values == NULL){
		values = get_json_node(node, "values");
	}
	json_node* parameter_node = get_json_node(node, "parameter");
	json_node* dimension_node = get_json_node(node, "dimension");
	bool centered = get_json_node_value_bool(node, "centered", true);
	Simplex* simplex = NULL;
	double* x = NULL;
	size_t dimension = 0;
	

	if(values != NULL){
		dimension = values->child_count;
		x = dvector(dimension);
        double sum = 0;
		for (size_t i = 0; i < dimension; i++) {
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
            for(size_t i = 0; i < dimension; i++){
                 x[i] /= sum;
            }
        }
        else{
            fprintf(stderr, "simplex %s does not sum to 1.0 (+- 0.0001)\n", id);
            exit(1);
        }

		simplex = new_Simplex_with_values(id, x, dimension);
	}
	else if(dimension_node != NULL){
		dimension = get_json_node_value_size_t(node, "dimension", 0);
		if(dimension < 2){
            fprintf(stderr, "the dimension (%lu) of simplex %s should be greater than 1\n", dimension, id);
            exit(1);
		}
		simplex = new_Simplex(id, dimension);
	}

	if(parameter_node != NULL){
		Parameter* p = NULL;
		if(parameter_node->node_type == MJSON_STRING){
			char* ref = parameter_node->value;
			// it was defined soemwhere else
			if(safe_is_reference(ref, id)){
				p = safe_get_reference_parameter(ref, hash, id);
				p = Hashtable_get(hash, ref+1);
				p->refCount++;
				// p should have no constraints (-inf, inf)
				simplex = new_Simplex_with_parameter(id, p);
			}
			// make it available
			else{
				p = simplex->parameter;
				Parameter_set_name(simplex->parameter, ref);
				Hashtable_add(hash, Parameter_name(p), p);
			}
		}
		else{
			p = new_Parameter_from_json(parameter_node, hash);
			Hashtable_add(hash, Parameter_name(p), p);
			// p should have no constraints (-inf, inf)
			simplex = new_Simplex_with_parameter(id, p);
		}
	}
	// else{
	// 	fprintf(stderr, "simplex %s requires a `dimension' or `values' attribute\n", id);
    //     exit(1);
	// }
	if(!centered){
		Simplex_use_stan_transform(simplex, false);
	}
	if(values != NULL){
		free(x);
	}

	Model* model = new_SimplexModel(id, simplex);
	return model;
}

