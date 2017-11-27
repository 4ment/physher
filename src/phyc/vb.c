//
//  vb.c
//  physher
//
//  Created by Mathieu Fourment on 24/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "vb.h"

#include <string.h>

#include "optimizer.h"
#include "matrix.h"
#include "gaussian.h"
#include "transforms.h"

#define MY_PI acos(-1.0)
#define LOG_TWO_PI (log(2.0)+log(MY_PI))

// Mean-field
double elbo_meanfield(struct variational_t* var){
	//const gsl_rng* r = var->rng;
	double elbo = 0;
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
	
	// Entropy: 0.5(1 + log(2 \pi s^2))
	double entropy = 0.5 * dim * (1.0 + LOG_TWO_PI);
	for (size_t i = 0; i < dim; i++) {
		entropy += Parameters_value(var->var_parameters, dim+i);
	}
	
	for (int i = 0; i < var->elbo_samples; i++) {
		double jacobian = 0.0;
		
		for (int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter* var_mean = Parameters_at(var->var_parameters, j);
			Parameter* var_sd = Parameters_at(var->var_parameters, dim+j);
			double zeta = rnorm() * exp(Parameter_value(var_sd)) + Parameter_value(var_mean);
			double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), &jacobian);
			Parameter_set_value(p, theta);
		}
		
		elbo += posterior->logP(posterior) + jacobian;
		
		if(isinf(elbo)){
			return elbo;
		}
	}
	return elbo/var->elbo_samples + entropy;
}

void grad_elbo_meanfield(struct variational_t* var, double* grads){
	//const gsl_rng* r = var->rng;
	Model* posterior = var->posterior;
	size_t dim = Parameters_count(var->parameters);
	double* eta = dvector(dim);
	double* zeta = dvector(dim);
	memset(grads, 0, sizeof(double)*dim*2);
	
	for (int i = 0; i < var->grad_samples; i++) {
		for ( int j = 0; j < dim; j++) {
			Parameter* var_mean = Parameters_at(var->var_parameters, j);
			Parameter* var_sd = Parameters_at(var->var_parameters, dim+j);
			eta[j] = rnorm();
			zeta[j] = eta[j] * exp(Parameter_value(var_sd)) + Parameter_value(var_mean);
			
			Parameter* p = Parameters_at(var->parameters, j);
			double theta = inverse_transform2(zeta[j], Parameter_lower(p), Parameter_upper(p));
			Parameter_set_value(p, theta);
		}
		
		for ( int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			double dlogP = Model_first_derivative(posterior, p, 0.001);
			const double gldit = grad_log_det_inverse_transform(zeta[j], Parameter_lower(p), Parameter_upper(p));
			double grad_mu = dlogP * grad_inverse_transform(zeta[j], Parameter_lower(p), Parameter_upper(p)) + gldit;
			grads[j] += grad_mu;
			grads[dim+j] += grad_mu * eta[j] * exp(Parameters_value(var->var_parameters, dim+j));
		}
	}
	
	for (int j = 0; j < dim; j++) {
		grads[j] /= var->grad_samples;
		grads[dim+j] = grads[dim+j]/var->grad_samples + 1.0;
		//printf("gradient %f %f\n", grads[i] ,grads[dim+i] );
	}
	//printf("\n\n");
	free(eta);
	free(zeta);
}

double elbo( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	return model->logP(model);
}

void grad_elbo( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	model->gradient(model, grad);
}

Model* new_Variational_from_json(json_node* node, Hashtable* hash){
	const char* posterior_string = get_json_node_value_string(node, "posterior");
	const char* var_string = get_json_node_value_string(node, "var");
	const char* dist_string = get_json_node_value_string(node, "distribution");
	const char* id = get_json_node_value_string(node, "id");
	
	struct variational_t* var = malloc(sizeof(struct variational_t));
	var->parameters = new_Parameters(1);
	get_parameters_references(node, hash, var->parameters);
	size_t dim = Parameters_count(var->parameters);
	
	var->posterior = Hashtable_get(hash, posterior_string+1);
	var->posterior->ref_count++;
	
	var->elbo_samples = get_json_node_value_size_t(node, "elbosamples", 100);
	var->grad_samples = get_json_node_value_size_t(node, "gradsamples", 1);
	
	StringBuffer* buffer = new_StringBuffer(10);
	if(strcasecmp(var_string, "meanfield") == 0){
		var->var_parameters = new_Parameters(dim*2); // mean + sd
		for(int i = 0; i < dim; i++){
			StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
			StringBuffer_append_string(buffer, ".mean");
			Parameters_move(var->var_parameters, new_Parameter(buffer->c, 1, NULL));
			StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
			StringBuffer_append_string(buffer, ".sd");
			Parameters_move(var->var_parameters, new_Parameter(buffer->c, 1, NULL));
		}
		var->elbofn = elbo_meanfield;
		var->grad_elbofn = grad_elbo_meanfield;
		var->f = elbo;
		var->grad_f = grad_elbo;
	}
	else if(strcasecmp(var_string, "fullrank") == 0){
		size_t n = dim+(dim*dim-dim)/2;
		var->var_parameters = new_Parameters(n);
		for(int i = 0; i < n; i++){
			Parameters_add(var->var_parameters, new_Parameter("", 1, NULL));
		}
	}
	free_StringBuffer(buffer);
	return new_VariationalModel(id, var);
}

static void _variational_model_free( Model *self ){
	if(self->ref_count == 1){
		printf("Free variational model %s\n", self->name);
		struct variational_t* var = (struct variational_t*)self->obj;
		var->posterior->free(var->posterior);
		free_Parameters(var->var_parameters);
		free_Parameters(var->parameters);
		free(var);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _variational_model_clone(Model* self, Hashtable *hash){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	return NULL;
}

static double _variational_model_logP(Model *self){
	struct variational_t* var = (struct variational_t*)self->obj;
	return var->elbofn(var);
}

static void _variational_model_gradient(Model *self, double* grad){
	struct variational_t* var = (struct variational_t*)self->obj;
	return var->grad_elbofn(var, grad);
}

static void _variational_model_get_free_parameters(Model* self, Parameters* parameters){
	struct variational_t* var = (struct variational_t*)self->obj;
	Parameters_add_free_parameters(parameters, var->var_parameters);
}

Model* new_VariationalModel(const char* name, struct variational_t* var){
	Model *model = new_Model("variational", name, var);
	model->free = _variational_model_free;
	model->clone = _variational_model_clone;
	model->logP = _variational_model_logP;
	model->gradient = _variational_model_gradient;
	model->get_free_parameters = _variational_model_get_free_parameters;
	return model;
}

void free_Variational(struct variational_t* var){
//	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
//		printf("%s %f\n", Parameters_name(var->var_parameters, i), Parameters_value(var->var_parameters, i));
//	}
	free_Parameters(var->var_parameters);
    free_Model(var->posterior);
    free_Parameters(var->parameters);
	free(var);
}
