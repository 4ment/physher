//
//  vb.c
//  physher
//
//  Created by Mathieu Fourment on 24/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "vb.h"

#include <string.h>
#include <strings.h>

#include "optimizer.h"
#include "tree.h"
#include "utils.h"
#include "transforms.h"

#include "gamvi.h"

#include "klpq.h"
#include "klqp.h"


Model* new_Variational_from_json(json_node* node, Hashtable* hash){
	const char* posterior_string = get_json_node_value_string(node, "posterior");
	const char* var_string = get_json_node_value_string(node, "var");
	const char* dist_string = get_json_node_value_string(node, "distribution");
	const char* id = get_json_node_value_string(node, "id");
	const char* filename = get_json_node_value_string(node, "log");
	json_node* parameters_node = get_json_node(node, "parameters");
	json_node* tree_node = get_json_node(node, "tree");
	
	variational_t* var = malloc(sizeof(variational_t));
	var->file = NULL;
	var->sample = NULL;
	var->sample_some = NULL;
	var->parameters = new_Parameters(1);
	var->var_parameters = NULL;
	
	if(parameters_node){
		get_parameters_references(node, hash, var->parameters);
	}
	if (tree_node != NULL) {
		char* tree_string = get_json_node_value_string(node, "tree");
		Model* mtree = Hashtable_get(hash, tree_string+1);
		Tree* tree = mtree->obj;
		Node** nodes = Tree_get_nodes(tree, POSTORDER);
	
		for(int i = 0; i < Tree_node_count(tree)-2; i++){
			Parameters_add(var->parameters, nodes[i]->distance);
		}
	}
	
	size_t dim = Parameters_count(var->parameters);
	
	var->posterior = Hashtable_get(hash, posterior_string+1);
	var->posterior->ref_count++;
	
	var->elbo_samples = get_json_node_value_size_t(node, "elbosamples", 100);
	var->grad_samples = get_json_node_value_size_t(node, "gradsamples", 1);
	
	json_node* var_parameters_node = get_json_node(node, "var_parameters");
	// Variational model parameters are provided
	if (var_parameters_node->node_type == MJSON_ARRAY) {
		var->var_parameters = new_Parameters(var_parameters_node->child_count);
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			json_node* var_p_node = var_parameters_node->children[i];
			Parameter* var_p = new_Parameter_from_json(var_p_node, hash);
			Parameters_move(var->var_parameters, var_p);
		}
	}
	
	StringBuffer* buffer = new_StringBuffer(10);
	if(strcasecmp(var_string, "meanfield") == 0){
//		printf("Creating meanfield variational model\n");
		if (strcasecmp(dist_string, "normal") == 0) {
			if(var->var_parameters == NULL){
				var->var_parameters = new_Parameters(dim*2); // mean + sd
				for(int i = 0; i < dim; i++){
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".mean");
					Parameter* p = new_Parameter(buffer->c, 0, NULL);
					p->id = i;
					Parameters_move(var->var_parameters, p);
				}
				for(int i = 0; i < dim; i++){
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".sd");
					Parameter* p = new_Parameter(buffer->c, 0, NULL);
					p->id = i+dim;
					Parameters_move(var->var_parameters, p);
				}
			}
			bool backward = get_json_node_value_bool(node, "backward", true); // true is KL(Q||P)
			if(backward){
				var->elbofn = klqp_meanfield_normal_elbo;
				var->grad_elbofn = klqp_meanfield_normal_grad_elbo;
			}
			else{
				var->elbofn = klpq_normal_meanfield;
				var->grad_elbofn = grad_klpq_normal_meanfield;
			}
			var->sample = klqp_meanfield_normal_sample;
			var->sample_some = klqp_meanfield_normal_sample_some;
			var->finalize = klqp_meanfield_normal_finalize;
			var->logP = klqp_meanfield_normal_logP;
			var->parameters_logP = klqp_meanfield_normal_logP_parameters;
		}
		else if (strcasecmp(dist_string, "gamma") == 0) {
			if(var->var_parameters == NULL){
				var->var_parameters = new_Parameters(dim*2); // alpha + beta
				for(int i = 0; i < dim; i++){
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".alpha");
					Parameter* p = new_Parameter(buffer->c, 0, NULL);
					p->id = i;
					Parameters_move(var->var_parameters, p);
				}
				for(int i = 0; i < dim; i++){
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".beta");
					Parameter* p = new_Parameter(buffer->c, 0, NULL);
					p->id = i+dim;
					Parameters_move(var->var_parameters, p);
				}
			}
			var->elbofn = elbo_gamma_meanfield;
			var->grad_elbofn = grad_elbo_gamma_meanfield;
			
			var->sample = variational_sample_gamma_meanfield;
			var->sample_some = variational_sample_some_gamma_meanfield;
			var->finalize = NULL; //TODO: implement function
			var->logP = variational_gamma_meanfield_logP;
			var->parameters_logP = variational_gamma_meanfield_parameters_logP;
		}
		else{
			fprintf(stderr, "distribution %s not recognized\n", dist_string);
			exit(1);
		}
	}
	else if(strcasecmp(var_string, "fullrank") == 0){
//		printf("Creating fullrank variational model\n");
		if(var->var_parameters == NULL){
			size_t n = 2*dim+(dim*dim-dim)/2; // mus, vars, covs
			var->var_parameters = new_Parameters(n);
			
			// means
			for(int i = 0; i < dim; i++){
				StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
				StringBuffer_append_string(buffer, ".mean");
				Parameter* p = new_Parameter(buffer->c, 1, NULL);
				Parameters_move(var->var_parameters, p);
			}
			
			// sparse covariance matrix
			if(tree_node != NULL){
				char* tree_string = get_json_node_value_string(node, "tree");
				Model* mtree = Hashtable_get(hash, tree_string+1);
				Tree* tree = mtree->obj;
				Node** nodes = Tree_get_nodes(tree, POSTORDER);
				for (int i = 0; i < Tree_node_count(tree)-2; i++) {
					Node* nodei = nodes[i];
					for (int j = 0; j < i; j++) {
						Node* nodej = nodes[j];
						StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
						StringBuffer_append_strings(buffer, 2, ".", Parameters_name(var->parameters, j));
						StringBuffer_append_string(buffer, ".cov");
						Parameter* p = new_Parameter(buffer->c, 1, NULL);
						if (Node_parent(nodei) != nodej && Node_parent(nodej) != nodei && Node_sibling(nodei) != nodej) {
							p->estimate = false;
							p->value = 0;
						}
	//					p->estimate = false;
	//					p->value = 0;
						Parameters_move(var->var_parameters, p);
					}
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".var");
					Parameter* p = new_Parameter(buffer->c, 1, NULL);
					Parameters_move(var->var_parameters, p);
				}
			}
			// full covariance matrix
			else{
				size_t row = dim;
				for(int i = 0; i < dim; i++){
					for(int j = 0; j < i; j++){
						StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
						StringBuffer_append_strings(buffer, 2, ".", Parameters_name(var->parameters, j));
						StringBuffer_append_string(buffer, ".cov");
						Parameter* p = new_Parameter(buffer->c, 1, NULL);
						Parameters_move(var->var_parameters, p);
						row++;
					}
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".var");
					Parameter* p = new_Parameter(buffer->c, 1, NULL);
					Parameters_move(var->var_parameters, p);
		//			p->estimate = false;
					row++;
				}
			}
		}
		
		var->elbofn = klqp_fullrank_normal_elbo;
		var->grad_elbofn = klqp_fullrank_normal_grad_elbo;
		var->sample = klqp_fullrank_normal_sample;
		var->logP = klqp_fullrank_normal_logP;
        var->parameters_logP = klqp_fullrank_normal_logP_parameters;
	}
	
	// variational parameters were not provided 
	if (var_parameters_node->node_type == MJSON_STRING) {
		Parameters_set_name2(var->var_parameters, get_json_node_value_string(node, "var_parameters"));
		Hashtable_add(hash, Parameters_name2(var->var_parameters), var->var_parameters);
	}
	
	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
		Hashtable_add(hash, Parameters_name(var->var_parameters, i), Parameters_at(var->var_parameters, i));
	}
	
	var->initialized = false;
	var->ready_to_sample = false;
	var->iter = 0;
	if (filename != NULL) {
		var->file = fopen(filename, "w");
		fprintf(var->file, "iteration");
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%s", Parameters_name(var->var_parameters, i));
		}
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%s.grad", Parameters_name(var->var_parameters, i));
		}
		fprintf(var->file, "\n");
	}
	
	free_StringBuffer(buffer);
	var->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return new_VariationalModel(id, var);
}


void free_Variational(variational_t* var){
	var->posterior->free(var->posterior);
	free_Parameters(var->var_parameters);
	free_Parameters(var->parameters);
	if (var->file != NULL) fclose(var->file);
	free(var);
}

static void _variational_model_free( Model *self ){
	if(self->ref_count == 1){
//		printf("Free variational model %s\n", self->name);
		variational_t* var = (variational_t*)self->obj;
		free_Variational(var);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static void _variational_model_handle_change( Model *self, Model *model, int index ){
	variational_t* var = (variational_t*)self->obj;
	size_t i = index;
	if(index >= Parameters_count(var->parameters)){
		index -= Parameters_count(var->parameters);
	}
//	self->listeners->fire( self->listeners, self, index );
}

static Model* _variational_model_clone(Model* self, Hashtable *hash){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	variational_t* var = self->obj;
	
	Model* posterior = var->posterior->clone(var->posterior, hash);
	Hashtable_add(hash, posterior->name, posterior);

	Parameters* parameters = new_Parameters(Parameters_count(var->parameters));
	for (int i = 0; i < Parameters_count(var->parameters); i++) {
		char* name = Parameters_name(var->parameters, i);
		if (Hashtable_exists(hash, name)) {
			Parameters_add(parameters, Hashtable_get(hash, name));
		}
		else{
			Parameter* p = clone_Parameter(Parameters_at(var->parameters, i));
			Parameters_move(parameters, p);
			Hashtable_add(hash, name, p);
		}
	}
	
	Parameters* var_parameters = new_Parameters(Parameters_count(var->var_parameters));
	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
		char* name = Parameters_name(var->var_parameters, i);
		if (Hashtable_exists(hash, name)) {
			Parameters_add(var_parameters, Hashtable_get(hash, name));
		}
		else{
			Parameter* p = clone_Parameter(Parameters_at(var->var_parameters, i));
			Parameters_move(var_parameters, p);
			Hashtable_add(hash, name, p);
		}
	}

	variational_t* clone = malloc(sizeof(variational_t));
	clone->posterior = posterior;
	clone->parameters = parameters;
	clone->var_parameters = var_parameters;
	clone->file = NULL;
	clone->sample = var->sample;
	clone->sample_some = var->sample_some;
	clone->elbofn = var->elbofn;
	clone->grad_elbofn = var->grad_elbofn;
	clone->elbo_samples = var->elbo_samples;
	clone->grad_samples = var->grad_samples;
	clone->logP = var->logP;
	clone->parameters_logP = var->parameters_logP;
	clone->initialized = var->initialized;
	clone->finalize = var->finalize;
	clone->iter = var->iter;
	clone->ready_to_sample = var->ready_to_sample;
	clone->rng = var->rng;
	Model* mclone = new_VariationalModel(self->name, clone);

	return mclone;
}

static double _variational_model_logP(Model *self){
	variational_t* var = (variational_t*)self->obj;
	return var->elbofn(var);
}

static void _variational_model_gradient(Model *self, double* grad){
	variational_t* var = (variational_t*)self->obj;
	return var->grad_elbofn(var, grad);
}

static void _variational_model_get_free_parameters(Model* self, Parameters* parameters){
	variational_t* var = (variational_t*)self->obj;
//	Parameters_add_free_parameters(parameters, var->var_parameters);
	Parameters_add_parameters(parameters, var->var_parameters);

}

void _variational_model_reset(Model* self){
	variational_t* var = (variational_t*)self->obj;
	if(var->initialized){
		// Restore the parameters of the posterior
		for (int i = 0; i < Parameters_count(var->parameters); i++) {
			Parameter_restore(Parameters_at(var->parameters, i));
		}
		var->ready_to_sample = false;
		var->initialized = false;
//		printf("reset restore\n");
	}
//	else{
//		
//		printf("reset store\n");
//	}

}

void _variational_model_sample(Model* self, double* samples, double* logP){
	variational_t* var = (variational_t*)self->obj;
	var->sample(var, samples);
	if(logP != NULL){
		*logP = var->logP(var, samples);
	}
}

Model* new_VariationalModel(const char* name, variational_t* var){
	Model *model = new_Model("variational", name, var);
	model->free = _variational_model_free;
	model->clone = _variational_model_clone;
	model->logP = _variational_model_logP;
	model->gradient = _variational_model_gradient;
	model->get_free_parameters = _variational_model_get_free_parameters;
	model->reset = _variational_model_reset;
	model->sample = _variational_model_sample;
	model->samplable = true;
	
	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
		Parameters_at(var->var_parameters, i)->listeners->add( Parameters_at(var->var_parameters, i)->listeners, model );
	}
	return model;
}
