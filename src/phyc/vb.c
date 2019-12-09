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
#include "matrix.h"

#include "gamvi.h"

#include "klpq.h"
#include "klqp.h"


Model* new_Variational_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"backward",
		"distribution",
		"elbosamples",
		"elbomulti",
		"gradsamples",
		"log",
		"log_samples",
		"parameters",
		"posterior",
		"simplices",
		"tree",
		"var",
		"var_parameters"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	const char* posterior_string = get_json_node_value_string(node, "posterior");
	const char* var_string = get_json_node_value_string(node, "var");
	const char* dist_string = get_json_node_value_string(node, "distribution");
	const char* id = get_json_node_value_string(node, "id");
	const char* filename = get_json_node_value_string(node, "log");
	json_node* parameters_node = get_json_node(node, "parameters");
	json_node* tree_node = get_json_node(node, "tree");
	json_node* simplices_node = get_json_node(node, "simplices");
	
	variational_t* var = malloc(sizeof(variational_t));
	var->file = NULL;
	var->sample = NULL;
	var->sample_some = NULL;
	var->parameters = new_Parameters(1);
	var->var_parameters = NULL;
	
	var->simplex_count = 0;
	var->simplices = NULL;
	if(simplices_node != NULL){
		if (simplices_node->node_type == MJSON_ARRAY) {
			var->simplices = malloc(simplices_node->child_count*sizeof(Model*));
			for (int i = 0; i < simplices_node->child_count; i++) {
				json_node* simplex_node = simplices_node->children[i];
				char* ref = simplex_node->value;
				Model* msimplex = Hashtable_get(hash, ref+1);
				Simplex* simplex = msimplex->obj;
				Parameters_add_parameters(var->parameters, simplex->parameters);
				var->simplices[i] = msimplex;
				msimplex->ref_count++;
				var->simplex_count++;
			}
		}
	}
	
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
	var->elbo_multi = get_json_node_value_int(node, "elbomulti", 1);
	var->log_samples = get_json_node_value_int(node, "log_samples", 0);
	
	json_node* var_parameters_node = get_json_node(node, "var_parameters");
	// Variational model parameters are provided
	if (var_parameters_node->node_type == MJSON_ARRAY) {
		var->var_parameters = new_Parameters(var_parameters_node->child_count);
		for(size_t i = 0; i < Parameters_count(var->var_parameters); i++){
			json_node* var_p_node = var_parameters_node->children[i];
			Parameter* var_p = new_Parameter_from_json(var_p_node, hash);
			Parameters_move(var->var_parameters, var_p);
		}
	}
	
	StringBuffer* buffer = new_StringBuffer(10);
	if(strcasecmp(var_string, "meanfield") == 0){
//		printf("Creating meanfield variational model\n");
        size_t param_count = 2;
        char** params = malloc(param_count*sizeof(char*));
        
        if (strcasecmp(dist_string, "normal") == 0 || strcasecmp(dist_string, "lognormal") == 0) {
            params[0] = String_clone(".mu");
            params[1] = String_clone(".sigma");
        }
        else if(strcasecmp(dist_string, "gamma") == 0){
            params[0] = String_clone(".alpha");
            params[1] = String_clone(".beta");
        }
        else{
            fprintf(stderr, "meanfield distribution should be: normal, lognormal, or gamma");
            exit(1);
        }
        if(var->var_parameters == NULL){
            var->var_parameters = new_Parameters(dim*2); // mean + sd
            for(size_t i = 0; i < dim; i++){
                StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
                StringBuffer_append_string(buffer, params[0]);
                Parameter* p = new_Parameter(buffer->c, 0, NULL);
                p->id = (int)i;
                Parameters_move(var->var_parameters, p);
            }
            for(size_t i = 0; i < dim; i++){
                StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
                StringBuffer_append_string(buffer, params[1]);
                Parameter* p = new_Parameter(buffer->c, 0, NULL);
                p->id = (int)(i+dim);
                Parameters_move(var->var_parameters, p);
            }
        }
        free_cmatrix(params, param_count);
        
		if (strcasecmp(dist_string, "normal") == 0) {
			bool backward = get_json_node_value_bool(node, "backward", true); // true is KL(Q||P)
			if(backward){
				var->elbofn = klqp_meanfield_normal_elbo;
				if (var->elbo_multi > 1) {
					var->elbofn = klqp_meanfield_normal_elbo_multi;
				}
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
            var->print = klqp_meanfield_normal_log_samples;
		}
        else if (strcasecmp(dist_string, "gamma") == 0) {
			var->elbofn = elbo_gamma_meanfield;
			var->grad_elbofn = grad_elbo_gamma_meanfield;
			
			var->sample = variational_sample_gamma_meanfield;
			var->sample_some = variational_sample_some_gamma_meanfield;
			var->finalize = NULL; //TODO: implement function
			var->logP = variational_gamma_meanfield_logP;
			var->parameters_logP = variational_gamma_meanfield_parameters_logP;
            var->print = meanfield_gamma_log_samples;
		}
        else if (strcasecmp(dist_string, "lognormal") == 0) {
            var->elbofn = klqp_meanfield_lognormal_elbo;
            var->grad_elbofn = klqp_meanfield_lognormal_grad_elbo;
            var->sample = klqp_meanfield_lognormal_sample;
            var->sample_some = klqp_meanfield_lognormal_sample_some;
            //TODO: implement lognormal finalize
            var->finalize = klqp_meanfield_normal_finalize;
            var->logP = klqp_meanfield_lognormal_logP;
            var->parameters_logP = klqp_meanfield_lognormal_logP_parameters;
            var->print = klqp_meanfield_lognormal_log_samples;
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
						Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
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
					Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
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
						Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
						Parameters_move(var->var_parameters, p);
						row++;
					}
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".var");
					Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
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
        var->print = klqp_fullrank_log_samples;
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
	if(var->simplices != NULL){
		for(int i = 0; i < var->simplex_count; i++){
			var->simplices[i]->free(var->simplices[i]);
		}
		free(var->simplices);
	}
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
	if (Parameters_name2(var->var_parameters) != NULL) {
		Parameters_set_name2(var_parameters, Parameters_name2(var->var_parameters));
	}
	Hashtable_add(hash, Parameters_name2(var_parameters), var_parameters);
	
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
	var->grad_elbofn(var, grad);
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

static void _variational_print(Model* self, FILE* file){
	variational_t* var = (variational_t*)self->obj;
	if(var->log_samples > 0){
        var->print(var, file);
		return;
	}

	fprintf(file, "muraw=c(");
	size_t i = 0;
	for (; i < Parameters_count(var->parameters); i++) {
		double value = Parameters_value(var->var_parameters, i);
		fprintf(file,"%f%s", value, i ==Parameters_count(var->parameters)-1 ? "":",");
	}
	
	fprintf(file, ")\n\n");
	fprintf(file, "sigma=c(");
	for (; i < Parameters_count(var->var_parameters); i++) {
		fprintf(file, "%f%s", Parameters_value(var->var_parameters, i), i ==Parameters_count(var->var_parameters)-1 ? "":",");
	}
	fprintf(file, ")\n");
	
	fprintf(file, "mu=c(");
	size_t dim = Parameters_count(var->parameters);
	for (i = 0; i < Parameters_count(var->parameters); i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		double sd = Parameters_value(var->var_parameters, i+dim);
		double value = Parameters_value(var->var_parameters, i);
		if(Parameter_lower(p) == 0){
			value = exp(value + sd*sd/2.0);
		}
		fprintf(file,"%f%s", value, i ==Parameters_count(var->parameters)-1 ? "":",");
	}
	fprintf(file, ")\n");
	
//	size_t row = dim;
//	for (int j = 0; j < dim; j++) {
//		double temp = 0;
//		// multiply L_j and eta
//		for (int k = 0; k < j+1; k++) {
//			if(Parameters_estimate(var->var_parameters, row)){
//				temp += Parameters_value(var->var_parameters, row)*eta[k];
//			}
//			row++;
//		}
//		zeta[j] = temp + Parameters_value(var->var_parameters, j); // add mu
//
//		Parameter* p = Parameters_at(var->parameters, j);
//		double theta = inverse_transform2(zeta[j], Parameter_lower(p), Parameter_upper(p));
//		Parameter_set_value(p, theta);
//	}
}

Model* new_VariationalModel(const char* name, variational_t* var){
	Model *model = new_Model(MODEL_VARIATIONAL, name, var);
	model->free = _variational_model_free;
	model->clone = _variational_model_clone;
	model->logP = _variational_model_logP;
	model->gradient = _variational_model_gradient;
	model->get_free_parameters = _variational_model_get_free_parameters;
	model->reset = _variational_model_reset;
	model->print = _variational_print;
	model->sample = _variational_model_sample;
	model->samplable = true;
	
	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
		Parameters_at(var->var_parameters, i)->listeners->add( Parameters_at(var->var_parameters, i)->listeners, model );
	}
	return model;
}
