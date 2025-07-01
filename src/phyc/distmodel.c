//
//  distmodel.c
//  physher
//
//  Created by Mathieu Fourment on 3/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "distmodel.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "matrix.h"
// #include "tree.h"

// #include "exponential.h"
// #include "gamma.h"

// #include "filereader.h"
// #include "statistics.h"
// #include "parametersio.h"
// #include "utilsio.h"


double DistributionModel_dlog_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_log_prob_grad_0(DistributionModel* dm, const Parameters* p){
	return 0;
}

double DistributionModel_d2log_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_ddlog_0(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	return 0.0;
}

static void _DistributionModel_error_sample(DistributionModel* dm){
	fprintf(stderr, "sample method not implemented in %s\n", DISTRIBUTION_NAME[dm->type]);
	exit(1);
}
static void _DistributionModel_error_rsample(DistributionModel* dm){
	fprintf(stderr, "rsample method not implemented in %s\n", DISTRIBUTION_NAME[dm->type]);
	exit(1);
}

// do not free tree and simplex since they are managed by their model
static void _free_partial_distribution(DistributionModel*dm){
	free_Parameters(dm->x);
    free_Parameters(dm->parameters);
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	if(dm->gradient != NULL) free(dm->gradient);
	// freeing data is left to the user
	free(dm);
}

//static void _free_full_distribution(DistributionModel*dm){
//    _free_partial_distribution(dm);
//	if(dm->simplex != NULL) free_Simplex(dm->simplex);
//	if(dm->tree != NULL) free_Tree(dm->tree);
//	// freeing data is left to the user
//	free(dm);
//}

static DistributionModel* DistributionMode_clone(DistributionModel* dm){
    Parameters* parameters = NULL;
    Parameters* x = NULL;
    // Simplex* simplex = NULL;
    
    if(dm->parameters != NULL){
		parameters = new_Parameters(Parameters_count(dm->parameters));
		for (size_t j = 0; j < Parameters_count(dm->parameters); j++) {
			Parameters_move(parameters, clone_Parameter(Parameters_at(dm->parameters, j)));
		}
    }
	if(dm->x != NULL){
		x = new_Parameters(Parameters_count(dm->x));
		for (size_t j = 0; j < Parameters_count(dm->x); j++) {
			Parameters_move(parameters, clone_Parameter(Parameters_at(dm->x, j)));
		}
	}
    
	return clone_DistributionModel_with_parameters(dm, parameters, x);
}

DistributionModel* clone_DistributionModel_with_parameters(DistributionModel* dm, Parameters* params, Parameters* x){
	DistributionModel* clone = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(clone);
	clone->type = dm->type;
	clone->parameters = NULL;
    size_t sizeP = 0;
	if(dm->parameters != NULL){
		clone->parameters = new_Parameters(Parameters_count(dm->parameters));
		for (size_t j = 0; j < Parameters_count(dm->parameters); j++) {
			Parameters_add(clone->parameters, Parameters_at(params, j));
			sizeP += Parameter_size(Parameters_at(params, j));
		}
	}
	clone->x = NULL;
	size_t sizeX = 0;
	if(x != NULL){
		clone->x = new_Parameters(Parameters_count(dm->x));
		for (size_t j = 0; j < Parameters_count(dm->x); j++) {
			Parameters_add(clone->x, Parameters_at(x, j));
			sizeX += Parameter_size(Parameters_at(x, j));
		}
	}

	clone->log_prob = dm->log_prob;
	clone->log_prob_grad = dm->log_prob_grad;
	clone->log_prob_hessian_diag = dm->log_prob_hessian_diag;
	clone->reparam_backprop = dm->reparam_backprop;
	clone->entropy = dm->entropy;
	clone->entropy_grad = dm->entropy_grad;
	clone->dlogP = dm->dlogP;
	clone->d2logP = dm->d2logP;
	clone->ddlogP = dm->ddlogP;
	clone->sample = dm->sample;
	clone->rsample = dm->rsample;
	clone->clone = dm->clone;
	clone->free = dm->free;
	clone->tempp = NULL;
	clone->tempx = NULL;
    if(dm->tempp != NULL && sizeP != 0){
        clone->tempp = clone_dvector(dm->tempp, sizeP);
    }
	if(dm->tempx != NULL && sizeX != 0){
		clone->tempx = clone_dvector(dm->tempx, sizeX);
	}
	clone->lp = dm->lp;
    clone->stored_lp = dm->stored_lp;
	clone->shift = dm->shift;
	clone->type = dm->type;
	clone->parameterization = dm->parameterization;
#ifndef GSL_DISABLED	
	clone->rng = dm->rng;
#endif
    clone->data = NULL;
	clone->gradient = NULL;
	clone->gradient_length = dm->gradient_length;
	clone->need_update_gradient = dm->need_update_gradient;
	clone->prepared_gradient = dm->prepared_gradient;
	if(dm->gradient != NULL){
		clone->gradient = clone_dvector(dm->gradient, dm->gradient_length);
	}
	clone->support[0] = dm->support[0];
	clone->support[1] = dm->support[1];
	return clone;
}

DistributionModel* new_DistributionModel(Parameters* p, Parameters* x){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
    dm->x = NULL;
	dm->parameters = NULL;
	
    if(p != NULL){
        dm->parameters = new_Parameters(Parameters_count(p));
        Parameters_add_parameters(dm->parameters, p);
    }

	size_t dimX = 0;
    if(x != NULL){
        dm->x = new_Parameters(Parameters_count(x));
        Parameters_add_parameters(dm->x, x);
	    for(size_t i = 0; i < Parameters_count(x); i++){
		    Parameter* p = Parameters_at(x, i);
		    dimX += Parameter_size(p);
	    }
    }
	dm->tree = NULL;
	dm->log_prob = NULL;
	dm->log_prob_grad = NULL;
	dm->log_prob_hessian_diag = NULL;
	dm->dlogP = NULL;
	dm->d2logP = NULL;
	dm->ddlogP = NULL;
	dm->sample = _DistributionModel_error_sample;
	dm->rsample = _DistributionModel_error_sample;
	dm->entropy = NULL;
	dm->entropy_grad = NULL;
	dm->reparam_backprop = NULL;
	dm->free = _free_partial_distribution;
	dm->clone = DistributionMode_clone;
	dm->data = NULL;
	dm->tempx = dvector(dimX);
	dm->tempp = dvector(dimX);
	dm->need_update = true;
	
	dm->prepared_gradient = 0;;
	dm->gradient = NULL;
	dm->gradient_length = 0;
	dm->need_update_gradient = true;

	dm->support[0] = -INFINITY;
	dm->support[1] = INFINITY;
	return dm;
}

//MARK: tree prior

double DistributionModel_log_uniform_tree(DistributionModel* dm){
	if(dm->need_update){
		int n = Tree_tip_count(dm->tree);
		dm->lp = -logDoubleFactorial(n*2-5);
		dm->need_update = false;
	}
	return dm->lp;
}

DistributionModel* new_UniformTreeDistribution(Tree* tree){
    DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
    assert(dm);
	dm->type = DISTRIBUTION_UNIFORM;
    dm->parameters = NULL;
    dm->x = NULL;
    dm->free = _free_partial_distribution;
    dm->clone = DistributionMode_clone;
    dm->tree = tree;
    dm->tempx = NULL;
    dm->tempp = NULL;
	dm->log_prob = DistributionModel_log_uniform_tree;
	dm->log_prob_grad = DistributionModel_log_prob_grad_0;
	dm->log_prob_hessian_diag = NULL; // no hessian for uniform tree
	dm->reparam_backprop = NULL; // no reparametrization for uniform tree
	dm->sample = _DistributionModel_error_sample;
	dm->rsample = _DistributionModel_error_rsample;
	dm->dlogP = DistributionModel_dlog_0;
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->entropy = NULL;
	dm->entropy_grad = NULL;
	dm->need_update = true;
	dm->support[0] = 0;
	dm->support[1] = INFINITY;
	dm->data = NULL;
    return dm;
}

//MARK: Model functions

void _dist_model_handle_change( Model *self, Model *model, Parameter* parameter, int index ){
	DistributionModel* dm = self->obj;
	dm->need_update = true;
	dm->need_update_gradient = true;
	self->listeners->fire( self->listeners, self, parameter, index );
}

void _dist_model_handle_restore( Model *self, Model *model, int index ){
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _dist_model_store(Model* self){
	self->storedLogP = self->lp;
	DistributionModel* dm = self->obj;
	dm->stored_lp = dm->lp;
	if(dm->parameters != NULL){
    	Parameters_store(dm->parameters);
	}
	if(dm->x != NULL){
		Parameters_store(dm->x);
	}
}

static void _dist_model_restore(Model* self){
	self->lp = self->storedLogP;
	DistributionModel* dm = self->obj;
	dm->lp = dm->stored_lp;
	Parameters_restore(dm->parameters);
	Parameters_restore(dm->x);
	//TODO: think about that
	// bool changed = false;
	// Parameter*p = NULL;
	// // restore the parameters of the model
	// for (size_t j = 0; j < Parameters_count(dm->parameters); j++) {
	// 	p = Parameters_at(dm->parameters, j);
	// 	if (Parameter_changed(p)) {
	// 		changed = true;
	// 		Parameter_restore_quietly(p);
	// 	}
	// }
	// if (changed) {
	// 	p->listeners->fire_restore(p->listeners, NULL, p->id);
	// }
	
	// // restore the domain
	// changed = false;
	// for (int i = 0; i < Parameters_count(dm->x); i++) {
	// 	p = Parameters_at(dm->x, i);
	// 	if (Parameter_changed(p)) {
	// 		changed = true;
	// 		Parameter_restore_quietly(p);
	// 	}
	// }
	// if (changed) {
	// 	p->listeners->fire_restore(p->listeners, NULL, p->id);
	// }
}

static double _dist_model_logP(Model *self){
	DistributionModel* dm = (DistributionModel*)self->obj;
	self->lp = dm->log_prob(dm);
	return self->lp;
}

static double _dist_model_gradient(Model *self, const Parameters* parameters){
	DistributionModel* dm = (DistributionModel*)self->obj;
	return dm->log_prob_grad(dm, parameters);
}

static double _dist_model_hessian_diag(Model *self, const Parameters* parameters){
	DistributionModel* dm = (DistributionModel*)self->obj;
	dm->log_prob_hessian_diag(dm, parameters);
	return 0;
}


static double _dist_model_dlogP(Model *self, const Parameter* p){
	DistributionModel* dm = (DistributionModel*)self->obj;
	return dm->dlogP(dm, p);
}

static double _dist_model_d2logP(Model *self, const Parameter* p){
	DistributionModel* dm = (DistributionModel*)self->obj;
	return dm->d2logP(dm, p);
}

static double _dist_model_ddlogP(Model *self, const Parameter* p1, const Parameter* p2){
	//TODO: implement
	// DistributionModel* dm = (DistributionModel*)self->obj;
	// bool found1 = false;
	// bool found2 = false;
	
	// for (int i = 0; i < Parameters_count(dm->x); i++) {
	// 	if(Parameters_at(dm->x, i) == p1){
	// 		found1 = true;
	// 	}
	// 	else if(Parameters_at(cm->x, i) == p2){
	// 		found2 = true;
	// 	}
	// }
	// if(found1 && found2) return dm->ddlogP(cm, p1, p2);
	
	// found1 = false;
	// found2 = false;
    // for (int i = 0; i < cm->parameter_count; i++) {
    //     for (int j = 0; j < Parameters_count(cm->parameters[i]); j++) {
    //         if(Parameters_at(cm->parameters[i], j) == p1){
    //             found1 = true;
    //         }
    //         else if(Parameters_at(cm->parameters[i], j) == p2){
    //             found2 = true;
    //         }
    //     }
    // }
	// if(found1 && found2) return cm->ddlogP(cm, p1, p2);
	return 0;
}

static void _dist_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free distribution model %s\n", self->name);
		DistributionModel* cm = (DistributionModel*)self->obj;
		// tree
		if(self->data != NULL){
			Model* model = self->data;
			model->free(model);
		}
        cm->free(cm);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _dist_model_clone( Model *self, Hashtable* hash ){
	if(Hashtable_exists(hash, self->name)){
		return Hashtable_get(hash, self->name);
	}
	DistributionModel* dm = (DistributionModel*)self->obj;
	
	// Model* msimplex = (Model*)self->data;
	// Model* msimplexclone = NULL;
	
	Parameters* x = new_Parameters(Parameters_count(dm->x));
	if(dm->x != NULL){
		for (size_t i = 0; i < Parameters_count(dm->x); i++) {
			char* name = Parameters_name(dm->x, i);
			if(Hashtable_exists(hash, name)){
				Parameters_add(x, Hashtable_get(hash, name));
			}
			else{
				Parameter* p = clone_Parameter(Parameters_at(dm->x, i));
				Parameters_move(x, p);
				Hashtable_add(hash, name, p);
			}
		}
	}

	Parameters* params = new_Parameters(Parameters_count(dm->parameters));
	
	// Flat dirichlet does not have parameters
	if(dm->parameters != NULL){
		for (size_t i = 0; i < Parameters_count(dm->parameters); i++) {
			char* name = Parameters_name(dm->parameters, i);
			if(Hashtable_exists(hash, name)){
				Parameters_add(params, Hashtable_get(hash, name));
			}
			else{
				Parameter* p = clone_Parameter(Parameters_at(dm->parameters, i));
				Parameters_move(params, p);
				Hashtable_add(hash, name, p);
			}
		}
	}
	// Simplex* s = NULL;
	// if(msimplexclone != NULL){
	// 	s = (Simplex*)msimplexclone->obj;
	// }
	DistributionModel* dmclone = clone_DistributionModel_with_parameters(dm, params, x);

    free_Parameters(x);
	free_Parameters(params);

	Model* clone = new_DistributionModel2(self->name, dmclone);
	Hashtable_add(hash, clone->name, clone);
	// if(msimplexclone != NULL){
	// 	msimplexclone->free(msimplexclone);
	// }
	clone->store = self->store;
	clone->restore = self->restore;
	clone->storedLogP = self->storedLogP;
	clone->lp = self->lp;
	clone->sample = self->sample;
	clone->rsample = self->rsample;
	clone->samplable = self->samplable;
    // clone->need_update = self->need_update;
	return clone;
}

static void _dist_model_sample(Model* model){
	DistributionModel* dm = (DistributionModel*)model->obj;
	dm->sample(dm);
}

static void _dist_model_rsample(Model* model){
	DistributionModel* dm = (DistributionModel*)model->obj;
	dm->rsample(dm);
}

Model* new_DistributionModel2(const char* name, DistributionModel* dm){
	Model *model = new_Model(MODEL_DISTRIBUTION,name, dm);
	model->logP = _dist_model_logP;
	model->gradient = _dist_model_gradient;
	model->hessian_diag = _dist_model_hessian_diag;
	model->dlogP = _dist_model_dlogP;
	// model->d2logP = _dist_model_d2logP;
	model->ddlogP = _dist_model_ddlogP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
	model->update = _dist_model_handle_change;
	model->handle_restore = _dist_model_handle_restore;
	model->sample = _dist_model_sample;
	model->rsample = _dist_model_rsample;
	model->samplable = false;

	Parameters_add_listener(dm->x, model);
	Parameters_add_parameters_recursively(model->parameters, dm->x);
	if(dm->parameters != NULL){
		Parameters_add_listener(dm->parameters, model);
		Parameters_add_parameters_recursively(model->parameters, dm->parameters);
	}

	return model;
}

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* amodel){
	Model *model = new_DistributionModel2(name, dm);
	model->data = amodel;
    if(amodel != NULL){
        amodel->ref_count++;
        amodel->listeners->add(amodel->listeners, model);
    }
	return model;
}

void distmodel_get_parameters(json_node* parameters_node, Hashtable* hash, Parameters* parameters){
    for (size_t i = 0; i < parameters_node->child_count; i++) {
        json_node* p_node = parameters_node->children[i];
        // it is a reference
        if(p_node->node_type == MJSON_STRING){
            char* ref = p_node->value;
            if(ref[0] == '&'){
                Parameter* p = Hashtable_get(hash, ref+1);
                Parameters_add(parameters, p);
            }
            else if (ref[0] == '%'){
                Parameters* ps = Hashtable_get(hash, ref+1);
                Parameters_add_parameters(parameters, ps);
            }
            else{
                fprintf(stderr, "Distribution only accepts Parameter, Parameters or a list of them\n");
                exit(1);
            }
        }
		else if(p_node->node_type == MJSON_ARRAY){
			for (size_t i = 0; i < p_node->child_count; i++) {
        		json_node* c_node = p_node->children[i];
				distmodel_get_parameters(c_node, hash, parameters);
			}
		}
        else if(p_node->node_type == MJSON_OBJECT &&  strcmp(get_json_node_value_string(p_node, "type"), "parameter") == 0){
            Parameter* parameter = new_Parameter_from_json(p_node, hash);
            Parameters_move(parameters, parameter);
            Hashtable_add(hash, Parameter_name(parameter), parameter);
        }
        else{
            fprintf(stderr, "Distribution only accepts Parameter, Parameters or a list of them\n");
            exit(1);
        }
    }
}

void distmodel_get_ref(const char* ref, Hashtable* hash, Parameters* parameters){
	if (ref[0] == '&') {
		Parameter* p = Hashtable_get(hash, ref+1);
		Parameters_add(parameters, p);
	}
	else if (ref[0] == '%') {
		// slicing
		if (ref[strlen(ref)-1] == ']') {
			get_parameters_slice(ref+1, parameters, hash);
		}
		else{
			Parameters* ps = Hashtable_get(hash, ref+1);
			Parameters_add_parameters(parameters, ps);
		}
	}
}

Parameters* distmodel_get_x(const char* who, json_node* node, Hashtable* hash){
//void get_parameters_from_node(json_node* node, Hashtable* hash, Parameters* parameters){
    Parameters* parameters = new_Parameters(1);
    // it's a ref
    if(node->node_type == MJSON_STRING){
        char* ref = (char*)node->value;
		distmodel_get_ref(ref, hash, parameters);
        Parameters_set_name2(parameters, ref+1);
    }
	else if(node->node_type == MJSON_ARRAY){
		if(node->children[0]->node_type == MJSON_PRIMITIVE){
			double* values = dvector(node->child_count);
			for (size_t i = 0; i < node->child_count; i++) {
				values[i] = atof(node->children[i]->value);
			}
			Parameter* parameter = new_Parameter_with_postfix2("","", values, node->child_count, new_Constraint(-INFINITY, INFINITY));
			Parameters_move(parameters, parameter);
			free(values);
		}
		else{
			for (int i = 0; i < node->child_count; i++) {
				json_node* child = node->children[i];
				char* ref = (char*)child->value;
				distmodel_get_ref(ref, hash, parameters);
			}
		}
	}
    else if(node->node_type == MJSON_OBJECT){
        json_node* p_node_dimension = get_json_node(node, "dimension");
        if(p_node_dimension != NULL ){
            Parameters* multi_parameter = new_MultiParameter_from_json(node, hash);
            Parameters_add_parameters(parameters, multi_parameter);
            Parameters_set_name2(parameters, Parameters_name2(multi_parameter));

            for (int i = 0; i < Parameters_count(multi_parameter); i++) {
                Hashtable_add(hash, Parameters_name(multi_parameter, i), Parameters_at(multi_parameter, i));
            }
            free_Parameters(multi_parameter);
            
            Hashtable_add(hash, Parameters_name2(parameters), parameters);
        }
        else{
            Parameter* p = new_Parameter_from_json(node, hash);
            Parameters_move(parameters, p);
            Parameters_set_name2(parameters, Parameter_name(p));

            Hashtable_add(hash, Parameter_name(p), p);
        }
    }
    else{
        fprintf(stderr, "Do not recognize node type of %s", node->key);
        exit(1);
    }
    return parameters;
}

Parameter* distmodel_parse_parameter(json_node* parameter_node, Hashtable* hash, const char* id, double lower, double upper){
	Parameter* parameter = NULL;
	if(parameter_node->node_type == MJSON_PRIMITIVE){
		parameter = new_Parameter(id, atof(parameter_node->value), new_Constraint(lower, upper));
	}
	else if(parameter_node->node_type == MJSON_ARRAY){
		double* values = dvector(parameter_node->child_count);
		for(size_t i = 0; i < parameter_node->child_count; i++){
			values[i] = atof((char*)parameter_node->children[i]->value);
		}
		parameter = new_Parameter2(id, values, parameter_node->child_count, new_Constraint(lower, upper));
		free(values);
	}
	else if(parameter_node->node_type == MJSON_OBJECT){
		parameter = new_Parameter_from_json(parameter_node, hash);
		Hashtable_add(hash, Parameter_name(parameter), parameter);
	}
	else{
		fprintf(stderr, "Do not recognize node type of %s", parameter_node->key);
		exit(1);
	}
	return parameter;
}
