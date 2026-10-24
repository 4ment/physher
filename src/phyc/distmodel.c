// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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


void DistributionModel_gradient2_0(DistributionModel* dm, Parameters* p){
}

static void _DistributionModel_error_sample(DistributionModel* dm){
	// fprintf(stderr, "sample method not implemented in %s\n", DISTRIBUTION_NAME[dm->parameterization]);
	exit(1);
}
static void _DistributionModel_error_rsample(DistributionModel* dm){
	// fprintf(stderr, "rsample method not implemented in %s\n", DISTRIBUTION_NAME[dm->parameterization]);
	exit(1);
}

// do not free tree and simplex since they are managed by their model
static void _free_partial_distribution(DistributionModel*dm){
	free_Parameters(dm->x);
    free_Parameters(dm->parameters);
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
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

// Copy of `dm` reading `params` and `x` (its own references are taken), with
// scratch buffers of its own. Whatever `data` and `tree` point to is left to the
// caller.
static DistributionModel* _clone_DistributionModel(const DistributionModel* dm,
                                                   Parameters* params, Parameters* x){
	DistributionModel* clone = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(clone);
	// Scalars and the distribution's function table carry over as they are.
	*clone = *dm;
	clone->parameters = NULL;
	if(params != NULL){
		clone->parameters = new_Parameters(Parameters_count(params));
		Parameters_add_parameters(clone->parameters, params);
	}
	clone->x = NULL;
	size_t sizeX = 0;
	if(x != NULL){
		clone->x = new_Parameters(Parameters_count(x));
		Parameters_add_parameters(clone->x, x);
		sizeX = Parameters_size(x);
	}
	// Both are sized by the random variable, as new_DistributionModel does.
	clone->tempx = dm->tempx == NULL ? NULL : clone_dvector(dm->tempx, sizeX);
	clone->tempp = dm->tempp == NULL ? NULL : clone_dvector(dm->tempp, sizeX);
	return clone;
}

static DistributionModel* DistributionMode_clone(DistributionModel* dm){
	Parameters* parameters = dm->parameters == NULL ? NULL : clone_Parameters(dm->parameters);
	Parameters* x = dm->x == NULL ? NULL : clone_Parameters(dm->x);
	DistributionModel* clone = _clone_DistributionModel(dm, parameters, x);
	free_Parameters(parameters);
	free_Parameters(x);
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
        // Parameters_set_name2(dm->parameters, Parameters_name2(p));
    }

	size_t dimX = 0;
    if(x != NULL){
        dm->x = new_Parameters(Parameters_count(x));
        Parameters_add_parameters(dm->x, x);
        // Parameters_set_name2(dm->x, Parameters_name2(x));
	    for(size_t i = 0; i < Parameters_count(x); i++){
		    Parameter* p = Parameters_at(x, i);
		    dimX += Parameter_size(p);
	    }
    }
	// dm->simplex = NULL;
	dm->tree = NULL;
	dm->logP = NULL;
	dm->gradient = NULL;
	dm->sample = _DistributionModel_error_sample;
	dm->rsample = _DistributionModel_error_sample;
	dm->entropy = NULL;
	dm->gradient_entropy = NULL;
	dm->free = _free_partial_distribution;
	dm->clone = DistributionMode_clone;
	dm->data = NULL;
	dm->tempx = dvector(dimX);
	dm->tempp = dvector(dimX);
	dm->need_update = true;

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
    // dm->simplex = NULL;
    dm->free = _free_partial_distribution;
    dm->clone = DistributionMode_clone;
    dm->tree = tree;
    dm->tempx = NULL;
    dm->tempp = NULL;
	dm->logP = DistributionModel_log_uniform_tree;
	dm->gradient = DistributionModel_gradient2_0;
	dm->entropy = NULL;
	dm->gradient_entropy = NULL;
	dm->need_update = true;
	dm->support[0] = 0;
	dm->support[1] = INFINITY;
    return dm;
}

//MARK: Model functions

void _dist_model_handle_change( Model *self, Model *model, Parameter* parameter, int index ){
	DistributionModel* dm = self->obj;
	dm->need_update = true;
	self->listeners->fire( self->listeners, self, parameter, index );
}

static void _dist_model_store(Model* self){
	if(!self->stored){
		self->storedLogP = self->lp;
		DistributionModel* dm = self->obj;
		dm->stored_lp = dm->lp;
		if(dm->parameters != NULL){
			Parameters_store(dm->parameters);
		}
		if(dm->x != NULL){
			Parameters_store(dm->x);
		}
		self->stored = true;
	}
}

static void _dist_model_restore(Model* self){
	if(self->stored){
		self->lp = self->storedLogP;
		DistributionModel* dm = self->obj;
		dm->lp = dm->stored_lp;
		if(dm->parameters != NULL){
			Parameters_restore(dm->parameters);
		}
		if(dm->x != NULL){
			Parameters_restore(dm->x);
		}
		self->stored = false;
	}
}

static void _dist_model_accept(Model* self){
	if(self->stored){
		DistributionModel* dm = self->obj;
		if(dm->parameters != NULL){
			Parameters_accept(dm->parameters);
		}
		if(dm->x != NULL){
			Parameters_accept(dm->x);
		}
		self->stored = false;
	}
}

static double _dist_model_logP(Model *self){
	DistributionModel* dm = (DistributionModel*)self->obj;
	self->lp = dm->logP(dm);
	return self->lp;
}

static void _dist_model_gradient(Model *self, Parameters* parameters){
	DistributionModel* dm = (DistributionModel*)self->obj;
	dm->gradient(dm, parameters);
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

Model* DistributionModel_clone_model(Model* self, Hashtable* hash){
	DistributionModel* dm = (DistributionModel*)self->obj;
	Parameters* x = clone_Parameters_into(dm->x, hash);
	Parameters* params = clone_Parameters_into(dm->parameters, hash);
	DistributionModel* dmclone = _clone_DistributionModel(dm, params, x);
	free_Parameters(x);
	free_Parameters(params);

	// The model the density reads (a tree, a coalescent), held both as the
	// submodel and through dm->tree or dm->data.
	Model* amodel = self->data;
	Model* amodelclone = Model_clone_into(amodel, hash);
	if(dm->tree != NULL){
		if(amodel == NULL || dm->tree != amodel->obj){
			fprintf(stderr, "clone: distribution %s reads a tree it does not hold\n", self->name);
			exit(1);
		}
		dmclone->tree = amodelclone->obj;
	}
	dmclone->data = NULL;
	if(dm->data != NULL && amodel != NULL && dm->data == amodel->obj){
		dmclone->data = amodelclone->obj;
	}

	Model* clone = new_DistributionModel3(self->name, dmclone, amodelclone);
	if(amodelclone != NULL) amodelclone->free(amodelclone);
	Model_copy_behaviour(clone, self);
	return clone;
}

static Model* _dist_model_clone( Model *self, Hashtable* hash ){
	Model* clone = DistributionModel_clone_model(self, hash);
	DistributionModel* dm = self->obj;
	if(dm->data != NULL && ((DistributionModel*)clone->obj)->data == NULL){
		fprintf(stderr, "clone: distribution %s holds data it cannot clone\n", self->name);
		exit(1);
	}
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
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
	model->accept = _dist_model_accept;
	model->update = _dist_model_handle_change;
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
            if(ref[0] == '&' || ref[0] == '%'){
                get_parameter_reference(ref, hash, parameters);
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
	get_parameter_reference(ref, hash, parameters);
	// simplex
	// else if (ref[0] == '$') {
	// 	Model* msimplex = Hashtable_get(hash, ref+1);
	// 	Simplex* simplex = msimplex->obj;
	// 	Parameters_add_parameters(parameters, simplex->parameters);
	// }
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
