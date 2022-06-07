/*
 *  demographicmodels.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/2/10.
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

#include "demographicmodels.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <strings.h>

#include "tree.h"
#include "mstring.h"
#include "matrix.h"
#include "combinatorics.h"
#include "mathconstant.h"
#include "gradient.h"

#define CHOOSE2(n) 0.5*n*(n-1)
#pragma mark Coalescent

static void _constant_calculate_gradient( Coalescent* coal );
static void _skyride_calculate_gradient( Coalescent* coal );
static void _skyride_calculate_gradient_real( Coalescent* coal );
static void _skygrid_calculate_gradient( Coalescent* coal );
static void _skygrid_calculate_gradient_real( Coalescent* coal );

Coalescent * create_Coalescent(demography type, size_t size);

void coalescentToLineage(Coalescent* c){
	int lineages = 0;
	for(int i = 0; i < c->n; i++){
		c->lineages[i] = lineages;
		if(c->iscoalescent[i]) lineages--;
		else lineages++;
	}
}

static void _coalescent_model_store(Model* self){
	Coalescent* coalescent = (Coalescent*)self->obj;
	for (int i = 0; i < Parameters_count(coalescent->p); i++) {
		Parameter_store(Parameters_at(coalescent->p, i));
	}
	memcpy(coalescent->stored_iscoalescent, coalescent->iscoalescent, coalescent->n*sizeof(bool));
	memcpy(coalescent->stored_times, coalescent->times, coalescent->n*sizeof(double));
	memcpy(coalescent->stored_lineages, coalescent->lineages, coalescent->n*sizeof(int));
	
	if(self->data != NULL){
		Model** models = self->data;
		if(models[1] != NULL) models[1]->store(models[1]);
	}
	self->storedLogP = self->lp;
	coalescent->stored_logP = coalescent->logP;
}

static void _coalescent_model_restore(Model* self){
	bool changed = false;
	Parameter*p = NULL;
	Coalescent* coalescent = (Coalescent*)self->obj;
	for (int i = 0; i < Parameters_count(coalescent->p); i++) {
		p = Parameters_at(coalescent->p, i);
		if (Parameter_changed(p)) {
			changed = true;
			Parameter_restore_quietly(p);
		}
	}
	// fire only once
	if (changed) {
		p->listeners->fire_restore(p->listeners, NULL, p->id);
	}
	if(self->data != NULL){
		Model** models = self->data;
		if(models[1] != NULL) models[1]->restore(models[1]);
	}
	self->lp = self->storedLogP;
}

static double _coalescent_model_logP(Model *self){
	Coalescent* mc = (Coalescent*)self->obj;
	self->lp = mc->calculate(mc);
	return self->lp;
}

static double _coalescent_model_full_logP(Model *self){
	Coalescent* mc = (Coalescent*)self->obj;
	mc->need_update = true;
	if(self->data != NULL) mc->need_update_intervals = true;
	self->lp = mc->calculate(mc);
	return self->lp;
}

static double _coalescent_model_dlogP(Model *self, const Parameter* p){
	Coalescent* mc = (Coalescent*)self->obj;
	return mc->dlogP(mc, p);
}

// Set flags for gradient calculation
void _coalescent_model_prepare_gradient(Model* self, const Parameters* ps){
	Coalescent* coal = (Coalescent*)self->obj;
	size_t paramCount = Parameters_count(ps);
	bool prepare_tree = false;
	bool prepare_theta = false;
	coal->prepared_gradient = 0;
	size_t gradient_length = 0;
	for (size_t i = 0; i < paramCount; i++) {
		Parameter* p = Parameters_at(ps, i);
		if (p->model == MODEL_COALESCENT && !prepare_theta) {
			prepare_theta = true;
			coal->prepared_gradient |= GRADIENT_FLAG_COALESCENT_THETA;
			gradient_length += Parameters_count(coal->p);
		}
		else if (p->model == MODEL_TREE && !prepare_tree) {
			prepare_tree = true;
			coal->prepared_gradient |= GRADIENT_FLAG_TREE_HEIGHTS;
			gradient_length += Tree_tip_count(coal->tree) - 1;
		}
		else if (p->model == MODEL_TREE_TRANSFORM && !prepare_tree) {
			prepare_tree = true;
			coal->prepared_gradient |= GRADIENT_FLAG_TREE_RATIOS;
			gradient_length += Tree_tip_count(coal->tree) - 1;
		}
	}
	if(coal->gradient == NULL){
		coal->gradient = calloc(gradient_length, sizeof(double));
		coal->gradient_length = gradient_length;
	}
	else if (coal->gradient_length < gradient_length) {
		coal->gradient = realloc(coal->gradient, sizeof(double)* gradient_length);
		coal->gradient_length = gradient_length;
	}
}

// by default we get the derivative wrt reparameterized node height
size_t Coalescent_initialize_gradient(Model *self, int flags){
	Coalescent* coal = (Coalescent*)self->obj;
	coal->prepared_gradient = flags;
	size_t gradient_length = 0;
	
	int prepare_tree_height = coal->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS;
	int prepare_tree_ratio = coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS;
	int prepare_theta = coal->prepared_gradient & GRADIENT_FLAG_COALESCENT_THETA;
	
	if(flags == 0){
		prepare_tree_ratio = true;
		prepare_tree_height = false;
		prepare_theta = true;
		coal->prepared_gradient = GRADIENT_FLAG_TREE_RATIOS | GRADIENT_FLAG_COALESCENT_THETA;
	}
	if(prepare_tree_ratio || prepare_tree_height){
		gradient_length += Tree_tip_count(coal->tree) - 1;
	}
	if (prepare_theta) {
		gradient_length += Parameters_count(coal->p);
		
	}

	if(coal->gradient == NULL){
		coal->gradient = calloc(gradient_length, sizeof(double));
		coal->gradient_length = gradient_length;
	}
	else if (coal->gradient_length < gradient_length) {
		coal->gradient = realloc(coal->gradient, sizeof(double)* gradient_length);
		coal->gradient_length = gradient_length;
	}
	return gradient_length;
}

double* Coalescent_gradient(Model *self){
	Coalescent* coal = (Coalescent*)self->obj;
	
//	if(coal->gradient_length == 0) return 0.0;
	if ( coal->need_update_intervals){
		coal->update_intervals(coal);
	}
	// coal->need_update is set to false if we call logP first
	if (coal->need_update_gradient) {
		coal->need_update_gradient = false;
		if(coal->type == COALESCENT_CONSTANT){
			_constant_calculate_gradient(coal);
		}
		else if(coal->type == COALESCENT_SKYRIDE){
			if(coal->parameterization == COALESCENT_THETA_LOG){
				_skyride_calculate_gradient(coal);
			}
			else if(coal->parameterization == COALESCENT_THETA){
				_skyride_calculate_gradient_real(coal);
			}
		}
		else if(coal->type == COALESCENT_SKYGRID){
			if(coal->parameterization == COALESCENT_THETA_LOG){
				_skygrid_calculate_gradient(coal);
			}
			else if(coal->parameterization == COALESCENT_THETA){
				_skygrid_calculate_gradient_real(coal);
			}
		}
	}
	return coal->gradient;
}

double _coalescent_model_dlogP_prepared(Model *self, const Parameter* p){
	Coalescent* coal = (Coalescent*)self->obj;
	
	if(coal->gradient_length == 0) return 0.0;
	Coalescent_gradient(self);
	
	int prepare_tree_height = coal->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS;
	int prepare_tree_ratio = coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS;
	int prepare_theta = coal->prepared_gradient & GRADIENT_FLAG_COALESCENT_THETA;
	size_t offset = 0;
	
	if (prepare_theta) {
		if (p->model == MODEL_COALESCENT) {
			return coal->gradient[p->id];
		}
		offset += Parameters_count(coal->p);
	}
	
	if(prepare_tree_height && p->model == MODEL_TREE){
		return coal->gradient[offset + Tree_node(coal->tree, p->id)->class_id];
	}
	else if(prepare_tree_ratio && p->model == MODEL_TREE_TRANSFORM){
		return coal->gradient[offset + p->id];
	}
	return 0;
}

static double _coalescent_model_d2logP(Model *self, const Parameter* p){
	Coalescent* mc = (Coalescent*)self->obj;
	return mc->d2logP(mc, p);
}

static double _coalescent_model_ddlogP(Model *self, const Parameter* p1, const Parameter* p2){
	Coalescent* mc = (Coalescent*)self->obj;
	return mc->ddlogP(mc, p1, p2);
}

static void _coalescent_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free coalescent model %s\n", self->name);
		Coalescent* mc = (Coalescent*)self->obj;
		free_Parameters(mc->p);
		free(mc->lineages);
		if(mc->nodes){
			for (int i = 0; i < mc->n; i++) {
				free(mc->nodes[i]);
			}
			free(mc->nodes);
		}
		free(mc->times);
		free(mc->iscoalescent);
		free(mc->stored_lineages);
		free(mc->stored_times);
		free(mc->stored_iscoalescent);
		if(mc->grid != NULL) free(mc->grid);
		if(mc->gradient != NULL) free(mc->gradient);
		free(mc);
		// tree or simplex
		if(self->data != NULL){
			Model** models = self->data;
			models[0]->free(models[0]); // tree
			if(models[1] != NULL) models[1]->free(models[1]); // groups
			free(self->data);
		}
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

//MARK: TODO
static Model* _coalescent_model_clone( Model *self, Hashtable* hash ){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	
	Model** models = (Model**)self->data;
	Model *mtree = models[0];
	Model *mgroup = models[1];
	Model* mtreeclone = NULL;
	Model* mgroupclone = NULL;
	DiscreteParameter* groupclone = NULL;
    
	if(models != NULL){
		if (Hashtable_exists(hash, mtree->name)) {
			mtreeclone = Hashtable_get(hash, mtree->name);
			mtreeclone->ref_count++; // it is decremented at the end using free
		}
		else{
			mtreeclone = mtree->clone(mtree, hash);
			Hashtable_add(hash, mtreeclone->name, mtreeclone);
		}
		
		if(mgroup != NULL){
			if (Hashtable_exists(hash, mgroup->name)) {
				mgroupclone = Hashtable_get(hash, mgroup->name);
				mgroupclone->ref_count++; // it is decremented at the end using free
			}
			else{
				mgroupclone = mgroup->clone(mgroup, hash);
				Hashtable_add(hash, mgroupclone->name, mgroupclone);
			}
			groupclone = mgroupclone->obj;
		}
	}
	Coalescent* coalescent = self->obj;
	size_t size = coalescent->n;
	Coalescent* clone = create_Coalescent(coalescent->type, size);
	clone->calculate = coalescent->calculate;
	clone->dlogP = coalescent->dlogP;
	clone->d2logP = coalescent->d2logP;
	clone->ddlogP = coalescent->ddlogP;
	memcpy(clone->lineages, coalescent->lineages, sizeof(int)*size);
	memcpy(clone->times, coalescent->times, sizeof(double)*size);
	memcpy(clone->stored_lineages, coalescent->stored_lineages, sizeof(int)*size);
	memcpy(clone->stored_times, coalescent->stored_times, sizeof(double)*size);
	memcpy(clone->iscoalescent, coalescent->iscoalescent, sizeof(bool)*size);
	memcpy(clone->stored_iscoalescent, coalescent->stored_iscoalescent, sizeof(bool)*size);
	
	clone->p = new_Parameters(Parameters_count(coalescent->p));
	if(Parameters_name2(coalescent->p) != NULL){
		Parameters_set_name2(clone->p, Parameters_name2(coalescent->p));
		Hashtable_add(hash, Parameters_name2(clone->p), clone->p);
	}
	
	for (int i = 0; i < Parameters_count(coalescent->p); i++) {
		char* name = Parameters_name(coalescent->p, i);
		if (Hashtable_exists(hash, name)) {
			Parameters_add(clone->p, Hashtable_get(hash, name));
		}
		else{
			Parameter* p = clone_Parameter(Parameters_at(coalescent->p, i));
			Parameters_move(clone->p, p);
			Hashtable_add(hash, name, p);
		}
	}
	
	clone->need_update_intervals = coalescent->need_update_intervals;
	clone->need_update = coalescent->need_update;
	
	if(mtreeclone != NULL){
		clone->tree = mtreeclone->obj;
		clone->nodes = malloc(size*sizeof(double_int_pair_t*));
		for (int i = 0; i < size; i++) {
			clone->nodes[i] = malloc(sizeof(double_int_pair_t));
		}
		clone->need_update_intervals = true;
		clone->need_update = true;
		clone->need_update_gradient = true;
	}
	
	clone->gridCount = coalescent->gridCount;
	clone->grid = NULL;
	if(coalescent->grid != NULL){
		clone->grid = clone_dvector(coalescent->grid, coalescent->gridCount-1); // not 0 at grid[0]
	}
	clone->groups = groupclone;
    clone->gradient = NULL;
	if (coalescent->gradient != NULL) {
		clone->gradient = clone_dvector(coalescent->gradient, size);
	}
	Model* mclone = new_CoalescentModel2(self->name, clone, mtreeclone, mgroupclone);
	if(mtreeclone!=NULL) mtreeclone->free(mtreeclone);
	if(mgroupclone != NULL)mgroupclone->free(mgroupclone);
	
	return mclone;
}

//MARK: TODO
static void _coalescent_model_sample(Model* model, double* samples, double* logP){
	error("todo _coalescent_model_sample\n");
}

//MARK: TODO
static double _coalescent_model_sample_evaluate(Model* model){
	error("todo _coalescent_model_sample_evaluate\n");
	return 0;
}

static void _coalescent_model_handle_change( Model *self, Model *model, int index ){
	Coalescent *c = (Coalescent*)self->obj;
	// theta has changed
	if(model == NULL){
		c->need_update = true;
		c->need_update_gradient = true;
	}
	else if ( model->type == MODEL_TREE ) {
		c->need_update_intervals = true;
		c->need_update = true;
		c->need_update_gradient = true;
	}
	else if ( model->type == MODEL_DISCRETE_PARAMETER ) {
		c->need_update = true;
		c->need_update_gradient = true;
	}
	else{
		fprintf(stderr, "_coalescent_model_handle_change %d %s\n", model->type, model->name);
		exit(1);
	}
	self->listeners->fire( self->listeners, self, index );
}

static void _coalescent_model_handle_restore( Model *self, Model *model, int index ){
	Coalescent *c = (Coalescent*)self->obj;
	memcpy(c->iscoalescent, c->stored_iscoalescent, c->n*sizeof(bool));
	memcpy(c->times, c->stored_times, c->n*sizeof(double));
	memcpy(c->lineages, c->stored_lineages, c->n*sizeof(int));
	c->logP = c->stored_logP;
	self->listeners->fire_restore( self->listeners, self, index );
}

Model* new_CoalescentModel2(const char* name, Coalescent* coalescent, Model* tree, Model* groups){
	Model* model = new_Model(MODEL_COALESCENT, name, coalescent);
	for ( int i = 0; i < Parameters_count(coalescent->p); i++ ) {
		Parameters_at(coalescent->p, i)->listeners->add( Parameters_at(coalescent->p, i)->listeners, model );
	}
	if(tree != NULL)tree->listeners->add( tree->listeners, model );
	
	model->logP = _coalescent_model_logP;
	model->full_logP = _coalescent_model_full_logP;
	model->dlogP = _coalescent_model_dlogP_prepared;
	model->d2logP = _coalescent_model_d2logP;
	model->ddlogP = _coalescent_model_ddlogP;
	model->free = _coalescent_model_free;
	model->clone = _coalescent_model_clone;
	model->store = _coalescent_model_store;
	model->restore = _coalescent_model_restore;
	model->sample = _coalescent_model_sample;
	model->sample_evaluate = _coalescent_model_sample_evaluate;
	
	model->update = _coalescent_model_handle_change;
	model->handle_restore = _coalescent_model_handle_restore;
	model->prepare_gradient = _coalescent_model_prepare_gradient;
	
	model->data = NULL;
	if(tree != NULL || groups != NULL){
		Model** models = calloc(2, sizeof(Model*));
		if(tree != NULL){
			models[0] = tree;
			tree->ref_count++;
		}
	
		if(groups != NULL){
			models[1] = groups;
			groups->ref_count++;
			groups->listeners->add(groups->listeners, model);
		}
		model->data = models;
	}
	return model;
}

Model* new_CoalescentModel(const char* name, Coalescent* coalescent, Model* tree){
	return new_CoalescentModel2(name, coalescent, tree, NULL);
}

Model* new_CoalescentModel_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"cutoff",
		"data",
		"groups",
		"model",
		"parameterization",
		"parameters",
		"tree"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* model = get_json_node_value_string(node, "model");
	Coalescent* c = NULL;
	
	json_node* data_node = get_json_node(node, "data");
	json_node* tree_node = get_json_node(node, "tree");
	char* parameterization_str = get_json_node_value_string(node, "parameterization");
	Model* mtree = NULL;
	Tree* tree = NULL;
	Model* mdp = NULL;
	
	// only for skyride and skygrid
	parameterization_t parameterization = COALESCENT_THETA_LOG;
	if(parameterization_str != NULL){
		if(strcasecmp(parameterization_str, "real")){
			parameterization = COALESCENT_THETA;
		}
		else if(strcasecmp(parameterization_str, "delta")){
			parameterization = COALESCENT_THETA_DELTA;
		}
	}
		
	double* times = NULL;
	bool* coalescent = NULL;
	int intervalCount = -1;
	if(tree_node != NULL){
		if(tree_node->node_type == MJSON_OBJECT){
			char* id = get_json_node_value_string(tree_node, "id");
			mtree = new_TreeModel_from_json(tree_node, hash);
			Hashtable_add(hash, id, mtree);
		}
		else{
			char* ref = (char*)tree_node->value;
			// check it starts with a &
			mtree = Hashtable_get(hash, ref+1);
			mtree->ref_count++;
		}
		tree = mtree->obj;
	}
	else{
		json_node* times_node = get_json_node(data_node, "times");
		json_node* coal_node = get_json_node(data_node, "coalescent");
		times = dvector(times_node->child_count);
		coalescent = bvector(times_node->child_count);
		intervalCount = times_node->child_count;
		for(int i = 0; i < times_node->child_count; i++){
			times[i] = atof(times_node->children[i]->value);
			coalescent[i] = atoi(coal_node->children[i]->value);
		}
	}
	
	json_node* parameters_node = get_json_node(node, "parameters");
    Parameters* ps = new_Parameters(parameters_node->child_count);
    
    for (int i = 0; i < parameters_node->child_count; i++) {
        json_node* p_node = parameters_node->children[i];
        json_node* dim_node = get_json_node(p_node, "dimension");
        if (dim_node != NULL) {
            Parameters* multi_parameter = new_MultiParameter_from_json(p_node, hash);
            Parameters_add_parameters(ps, multi_parameter);
            Hashtable_add(hash, Parameters_name2(multi_parameter), multi_parameter);
        }
        else{
            Parameter* parameter = new_Parameter_from_json(p_node, hash);
            Parameters_move(ps, parameter);
            Hashtable_add(hash, Parameter_name(parameter), parameter);
        }
        
    }
	
	if(strcasecmp(model, "constant") == 0){
		if(tree != NULL)
			c = new_ConstantCoalescent(tree, Parameters_at(ps, 0));
		else
			c = new_ConstantCoalescent_with_data(Parameters_at(ps, 0), times, coalescent, intervalCount);
		Hashtable_add(hash, Parameters_name(ps, 0), Parameters_at(ps, 0));
	}
	else if(strcasecmp(model, "exponential") == 0){
		if (strcasecmp(parameters_node->children[0]->key, "n0") != 0) {
			Parameters_swap_index(ps, 0, 1);
		}
		if(tree != NULL)
			c = new_ExponentialCoalescent(tree, ps);
		else
			c = new_ExponentialCoalescent_with_data(ps, times, coalescent, intervalCount);
		Hashtable_add(hash, Parameters_name(ps, 0), Parameters_at(ps, 0));
		Hashtable_add(hash, Parameters_name(ps, 1), Parameters_at(ps, 1));
	}
	else if(strcasecmp(model, "skygrid") == 0){
		int gridCount = Parameters_count(ps);
		
		double cutoff = get_json_node_value_double(node, "cutoff", -1);
		if(cutoff <= 0){
			fprintf(stderr, "cutoff is not optional and should be greater than 0\n");
			exit(3);
		}
		if(tree != NULL)
			c = new_GridCoalescent(tree, ps, gridCount, cutoff, parameterization);
		else
			c = new_GridCoalescent_with_data(ps, times, coalescent, intervalCount, gridCount, cutoff, parameterization);

	}
	else if(strcasecmp(model, "skyline") == 0){
		json_node* groups_node = get_json_node(node, "groups");
		if (groups_node->node_type == MJSON_OBJECT) {
			mdp = new_DiscreteParameterModel_from_json(groups_node, hash);
			json_node* id_node = get_json_node(groups_node, "id");
			Hashtable_add(hash, (char*)id_node->value, mdp);
		}
		else if(groups_node->node_type == MJSON_STRING){
			char* ref = (char*)groups_node->value;
			// check it starts with a &
			mdp = Hashtable_get(hash, ref+1);
			mdp->ref_count++;
		}
		DiscreteParameter* dp = mdp->obj;
		if(tree != NULL)
			c = new_SkylineCoalescent(tree, ps, dp);
		else
			c = new_SkylineCoalescent_with_data(ps, times, coalescent, intervalCount, dp);
	}
	else if(strcasecmp(model, "classical") == 0){
		c = new_ClassicalSkylineCoalescent_with_parameters(tree, ps);
	}
	else if(strcasecmp(model, "skyride") == 0){
		if(tree != NULL)
			c = new_SkyrideCoalescent(tree, ps, parameterization);
		else
			c = new_SkyrideCoalescent_with_data(ps, times, coalescent, intervalCount, parameterization);
	}
	else{
		fprintf(stderr, "Coalescent model unknown %s\n", model);
		exit(1);
	}
	
	free_Parameters(ps);
	
	// add Paramters to hash
	char* id_ps = get_json_node_value_string(parameters_node, "id");
	if(id_ps != NULL){
		Parameters_set_name2(c->p, id_ps);
		Hashtable_add(hash, id_ps, c->p);
	}
	
	char* id = get_json_node_value_string(node, "id");
	Model* mc = new_CoalescentModel2(id, c, mtree, mdp);
	if(mtree!=NULL) mtree->free(mtree);
	if(mdp != NULL)mdp->free(mdp);
	
	if(times != NULL){
		free(times);
		free(coalescent);
		c->need_update_intervals = false;
	}
	
	return mc;
}

#pragma mark -

void free_Coalescent( Coalescent *coalescent ){
	free_Parameters(coalescent->p);
	free(coalescent->lineages);
	if(coalescent->nodes != NULL){
		for (int i = 0; i < coalescent->n; i++) {
			free(coalescent->nodes[i]);
		}
		free(coalescent->nodes);
	}
	free(coalescent->times);
	free(coalescent->iscoalescent);
	free(coalescent->stored_lineages);
	free(coalescent->stored_times);
	free(coalescent->stored_iscoalescent);
	if(coalescent->grid != NULL) free(coalescent->grid);
	if(coalescent->groups != NULL) coalescent->groups->free(coalescent->groups);
	if(coalescent->gradient != NULL) free(coalescent->gradient);
	free(coalescent);
}


void _update_intervals(Coalescent* coal){
	Node **nodes = Tree_get_nodes(coal->tree, POSTORDER);
	size_t nodeCount = Tree_node_count(coal->tree);
	memset(coal->times, 0, nodeCount*sizeof(double));
	memset(coal->iscoalescent, 0, nodeCount*sizeof(bool));
	memset(coal->lineages, 0,nodeCount*sizeof(int));
	Tree_update_heights(coal->tree);
	size_t tipCount = Tree_tip_count(coal->tree);
	for (int i = 0; i < nodeCount; i++) {
		coal->nodes[i]->index = Node_id(nodes[i]);
		coal->nodes[i]->value = Node_height(nodes[i]);
	}
	qsort(coal->nodes, nodeCount, sizeof(struct double_int_pair_t*), cmp_double_int_pair_asc);
	double start = coal->nodes[0]->value;
	int lineageCount = 0;
	int nodeIndex = 0;
	int intervalCount = 0;
	double finish;
	
	while (nodeIndex < nodeCount) {
		bool isSampling = coal->nodes[nodeIndex]->index < tipCount;
		finish = coal->nodes[nodeIndex]->value;
		nodeIndex++;
		
		// sampling event
		if(isSampling){
			coal->times[intervalCount] = finish - start;
			coal->lineages[intervalCount] = lineageCount;
			intervalCount++;
			start = finish;
			
			lineageCount++;
		}
		// coalescent event
		else {
			coal->times[intervalCount] = finish - start;
			coal->lineages[intervalCount] = lineageCount;
			coal->iscoalescent[intervalCount] = true;
			start = finish;
			intervalCount++;
			
			lineageCount--;
		}
	}
	
	coal->n = intervalCount;
	coal->need_update_intervals = false;
	coal->need_update = true;
	coal->need_update_gradient = true;
}

void _update_intervals_grid(Coalescent* coal){
	Node **nodes = Tree_get_nodes(coal->tree, POSTORDER);
	size_t nodeCount = Tree_node_count(coal->tree);
	size_t points_count = nodeCount + coal->gridCount;
	memset(coal->times, 0, points_count*sizeof(double));
	memset(coal->iscoalescent, 0, points_count*sizeof(bool));
	memset(coal->lineages, 0, points_count*sizeof(int));
	Tree_update_heights(coal->tree);

	size_t tipCount = Tree_tip_count(coal->tree);
	size_t i = 0;
	for ( ; i < nodeCount; i++) {
		coal->nodes[i]->index = Node_id(nodes[i]);
		coal->nodes[i]->value = Node_height(nodes[i]);
	}
	for (size_t j = 0; j < coal->gridCount; j++, i++) {
		coal->nodes[i]->index = -1;
		coal->nodes[i]->value = coal->grid[j];
	}
	qsort(coal->nodes, points_count, sizeof(struct double_int_pair_t*), cmp_double_int_pair_asc);
	
	// never starts with a grid point
	double start = coal->nodes[0]->value;
	int lineageCount = 0;
	double finish;
	
	for (size_t i = 0; i < points_count; i++) {
		finish = coal->nodes[i]->value;
		coal->lineages[i] = lineageCount;
		
		// not a grid point
		if (coal->nodes[i]->index >= 0) {
			// sampling event
			if (coal->nodes[i]->index < tipCount) {
				lineageCount++;
			}
			// coalescent event
			else {
				coal->iscoalescent[i] = true;
				lineageCount--;
			}
		}
		coal->times[i] = finish - start;
		start = finish;
	}
	
	coal->n = points_count;
	coal->need_update_intervals = false;
	coal->need_update = true;
	coal->need_update_gradient = true;
}


static void _premultiply_proportions(Node* node, double* descendant, Parameters* reparams){
    if(!Node_isleaf(node)){
		size_t node_id = Node_id(node);
        descendant[node_id] = descendant[Node_id(node->parent)]*Parameters_value(reparams, Node_class_id(node));
        _premultiply_proportions(node->left, descendant, reparams);
        _premultiply_proportions(node->right, descendant, reparams);
    }
}

void height_gradient_from_interval_gradient(Coalescent* coal, const double* interval_gradient, double* height_gradient){
	Node** nodes = Tree_nodes(coal->tree);
	size_t i = 1;
	for( ; i < coal->n - 1; i++){
		if(coal->iscoalescent[i]){
			size_t node_class_id = Node_class_id(nodes[coal->nodes[i]->index]);
			height_gradient[node_class_id] += interval_gradient[i];
			height_gradient[node_class_id] -= interval_gradient[i+1];
		}
	}
	height_gradient[Node_class_id(Tree_root(coal->tree))] += interval_gradient[i];
}

#pragma mark -
#pragma mark Constant coalescent


double _constant_calculate( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	if ( coal->need_update ) {
		coal->logP = 0;
		double theta = Parameters_value(coal->p, 0);
		double logTheta = log(theta);
	 
		for( int i = 0; i< coal->n; i++  ){
			double lambda = CHOOSE2(coal->lineages[i])/theta;
			coal->logP -= lambda * coal->times[i];
			
			if( coal->iscoalescent[i] ){
				coal->logP -= logTheta;
			}
			
		}
		coal->need_update = false;
	}
	return coal->logP;
}

void _constant_calculate_gradient( Coalescent* coal ){
	size_t offset = 0;
	double* chooses = dvector(coal->n);
	for( int i = 0; i< coal->n; i++  ){
		chooses [i] = CHOOSE2(coal->lineages[i]);
	}
	double theta = Parameters_value(coal->p, 0);
	double theta_inv = 1.0/theta;

    if(coal->prepared_gradient & GRADIENT_FLAG_COALESCENT_THETA){
        double theta2 = theta*theta;
		double dlogP = 0;
        for( int i = 0; i< coal->n; i++  ){
            dlogP += chooses[i] * coal->times[i] /theta2;
			if (coal->iscoalescent[i]) {
				dlogP -= theta_inv;
			}
        }
		coal->gradient[0] = dlogP;
		offset++;
    }
	
	if(coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS
	   || coal->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS){
		double* interval_gradient = dvector(coal->n);
		for(size_t i = 0; i < coal->n; i++){
			interval_gradient[i] = -chooses[i]*theta_inv;
		}
		// derivatives wrt to reparameterization
		if(get_reparams(coal->tree) != NULL && coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS){
			double* height_gradient = dvector(Tree_tip_count(coal->tree)-1);
			height_gradient_from_interval_gradient(coal, interval_gradient, height_gradient);
			Tree_node_transform_jvp(coal->tree, height_gradient, coal->gradient+offset);
			free(height_gradient);
		}
		// derivatives wrt to node heights
		else{
			height_gradient_from_interval_gradient(coal, interval_gradient, coal->gradient+offset);
		}
		free(interval_gradient);
	}
	free(chooses);
}

double _constant_calculate_dlogP( Coalescent* coal, const Parameter* p ){
	if(p != Parameters_at(coal->p, 0) && p->model != MODEL_TREE && p->model != MODEL_TREE_TRANSFORM) return 0;
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}

    double theta = Parameters_value(coal->p, 0);
    double dlogP = 0.0;
    if(p == Parameters_at(coal->p, 0)){
        double theta2 = theta*theta;
        for( int i = 0; i< coal->n; i++  ){
            dlogP += CHOOSE2(coal->lineages[i]) * coal->times[i] /theta2;
            
            if( coal->iscoalescent[i] ){
                dlogP -= 1.0/theta;
            }
            
        }
    }
    else{
        Node* node = Tree_node_from_parameter(coal->tree, p);
        size_t node_id = Node_id(node);
        size_t index = 0;
        double* proportion_derivatives = dvector(coal->n);
        for(; index < coal->n; index++){
            if(node_id == coal->nodes[index]->index){
                break;
            }
        }
        Parameters* reparams = get_reparams(coal->tree);
        double lower = Tree_lowers(coal->tree)[Node_id(node)];
        
        if (Node_isroot(node)) {
            proportion_derivatives[Node_id(node)] = 1;
        }
        else{
            proportion_derivatives[Node_id(node)] = Node_height(Node_parent(node)) - lower;
            dlogP = proportion_derivatives[Node_id(node)]*CHOOSE2(coal->lineages[index+1]);
        }
        _premultiply_proportions(node->left, proportion_derivatives, reparams);
        _premultiply_proportions(node->right, proportion_derivatives, reparams);
        
        index--;
        Node* parent = node;
        Node* n = Tree_node(coal->tree, coal->nodes[index]->index);
        while(index != 0){
            dlogP -= (proportion_derivatives[Node_id(parent)] - proportion_derivatives[Node_id(n)])*CHOOSE2(coal->lineages[index+1]);
            parent = n;
            index--;
            n = Tree_node(coal->tree, coal->nodes[index]->index);
        }
        free(proportion_derivatives);
        dlogP /= theta;
    }
	return dlogP;
}

double _constant_calculate_d2logP( Coalescent* coal, const Parameter* p ){
	if(p != Parameters_at(coal->p, 0)) return 0;
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	double theta = Parameters_value(coal->p, 0);
	double theta2 = theta*theta;
	double theta3 = theta2*theta;
	double d2logP = 0.0;
	for( int i = 0; i< coal->n; i++  ){
		d2logP -= CHOOSE2(coal->lineages[i]) * coal->times[i] *2.0/theta3;
		
		if( coal->iscoalescent[i] ){
			d2logP += 1.0/theta2;
		}
		
	}
	return d2logP;
}

double _constant_calculate_ddlogP( Coalescent* coal, const Parameter* p1, const Parameter* p2 ){
	if( p1 == p2 && p1 == Parameters_at(coal->p, 0) ) return _constant_calculate_d2logP(coal, p1);
	return 0.0;
}

Coalescent * create_GridCoalescent(demography type, size_t size, size_t grid_count){
	Coalescent *coal = (Coalescent*)malloc(sizeof(Coalescent));
	assert(coal);
	coal->p = NULL;
	coal->tree = NULL;
	coal->type = type;
	coal->calculate = NULL;
	coal->dlogP = NULL;
	coal->d2logP = NULL;
	coal->ddlogP = NULL;
	if(grid_count == 0){
		coal->update_intervals = _update_intervals;
	}
	else{
		coal->update_intervals = _update_intervals_grid;
	}
	coal->lineages = ivector(size + grid_count);
	coal->times = dvector(size + grid_count);
	coal->nodes = NULL;
	coal->stored_lineages = ivector(size + grid_count);
	coal->stored_times = dvector(size + grid_count);
	coal->iscoalescent = bvector(size + grid_count);
	coal->stored_iscoalescent = bvector(size + grid_count);
	coal->n = size + grid_count;
	coal->need_update_intervals = true;
	coal->need_update_gradient = true;
	coal->need_update = true;
	coal->grid = NULL;
	coal->gridCount = 0;
	coal->groups = NULL;
    coal->gradient = NULL;
	coal->parameterization = COALESCENT_THETA;
	
	return coal;
}

Coalescent * create_Coalescent(demography type, size_t size){
	return create_GridCoalescent(type, size, 0);
}

Coalescent * new_ConstantCoalescent_with_parameter( Parameter* theta, int size ){
	Coalescent *coal = create_Coalescent(COALESCENT_CONSTANT, size);
	coal->p = new_Parameters(1);
	Parameters_add(coal->p, theta);
	theta->model = MODEL_COALESCENT;
	coal->calculate = _constant_calculate;
	coal->dlogP = _constant_calculate_dlogP;
	coal->d2logP = _constant_calculate_d2logP;
	coal->ddlogP = _constant_calculate_ddlogP;
	return coal;
}

Coalescent * new_ConstantCoalescent( Tree* tree, Parameter* theta ){
	Coalescent *coal = new_ConstantCoalescent_with_parameter(theta, Tree_node_count(tree));
	coal->tree = tree;
	coal->nodes = malloc(Tree_node_count(tree)*sizeof(double_int_pair_t*));
	for (int i = 0; i < Tree_node_count(tree); i++) {
		coal->nodes[i] = malloc(sizeof(double_int_pair_t));
	}
	return coal;
}

Coalescent * new_ConstantCoalescent_with_data( Parameter* theta, double* times, bool* coalescent, int size ){
	Coalescent *coal = new_ConstantCoalescent_with_parameter(theta, size);
	memcpy(coal->times, times, sizeof(double)*size);
	memcpy(coal->iscoalescent, coalescent, sizeof(bool)*size);
	coalescentToLineage(coal);
	coal->need_update_intervals = false;
	return coal;
}

#pragma mark -
#pragma mark Exponential coalescent

double _coalescent_exponential_calculate( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	if ( coal->need_update ) {
		coal->logP = 0;
		double n0 = Parameters_value(coal->p, 0);
		double rate = Parameters_value(coal->p, 1);
		double start = 0;
		double logN0 = log(n0);
		double n0rate = n0*rate;
		
		for( int i = 0; i< coal->n; i++  ){
			double finish = start + coal->times[i];
			double integral;
			if(rate == 0.0){
				integral = (finish - start)/n0;
			}
			else{
				integral = (exp(finish*rate) - exp(start*rate))/n0rate;
			}
			
			coal->logP -= CHOOSE2(coal->lineages[i])*integral;
			
			if( coal->iscoalescent[i] ){
				// coal->logP -= log(n0*exp( -finish * rate));
				coal->logP -= logN0;
				if (rate != 0.0) {
					coal->logP += finish*rate;
				}
			}
			start = finish;
		}
		coal->need_update = false;
	}
	return coal->logP;
}

double _coalescent_exponential_dlogP( Coalescent* coal, const Parameter* p ){
	if(p != Parameters_at(coal->p, 0) && p != Parameters_at(coal->p, 1)) return 0;
	
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	
	double dlogP = 0;
	double n0 = Parameters_value(coal->p, 0);
	double rate = Parameters_value(coal->p, 1);
	double start = 0;
	double n02 = n0*n0;
	double rate2 = rate*rate;
	
	for( int i = 0; i< coal->n; i++  ){
		double finish = start + coal->times[i];
		double integral;
		if(rate == 0.0){
			if(Parameters_at(coal->p, 0) == p){
				integral = -(finish - start)/n02;
				dlogP -= CHOOSE2(coal->lineages[i])*integral;
			}
		}
		else{
			if(Parameters_at(coal->p, 0) == p){
				integral = -(exp(finish*rate) - exp(start*rate))/n02/rate;
			}
			else{
				integral = (exp(finish*rate)*(finish*rate -1.0) + exp(start*rate)*(1.0 - start*rate))/rate2/n0;
			}
			dlogP -= CHOOSE2(coal->lineages[i])*integral;
		}
		
		
		if( coal->iscoalescent[i] ){
			if(Parameters_at(coal->p, 0) == p){
				dlogP -= 1.0/n0;
			}
			else{
				dlogP += finish;
			}
		}
		start = finish;
	}
	
	return dlogP;
}

double _coalescent_exponential_d2logP( Coalescent* coal, const Parameter* p ){
	if(p != Parameters_at(coal->p, 0) && p != Parameters_at(coal->p, 1)) return 0;
	
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	
	double d2logP = 0;
	double n0 = Parameters_value(coal->p, 0);
	double rate = Parameters_value(coal->p, 1);
	double start = 0;
	double n02 = n0*n0;
	double n03 = n02*n0;
	double rate2 = rate*rate;
	double rate3 = rate2*rate;
	
	for( int i = 0; i< coal->n; i++  ){
		double finish = start + coal->times[i];
		double integral;
		if(rate == 0.0){
			if(Parameters_at(coal->p, 0) == p){
				integral = (finish - start)*2.0/n03;
				d2logP -= CHOOSE2(coal->lineages[i])*integral;
			}
		}
		else{
			if(Parameters_at(coal->p, 0) == p){
				integral = (exp(finish*rate) - exp(start*rate))*2.0/n03/rate;
			}
			else{
				integral = (finish*finish*rate2*exp(finish*rate) - 2.0*finish*rate*exp(finish*rate) + 2.0*exp(finish*rate) - start*start*rate2*exp(start*rate) + 2.0*start*rate*exp(start*rate) - (2.0*exp(start*rate)))/rate3/n0;
			}
			d2logP -= CHOOSE2(coal->lineages[i])*integral;
		}
		
		
		if( coal->iscoalescent[i] ){
			if(Parameters_at(coal->p, 0) == p){
				d2logP += -2.0/n03;
			}
		}
		start = finish;
	}
	
	return d2logP;
}

double _coalescent_exponential_ddlogP( Coalescent* coal, const Parameter* p1, const Parameter* p2 ){
	if(p1 == p2 && p1 == Parameters_at(coal->p, 0)) return _coalescent_exponential_d2logP(coal, p1);
	if(Parameters_value(coal->p, 1) == 0.0 ||
	   !((p1 == Parameters_at(coal->p, 0) && p2 ==Parameters_at(coal->p, 1)) ||
		 (p1 == Parameters_at(coal->p, 1) && p2 ==Parameters_at(coal->p, 0))) ) return 0;
	
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	
	double ddlogP = 0;
	double n0 = Parameters_value(coal->p, 0);
	double rate = Parameters_value(coal->p, 1);
	double start = 0;
	double n02 = n0*n0;
	double rate2 = rate*rate;
	
	for( int i = 0; i< coal->n; i++  ){
		double finish = start + coal->times[i];
		double integral = -(exp(finish*rate)*(finish*rate - 1.0) + exp(start*rate)*(1.0 - start*rate))/rate2/n02;
		
		ddlogP -= CHOOSE2(coal->lineages[i])*integral;
		start = finish;
	}
	
	return ddlogP;
}

Coalescent * new_ExponentialCoalescent_with_parameters( Parameters* parameters, int size ){
	Coalescent *coal = create_Coalescent(COALESCENT_EXPONENTIAL, size);
	coal->p = new_Parameters(2);
	Parameters_add(coal->p, Parameters_at(parameters, 0));
	Parameters_add(coal->p, Parameters_at(parameters, 1));
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameters_at(parameters, i)->model = MODEL_COALESCENT;
		Parameters_at(parameters, i)->id = i;
	}
	coal->calculate = _coalescent_exponential_calculate;
	coal->dlogP = _coalescent_exponential_dlogP;
	coal->d2logP = _coalescent_exponential_d2logP;
	coal->ddlogP = _coalescent_exponential_ddlogP;
	return coal;
}

Coalescent * new_ExponentialCoalescent( Tree *tree, Parameters* parameters ){
	Coalescent *coal = new_ExponentialCoalescent_with_parameters(parameters, Tree_node_count(tree));
	coal->tree = tree;
	coal->nodes = malloc(Tree_node_count(tree)*sizeof(double_int_pair_t*));
	for (int i = 0; i < Tree_node_count(tree); i++) {
		coal->nodes[i] = malloc(sizeof(double_int_pair_t));
	}
	return coal;
}

Coalescent * new_ExponentialCoalescent_with_data( Parameters* parameters, double* times, bool* coalescent, int size ){
	Coalescent *coal = new_ExponentialCoalescent_with_parameters(parameters, size);
	memcpy(coal->times, times, sizeof(double)*size);
	memcpy(coal->iscoalescent, coalescent, sizeof(bool)*size);
	coalescentToLineage(coal);
	coal->need_update_intervals = false;
	return coal;
}

#pragma mark -
#pragma mark classic skyline

double _coalescent_classical_skyline_calculate( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	if ( coal->need_update ) {
		coal->logP = 0;
		double popSize = -1;
		int index = 0;
		for( int i = 0; i< coal->n; i++  ){
			// t==0 for consecutive samling events
			if(coal->times[i] != 0.0){
				popSize = CHOOSE2(coal->lineages[i])/coal->times[i];
				Parameters_set_value(coal->p, index++, popSize);
				coal->logP -= coal->times[i]*coal->times[i];
			}
			
			if( coal->iscoalescent[i] ){
				coal->logP -= log(popSize);
			}
		}
		
		coal->need_update = false;
	}
	return coal->logP;
}

Coalescent * new_ClassicalSkylineCoalescent_with_parameters( Tree *tree, Parameters* parameters){
	Coalescent *coal = create_Coalescent(COALESCENT_SKYLINE_CLASSIC, Tree_node_count(tree));
	coal->p = new_Parameters(Parameters_count(parameters));
	Parameters_add_parameters(coal->p, parameters);
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameters_at(parameters, i)->model = MODEL_COALESCENT;
		Parameters_at(parameters, i)->id = i;
	}
	coal->tree = tree;
	coal->calculate = _coalescent_classical_skyline_calculate;
	coal->dlogP = NULL;
	coal->d2logP = NULL;
	coal->ddlogP = NULL;
	for (int i = 0; i < Tree_node_count(tree); i++) {
		coal->nodes[i] = malloc(sizeof(double_int_pair_t));
	}
	return coal;
}

#pragma mark -
#pragma mark skyride

double _coalescent_skyride_calculate_real_space( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	if ( coal->need_update ) {
		coal->logP = 0;
		double popSize = INFINITY;
		double logPopSize = 0;
		int index = 0;
		for( int i = 0; i< coal->n; i++  ){
			// t==0 for consecutive samling events
			if(coal->times[i] != 0.0){
				popSize = Parameters_value(coal->p, index);
				logPopSize = log(popSize);
				coal->logP -= coal->times[i]*CHOOSE2(coal->lineages[i])/popSize;
			}
			
			if( coal->iscoalescent[i] ){
				coal->logP -= logPopSize;
				index++;
			}
		}
		
		coal->need_update = false;
	}
	return coal->logP;
}

double _coalescent_skyride_calculate( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	if ( coal->need_update ) {
		coal->logP = 0;
		double logPopSize;
		double popSize;
		int index = 0;
		for( int i = 0; i< coal->n; i++  ){
			// t==0 for consecutive samling events
			if(coal->times[i] != 0.0){
				logPopSize = Parameters_value(coal->p, index);
				popSize = exp(logPopSize);
				if( coal->iscoalescent[i] ){
					coal->logP -= logPopSize;
					index++;
				}
				coal->logP -= coal->times[i]*CHOOSE2(coal->lineages[i])/popSize;
			}
			
		}
		coal->need_update = false;
	}
	return coal->logP;
}

// parameterized as \Delta_i = t_{i+1} - t_i
double _coalescent_skyride_calculate_deltas( Coalescent* coal ){
    if ( coal->need_update_intervals ) {
        _update_intervals(coal);
    }
    size_t parameter_count = Parameters_count(coal->p);
    double zgam = Parameters_value(coal->p, parameter_count-2);
    double tau = Parameters_value(coal->p, parameter_count-1);
    double zeta = 0.015;
    double gam = zgam/tau;
    
    if ( coal->need_update ) {
        coal->logP = 0;
        double logPopSize = NAN;
        double popSize = NAN;
        int index = 0;
        for( int i = 0; i< coal->n; i++  ){
            // t==0 for consecutive samling events
            if(coal->times[i] != 0.0){
                if (index == 0) {
                    popSize = Parameters_value(coal->p, 0);
                }
                else{
                    popSize = exp(zeta*gam*Parameters_value(coal->p, index) + logPopSize);
                }
                logPopSize = log(popSize);
                if( coal->iscoalescent[i] ){
                    coal->logP -= logPopSize;
                    index++;
                }
                coal->logP -= coal->times[i]*CHOOSE2(coal->lineages[i])/popSize;
            }
            
        }
        coal->need_update = false;
    }
    return coal->logP;
}

void _skyride_calculate_gradient_real( Coalescent* coal ){
	size_t offset = 0;
	double* chooses = dvector(coal->n);
	double invPop = 1.0/Parameters_value(coal->p, 0);
	size_t index = 0;
	size_t i = 1;
	for( ; i < coal->n - 1; i++){
		chooses[i] = CHOOSE2(coal->lineages[i])*invPop;
		if (coal->iscoalescent[i]) {
			invPop = 1.0/Parameters_value(coal->p, ++index);
		}
	}
	chooses [i] = CHOOSE2(coal->lineages[i])*invPop;

	if(coal->prepared_gradient & GRADIENT_FLAG_COALESCENT_THETA){
		memset(coal->gradient, 0, Parameters_count(coal->p)*sizeof(double));
		size_t i = 1;
		index = 0;
		double invPop = 1.0/Parameters_value(coal->p, 0);
		for( ; i < coal->n - 1; i++){
			coal->gradient[index] += coal->times[i]*chooses[i]*invPop;
			if (coal->iscoalescent[i]) {
				coal->gradient[index] -= invPop;
				invPop = 1.0/Parameters_value(coal->p, ++index);
			}
		}
		coal->gradient[index] += invPop*(coal->times[i]*chooses[i] - 1.0);
		offset += Parameters_count(coal->p);
	}

	if(coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS
	   || coal->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS){
		double* interval_gradient = dvector(coal->n);
		
		for(size_t i = 1; i < coal->n; i++){
			interval_gradient[i] = -chooses[i];
		}
		// derivatives wrt to reparameterization
		if(get_reparams(coal->tree) != NULL && coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS){
			double* height_gradient = dvector(Tree_tip_count(coal->tree)-1);
			height_gradient_from_interval_gradient(coal, interval_gradient, height_gradient);
			Tree_node_transform_jvp(coal->tree, height_gradient, coal->gradient+offset);
			free(height_gradient);
		}
		// derivatives wrt to node heights
		else{
			height_gradient_from_interval_gradient(coal, interval_gradient, coal->gradient+offset);
		}
		free(interval_gradient);
		
	}
	free(chooses);
}

void _skyride_calculate_gradient( Coalescent* coal ){
	size_t offset = 0;
	double* chooses = dvector(coal->n);
	double mexpPop = exp(-Parameters_value(coal->p, 0));
	size_t index = 0;
	size_t i = 1;
	for( ; i < coal->n - 1; i++){
		chooses[i] = CHOOSE2(coal->lineages[i])*mexpPop;
		if (coal->iscoalescent[i]) {
			mexpPop = exp(-Parameters_value(coal->p, ++index));
		}
	}
	chooses [i] = CHOOSE2(coal->lineages[i])*exp(-Parameters_value(coal->p, index));

    if(coal->prepared_gradient & GRADIENT_FLAG_COALESCENT_THETA){
		memset(coal->gradient, 0, Parameters_count(coal->p)*sizeof(double));
		size_t i = 1;
		index = 0;
        for( ; i < coal->n - 1; i++){
			coal->gradient[index] += coal->times[i]*chooses[i];
			if (coal->iscoalescent[i]) {
				coal->gradient[index++] -= 1.0;
			}
        }
		coal->gradient[index] += coal->times[i]*chooses[i] - 1.0;
		offset += Parameters_count(coal->p);
    }
	// -log(theta) - k/theta t
	// -ltheta - k/exp(ltheta) t
	if(coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS
	   || coal->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS){
		double* interval_gradient = dvector(coal->n);
		for(size_t i = 1; i < coal->n; i++){
			interval_gradient[i] = -chooses[i];
		}
		// derivatives wrt to reparameterization
		if(get_reparams(coal->tree) != NULL && coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS){
			double* height_gradient = dvector(Tree_tip_count(coal->tree)-1);
			height_gradient_from_interval_gradient(coal, interval_gradient, height_gradient);
			Tree_node_transform_jvp(coal->tree, height_gradient, coal->gradient+offset);
			free(height_gradient);
		}
		// derivatives wrt to node heights
		else{
			height_gradient_from_interval_gradient(coal, interval_gradient, coal->gradient+offset);
		}
		free(interval_gradient);
		
	}
	free(chooses);
}

double _skyride_calculate_dlogP( Coalescent* coal, const Parameter* p ){
	int ii = 0;
    if(p->model == MODEL_TREE){
        
    }
    else{
        for(; ii < Parameters_count(coal->p); ii++){
            if(Parameters_at(coal->p, ii) == p){
                break;
            }
        }
        if(ii == Parameters_count(coal->p)) return 0;
    }
    
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	double dlogP = 0.0;
    if(p->model != MODEL_TREE){
        int index = 0;
        double mexpPop = exp(-Parameters_value(coal->p, ii));
        for( int i = 0; i< coal->n; i++  ){
            // t==0 for consecutive samling events
            if(coal->times[i] != 0.0){
                if(index == ii){
                    dlogP += coal->times[i]*CHOOSE2(coal->lineages[i])*mexpPop;
                }
            }
            
            if(coal->iscoalescent[i]){
                if( index == ii){
                    dlogP -= 1.0;
                    break;
                }
                index++;
            }
        }
    }
    else{
        Node* node = Tree_node_from_parameter(coal->tree, p);
        size_t node_id = Node_id(node);
        size_t index = 0;
        double lower = 0;
        for(; index < coal->n; index++){
            if(node_id == coal->nodes[index]->index){
                break;
            }
        }
        // another tree
        if(index == coal->n) return 0;
        
        double* proportion_derivatives = dvector(coal->n);
        Parameters* reparams = get_reparams(coal->tree);
        lower = Tree_lowers(coal->tree)[Node_id(node)];
        size_t theta_index = 0;
        for(int i = 0; i < index; i++){
            if(coal->iscoalescent[i]){
                theta_index++;
            }
        }
        
        if (Node_isroot(node)) {
            proportion_derivatives[Node_id(node)] = 1;
        }
        else{
            proportion_derivatives[Node_id(node)] = Node_height(Node_parent(node)) - lower;
            dlogP = proportion_derivatives[Node_id(node)]*choose(coal->lineages[index+1], 2)/exp(Parameters_value(coal->p, theta_index+1));
        }
        _premultiply_proportions(node->left, proportion_derivatives, reparams);
        _premultiply_proportions(node->right, proportion_derivatives, reparams);
        index--;
        Node* parent = node;
        Node* n = Tree_node(coal->tree, coal->nodes[index]->index);
        while(index != 0){
            dlogP -= (proportion_derivatives[Node_id(parent)] - proportion_derivatives[Node_id(n)])*choose(coal->lineages[index+1], 2)/exp(Parameters_value(coal->p, theta_index));
            parent = n;
            index--;
            n = Tree_node(coal->tree, coal->nodes[index]->index);
            // t==0 for consecutive samling events
            if(coal->times[index+1] != 0.0){
                theta_index--;
            }
        }
        free(proportion_derivatives);
    }
	
	return dlogP;
}


double _skyride_calculate_d2logP( Coalescent* coal, const Parameter* p ){
	int ii = 0;
	for(; ii < Parameters_count(coal->p); ii++){
		if(Parameters_at(coal->p, ii) == p){
			break;
		}
	}
	if(ii == Parameters_count(coal->p)) return 0;
	
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	
	double popSize = exp(Parameters_value(coal->p, ii));
	int index = 0;
	double dlogP = 0;
	for( int i = 0; i< coal->n; i++  ){
		// t==0 for consecutive samling events
		if(coal->times[i] != 0.0){
			if(index == ii){
				dlogP -= coal->times[i]*CHOOSE2(coal->lineages[i])/popSize;
			}
		}
		
		if(coal->iscoalescent[i]){
			if( index == ii){
				break;
			}
			index++;
		}
	}
	return dlogP;
}

double _coalescent_skyride_ddlogP( Coalescent* coal, const Parameter* p1, const Parameter* p2 ){
if(p1 == p2) return _skyride_calculate_d2logP(coal, p1);
	return 0;
}

Coalescent * create_SkyrideCoalescent_with_parameters(Parameters* parameters, int size, parameterization_t parameterization){
	Coalescent *coal = create_Coalescent(COALESCENT_SKYRIDE, size);
	coal->p = new_Parameters(Parameters_count(parameters));
	Parameters_add_parameters(coal->p, parameters);
	
	coal->parameterization = parameterization;
	if(parameterization == COALESCENT_THETA_LOG){
		coal->calculate = _coalescent_skyride_calculate;
	}
	else if(parameterization == COALESCENT_THETA){
		coal->calculate = _coalescent_skyride_calculate_real_space;
	}
	else if(parameterization == COALESCENT_THETA_DELTA){
		coal->calculate = _coalescent_skyride_calculate_deltas;
	}
	
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameters_at(parameters, i)->model = MODEL_COALESCENT;
		Parameters_at(parameters, i)->id = i;
	}
	//TODO: size is never equal to param count
    if (Parameters_count(parameters) != size) {
//        coal->calculate = _coalescent_skyride_calculate_deltas;
    }
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameters_at(parameters, i)->model = MODEL_COALESCENT;
	}
	coal->dlogP = _skyride_calculate_dlogP;
	coal->d2logP = _skyride_calculate_d2logP;
	coal->ddlogP = _coalescent_skyride_ddlogP;
	return coal;
}

Coalescent * new_SkyrideCoalescent( Tree *tree, Parameters* parameters, parameterization_t parameterization){
	Coalescent *coal =create_SkyrideCoalescent_with_parameters(parameters, Tree_node_count(tree), parameterization);
	coal->tree = tree;
	coal->nodes = malloc(Tree_node_count(tree)*sizeof(double_int_pair_t*));
	for (int i = 0; i < Tree_node_count(tree); i++) {
		coal->nodes[i] = malloc(sizeof(double_int_pair_t));
	}
	return coal;
}

Coalescent * new_SkyrideCoalescent_with_data(Parameters* parameters, double* times, bool* coalescent, int size, parameterization_t parameterization){
	Coalescent *coal = create_SkyrideCoalescent_with_parameters(parameters, size, parameterization);
	memcpy(coal->times, times, sizeof(double)*size);
	memcpy(coal->iscoalescent, coalescent, sizeof(bool)*size);
	coalescentToLineage(coal);
	coal->need_update_intervals = false;
	return coal;
}

#pragma mark -
#pragma mark skygrid

double _coalescent_grid_calculate( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		coal->update_intervals(coal);
	}
	
	if ( coal->need_update ) {
		coal->logP = 0;
		size_t currentGridIndex = 0;
		double popSize = Parameters_value(coal->p, 0);
		double logPopSize = log(popSize);
		double lchoose2;
		
		for(size_t i = 0; i < coal->n; i++){
			if(coal->times[i] != 0.0){
				lchoose2 = CHOOSE2(coal->lineages[i]);
				coal->logP -= coal->times[i]*lchoose2/popSize;
				if( coal->iscoalescent[i] && coal->nodes[i]->index >= 0 ){
					coal->logP -= logPopSize;
				}
				else if (coal->nodes[i]->index < 0 && currentGridIndex < Parameters_count(coal->p) - 1) {
					popSize = Parameters_value(coal->p, ++currentGridIndex);
					logPopSize = log(popSize);
				}
			}
		}
		coal->need_update = false;
	}
	return coal->logP;
}

double _coalescent_grid_calculate_log_space( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		coal->update_intervals(coal);
	}
	
	if ( coal->need_update ) {
		coal->logP = 0;
		size_t currentGridIndex = 0;
		double logPopSize = Parameters_value(coal->p, 0);
		double popSize = exp(logPopSize);
		double lchoose2;
		
		for(size_t i = 0; i < coal->n; i++){
			if(coal->times[i] != 0.0){
				lchoose2 = CHOOSE2(coal->lineages[i]);
				coal->logP -= coal->times[i]*lchoose2/popSize;
				if( coal->iscoalescent[i] && coal->nodes[i]->index >= 0 ){
					coal->logP -= logPopSize;
				}
				else if (coal->nodes[i]->index < 0 &&currentGridIndex < Parameters_count(coal->p) - 1) {
					logPopSize = Parameters_value(coal->p, ++currentGridIndex);
					popSize = exp(logPopSize);
				}
			}
		}
		coal->need_update = false;
	}
	return coal->logP;
}

double _coalescent_grid_calculate_dlogP_space( Coalescent* coal, const Parameter* p ){
	int ii = 0;
	for(; ii < Parameters_count(coal->p); ii++){
		if(Parameters_at(coal->p, ii) == p){
			break;
		}
	}
	if(ii == Parameters_count(coal->p)) return 0;
	
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	
	double start = 0;
	size_t currentGridIndex = 0;
	
	double finish;
	double lchoose2;
	double dlogP = 0;
	
	for( int i = 0; i< coal->n; i++  ){
		finish = start + coal->times[i];
		
		if(coal->times[i] != 0.0){
			lchoose2 = CHOOSE2(coal->lineages[i]);
			
			// grid splits an interval
			while(currentGridIndex < coal->gridCount-1 && finish > coal->grid[currentGridIndex]){
				double end = fmin(coal->grid[currentGridIndex], finish);
				if(currentGridIndex == ii){
					coal->need_update = false;
					return dlogP + (end - start)*lchoose2*Parameters_value(coal->p, currentGridIndex);
				}
				
				start = end;
				
				if(currentGridIndex < coal->gridCount-1){
					currentGridIndex++;
				}
			}
			if(currentGridIndex == ii){
				dlogP += (finish - start)*lchoose2*Parameters_value(coal->p, ii);
				if(coal->iscoalescent[i]){
					dlogP -= 1.0;
				}
			}
		}
		
		start = finish;
	}
	coal->need_update = false;
	
	return dlogP;
}

void _skygrid_calculate_gradient( Coalescent* coal ){
	size_t offset = 0;
	double* chooses = dvector(coal->n);
	double mexpPop = exp(-Parameters_value(coal->p, 0));
	size_t index = 0;
	for(size_t i = 0; i < coal->n; i++){
		chooses [i] = CHOOSE2(coal->lineages[i])*mexpPop;
		if ( coal->nodes[i]->index < 0 && index < Parameters_count(coal->p)-1) {
			mexpPop = exp(-Parameters_value(coal->p, ++index));
		}
	}

    if(coal->prepared_gradient & GRADIENT_FLAG_COALESCENT_THETA){
		memset(coal->gradient, 0, Parameters_count(coal->p)*sizeof(double));
		size_t index = 0;
		for(size_t i = 0; i < coal->n; i++){
			coal->gradient[index] += coal->times[i]*chooses[i];
			if (coal->iscoalescent[i] && coal->nodes[i]->index >= 0) {
				coal->gradient[index] -= 1.0;
			}
			else if (coal->nodes[i]->index < 0) {
				index++;
			}
		}
		offset += Parameters_count(coal->p);
    }
	
	if(coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS
	   || coal->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS){
		double* interval_gradient = dvector(coal->n);
		for(size_t i = 1; i < coal->n; i++){
			interval_gradient[i] = -chooses[i];
		}

		// derivatives wrt to reparameterization
		if(get_reparams(coal->tree) != NULL && coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS){
			double* height_gradient = dvector(Tree_tip_count(coal->tree)-1);
			height_gradient_from_interval_gradient(coal, interval_gradient, height_gradient);
			Tree_node_transform_jvp(coal->tree, height_gradient, coal->gradient+offset);
			free(height_gradient);
		}
		else{
			height_gradient_from_interval_gradient(coal, interval_gradient, coal->gradient+offset);
		}
		free(interval_gradient);
	}
	free(chooses);
}

void _skygrid_calculate_gradient_real( Coalescent* coal ){
	size_t offset = 0;
	double* chooses = dvector(coal->n);
	double invPop = 1.0/Parameters_value(coal->p, 0);
	size_t index = 0;
	for(size_t i = 0; i < coal->n; i++){
		chooses [i] = CHOOSE2(coal->lineages[i])*invPop;
		if ( coal->nodes[i]->index < 0 && index < Parameters_count(coal->p)-1) {
			invPop = 1.0/Parameters_value(coal->p, ++index);
		}
	}

	if(coal->prepared_gradient & GRADIENT_FLAG_COALESCENT_THETA){
		memset(coal->gradient, 0, Parameters_count(coal->p)*sizeof(double));
		size_t index = 0;
		size_t i = 0;
		double invPop = 1.0/Parameters_value(coal->p, 0);
		for(; i < coal->n; i++){
			coal->gradient[index] += coal->times[i]*chooses[i]*invPop;
			if (coal->iscoalescent[i] && coal->nodes[i]->index >= 0) {
				coal->gradient[index] -= invPop;
			}
			else if (coal->nodes[i]->index < 0) {
				invPop = 1.0/Parameters_value(coal->p, ++index);
			}
		}
		// the last grid is empty (i.e. cutoff older than the root)
//		if(chooses[i] != 0.0){
//			coal->gradient[index] += invPop*(coal->times[i]*chooses[i] - 1.0);
//		}
		offset += Parameters_count(coal->p);
	}
	
	if(coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS
	   || coal->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS){
		double* interval_gradient = dvector(coal->n);
		for(size_t i = 1; i < coal->n; i++){
			interval_gradient[i] = -chooses[i];
		}
		// derivatives wrt to reparameterization
		if(get_reparams(coal->tree) != NULL && coal->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS){
			double* height_gradient = dvector(Tree_tip_count(coal->tree)-1);
			height_gradient_from_interval_gradient(coal, interval_gradient, height_gradient);
			Tree_node_transform_jvp(coal->tree, height_gradient, coal->gradient+offset);
			free(height_gradient);
		}
		// derivatives wrt to node heights
		else{
			height_gradient_from_interval_gradient(coal, interval_gradient, coal->gradient+offset);
		}
		free(interval_gradient);
	}
	free(chooses);
}

double _coalescent_grid_calculate_dlogP( Coalescent* coal, const Parameter* p ){
	int ii = 0;
	for(; ii < Parameters_count(coal->p); ii++){
		if(Parameters_at(coal->p, ii) == p){
			break;
		}
	}
	if(ii == Parameters_count(coal->p)) return 0;
	
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	
	double start = 0;
	size_t currentGridIndex = 0;
	
	double finish;
	double lchoose2;
	double dlogP = 0;
	
	for( int i = 0; i< coal->n; i++  ){
		finish = start + coal->times[i];
		
		if(coal->times[i] != 0.0){
			lchoose2 = CHOOSE2(coal->lineages[i]);
			
			// grid splits an interval
			while(currentGridIndex < coal->gridCount-1 && finish > coal->grid[currentGridIndex]){
				double end = fmin(coal->grid[currentGridIndex], finish);
				//if(currentGridIndex == ii) dlogP += (end - start)*lchoose2*exp(-Parameters_value(coal->p, currentGridIndex));
				if(currentGridIndex == ii){
					coal->need_update = false;
					return dlogP + (end - start)*lchoose2*exp(-Parameters_value(coal->p, currentGridIndex));
				}
				
				start = end;
				
				if(currentGridIndex < coal->gridCount-1){
					currentGridIndex++;
				}
			}
			if(currentGridIndex == ii){
				dlogP += (finish - start)*lchoose2*exp(-Parameters_value(coal->p, ii));
				if(coal->iscoalescent[i]){
					dlogP -= 1.0;
				}
			}
		}
		
		start = finish;
	}
	coal->need_update = false;
	
	return dlogP;
}

Coalescent * create_GridCoalescent_with_parameters(Parameters* parameters, int size, int grid, double cutoff, parameterization_t parameterization){
	Coalescent *coal = create_GridCoalescent(COALESCENT_SKYGRID, size, grid - 1);
	assert(coal);
	coal->grid = dvector(grid-1); // not 0 at grid[0]
	coal->gridCount = grid-1; // number of lines-1
	coal->p = new_Parameters(Parameters_count(parameters));
	Parameters_add_parameters(coal->p, parameters);
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameters_at(parameters, i)->model = MODEL_COALESCENT;
		Parameters_at(parameters, i)->id = i;
	}
	coal->parameterization = parameterization;
	if(parameterization == COALESCENT_THETA_LOG){
		coal->calculate = _coalescent_grid_calculate_log_space;
	}
	else if(parameterization == COALESCENT_THETA){
		coal->calculate = _coalescent_grid_calculate;
	}
	coal->dlogP = _coalescent_grid_calculate_dlogP_space;
	coal->d2logP = NULL;
	coal->ddlogP = NULL;
	for(int i = 1; i <= coal->gridCount; i++){
		coal->grid[i-1] = cutoff*i/coal->gridCount;
	}
	return coal;
}


Coalescent * new_GridCoalescent( Tree *tree, Parameters* parameters, int grid, double cutoff, parameterization_t parameterization ){
	Coalescent *coal = create_GridCoalescent_with_parameters(parameters, Tree_node_count(tree), grid, cutoff, parameterization);
	coal->tree = tree;
	coal->nodes = malloc(coal->n*sizeof(double_int_pair_t*));
	for (int i = 0; i < coal->n; i++) {
		coal->nodes[i] = malloc(sizeof(double_int_pair_t));
	}
	return coal;
}

Coalescent * new_GridCoalescent_with_data(Parameters* parameters, double* times, bool* coalescent, int size, int grid, double cutoff, parameterization_t parameterization){
	Coalescent *coal = create_GridCoalescent_with_parameters(parameters, size, grid, cutoff, parameterization);
	memcpy(coal->times, times, sizeof(double)*size);
	memcpy(coal->iscoalescent, coalescent, sizeof(bool)*size);
	coalescentToLineage(coal);
	coal->need_update_intervals = false;
	return coal;
}

#pragma mark -
#pragma mark skyline

double _coalescent_skyline_calculate( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		_update_intervals(coal);
	}
	if ( coal->need_update ) {
		coal->logP = 0;
		size_t currentGroupIndex = 0;
		
		size_t cum = coal->groups->values[currentGroupIndex];
		size_t coalescentCount = 0;
		double popSize = Parameters_value(coal->p, currentGroupIndex);
		double logPopSize = log(popSize);
		
		for( int i = 0; i< coal->n; i++  ){
			if (coalescentCount == cum) {
				currentGroupIndex++;
				cum += coal->groups->values[currentGroupIndex];
				popSize = Parameters_value(coal->p, currentGroupIndex);
				logPopSize = log(popSize);
			}
//			printf("%d %d %f\n", i, currentGroupIndex, Parameters_value(coal->p, currentGroupIndex));
			if(coal->times[i] != 0.0){
				coal->logP -= coal->times[i]*CHOOSE2(coal->lineages[i])/popSize;
			}
			
			if( coal->iscoalescent[i] ){
				coal->logP -= logPopSize;
				coalescentCount++;
			}
		}

		coal->need_update = false;
	}
//	printf("%f\n",coal->logP);
	return coal->logP;
}

Coalescent * new_SkylineCoalescent_with_parameters( Parameters* parameters, int size, DiscreteParameter* groups ){
	Coalescent *coal = create_Coalescent(COALESCENT_SKYLINE, size);
	coal->groups = groups;
	coal->p = new_Parameters(Parameters_count(parameters));
	Parameters_add_parameters(coal->p, parameters);
	for (int i = 0; i < Parameters_count(parameters); i++) {
		Parameters_at(parameters, i)->model = MODEL_COALESCENT;
		Parameters_at(parameters, i)->id = i;
	}
	coal->calculate = _coalescent_skyline_calculate;
	coal->dlogP = _coalescent_exponential_dlogP;
	coal->d2logP = _coalescent_exponential_d2logP;
	coal->ddlogP = _coalescent_exponential_ddlogP;
	return coal;
}

Coalescent * new_SkylineCoalescent( Tree *tree, Parameters* parameters, DiscreteParameter* groups ){
	Coalescent *coal = new_SkylineCoalescent_with_parameters(parameters, Tree_node_count(tree), groups);
	coal->tree = tree;
	coal->nodes = malloc(Tree_node_count(tree)*sizeof(double_int_pair_t*));
	for (int i = 0; i < Tree_node_count(tree); i++) {
		coal->nodes[i] = malloc(sizeof(double_int_pair_t));
	}
	return coal;
}

Coalescent * new_SkylineCoalescent_with_data( Parameters* parameters, double* times, bool* coalescent, int size, DiscreteParameter* groups ){
	Coalescent *coal = new_SkylineCoalescent_with_parameters(parameters, size, groups);
	memcpy(coal->times, times, sizeof(double)*size);
	memcpy(coal->iscoalescent, coalescent, sizeof(bool)*size);
	coalescentToLineage(coal);
	coal->need_update_intervals = false;
	return coal;
}
