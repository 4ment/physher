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


#pragma mark Coalescent

static void _coalescent_model_store(Model* self){
	Coalescent* coalescent = (Coalescent*)self->obj;
	for (int i = 0; i < Parameters_count(coalescent->p); i++) {
		Parameter_store(Parameters_at(coalescent->p, i));
	}
	memcpy(coalescent->stored_iscoalescent, coalescent->iscoalescent, coalescent->n*sizeof(bool));
	memcpy(coalescent->stored_times, coalescent->times, coalescent->n*sizeof(double));
	memcpy(coalescent->stored_lineages, coalescent->lineages, coalescent->n*sizeof(int));
	self->storedLogP = self->lp;
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
	memcpy(coalescent->iscoalescent, coalescent->stored_iscoalescent, coalescent->n*sizeof(bool));
	memcpy(coalescent->times, coalescent->stored_times, coalescent->n*sizeof(double));
	memcpy(coalescent->lineages, coalescent->stored_lineages, coalescent->n*sizeof(int));
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
	mc->need_update_intervals = true;
	self->lp = mc->calculate(mc);
	return self->lp;
}

static double _coalescent_model_dlogP(Model *self, const Parameter* p){
	Coalescent* mc = (Coalescent*)self->obj;
	return mc->dlogP(mc, p);
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
		// tree or simplex
		if(self->data != NULL){
			Model* model = self->data;
			model->free(model);
		}
		free_Coalescent(mc);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

//MARK: TODO
static Model* _coalescent_model_clone( Model *self, Hashtable* hash ){
	error("todo _coalescent_model_clone\n");
	return NULL;
}

//MARK: TODO
static void _coalescent_model_get_free_parameters(Model* model, Parameters* parameters){
	error("todo _coalescent_model_get_free_parameters\n");
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
	}
	else if ( model->type == MODEL_TREE ) {
		c->need_update_intervals = true;
		c->need_update = true;
	}
	else{
		fprintf(stderr, "_coalescent_model_handle_change\n");
		exit(1);
	}
	self->listeners->fire( self->listeners, self, index );
}

static void _coalescent_model_handle_restore( Model *self, Model *model, int index ){
	Coalescent *c = (Coalescent*)self->obj;
	c->need_update = true;
	c->need_update_intervals = true;
	self->listeners->fire_restore( self->listeners, self, index );
}

Model* new_CoalescentModel(const char* name, Coalescent* coalescent, Model* tree){
	Model* model = new_Model(MODEL_COALESCENT, name, coalescent);
	for ( int i = 0; i < Parameters_count(coalescent->p); i++ ) {
		Parameters_at(coalescent->p, i)->listeners->add( Parameters_at(coalescent->p, i)->listeners, model );
	}
	tree->listeners->add( tree->listeners, model );
	
	model->logP = _coalescent_model_logP;
	model->full_logP = _coalescent_model_full_logP;
	model->dlogP = _coalescent_model_dlogP;
	model->d2logP = _coalescent_model_d2logP;
	model->ddlogP = _coalescent_model_ddlogP;
	model->free = _coalescent_model_free;
	model->clone = _coalescent_model_clone;
	model->get_free_parameters = _coalescent_model_get_free_parameters;
	model->store = _coalescent_model_store;
	model->restore = _coalescent_model_restore;
	model->sample = _coalescent_model_sample;
	model->sample_evaluate = _coalescent_model_sample_evaluate;
	
	model->update = _coalescent_model_handle_change;
	model->handle_restore = _coalescent_model_handle_restore;
	
	model->data = tree;
	tree->ref_count++;
	return model;
}

Model* new_CoalescentModel_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"model",
		"parameters",
		"tree"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* model = get_json_node_value_string(node, "model");
	Coalescent* c = NULL;
	
	json_node* tree_node = get_json_node(node, "tree");
	char* ref = (char*)tree_node->value;
	// check it starts with a &
	Model* mtree = Hashtable_get(hash, ref+1);
	mtree->ref_count++;
	
	if(strcasecmp(model, "constant") == 0){
		Parameters* ps = new_Parameters(1);
		get_parameters_references(node, hash, ps);
		c = new_ConstantCoalescent_with_parameter(mtree->obj, Parameters_at(ps, 0));
		Hashtable_add(hash, Parameters_name(ps, 0), Parameters_at(ps, 0));
		free_Parameters(ps);
	}
	else if(strcasecmp(model, "exponential") == 0){
		Parameters* ps = new_Parameters(2);
		get_parameters_references(node, hash, ps);
		if (strcasecmp(Parameters_name(ps, 0), "n0") != 0) {
			Parameters_swap_index(ps, 0, 1);
		}
		c = new_ExponentialCoalescent_with_parameters(mtree->obj, ps);
		Hashtable_add(hash, Parameters_name(ps, 0), Parameters_at(ps, 0));
		Hashtable_add(hash, Parameters_name(ps, 1), Parameters_at(ps, 1));
		free_Parameters(ps);
	}
	else{
		fprintf(stderr, "BranchModel type unknown %s\n", model);
		exit(1);
	}
	
	char* id = get_json_node_value_string(node, "id");
	Model* mc = new_CoalescentModel(id, c, mtree);
	mtree->free(mtree);
	return mc;
}

void free_Coalescent( Coalescent *coalescent ){
	free_Parameters(coalescent->p);
	free(coalescent->lineages);
	free(coalescent->nodes);
	free(coalescent->times);
	free(coalescent->iscoalescent);
	free(coalescent->stored_lineages);
	free(coalescent->stored_times);
	free(coalescent->stored_iscoalescent);
	free(coalescent);
}

#pragma mark -
#pragma mark Constant coalescent

static double _constant_calculate( Coalescent* coal );
static void _update_intervals( Coalescent* coal );
static void _update_intervals_heterochronous( Coalescent* coal );

double _constant_calculate_dlogP( Coalescent* coal, const Parameter* p ){
	if(p != Parameters_at(coal->p, 0)) return 0;
	if ( coal->need_update_intervals ) {
		if(Tree_homochronous(coal->tree)){
			_update_intervals(coal);
		}
		else{
			_update_intervals_heterochronous(coal);
		}
		coal->need_update = true;
	}
	double theta = Parameters_value(coal->p, 0);
	double theta2 = theta*theta;
	double dlogP = 0.0;
	for( int i = 0; i< coal->n; i++  ){
		dlogP += choose(coal->lineages[i], 2) * coal->times[i] /theta2;
		
		if( coal->iscoalescent[i] ){
			dlogP -= 1.0/theta;
		}
		
	}
	return dlogP;
}

double _constant_calculate_d2logP( Coalescent* coal, const Parameter* p ){
	if(p != Parameters_at(coal->p, 0)) return 0;
	if ( coal->need_update_intervals ) {
		if(Tree_homochronous(coal->tree)){
			_update_intervals(coal);
		}
		else{
			_update_intervals_heterochronous(coal);
		}
		coal->need_update = true;
	}
	double theta = Parameters_value(coal->p, 0);
	double theta2 = theta*theta;
	double theta3 = theta2*theta;
	double d2logP = 0.0;
	for( int i = 0; i< coal->n; i++  ){
		d2logP -= choose(coal->lineages[i], 2) * coal->times[i] *2.0/theta3;
		
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

Coalescent * new_ConstantCoalescent_with_parameter( Tree *tree, Parameter* theta ){
	Coalescent *coal = (Coalescent*)malloc(sizeof(Coalescent));
	assert(coal);
	coal->p = new_Parameters(1);
	Parameters_add(coal->p, theta);
	coal->tree = tree;
	coal->type = CONSTANT_DEMOGRAPHY;
	coal->calculate = _constant_calculate;
	coal->dlogP = _constant_calculate_dlogP;
	coal->d2logP = _constant_calculate_d2logP;
	coal->ddlogP = _constant_calculate_ddlogP;
	coal->lineages = ivector(Tree_node_count(tree));
	coal->times = dvector(Tree_node_count(tree));
	coal->nodes = ivector(Tree_node_count(tree));
	coal->stored_lineages = ivector(Tree_node_count(tree));
	coal->stored_times = dvector(Tree_node_count(tree));
	coal->iscoalescent = bvector(Tree_node_count(tree));
	coal->stored_iscoalescent = bvector(Tree_node_count(tree));
	coal->n = Tree_node_count(tree);
	coal->need_update_intervals = true;
	coal->need_update = true;
	return coal;
}

static void _sort_ascending_time( double *times, int *map, int size ){
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( times[i] > times[i+1] ) {
				done = false;
				dswap(&times[i], &times[i+1]);
				swap_int(&map[i], &map[i+1]);
			}
		}
		size--;
	}
}

void _update_intervals( Coalescent* coal ){
	int nodeCount = Tree_node_count(coal->tree);
	int internalCount = nodeCount - Tree_tip_count(coal->tree);
	int n = 0;
	for(int i = 0; i < nodeCount; i++){
		Node* node = Tree_node(coal->tree, i);
		if(!Node_isleaf(node)){
			coal->nodes[n] = Node_id(node);
			coal->times[n] = Node_height(node);
			n++;
		}
	}
	_sort_ascending_time(coal->times, coal->nodes, n);
	
	for ( int i = n-1; i > 0; i-- ) {
		coal->times[i] = coal->times[i] - coal->times[i-1];
		coal->iscoalescent[i] = true;
		coal->lineages[i] = n+1-i;
	}
	coal->lineages[0] = n+1;
	coal->iscoalescent[0] = true;
	coal->n = internalCount;
	coal->need_update_intervals = false;
}

// do not change order of nodes
static void _sort_asc_node_height( Node** nodes, int *order, size_t size ){
	
	for ( int i = 0; i < size; i++ ) {
		order[i] = i;
	}
	
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( Node_height(nodes[order[i]]) > Node_height(nodes[order[i+1]]) ) {
				done = false;
				swap_int(&order[i], &order[i+1]);
			}
		}
		size--;
	}
}

void _update_intervals_heterochronous(Coalescent* coal){
	Node **nodes = Tree_get_nodes(coal->tree, POSTORDER);
	memset(coal->times, 0, Tree_node_count(coal->tree)*sizeof(double));
	memset(coal->iscoalescent, 0, Tree_node_count(coal->tree)*sizeof(bool));
	memset(coal->lineages, 0, Tree_node_count(coal->tree)*sizeof(int));
	
	size_t nodeCount = Tree_node_count(coal->tree);
	
	_sort_asc_node_height(nodes, coal->nodes, nodeCount);

	Node* node = nodes[coal->nodes[0]];
	double start = Node_height(node);
	int lineageCount = 0;
	int nodeIndex = 0;
	int intervalCount = 0;
	double finish;
	
	while (nodeIndex < nodeCount) {
		node = nodes[coal->nodes[nodeIndex]];
		
		finish = Node_height(node);
		nodeIndex++;
		
		// sampling event
		if (Node_isleaf(node)) {
			if (intervalCount > 0) {
				coal->times[intervalCount] = finish - start;
				coal->lineages[intervalCount] = lineageCount;
				intervalCount++;
			}
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
}

double _constant_calculate( Coalescent* coal ){
    if ( coal->need_update_intervals ) {
		if(Tree_homochronous(coal->tree)){
			_update_intervals(coal);
		}
		else{
			_update_intervals_heterochronous(coal);
		}
		coal->need_update = true;
    }
	if ( coal->need_update ) {
		coal->logP = 0;
		double theta = Parameters_value(coal->p, 0);
		double logTheta = log(theta);
	 
		for( int i = 0; i< coal->n; i++  ){
			double lambda = choose(coal->lineages[i], 2)/theta;
			coal->logP -= lambda * coal->times[i];
			
			if( coal->iscoalescent[i] ){
				coal->logP -= logTheta;
			}
			
		}
		coal->need_update = false;
	}
    return coal->logP;
}

#pragma mark -
#pragma mark Exponential coalescent

double _coalescent_exponential_calculate( Coalescent* coal ){
	if ( coal->need_update_intervals ) {
		if(Tree_homochronous(coal->tree)){
			_update_intervals(coal);
		}
		else{
			_update_intervals_heterochronous(coal);
		}
		coal->need_update = true;
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
			
			coal->logP -= choose(coal->lineages[i], 2)*integral;
			
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
		if(Tree_homochronous(coal->tree)){
			_update_intervals(coal);
		}
		else{
			_update_intervals_heterochronous(coal);
		}
		coal->need_update = true;
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
				dlogP -= choose(coal->lineages[i], 2)*integral;
			}
		}
		else{
			if(Parameters_at(coal->p, 0) == p){
				integral = -(exp(finish*rate) - exp(start*rate))/n02/rate;
			}
			else{
				integral = (exp(finish*rate)*(finish*rate -1.0) + exp(start*rate)*(1.0 - start*rate))/rate2/n0;
			}
			dlogP -= choose(coal->lineages[i], 2)*integral;
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
		if(Tree_homochronous(coal->tree)){
			_update_intervals(coal);
		}
		else{
			_update_intervals_heterochronous(coal);
		}
		coal->need_update = true;
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
				d2logP -= choose(coal->lineages[i], 2)*integral;
			}
		}
		else{
			if(Parameters_at(coal->p, 0) == p){
				integral = (exp(finish*rate) - exp(start*rate))*2.0/n03/rate;
			}
			else{
				integral = (finish*finish*rate2*exp(finish*rate) - 2.0*finish*rate*exp(finish*rate) + 2.0*exp(finish*rate) - start*start*rate2*exp(start*rate) + 2.0*start*rate*exp(start*rate) - (2.0*exp(start*rate)))/rate3/n0;
			}
			d2logP -= choose(coal->lineages[i], 2)*integral;
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
		if(Tree_homochronous(coal->tree)){
			_update_intervals(coal);
		}
		else{
			_update_intervals_heterochronous(coal);
		}
		coal->need_update = true;
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
		
		ddlogP -= choose(coal->lineages[i], 2)*integral;
		start = finish;
	}
	
	return ddlogP;
}

Coalescent * new_ExponentialCoalescent_with_parameters( Tree *tree, Parameters* parameters ){
	Coalescent *coal = (Coalescent*)malloc(sizeof(Coalescent));
	assert(coal);
	coal->p = new_Parameters(2);
	Parameters_add(coal->p, Parameters_at(parameters, 0));
	Parameters_add(coal->p, Parameters_at(parameters, 1));
	coal->tree = tree;
	coal->type = EXP_DEMOGRAPHY;
	coal->calculate = _coalescent_exponential_calculate;
	coal->dlogP = _coalescent_exponential_dlogP;
	coal->d2logP = _coalescent_exponential_d2logP;
	coal->ddlogP = _coalescent_exponential_ddlogP;
	coal->lineages = ivector(Tree_node_count(tree));
	coal->times = dvector(Tree_node_count(tree));
	coal->nodes = ivector(Tree_node_count(tree));
	coal->stored_lineages = ivector(Tree_node_count(tree));
	coal->stored_times = dvector(Tree_node_count(tree));
	coal->iscoalescent = bvector(Tree_node_count(tree));
	coal->stored_iscoalescent = bvector(Tree_node_count(tree));
	coal->n = Tree_node_count(tree);
	coal->need_update_intervals = true;
	coal->need_update = true;
	return coal;
}
