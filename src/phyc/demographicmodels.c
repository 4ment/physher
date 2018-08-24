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

#include "tree.h"
#include "mstring.h"
#include "matrix.h"
#include "combinatorics.h"


#pragma mark Coalescent

static void _coalescent_model_store(Model* self){
	self->storedLogP = self->lp;
}

static void _coalescent_model_restore(Model* self){
	self->lp = self->storedLogP;
}

static double _coalescent_model_logP(Model *self){
	Coalescent* mc = (Coalescent*)self->obj;
	self->lp = mc->calculate(mc);
	return self->lp;
}

static double _coalescent_model_dlogP(Model *self, const Parameter* p){
	error("can't do that _coalescent_model_dlogP\n");
	return 0;
}

static double _coalescent_model_d2logP(Model *self, const Parameter* p){
	error("can't do that _coalescent_model_d2logP\n");
	return 0;
}

static double _coalescent_model_ddlogP(Model *self, const Parameter* p1, const Parameter* p2){
	error("can't do that _coalescent_model_ddlogP\n");
	return 0;
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
	else if ( strcmp(model->type, "tree") == 0 ) {
		c->need_update_intervals = true;
		c->need_update = true;
	}
	else{
		fprintf(stderr, "_coalescent_model_handle_change\n");
		exit(1);
	}
	self->listeners->fire( self->listeners, self, index );
}

Model* new_CoalescentModel(const char* name, Coalescent* coalescent, Model* tree){
	Model* model = new_Model("distribution", name, coalescent);
	for ( int i = 0; i < Parameters_count(coalescent->p); i++ ) {
		Parameters_at(coalescent->p, i)->listeners->add( Parameters_at(coalescent->p, i)->listeners, model );
	}
	tree->listeners->add( tree->listeners, model );
	
	model->logP = _coalescent_model_logP;
	model->dlogP = _coalescent_model_dlogP;
	model->d2logP = _coalescent_model_d2logP;
	model->ddlogP = _coalescent_model_ddlogP;
	model->free = _coalescent_model_free;
	model->clone = _coalescent_model_clone;
	model->get_free_parameters = _coalescent_model_get_free_parameters;
	model->store = _coalescent_model_store;
	model->restore = _coalescent_model_restore;
	
	model->update = _coalescent_model_handle_change;
	
	model->data = tree;
	tree->ref_count++;
	return model;
}

Model* new_CoalescentModel_from_json(json_node* node, Hashtable* hash){
	char* model = get_json_node_value_string(node, "model");
	Coalescent* c = NULL;
	
	json_node* tree_node = get_json_node(node, "tree");
	char* ref = (char*)tree_node->value;
	// check it starts with a &
	Model* mtree = Hashtable_get(hash, ref+1);
	mtree->ref_count++;
	
	if(strcasecmp(model, "constant") == 0){
		json_node* p_node = get_json_node(node, "theta");
		Parameter* p = new_Parameter_from_json(p_node, hash);
		c = new_ConstantCoalescent_with_parameter(mtree->obj, p);
		Hashtable_add(hash, Parameter_name(p), p);
		free_Parameter(p);
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

    switch ( coalescent->type ) {
        case CONSTANT_DEMOGRAPHY:{
            free_ConstantCoalescent(coalescent);
            break;
        }
        default:
            assert(0);
    }
}

#pragma mark -
#pragma mark Constant coalescent

static double _constant_calculate( Coalescent* coal );
static void _update_intervals( Coalescent* coal );

Coalescent * new_ConstantCoalescent_with_parameter( Tree *tree, Parameter* theta ){
	Coalescent *coal = (Coalescent*)malloc(sizeof(Coalescent));
	assert(coal);
	coal->p = new_Parameters(1);
	Parameters_add(coal->p, theta);
	coal->tree = tree;
	coal->type = CONSTANT_DEMOGRAPHY;
	coal->calculate = _constant_calculate;
	coal->lineages = ivector(Tree_node_count(tree));
	coal->times = dvector(Tree_node_count(tree));
	coal->nodes = ivector(Tree_node_count(tree));
	coal->iscoalescent = bvector(Tree_node_count(tree));
	coal->n = Tree_node_count(tree);
	coal->need_update_intervals = true;
	coal->need_update = true;
	return coal;
}

void free_ConstantCoalescent( Coalescent *coal ){
	free_Parameters(coal->p);
    free(coal->lineages);
    free(coal->times);
    free(coal->iscoalescent);
	free(coal);
}

double get_demographic_constant( const Coalescent *coal, double time ){
	return Parameters_value(coal->p, 0);
}

double get_intensity_constant( const Coalescent *coal, double time ){
	return time/Parameters_value(coal->p, 0);
}

void _post_traversal( Node *node, double lower, double upper, int *lineages, bool *iscoalescent ){
    if( node == NULL ){
        return;
    }
    _post_traversal(node->left, lower, upper, lineages, iscoalescent);
    _post_traversal(node->right, lower, upper, lineages, iscoalescent);
    
    if( !Node_isleaf(node) && Node_height(node) >= upper ){
        if( Node_height( Node_left(node) ) <= lower ){
            //fprintf(stderr, "%s %f %f\n", Node_name(Node_left(node)), Node_height( Node_left(node) ) , lower);
            (*lineages)++;
        }
        if( Node_height( Node_right(node) ) <=lower ){
            //fprintf(stderr, "%s %f %f\n", Node_name(Node_right(node)), Node_height( Node_right(node) ) , lower);
            (*lineages)++;
        }
        if( Node_height(node) == upper ){
            *iscoalescent = true;
        }
    }
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
//		printf("%f %f\n", coal->times[i],coal->times[i-1]);
	}
	coal->lineages[0] = n+1;
	coal->iscoalescent[0] = true;
	coal->n = internalCount;
	coal->need_update_intervals = false;
}

void _update_intervals2( Coalescent* coal ){
    Node **nodes = Tree_get_nodes(coal->tree, POSTORDER);
    memset(coal->times, 0, Tree_node_count(coal->tree)*sizeof(double));
    memset(coal->iscoalescent, 0, Tree_node_count(coal->tree)*sizeof(bool));
    memset(coal->lineages, 0, Tree_node_count(coal->tree)*sizeof(int));
    
    int n = 0;
    for ( int i = 0; i < Tree_node_count(coal->tree); i++ ) {
        if( Node_isleaf(nodes[i]) && Node_height(nodes[i]) != 0.0 ){
            coal->times[n] = Node_height(nodes[i]);
            n++;
        }
    }
    qsort(coal->times, n, sizeof(double), qsort_asc_dvector);
	
    //print_dvector(coal->times, n);

    int len = Tree_tip_count(coal->tree);
    n = 1;
    for ( int i = 1; i < Tree_node_count(coal->tree)-1; i++ ) {
        if( Node_isleaf(nodes[i]) && Node_height(nodes[i]) != 0.0  ){
            if( coal->times[n] == coal->times[n-1] ){
                len--;
                memmove( coal->times+n, coal->times+n+1, len*sizeof(double) );
            }
            else {
                n++;
            }
        }
    }
   
    //print_dvector(coal->times, n);

    for ( int i = 0; i < Tree_node_count(coal->tree); i++ ) {
        if( !Node_isleaf(nodes[i]) ){
            coal->times[n] = Node_height(nodes[i]);
            n++;
        }
    }
    
    qsort(coal->times, n, sizeof(double), qsort_asc_dvector);
    //print_dvector(coal->times, n);
	
    _post_traversal(Tree_root(coal->tree), 0.0, coal->times[0], &coal->lineages[0], &coal->iscoalescent[0]);
    for ( int i = 1; i < n; i++ ) {
        _post_traversal(Tree_root(coal->tree), coal->times[i-1], coal->times[i], &coal->lineages[i], &coal->iscoalescent[i]);
    }
    //print_ivector(coal->lineages, n);

    for ( int i = n-1; i > 0; i-- ) {
        coal->times[i] = coal->times[i] - coal->times[i-1];
    }
    coal->n = n;
    coal->need_update_intervals = false;
}

double _constant_calculate( Coalescent* coal ){
    if ( coal->need_update_intervals ) {
        _update_intervals(coal);
    }
	if ( coal->need_update ) {
		coal->logP = 0;
		double theta = Parameters_value(coal->p, 0);
	 
		for( int i = 0; i< coal->n; i++  ){
			double lambda = choose(coal->lineages[i], 2)/theta;
			coal->logP -= lambda * coal->times[i];
			
			if( coal->iscoalescent[i] ){
				coal->logP -= log(theta);
			}
			
		}
		coal->need_update = false;
	}
    return coal->logP;
}


exponentialdemography * new_exponentialdemography( double n0, double r ){
	exponentialdemography *cd = (exponentialdemography*)malloc(sizeof(exponentialdemography));	
	assert(cd);
    cd->name = String_clone("exponential.demography");	
	cd->type = CONSTANT_DEMOGRAPHY;
	cd->n0 = n0;
	cd->r = r;
	cd->lower = 0;
	cd->upper = INFINITY;
	return cd;
}

void free_exponentialdemography( exponentialdemography *demo ){
	free(demo->name);
	free(demo);
}

double get_demographic_exp( exponentialdemography * demo, double time ){
	if( demo->r == 0. ) return demo->n0;
	else return demo->n0 * exp( -time * demo->r);
}

double get_intensity_exp( exponentialdemography * demo, double time ){
	if( demo->r == 0. ) return time/demo->n0;
	else return exp( time * demo->r)/demo->n0/demo->r;
}
