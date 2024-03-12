/*
 *  treelikelihood.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/2/10.
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

#include "treelikelihood.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "tree.h"
#include "node.h"
#include "substmodel.h"
#include "sitepattern.h"

#include "matrix.h"
#include "parameters.h"
#include "mstring.h"
#include "mathconstant.h"


#include "hessenberg.h" // for det
#include "solve.h"

#include "treelikelihood4.h"
#include "treelikelihood20.h"
#include "treelikelihoodX.h"
#include "treelikelihoodCodon.h"

#ifdef SSE3_ENABLED
#if defined(__aarch64__)
#include "neon2sse.h"
#else
#include <xmmintrin.h> // SSE
#include <pmmintrin.h> // SSE3
#endif
#endif

#ifdef AVX_ENABLED
#include <immintrin.h>
#endif


// MARK:  Private function declaration


static bool _calculate_partials( SingleTreeLikelihood *tlk, Node *n  );
static double _calculate( SingleTreeLikelihood *tlk );

double _calculate_uppper( SingleTreeLikelihood *tlk, Node *node );

void TreeLikelihood_calculate_gradient( Model *model, double* grads );

void allocate_storage(SingleTreeLikelihood* tlk, size_t index);

#pragma mark -
#pragma mark TreeLikelihoodModel

void _treelikelihood_handle_change( Model *self, Model *model, int index ){
	SingleTreeLikelihood *tlk = (SingleTreeLikelihood*)self->obj;

	if ( model->type == MODEL_TREE ) {
		if(index == -1){
			SingleTreeLikelihood_update_all_nodes(tlk);
			// a reparameterized node height has changed but we wait for the node heights
			// to be updated to know which partials need to be recomputed
			// tlk->update = true;
			// tlk->update_upper = true;
			// printf("MODEL_TREE1\n");
		}
		else{
			// printf("MODEL_TREE\n");
			tlk->update_nodes[index] = true;
			tlk->update = true;
			tlk->update_upper = true;
		}
	}
	else if ( model->type == MODEL_BRANCHMODEL ) {
		if(index == -1){
			SingleTreeLikelihood_update_all_nodes(tlk);
			// printf("MODEL_BRANCHMODEL2\n");
		}
		else{
			// printf("MODEL_BRANCHMODEL\n");
			tlk->update_nodes[index] = true;
			tlk->update = true;
			tlk->update_upper = true;
		}
	}
	else if ( model->type == MODEL_SITEMODEL ) {
		SingleTreeLikelihood_update_all_nodes(tlk);
	}
	else if ( model->type == MODEL_SUBSTITUTION ) {
		SingleTreeLikelihood_update_all_nodes(tlk);
	}
	else {
		fprintf(stderr, "%s of type %s\n", model->name, model_type_strings[model->type]);
		error("Unknown change in SingleLikelihood\n");
	}
}

void _treelikelihood_handle_restore( Model *self, Model *model, int index ){
// parameters are restored of evry model is restored but models can be dirty
//	((SingleTreeLikelihood*)self->obj)->update = true;
	SingleTreeLikelihood* tlk = self->obj;
	memcpy(tlk->current_matrices_indexes, tlk->stored_matrices_indexes, sizeof(unsigned)*Tree_node_count(tlk->tree)*2);
	memcpy(tlk->current_partials_indexes, tlk->stored_partials_indexes, sizeof(unsigned)*Tree_node_count(tlk->tree)*2);
	tlk->lk = tlk->stored_lk;
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _singleTreeLikelihood_store(Model* self){
	SingleTreeLikelihood* tlk = self->obj;
	if(tlk->partials[1] == NULL){
		allocate_storage(tlk, 1);
		if (!tlk->use_tip_states) {
			Node **nodes = Tree_get_nodes( tlk->tree, POSTORDER );
			for (size_t i = 0; i < Tree_node_count(tlk->tree); i++) {
				if(Node_isleaf(nodes[i])){
					memcpy(tlk->partials[1][Node_id(nodes[i])], tlk->partials[0][Node_id(nodes[i])], sizeof(double)*tlk->partials_size);
				}
			}
		}
	}
	memcpy(tlk->stored_matrices_indexes, tlk->current_matrices_indexes, sizeof(unsigned)*Tree_node_count(tlk->tree)*2);
	memcpy(tlk->stored_partials_indexes, tlk->current_partials_indexes, sizeof(unsigned)*Tree_node_count(tlk->tree)*2);
	Model** models = (Model**)self->data;
	models[0]->store(models[0]); // tree
	models[1]->store(models[1]); // substitutionmodel
	models[2]->store(models[2]); // sitemodel
	if(models[3] != NULL){
		models[3]->store(models[3]); // branchmodel
	}
	self->storedLogP = self->lp;
	tlk->stored_lk = tlk->lk;
}

static void _singleTreeLikelihood_restore(Model* self){
	Model** models = (Model**)self->data;
	models[0]->restore(models[0]); // tree
	models[1]->restore(models[1]); // substitutionmodel
	models[2]->restore(models[2]); // sitemodel
	if(models[3] != NULL){
		models[3]->restore(models[3]); // branchmodel
	}
	self->lp = self->storedLogP;
}

double _singleTreeLikelihood_logP(Model *self){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	self->lp = tlk->calculate(tlk);
	if(tlk->include_jacobian){
		Model** models = (Model**)self->data;
		Model* mtree = models[0];
		self->lp += mtree->logP(mtree);
	}
	return self->lp;
}

double _singleTreeLikelihood_full_logP(Model *self){
	SingleTreeLikelihood_update_all_nodes((SingleTreeLikelihood*)self->obj);
	return _singleTreeLikelihood_logP(self);
}

// Set flags for gradient calculation
void _model_prepare_gradient(Model* self, const Parameters* ps){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	
	size_t paramCount = Parameters_count(ps);
	bool prepare_tree = false;
	bool prepare_subsitution_model = false;
	bool prepare_site_model = false;
	bool prepare_branch_model = false;
	tlk->prepared_gradient = 0;
	size_t gradient_length = 0;
	tlk->include_root_freqs = true;
	for (size_t i = 0; i < paramCount; i++) {
		Parameter* p = Parameters_at(ps, i);
		if ((p->model == MODEL_TREE || p->model == MODEL_TREE_TRANSFORM) && !prepare_tree) {
			prepare_tree = true;
			tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_TREE_MODEL;
			if(!Tree_is_time_mode(tlk->tree)){
				gradient_length += Tree_node_count(tlk->tree);
			}
			else{
				gradient_length += Tree_tip_count(tlk->tree) - 1;
			}
		}
		else if (p->model == MODEL_SITEMODEL && !prepare_site_model) {
			prepare_site_model = true;
			tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_SITE_MODEL;
			gradient_length += Parameters_count(tlk->sm->rates);
			if(tlk->sm->proportions != NULL){
				gradient_length += Parameters_count(tlk->sm->proportions->parameters);
			}
			if(tlk->sm->mu != NULL){
				gradient_length++;
			}
		}
		else if (p->model == MODEL_BRANCHMODEL && !prepare_branch_model) {
			prepare_branch_model = true;
			tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_BRANCH_MODEL;
			gradient_length += Parameters_count(tlk->bm->rates);
		}
		else if (p->model == MODEL_SUBSTITUTION && !prepare_subsitution_model) {
			prepare_subsitution_model = true;
			tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_UNCONSTRAINED;
			gradient_length += tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K - 1;
			gradient_length += tlk->m->simplex->K - 1;
			tlk->include_root_freqs = false;
		}
	}
	if(tlk->gradient == NULL){
		tlk->gradient = calloc(gradient_length, sizeof(double));
		tlk->gradient_length = gradient_length;
	}
	else if (tlk->gradient_length < gradient_length) {
		tlk->gradient = realloc(tlk->gradient, sizeof(double)* gradient_length);
		tlk->gradient_length = gradient_length;
	}
}

size_t TreeLikelihood_initialize_gradient(Model *self, int flags){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	tlk->prepared_gradient = flags;
	size_t gradient_length = 0;
	tlk->include_root_freqs = true;
	
	int prepare_tree = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_TREE_MODEL;
	int prepare_site_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SITE_MODEL;
	int prepare_branch_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_BRANCH_MODEL;
	int prepare_substitution_model_unconstrained = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_UNCONSTRAINED;
	int prepare_substitution_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL;
	int prepare_substitution_model_rates = prepare_substitution_model | (tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_RATES);
	int prepare_substitution_model_frequencies = prepare_substitution_model | (tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_FREQUENCIES);
	
	if(prepare_substitution_model_unconstrained && prepare_substitution_model){
		fprintf(stderr, "Can only request unconstrained and constrained gradient at the same time\n");
		exit(2);
	}
	if(flags == 0){
		prepare_tree = true;
		prepare_site_model = tlk->sm->proportions != NULL || Parameters_count(tlk->sm->rates) > 0 || tlk->sm->mu != NULL;
		prepare_branch_model = tlk->bm!= NULL;
		prepare_substitution_model = tlk->m->dPdp != NULL;
		tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_TREE_MODEL;
		if(prepare_site_model){
			tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_SITE_MODEL;
		}
		if(prepare_branch_model){
			tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_BRANCH_MODEL;
		}
		if(prepare_substitution_model){
			tlk->prepared_gradient |= TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL;
		}
	}
	if(prepare_tree){
		if(!Tree_is_time_mode(tlk->tree)){
			gradient_length += Tree_node_count(tlk->tree);
		}
		else{
			gradient_length += Tree_tip_count(tlk->tree) - 1;
		}
	}
	if (prepare_site_model) {
		gradient_length += Parameters_count(tlk->sm->rates);
		if(tlk->sm->proportions != NULL){
			gradient_length += Parameters_count(tlk->sm->proportions->parameters);
		}
		if(tlk->sm->mu != NULL){
			gradient_length++;
		}
	}
	if (prepare_branch_model) {
		gradient_length += Parameters_count(tlk->bm->rates);
	}
	if (prepare_substitution_model_unconstrained) {
		gradient_length += tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K - 1;
		gradient_length += tlk->m->simplex->K - 1;
		tlk->include_root_freqs = false;
	}
	else{
		if (prepare_substitution_model_rates) {
			gradient_length += tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K;
			tlk->m->grad_wrt_reparam = false;
			tlk->include_root_freqs = false;
		}
		if (prepare_substitution_model_frequencies) {
			gradient_length += tlk->m->simplex->K;
			tlk->m->grad_wrt_reparam = false;
			tlk->include_root_freqs = false;
		}
	} 

	if(tlk->gradient == NULL){
		tlk->gradient = calloc(gradient_length, sizeof(double));
		tlk->gradient_length = gradient_length;
	}
	else if (tlk->gradient_length < gradient_length) {
		tlk->gradient = realloc(tlk->gradient, sizeof(double)* gradient_length);
		tlk->gradient_length = gradient_length;
	}
	return gradient_length;
}

double* TreeLikelihood_gradient(Model *self){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	
	if(tlk->update_upper){
		if(Tree_is_time_mode(tlk->tree)){
			Tree_update_heights(tlk->tree);
		}
		double logP = tlk->calculate(tlk); // make sure it is updated
		if (isnan(logP) || isinf(logP)) {
			for(size_t i = 0; i < tlk->gradient_length; i++){
				tlk->gradient[i] = NAN;
			}
		}
		else{
			update_upper_partials(tlk, Tree_root(tlk->tree), tlk->include_root_freqs);
			TreeLikelihood_calculate_gradient(self, tlk->gradient);
			tlk->update_upper = false;
		}
	}
	return tlk->gradient;
}

double _singleTreeLikelihood_dlogP_prepared(Model *self, const Parameter* p){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	
	if(tlk->gradient_length == 0) return 0.0;
	if(p->model != MODEL_TREE && p->model != MODEL_TREE_TRANSFORM && p->model != MODEL_SITEMODEL
	    && p->model != MODEL_SUBSTITUTION && p->model != MODEL_BRANCHMODEL) return 0.0;
			
	if(tlk->update_upper){
		if(Tree_is_time_mode(tlk->tree)){
			Tree_update_heights(tlk->tree);
		}
		double logP = tlk->calculate(tlk); // make sure it is updated
		if (isnan(logP) || isinf(logP)) {
				return logP;
		}
		update_upper_partials(tlk, Tree_root(tlk->tree), tlk->include_root_freqs);
		TreeLikelihood_calculate_gradient(self, tlk->gradient);
		tlk->update_upper = false;
	}
	
	int prepare_tree = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_TREE_MODEL;
	int prepare_site_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SITE_MODEL;
	int prepare_branch_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_BRANCH_MODEL;
	int prepare_substitution_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_UNCONSTRAINED;
	size_t offset = 0;
	
	if(prepare_tree){
		if ((p->model == MODEL_TREE || p->model == MODEL_TREE_TRANSFORM)) {
			return tlk->gradient[p->id];
		}
		if(!Tree_is_time_mode(tlk->tree)){
			offset += Tree_node_count(tlk->tree);
		}
		else{
			offset += Tree_tip_count(tlk->tree) - 1;
		}
	}
	if (prepare_site_model) {
		if (p->model == MODEL_SITEMODEL) {
			if (Parameters_count(tlk->sm->rates) != 0) {
				if(p == Parameters_at(tlk->sm->rates, 0)){
					return tlk->gradient[offset + p->id];
				}
				offset++;
			}
			if (tlk->sm->proportions != NULL) {
				if(p == Parameters_at(tlk->sm->proportions->parameters, 0)){
					return tlk->gradient[offset + p->id];
				}
				offset++;
			}
			if (tlk->sm->mu != NULL) {
				if(p == tlk->sm->mu){
					return tlk->gradient[offset + p->id];
				}
				offset++;
			}
		}
		else{
			offset += Parameters_count(tlk->sm->rates);
			if (tlk->sm->proportions != NULL) {
				offset++;
			}
			if (tlk->sm->mu != NULL) {
				offset++;
			}
		}
	}
	if (prepare_branch_model) {
		size_t rateCount = Parameters_count(tlk->bm->rates);
		if (p->model == MODEL_BRANCHMODEL) {
			if(rateCount == 1){
				return tlk->gradient[offset];
			}
			else{
				// gradient is indexed using class_id
				// node and rate parameter have the same id
				Node* node = Tree_node(tlk->tree, p->id);
				size_t index = Node_isleaf(node) ? node->class_id : node->class_id + Tree_tip_count(tlk->tree);
				return tlk->gradient[offset + p->id];
			}
		}
		offset += rateCount;
	}
	if (prepare_substitution_model) {
		if (p->model == MODEL_SUBSTITUTION) {
			if (Parameters_count(tlk->m->rates) != 0) {
				if(p->id < Parameters_count(tlk->m->rates) && p == Parameters_at(tlk->m->rates, p->id)){
					return tlk->gradient[offset + p->id];
				}
				offset += Parameters_count(tlk->m->rates);
			}
			else if (tlk->m->rates_simplex != NULL) {
				if(p->id < (tlk->m->rates_simplex->K - 1) && p == Parameters_at(tlk->m->rates_simplex->parameters, p->id)){
					return tlk->gradient[offset + p->id];
				}
				offset += tlk->m->rates_simplex->K - 1;
			}
			if(p == Parameters_at(tlk->m->simplex->parameters, p->id)){
				return tlk->gradient[offset + p->id];
			}
		}
	}
	
	return 0;
}

static void _calculate_dlog_jacobian(const Node* noderef, Node* node, double* dlogP, double* descendant, unsigned *map, Parameters* reparams, double* lowers){
    if(!Node_isleaf(node)){
        if(!Node_isroot(node) && node != noderef){
            descendant[Node_id(node)] = descendant[Node_id(node->parent)]*Parameters_value(reparams, map[Node_id(node)]);
        }
        else if(!Node_isroot(node)){
            descendant[Node_id(node)] = Node_height(Node_parent(node)) - lowers[Node_id(node)];
        }
        else{
            descendant[Node_id(node)] = 1;
        }
        _calculate_dlog_jacobian(noderef, node->left, dlogP, descendant, map, reparams, lowers);
        _calculate_dlog_jacobian(noderef, node->right, dlogP, descendant, map, reparams, lowers);
        
        if(!Node_isroot(node) && node != noderef){
            *dlogP += descendant[Node_id(node->parent)]/(Node_height(Node_parent(node)) - lowers[Node_id(node)]);
        }
    }
}

double _singleTreeLikelihood_d2logP(Model *self, const Parameter* p){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	
	
	Node* node = NULL;
	int i = 0;
	for (; i < Tree_node_count(tlk->tree); i++) {
		node = Tree_node(tlk->tree, i);
		if (strcmp(node->distance->name, Parameter_name(p)) == 0) {
			break;
		}
	}
	
	if(Tree_node_count(tlk->tree) == i){
		return Model_second_derivative(self, p, NULL, 0.001);
	}
	
	double* pattern_likelihoods = tlk->pattern_lk + tlk->sp->count;
	
	if(tlk->update_upper){
		//		for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
		//			if (tlk->update_nodes[i]) {
		//				double lk = tlk->calculate_upper(tlk, Tree_node(tlk->tree, i));
		//				tlk->node_upper = Tree_node(tlk->tree, i);
		//				tlk->use_upper = true;
		//				//					printf("+ %d %f\n", i, lk);
		//				tlk->update_nodes[i] = false;
		//				return lk;
		//			}
		//		}
//		SingleTreeLikelihood_update_all_nodes(tlk);
		
		double logP = tlk->calculate(tlk); // make sure it is updated
		if (isnan(logP) || isinf(logP)) {
			return logP;
		}
		update_upper_partials(tlk, Tree_root(tlk->tree), false);
		
		for (int i = 0; i < tlk->sp->count; i++) {
			pattern_likelihoods[i] = exp(tlk->pattern_lk[i]);
		}
		tlk->update_upper = false;
	}

	double* pattern_dlikelihoods = tlk->pattern_lk + 2*tlk->sp->count;
	calculate_dldt_uppper(tlk, node, pattern_dlikelihoods);
//	print_dvector(pattern_dlikelihoods, tlk->sp->count);
//	print_dvector(tlk->pattern_lk, tlk->sp->count);
//	exit(1);
	
	double d2logP = d2lnldt2_uppper(tlk, node, pattern_likelihoods, pattern_dlikelihoods);
	if(isnan(d2logP)){
		SingleTreeLikelihood_update_all_nodes(tlk);
		//		Tree_print_parameters(tlk->tree);
		//		printf("%f\n", tlk->calculate(tlk));
		//		printf("%f\n", Model_first_derivative(self, p, 0.001));
		//		printf("%s %f\n", p->name, p->value);
		//		exit(10);
	}
	//printf("%f\n", d2logP);exit(10);
	return d2logP;
}


//TODO: update_upper is always true since lower likelihoods are used to calculate mixed derivative of the likelihood
double _singleTreeLikelihood_ddlogP(Model *self, const Parameter* p1, const Parameter* p2){
    SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
    
    Node* node1 = NULL;
    Node* node2 = NULL;

    for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
        Node* node = Tree_node(tlk->tree, i);
        if (strcmp(node->distance->name, Parameter_name(p1)) == 0) {
            node1 = node;
        }
        else if (strcmp(node->distance->name, Parameter_name(p2)) == 0) {
            node2 = node;
        }
    }
    
    if(node1 == NULL || node2 == NULL){
        return Model_mixed_derivative(self, p1, p2);
    }
	
	double* pattern_likelihoods = tlk->pattern_lk + tlk->sp->count;
	
	if(tlk->update_upper){
		//		for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
		//			if (tlk->update_nodes[i]) {
		//				double lk = tlk->calculate_upper(tlk, Tree_node(tlk->tree, i));
		//				tlk->node_upper = Tree_node(tlk->tree, i);
		//				tlk->use_upper = true;
		//				//					printf("+ %d %f\n", i, lk);
		//				tlk->update_nodes[i] = false;
		//				return lk;
		//			}
		//		}
		//		SingleTreeLikelihood_update_all_nodes(tlk);
		
		double logP = tlk->calculate(tlk); // make sure it is updated
		if (isnan(logP) || isinf(logP)) {
			return logP;
		}
		update_upper_partials(tlk, Tree_root(tlk->tree), false);
		
		for (int i = 0; i < tlk->sp->count; i++) {
			pattern_likelihoods[i] = exp(tlk->pattern_lk[i]);
		}
		tlk->update_upper = false;
	}
	
	double* pattern_dlikelihoods1 = tlk->pattern_lk + tlk->sp->count*2;
	double* pattern_dlikelihoods2 = tlk->pattern_lk + tlk->sp->count*3;
	calculate_dldt_uppper(tlk, node1, pattern_dlikelihoods1);
	calculate_dldt_uppper(tlk, node2, pattern_dlikelihoods2);
	
    double *pmats1 = dvector(tlk->sm->cat_count * tlk->matrix_size);
    double *pmats2 = dvector(tlk->sm->cat_count * tlk->matrix_size);
    
    memcpy(pmats1, tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)], tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
    memcpy(pmats2, tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)], tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
    
    double bl1 = Node_distance(node1);
    double bl2 = Node_distance(node2);
	bool node1_transpose = tlk->partials[0][Node_id(node1)] == NULL;
	bool node2_transpose = tlk->partials[0][Node_id(node2)] == NULL;
	
    for (int i = 0; i < tlk->sm->cat_count; i++) {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
        if( tlk->use_SIMD ){
            if( node1_transpose ){
                tlk->m->dp_dt_transpose(tlk->m,
                                            bl1 * tlk->sm->get_rate(tlk->sm, i),
                                            &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
            }
            else {
                tlk->m->dp_dt(tlk->m,
                                  bl1 * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
            }
            if( node2_transpose ){
                tlk->m->dp_dt_transpose(tlk->m,
                                            bl2 * tlk->sm->get_rate(tlk->sm, i),
                                            &tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size]);
            }
            else {
                tlk->m->dp_dt(tlk->m,
                                  bl2 * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size]);
            }
        }
        else{
            tlk->m->dp_dt(tlk->m,
                              bl1 * tlk->sm->get_rate(tlk->sm, i),
                              &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
            
            tlk->m->dp_dt(tlk->m,
                              bl2 * tlk->sm->get_rate(tlk->sm, i),
                              &tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size]);
        }
#else
        tlk->m->dp_dt(tlk->m,
                          bl1 * tlk->sm->get_rate(tlk->sm, i),
                          &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
        
        tlk->m->dp_dt(tlk->m,
                          bl2 * tlk->sm->get_rate(tlk->sm, i),
                          &tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size]);
#endif
        
        for (int k = 0; k < tlk->matrix_size; k++) {
            tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
            tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
        }
        
    }
    
    SingleTreeLikelihood_update_one_node(tlk, node1);
    SingleTreeLikelihood_update_one_node(tlk, node2);
    
    Node* root = Tree_root(tlk->tree);
    
    _calculate_partials( tlk, root );
    
    tlk->integrate_partials(tlk, tlk->partials[tlk->current_matrices_indexes[Node_id(root)]][Node_id(root)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
    int v = 0;
    int i = 0;
    
    const int nstate   = tlk->m->nstate;
    const int patternCount = tlk->sp->count;
    const double* freqs = tlk->get_root_frequencies(tlk);
    
    for ( int k = 0; k < patternCount; k++ ) {
        
        tlk->pattern_lk[k] = 0;
        for ( i = 0; i < nstate; i++ ) {
            
            tlk->pattern_lk[k] += freqs[i] * tlk->partials[tlk->current_matrices_indexes[Node_id(root)]][Node_id(root)][v];
            v++;
        }
        //tlk->pattern_lk[k] = log(tlk->pattern_lk[k]);
        
        if ( tlk->scale ) {
            //printf("scaling\n");
            tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
        }
    }
	
    double ddlogP = 0;
    for ( int i = 0; i < tlk->sp->count; i++) {
        ddlogP += ((tlk->pattern_lk[i] * pattern_likelihoods[i] - (pattern_dlikelihoods1[i]*pattern_dlikelihoods2[i])) / (pattern_likelihoods[i]*pattern_likelihoods[i])) * tlk->sp->weights[i];
    }
    
    memcpy( tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)], pmats1, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
    SingleTreeLikelihood_update_one_node(tlk, node1);
    
    memcpy( tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)], pmats2, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
    SingleTreeLikelihood_update_one_node(tlk, node2);
    
    free(pmats1);
    free(pmats2);
    return ddlogP;
}

static void _treeLikelihood_model_free( Model *self ){
#ifdef DEBUG_REF_COUNTING
	printf("Free treelikelihood model: %d\n", self->ref_count);
#endif
	if(self->ref_count == 1){
		SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
		int count = (tlk->bm==NULL?3:4);
		Model** list = (Model**)self->data;
		for(int i = 0; i < count; i++){
			list[i]->free(list[i]);
		}
		free_SitePattern(tlk->sp);
		free_SingleTreeLikelihood_internals(tlk);
		free(self->data);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _treeLikelihood_model_clone(Model* self, Hashtable* hash){
	Model** list = (Model**)self->data;
	Model* mtree = list[0];
	Model* mm = list[1];
	Model* msm = list[2];
	Model* mbm = list[3];
	Model *mtreeclone = NULL;
	Model *mmclone = NULL;
	Model *msmclone = NULL;
	Model* mbmclone = NULL;
	
	if(Hashtable_exists(hash, mtree->name)){
		mtreeclone = Hashtable_get(hash, mtree->name);
		mtreeclone->ref_count++;
	}
	else{
		mtreeclone = mtree->clone(mtree, hash);
		Hashtable_add(hash, mtreeclone->name, mtreeclone);
	}
	
	if(Hashtable_exists(hash, mm->name)){
		mmclone = Hashtable_get(hash, mm->name);
		mmclone->ref_count++;
	}
	else{
		mmclone = mm->clone(mm, hash);
		Hashtable_add(hash, mmclone->name, mmclone);
	}
	
	if(Hashtable_exists(hash, msm->name)){
		msmclone = Hashtable_get(hash, msm->name);
		msmclone->ref_count++;
	}
	else{
		msmclone = msm->clone(msm, hash);
		Hashtable_add(hash, msmclone->name, msmclone);
	}
	
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	SingleTreeLikelihood* clonetlk;
	
	if(mbm != NULL){
		if(Hashtable_exists(hash, mbm->name)){
			mbmclone = Hashtable_get(hash, mbm->name);
			mbmclone->ref_count++;
		}
		else{
			mbmclone = mbm->clone(mbm, hash);
			Hashtable_add(hash, mbmclone->name, mbmclone);
		}
		clonetlk = clone_SingleTreeLikelihood_with(tlk, (Tree*)mtreeclone->obj, (SubstitutionModel*)mmclone->obj, (SiteModel*)msmclone->obj, tlk->sp, (BranchModel*)mbmclone->obj);
	}
	else{
		clonetlk = clone_SingleTreeLikelihood_with(tlk, (Tree*)mtreeclone->obj, (SubstitutionModel*)mmclone->obj, (SiteModel*)msmclone->obj, tlk->sp, NULL);
	}
	Model* clone = new_TreeLikelihoodModel(self->name, clonetlk, mtreeclone, mmclone, msmclone, mbmclone);
	Hashtable_add(hash, clone->name, clone);
	mtreeclone->free(mtreeclone);
	mmclone->free(mmclone);
	msmclone->free(msmclone);
	if(mbmclone) mbmclone->free(mbmclone);
	clone->store = self->store;
	clone->restore = self->restore;
	clone->storedLogP = self->storedLogP;
	clone->lp = self->lp;
	clone->full_logP = self->full_logP;
	return clone;
}

// TreeLikelihood listen to the TreeModel, SiteModel, BranchModel
Model * new_TreeLikelihoodModel( const char* name, SingleTreeLikelihood *tlk,  Model *tree, Model *m, Model *sm, Model *bm ){
	Model *model = new_Model(MODEL_TREELIKELIHOOD,name, tlk);

	tree->listeners->add( tree->listeners, model );
	m->listeners->add( m->listeners, model );
	if(bm != NULL)bm->listeners->add( bm->listeners, model );
	sm->listeners->add( sm->listeners, model );
	model->handle_restore = _treelikelihood_handle_restore;
	
	model->logP = _singleTreeLikelihood_logP;
	model->full_logP = _singleTreeLikelihood_full_logP;
	model->dlogP = _singleTreeLikelihood_dlogP_prepared;
	model->d2logP = _singleTreeLikelihood_d2logP;
	model->ddlogP = _singleTreeLikelihood_ddlogP;
	model->update = _treelikelihood_handle_change;
	model->free = _treeLikelihood_model_free;
	model->clone = _treeLikelihood_model_clone;
	model->data = (Model**)malloc(sizeof(Model*)*4);
	model->store = _singleTreeLikelihood_store;
	model->restore = _singleTreeLikelihood_restore;
	Model** list = (Model**)model->data;
	list[0] = tree;
	list[1] = m;
	list[2] = sm;
	list[3] = bm;
	tree->ref_count++;
	m->ref_count++;
	sm->ref_count++;
	if(bm != NULL) bm->ref_count++;
	
	model->prepare_gradient = _model_prepare_gradient;
	return model;
}

Model * new_TreeLikelihoodModel_from_json(json_node*node, Hashtable*hash){
	char* allowed[] = {
		"branchmodel",
		"include_jacobian",
		"reparameterized", // deprecated. use include_jacobian
		"root_frequencies",
		"sitemodel",
		"sitepattern",
		"sse",
		"substitutionmodel",
		"tipstates",
		"tree"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* patterns_node = get_json_node(node, "sitepattern");
	json_node* tree_node = get_json_node(node, "tree");
	json_node* m_node = get_json_node(node, "substitutionmodel");
	json_node* sm_node = get_json_node(node, "sitemodel");
	json_node* bm_node = get_json_node(node, "branchmodel");
	bool include_jacobian = get_json_node_value_bool(node, "reparameterized", false);
	include_jacobian |= get_json_node_value_bool(node, "include_jacobian", false);
	bool use_tip_states = get_json_node_value_bool(node, "tipstates", true);
	
	Model* mtree = NULL;
	Model* mm = NULL;
	Model* msm = NULL;
	SitePattern* patterns = NULL;
	Model* mbm = NULL;
	BranchModel* bm = NULL;
	
	if(patterns_node->node_type == MJSON_STRING){
		char* ref = (char*)patterns_node->value;
		// check it starts with a &
		patterns = Hashtable_get(hash, ref+1);
		patterns->ref_count++;
	}
	else{
		char* id = get_json_node_value_string(patterns_node, "id");
		patterns = new_SitePattern_from_json(patterns_node, hash);
		Hashtable_add(hash, id, patterns);
	}
	
	if (tree_node->node_type == MJSON_STRING) {
		char* ref = (char*)tree_node->value;
		// check it starts with a &
		mtree = Hashtable_get(hash, ref+1);
		mtree->ref_count++;
	}
	else{
		char* id = get_json_node_value_string(tree_node, "id");
		mtree = new_TreeModel_from_json(tree_node, hash);
		Hashtable_add(hash, id, mtree);
	}
	
	if (sm_node->node_type == MJSON_STRING) {
		char* ref = (char*)sm_node->value;
		// check it starts with a &
		msm = Hashtable_get(hash, ref+1);
		msm->ref_count++;
	}
	else{
		char* id = get_json_node_value_string(sm_node, "id");
		msm = new_SiteModel_from_json(sm_node, hash);
		Hashtable_add(hash, id, msm);
	}
	
	// the old way
	if (m_node == NULL) {
		m_node = get_json_node(sm_node, "substitutionmodel");
	}
	if (m_node->node_type == MJSON_STRING && ((char*)m_node->value)[0] == '&') {
		char* ref = (char*)m_node->value;
		mm = Hashtable_get(hash, ref+1);
		mm->ref_count++;
	}
	else{
		mm = new_SubstitutionModel_from_json(m_node, hash);
	}

	
	if(bm_node != NULL){
		if (bm_node->node_type == MJSON_STRING) {
			char* ref = (char*)bm_node->value;
			// check it starts with a &
			mbm = Hashtable_get(hash, ref+1);
			mbm->ref_count++;
		}
		else{
			char* id = get_json_node_value_string(bm_node, "id");
			mbm = new_BranchModel_from_json(bm_node, hash);
			Hashtable_add(hash, id, mbm);
		}
		bm = mbm->obj;
	}
	
	if(patterns->partials != NULL){
		use_tip_states = false;
	}
	
	SingleTreeLikelihood* tlk = new_SingleTreeLikelihood((Tree*)mtree->obj, (SubstitutionModel*)mm->obj, (SiteModel*)msm->obj, patterns, bm, use_tip_states);
	char* id = get_json_node_value_string(node, "id");
	Model* model = new_TreeLikelihoodModel(id, tlk, mtree, mm, msm, mbm);
	mm->free(mm);
	mtree->free(mtree);
	msm->free(msm);
	if(mbm != NULL) mbm->free(mbm);
	
	bool root_freqs_unknown = get_json_node_value_bool(node, "root_frequencies", false);
	if(root_freqs_unknown){
		tlk->get_root_frequencies = get_root_frequencies_fixed;
		tlk->root_frequencies = dvector(tlk->m->nstate);
		for (int i = 0; i < tlk->m->nstate; i++) {
			tlk->root_frequencies[i] = 1;
		}
	}
	
	bool use_sse = get_json_node_value_bool(node, "sse", true);
	if (!use_sse) {
		SingleTreeLikelihood_enable_SSE(tlk, false);
	}
	tlk->include_jacobian = include_jacobian;
	return model;
}

#pragma mark -
// MARK: SingleTreeLikelihood

void allocate_storage(SingleTreeLikelihood* tlk, size_t index){
    Tree* tree = tlk->tree;
    size_t nodeCount = Tree_node_count(tree);
    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
    
    if(index == 0){
        tlk->current_matrices_indexes = uivector(2*nodeCount);
        tlk->current_partials_indexes = uivector(2*nodeCount);
        tlk->stored_matrices_indexes = NULL;
        tlk->stored_partials_indexes = NULL;
        tlk->partials = (double***)malloc(2*sizeof(double**));
        tlk->partials[0] = (double**)malloc( tlk->partials_dim*sizeof(double*));
        tlk->partials[1] = NULL;
        
        int i = 0;
        for ( ; i < nodeCount; i++ ) {
            if( Node_isleaf( nodes[i] ) && tlk->use_tip_states ){
                tlk->partials[0][Node_id( nodes[i] )] = NULL;
            }
            else {
                tlk->partials[0][Node_id( nodes[i] )] = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
            }
        }
        for ( ; i < tlk->partials_dim; i++ ) {
            tlk->partials[0][i] = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
        }
        tlk->matrices = (double***)malloc( 2*sizeof(double**) );
        tlk->matrices[0] = (double**)malloc( tlk->matrix_dim*sizeof(double*) );
        tlk->matrices[1] = NULL;
        
        int mat_len = tlk->matrix_size*tlk->sm->cat_count;
        for ( int i = 0; i < tlk->matrix_dim; i++ ) {
            tlk->matrices[0][i] = aligned16_malloc( mat_len * sizeof(double) );
        }
    }
    else{
        tlk->stored_matrices_indexes = uivector(nodeCount*2);
        tlk->stored_partials_indexes = uivector(Tree_node_count(tree)*2);
        tlk->partials[1] = (double**)malloc( tlk->partials_dim*sizeof(double*));
        int i = 0;
        for ( ; i < Tree_node_count(tree); i++ ) {
            if( Node_isleaf( nodes[i] ) && tlk->use_tip_states ){
                tlk->partials[1][Node_id( nodes[i] )] = NULL;
            }
            else {
                tlk->partials[1][Node_id( nodes[i] )] = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
            }
        }
        for ( ; i < tlk->partials_dim; i++ ) {
            tlk->partials[1][i] = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
        }
        tlk->matrices[1] = (double**)malloc( tlk->matrix_dim*sizeof(double*) );
        
        int mat_len = tlk->matrix_size*tlk->sm->cat_count;
        for ( int i = 0; i < tlk->matrix_dim; i++ ) {
            tlk->matrices[1][i] = aligned16_malloc( mat_len * sizeof(double) );
        }
    }
}

SingleTreeLikelihood * new_SingleTreeLikelihood( Tree *tree, SubstitutionModel *m, SiteModel *sm, SitePattern *sp, BranchModel *bm, bool use_tip_states ){
	SingleTreeLikelihood *tlk = (SingleTreeLikelihood *)malloc( sizeof(SingleTreeLikelihood));
	assert(tlk);
	
	tlk->tree = tree;
	tlk->m = m;
	tlk->sm = sm;
	tlk->sp = sp;
	tlk->bm = bm;
	
	tlk->pattern_count = tlk->sp->count;
	if(tlk->sm->site_category == NULL){
		tlk->cat_count = sm->cat_count;
	}
	else{
		tlk->cat_count = 1;
	}
	tlk->use_tip_states = use_tip_states;
	
	Node **nodes = Tree_get_nodes( tree, POSTORDER );
	
	tlk->partials_size = sp->count * m->nstate * tlk->cat_count;
	
	
	tlk->matrix_dim = Tree_node_count(tree);
	tlk->matrix_size = m->nstate * m->nstate;
	
	tlk->upper_partial_indexes = ivector(Tree_node_count(tree));
	
	tlk->root_partials_size = sp->count*m->nstate*3;
	tlk->pattern_lk_size = sp->count*4;
	
	tlk->partials_dim = Tree_node_count(tree)*2; // allocate *2 for upper likelihoods
	
    // odd number of state
//    if( sm->nstate & 1 ){
//        tlk->partials_size += sp->count * sm->cat_count;
//        tlk->matrix_size   += sm->nstate;
//    }
	
    allocate_storage(tlk, 0);
	// This now allocated when store is called for the first time
    // allocate_storage(tlk, 1); // used by mcmc with store/restore
	
	tlk->root_partials = aligned16_malloc( tlk->root_partials_size * sizeof(double) );
    assert(tlk->root_partials);
	
	tlk->update_nodes = bvector(Tree_node_count(tree));
	for (int i = 0; i < Tree_node_count(tree); i++){
		tlk->update_nodes[i] = true;
	}
	tlk->update = true;
	tlk->update_upper = true;
	
	tlk->pattern_lk = dvector(tlk->pattern_lk_size);
	tlk->lk = 0.;
	
	tlk->calculate = _calculate;
	tlk->calculate_upper = _calculate_uppper;
	
	tlk->update_partials      = update_partials_general;
	tlk->integrate_partials   = integrate_partials_general;
	tlk->node_log_likelihoods = node_log_likelihoods_general;
	tlk->calculate_per_cat_partials = calculate_branch_partials;
	tlk->update_partials_flexible = NULL;
	
	if ( m->nstate == 4 ) {
		tlk->update_partials      = update_partials_4;
		tlk->update_partials_flexible = update_partials_flexible_4;
		tlk->integrate_partials   = integrate_partials_4;
		tlk->node_log_likelihoods = node_log_likelihoods_4;
		tlk->calculate_per_cat_partials = calculate_branch_partials_4;
	}
	else if( m->nstate == 20 ){
		tlk->update_partials      = update_partials_general;
		tlk->integrate_partials   = integrate_partials_general;
		tlk->node_log_likelihoods = node_log_likelihoods_general;
		tlk->calculate_per_cat_partials = calculate_branch_partials_20;
	}
	else if( m->nstate >= 60 ){
		tlk->update_partials      = update_partials_codon;
		tlk->integrate_partials   = integrate_partials_codon;
		tlk->node_log_likelihoods = node_log_likelihoods_codon;
	}
    
	tlk->mapping = ivector(Tree_node_count(tree));
	
	// map node names to sequence names
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
		if( Node_isleaf( nodes[i] ) ){
			tlk->mapping[Node_id(nodes[i])] = get_sequence_index(tlk->sp, nodes[i]->name);
			if(tlk->mapping[Node_id(nodes[i])] == -1){
				printf("Could not find taxon `%s` in alignment\n", nodes[i]->name);
				exit(1);
			}
		}
		else tlk->mapping[Node_id( nodes[i] )] = -1;
	}
	
	if (!use_tip_states) {
		for (size_t i = 0; i < Tree_node_count(tlk->tree); i++) {
			if(Node_isleaf(nodes[i])){
				sp->get_partials(sp, tlk->mapping[Node_id(nodes[i])], tlk->partials[0][Node_id(nodes[i])]);
				
				for(size_t k = 1; k < tlk->sm->cat_count; k++){
					size_t offset = k*m->nstate*tlk->pattern_count;
					memcpy(tlk->partials[0][Node_id(nodes[i])] + offset, tlk->partials[0][Node_id(nodes[i])], sizeof(double)*m->nstate*tlk->pattern_count);
				}
			}
		}
	}
	
	tlk->scale = false;
	tlk->scaling_factors = NULL;
	tlk->scaling_threshold = 1.E-40;
	
	tlk->node_id = -1;
    
    // Upper likelihood calculation
    // Variables are instanciated when we need them using SingleTreeLikelihood_use_upper
    tlk->use_upper = false;
    tlk->node_upper = NULL;
    
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	tlk->use_SIMD = false;
#endif
	
#ifdef SSE3_ENABLED
	if ( m->nstate == 4 ) {
		tlk->update_partials      = update_partials_4_SSE;
		tlk->update_partials_flexible = update_partials_flexible_4_SSE;
		tlk->integrate_partials   = integrate_partials_4_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_4_SSE;
		tlk->calculate_per_cat_partials = calculate_branch_partials_4_SSE;
        tlk->use_SIMD = true;
	}
	else if ( m->nstate == 20 ) {
		tlk->update_partials      = update_partials_20_SSE;
		tlk->integrate_partials   = integrate_partials_general;
		tlk->node_log_likelihoods = node_log_likelihoods_general;
		tlk->calculate_per_cat_partials = calculate_branch_partials_20_SSE;
        tlk->use_SIMD = true;
	}
    // Does not work for matrix with odd dimension.
    // If we want to load doubles at indexes 1 and 2 from matrices it will crash
    // same issue with partials, it will try to store at odd indexes
	else if( m->nstate >= 60 && !(m->nstate & 1) ){
		tlk->update_partials      = update_partials_codon_SSE;
		tlk->integrate_partials   = integrate_partials_codon_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_codon_SSE;
        tlk->use_SIMD = true;
	}
    else if(!(m->nstate & 1)){
        tlk->update_partials      = update_partials_general_even_SSE;
        tlk->integrate_partials   = integrate_partials_general_even_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_general_even_SSE;
		tlk->calculate_per_cat_partials = calculate_branch_partials_even_SSE;
        tlk->use_SIMD = true;
    }
#elif AVX_ENABLED
	if ( sm->nstate == 4 ) {
		tlk->update_partials      = update_partials_4_AVX;
		tlk->integrate_partials   = integrate_partials_4_AVX;
		tlk->node_log_likelihoods = node_log_likelihoods_4_AVX;
        tlk->use_SIMD = true;
	}
	
#endif

    tlk->nthreads = 1;
	tlk->root_frequencies = NULL;
	tlk->get_root_frequencies = get_root_frequencies;

	tlk->tripod = false;
	tlk->gradient = NULL;
	tlk->gradient_length = 0;
	tlk->prepared_gradient = 0;
	return tlk;
}

void free_SingleTreeLikelihood_internals( SingleTreeLikelihood *tlk ){
	if(tlk->scaling_factors != NULL ){
		free_dmatrix( tlk->scaling_factors[0], tlk->partials_dim);
        if(tlk->scaling_factors[1] != NULL){
            free_dmatrix( tlk->scaling_factors[1], tlk->partials_dim);
        }
		free(tlk->scaling_factors);
	}
	for ( int i = 0; i < tlk->partials_dim; i++ ) {
		if( tlk->partials[0][i] != NULL ){
			free(tlk->partials[0][i]);
		}
	}
    free(tlk->partials[0]);
	free_dmatrix(tlk->matrices[0], tlk->matrix_dim);
	free(tlk->current_partials_indexes);
	free(tlk->current_matrices_indexes);
    
    if(tlk->partials[1] != NULL){
        for ( int i = 0; i < tlk->partials_dim; i++ ) {
            if( tlk->partials[1][i] != NULL ){
                free(tlk->partials[1][i]);
            }
        }
        free(tlk->partials[1]);
        free_dmatrix(tlk->matrices[1], tlk->matrix_dim);
        free(tlk->stored_partials_indexes);
        free(tlk->stored_matrices_indexes);
    }
    
    free(tlk->partials);
    free(tlk->matrices);
	
	free(tlk->mapping);
	free(tlk->update_nodes);
	free(tlk->pattern_lk);
	free(tlk->root_partials);
	
	free(tlk->upper_partial_indexes);
	if(tlk->root_frequencies != NULL)free(tlk->root_frequencies);
	if(tlk->gradient != NULL) free(tlk->gradient);
	free(tlk);
}

void free_SingleTreeLikelihood( SingleTreeLikelihood *tlk ){
	if ( tlk->m != NULL ) tlk->m->free( tlk->m );
	if ( tlk->bm != NULL ) tlk->bm->free( tlk->bm, false );
	if ( tlk->sm != NULL ) free_SiteModel( tlk->sm );
	if ( tlk->sp != NULL ) free_SitePattern( tlk->sp );
	free_SingleTreeLikelihood_internals(tlk);
}

SingleTreeLikelihood * clone_SingleTreeLikelihood( SingleTreeLikelihood *tlk ){
	SitePattern *sp = clone_SitePattern(tlk->sp);
	SubstitutionModel *m   = clone_substitution_model(tlk->m);
	SiteModel *sm   = clone_SiteModel(tlk->sm);
	Tree* tree = clone_Tree(tlk->tree);
	BranchModel* bm = NULL;
	if(tlk->bm != NULL){
		DiscreteParameter* dp = NULL;
		if (tlk->bm->ssvs != NULL) {
			dp = tlk->bm->ssvs->clone(tlk->bm->ssvs);
		}
		bm = clone_BranchModel(tlk->bm, tree, dp);
	}
	return clone_SingleTreeLikelihood_with( tlk, tree, m, sm, sp, bm  );
}

static void _singletreellikelihood_copy(SingleTreeLikelihood *source, SingleTreeLikelihood *dest, Node *snode, Node *dnode){
    if( !Node_isleaf(snode) ){
        _singletreellikelihood_copy(source, dest, snode->left, dnode->left);
        _singletreellikelihood_copy(source, dest, snode->right, dnode->right);
    }
}

void SingleTreeLikelihood_copy(SingleTreeLikelihood *source, SingleTreeLikelihood *dest){
    _singletreellikelihood_copy(source, dest, Tree_root(source->tree), Tree_root(dest->tree));
}


SingleTreeLikelihood * clone_SingleTreeLikelihood_with( SingleTreeLikelihood *tlk, Tree *tree, SubstitutionModel* m, SiteModel *sm, SitePattern *sp, BranchModel *bm){
	SingleTreeLikelihood *newtlk = (SingleTreeLikelihood *)malloc( sizeof(SingleTreeLikelihood));
	assert(newtlk);
	
	newtlk->tree = tree;
	newtlk->m   = m;
	newtlk->sm   = sm;
	newtlk->sp   = sp;
    newtlk->bm = bm;
	
	newtlk->cat_count = tlk->cat_count;
	newtlk->pattern_count = tlk->pattern_count;
	newtlk->use_tip_states = tlk->use_tip_states;
	
	newtlk->mapping = clone_ivector(tlk->mapping, Tree_node_count(tlk->tree));
		
	newtlk->partials_size = tlk->partials_size;
	newtlk->pattern_lk_size = tlk->pattern_lk_size;
	newtlk->root_partials_size = tlk->root_partials_size;
	
	newtlk->scaling_factors = NULL;
	if ( tlk->scaling_factors != NULL ){
		newtlk->scaling_factors = (double***)malloc(2*sizeof(double**));
		newtlk->scaling_factors[0] = clone_dmatrix( tlk->scaling_factors[0], tlk->partials_dim, tlk->sp->count );
        newtlk->scaling_factors[1] = NULL;
        if(tlk->scaling_factors[1] != NULL){
            newtlk->scaling_factors[1] = clone_dmatrix( tlk->scaling_factors[1], tlk->partials_dim, tlk->sp->count );
        }
	}
	
	newtlk->scale = tlk->scale;
	newtlk->scaling_threshold = tlk->scaling_threshold;
	
	newtlk->matrix_size = tlk->matrix_size;
	
	newtlk->update_nodes = clone_bvector( tlk->update_nodes, Tree_node_count(newtlk->tree));
	newtlk->update = tlk->update;
	newtlk->update_upper = tlk->update_upper;
	
	
	newtlk->pattern_lk = clone_dvector( tlk->pattern_lk, tlk->pattern_lk_size );
	newtlk->lk = tlk->lk;
	
	newtlk->calculate            = tlk->calculate;
	newtlk->update_partials      = tlk->update_partials;
	newtlk->update_partials_flexible = tlk->update_partials_flexible;
	newtlk->integrate_partials   = tlk->integrate_partials;
	newtlk->node_log_likelihoods = tlk->node_log_likelihoods;
	newtlk->calculate_per_cat_partials = tlk->calculate_per_cat_partials;
	
	newtlk->upper_partial_indexes = clone_ivector(tlk->upper_partial_indexes, Tree_node_count(tlk->tree));
	
	newtlk->node_id = tlk->node_id;
    
    newtlk->nthreads = tlk->nthreads;
    
	newtlk->use_upper = false;
    newtlk->calculate_upper = NULL;
    newtlk->node_upper = NULL;
	
	newtlk->matrix_dim = tlk->matrix_dim;
	newtlk->partials_dim = tlk->partials_dim;
	
	Node **nodes = Tree_get_nodes(newtlk->tree, POSTORDER);
	
	newtlk->use_upper = tlk->use_upper;
	newtlk->calculate_upper = tlk->calculate_upper;
	
	newtlk->node_upper = NULL;
	newtlk->tripod = tlk->tripod;
	
    
    // partials and matrices need to be aligned when compiled with GCC >= 4.7
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
    newtlk->use_SIMD = tlk->use_SIMD;
#endif
    
    size_t nodeCount = Tree_node_count(tlk->tree);
    
    allocate_storage(newtlk, 0);
    memcpy(newtlk->current_matrices_indexes, tlk->current_matrices_indexes, 2*nodeCount*sizeof(unsigned));
    memcpy(newtlk->current_partials_indexes, tlk->current_partials_indexes, 2*nodeCount*sizeof(unsigned));
	int i = 0;
	for ( ; i < nodeCount; i++ ) {
		if( newtlk->partials[0][Node_id(nodes[i])] != NULL ){
			memcpy(newtlk->partials[0][Node_id( nodes[i] )], tlk->partials[0][Node_id( nodes[i] )], tlk->partials_size * sizeof(double));
		}
	}
	for ( ; i < tlk->partials_dim; i++ ) {
		memcpy(newtlk->partials[0][i], tlk->partials[0][i], tlk->partials_size * sizeof(double));
	}

	for ( int i = 0; i < tlk->matrix_dim; i++ ) {
		memcpy(newtlk->matrices[0][i], tlk->matrices[0][i], tlk->matrix_size*tlk->sm->cat_count * sizeof(double));
	}
    
   if(tlk->partials[1] != NULL){
        allocate_storage(newtlk, 1);
        memcpy(newtlk->stored_matrices_indexes, tlk->stored_matrices_indexes, 2*nodeCount*sizeof(unsigned));
        memcpy(newtlk->stored_partials_indexes, tlk->stored_partials_indexes, 2*nodeCount*sizeof(unsigned));
        i = 0;
        for ( ; i < nodeCount; i++ ) {
            if( newtlk->partials[1][Node_id(nodes[i])] != NULL ){
                memcpy(newtlk->partials[1][Node_id( nodes[i] )], tlk->partials[1][Node_id( nodes[i] )], tlk->partials_size * sizeof(double));
            }
        }
        for ( ; i < tlk->partials_dim; i++ ) {
            memcpy(newtlk->partials[1][i], tlk->partials[1][i], tlk->partials_size * sizeof(double));
        }

        for ( int i = 0; i < tlk->matrix_dim; i++ ) {
            memcpy(newtlk->matrices[1][i], tlk->matrices[1][i], tlk->matrix_size*tlk->sm->cat_count * sizeof(double));
        }
   }
	
	newtlk->root_partials = (double*)malloc( newtlk->root_partials_size*sizeof(double) );
    assert(newtlk->root_partials);
	memcpy(newtlk->root_partials, tlk->root_partials, newtlk->root_partials_size * sizeof(double));
    
	newtlk->sp->ref_count++;
	
	newtlk->root_frequencies = NULL;
	if(tlk->root_frequencies != NULL){
		newtlk->root_frequencies = clone_dvector(tlk->root_frequencies, tlk->m->nstate);
	}
	newtlk->get_root_frequencies = tlk->get_root_frequencies;

	newtlk->include_jacobian = tlk->include_jacobian;
	
	newtlk->gradient = NULL;
	if(tlk->gradient != NULL){
		newtlk->gradient = clone_dvector(tlk->gradient, tlk->gradient_length);
		newtlk->gradient_length = tlk->gradient_length;
		newtlk->prepared_gradient = tlk->prepared_gradient;
	}
	newtlk->include_root_freqs = tlk->include_root_freqs;
	return newtlk;
}

bool SingleTreeLikelihood_rescaling( SingleTreeLikelihood *tlk ){
	return tlk->scale;
}

void SingleTreeLikelihood_use_rescaling( SingleTreeLikelihood *tlk, bool use ){
	tlk->scale = use;
	if ( tlk->scaling_factors == NULL && use ) {
		tlk->scaling_factors = (double***)malloc(2*sizeof(double**));
		tlk->scaling_factors[0] = dmatrix(tlk->partials_dim, tlk->sp->count );
		for(size_t i = 0; i < tlk->partials_dim; i++){
			memset(tlk->scaling_factors[0][i], 0.0, sizeof(double)*tlk->sp->count);
		}
        tlk->scaling_factors[1] = NULL;
        if(tlk->partials[1] != NULL){
            tlk->scaling_factors[1] = clone_dmatrix( tlk->scaling_factors[0], tlk->partials_dim, tlk->sp->count );
        }
	}
}

//TODO: Only works for codon and nucleotide without SSE or AVX
void SingleTreeLikelihood_set_nthreads( SingleTreeLikelihood *tlk, int nthreads ){
    #if defined (_OPENMP) && !((SSE3_ENABLED) || (AVX_ENABLED))
    //#if defined _OPENMP

     if( tlk->m->nstate >= 60 || tlk->m->nstate == 4 ){
        if( nthreads > 1 ){
            tlk->nthreads = nthreads;
            if( tlk->m->nstate == 4 ){
                tlk->update_partials = update_partials_4_openmp;
            }
            else {
                tlk->update_partials = update_partials_codon_openmp;
            }            
        }
        else {
            tlk->nthreads = 1;
            if( tlk->m->nstate == 4 ){
                tlk->update_partials = update_partials_4;
            }
            else {
                tlk->update_partials = update_partials_codon;
            }
        }
     }
    #endif
    
}

double _calculate_simple( SingleTreeLikelihood *tlk ){
#ifdef UPPER_PARTIALS
	printf("simple\n");
#endif
	if( !tlk->update ){
		return tlk->lk;
	}

	bool success = tlk->sm->update(tlk->sm);
	if(!success){
		tlk->lk = NAN;
		return NAN;
	}

	if(Tree_is_time_mode(tlk->tree)){
		Tree_update_heights(tlk->tree);
	}
	_calculate_partials( tlk, Tree_root(tlk->tree) );
	int nodeID = Node_id(Tree_root(tlk->tree));
	if(tlk->sm->integrate){
		tlk->integrate_partials(tlk, tlk->partials[tlk->current_partials_indexes[nodeID]][nodeID], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
	}
	else{
		memcpy(tlk->root_partials, tlk->partials[tlk->current_partials_indexes[nodeID]][nodeID], sizeof(double)*tlk->sp->count*tlk->sp->nstate);
	}
	
	tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
	
	tlk->lk = 0;
	int i = 0;
	
	for ( ; i < tlk->sp->count; i++) {
		tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
	}

	if(isnan(tlk->lk)){
//		printf("NAN %f\n", tlk->lk);
//		Tree_print_parameters(tlk->tree);exit(10);
		for ( i = 0; i < Tree_node_count(tlk->tree); i++) tlk->update_nodes[i] = true;
		tlk->update = true;
		tlk->update_upper = true;
	}
	else if( isinf(tlk->lk) ){
		fprintf(stdout, "_calculate: rescaling %f\n", tlk->lk);
		SingleTreeLikelihood_use_rescaling(tlk, true );
		
		SingleTreeLikelihood_update_all_nodes( tlk );
		tlk->lk = 0;
		_calculate_partials( tlk, Tree_root(tlk->tree) );
		
		int nodeID = Node_id(Tree_root(tlk->tree));
		if(tlk->sm->integrate){
			tlk->integrate_partials(tlk, tlk->partials[tlk->current_partials_indexes[nodeID]][nodeID], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
		}
		else{
			memcpy(tlk->root_partials, tlk->partials[tlk->current_partials_indexes[nodeID]][nodeID], sizeof(double)*tlk->sp->count*tlk->sp->nstate);
		}
	
		tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
		
		
		for ( i = 0; i < tlk->sp->count; i++) {
			tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
		}
		
	}
	
	for ( i = 0; i < Tree_node_count(tlk->tree); i++) tlk->update_nodes[i] = false;
	tlk->update = false;
	tlk->update_upper = true;
	//printf("calculate lower %.11f\n", tlk->lk);
	return tlk->lk;
}

//#define UPPER_PARTIALS 1

void SingleTreeLikelihood_update_uppers(SingleTreeLikelihood *tlk){
	_calculate_simple(tlk);
#ifdef UPPER_PARTIALS
	printf("calculate_upper root update\n");
#endif
	update_upper_partials(tlk, Tree_root(tlk->tree), false);
	tlk->update_upper = false;
	tlk->use_upper = true;
}

void SingleTreeLikelihood_update_uppers2(SingleTreeLikelihood *tlk){
	update_upper_partials2(tlk, Tree_root(tlk->tree));
	memset(tlk->update_nodes, 0, sizeof(bool)*Tree_node_count(tlk->tree));
	tlk->update_upper = false;
	tlk->use_upper = true;
}

// use_upper==true usage:
// tlk->calculate(tlk) should be called before using upper
// In brent optimization the first call will have no updates since it is evaluated at the current value
// If only one branch has changed then it is part of the brent optimization
// Once we optimize the next branch there should be 2 update_nodes[i]==true
double _calculate( SingleTreeLikelihood *tlk ){
	if (tlk->use_upper && tlk->tripod) {
		tlk->lk = tlk->calculate_upper(tlk, tlk->node_upper);
		return tlk->lk;
	}

	if (tlk->use_upper) {
#ifdef UPPER_PARTIALS
		printf("-----------------\n");
#endif
		// In brent the initial point is calculated without changing the parameter
		if(tlk->update){
			int updateCount = 0;
			int indexUpdate = -1;
			int nodeCount = Tree_node_count(tlk->tree);
			for (int i = 0; i < nodeCount; i++) {
				if(tlk->update_nodes[i]){
					updateCount++;
					indexUpdate = i;
				}
			}
			// only one branch has changed
			if (updateCount == 1) {
				tlk->lk = tlk->calculate_upper(tlk, Tree_node(tlk->tree, indexUpdate));
				tlk->node_upper = Tree_node(tlk->tree, indexUpdate);
#ifdef UPPER_PARTIALS
				printf("fast %d\n", indexUpdate);
#endif
				return tlk->lk;
			}
			
		}
		else{
#ifdef UPPER_PARTIALS
			printf("no update\n");
#endif
			return  tlk->lk;
		}
		tlk->use_upper = false;
		int nodeCount = Tree_node_count(tlk->tree);
		int indexUpdate = -1;
		int previousID = Node_id(tlk->node_upper);
		for (int i = 0; i < nodeCount; i++) {
			if(tlk->update_nodes[i] && i != previousID){
				indexUpdate = i;
			}
		}
#ifdef UPPER_PARTIALS
		printf("more than 1 branch %d\n", indexUpdate);
#endif
		tlk->update_nodes[previousID] = false;
		tlk->lk = tlk->calculate_upper(tlk, Tree_node(tlk->tree, indexUpdate));
		tlk->node_upper = Tree_node(tlk->tree, indexUpdate);
		tlk->use_upper = true;
		return tlk->lk;
	}
#ifdef UPPER_PARTIALS
	printf("standard\n");
#endif
	return _calculate_simple(tlk);
}

void SingleTreeLikelihood_update_Q(SingleTreeLikelihood* tlk, Node* n){
	double bl = Node_distance(n);
	bool node_transpose = tlk->partials[0][Node_id(n)] == NULL;
	
	for (int i = 0; i < tlk->sm->cat_count; i++) {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		if( tlk->use_SIMD ){
			if( node_transpose ){
				tlk->m->p_t_transpose(tlk->m,
										  bl * tlk->sm->get_rate(tlk->sm, i),
										  &tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
			}
			else {
				tlk->m->p_t(tlk->m,
								bl * tlk->sm->get_rate(tlk->sm, i),
								&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
			}
		}
		else{
			tlk->m->p_t(tlk->m,
							bl * tlk->sm->get_rate(tlk->sm, i),
							&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
		}
#else
		tlk->m->p_t(tlk->m,
						bl * tlk->sm->get_rate(tlk->sm, i),
						&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
#endif
	}
}

bool _calculate_partials( SingleTreeLikelihood *tlk, Node *n  ){
	bool updated = false;
	
	if( tlk->update_nodes[ Node_id(n) ] && !Node_isroot(n) ){
		// a bit hackish here, for local optimization of clades
		if ( tlk->node_id == -1 || Node_id(n) != tlk->node_id ) {
			
			double bl = 0;
			if( tlk->bm == NULL || !Tree_is_time_mode(tlk->tree) ){
				bl = Node_distance(n);
			}
			else{
				bl = tlk->bm->get(tlk->bm, n) * Node_time_elapsed(n);
				
				if(bl < 0 ){
					fprintf(stderr, "calculate_partials: %s branch length = %E rate = %f height = %f - parent height [%s]= %f (%f)\n", n->name, bl, tlk->bm->get(tlk->bm, n), Node_height(n), n->parent->name, Node_height(Node_parent(n)), Node_distance(n));
					exit(1);
				}
			}
			
			if (!tlk->use_upper && tlk->partials[1] != NULL) {
				tlk->current_matrices_indexes[Node_id(n)] = 1 - tlk->current_matrices_indexes[Node_id(n)];
				tlk->current_matrices_indexes[Node_id(n)+Tree_node_count(tlk->tree)] = tlk->current_matrices_indexes[Node_id(n)];
			}

			bool node_transpose = tlk->partials[0][Node_id(n)] == NULL;

			for (int i = 0; i < tlk->sm->cat_count; i++) {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
				if( tlk->use_SIMD ){
					if( node_transpose ){
						tlk->m->p_t_transpose(tlk->m,
												  bl * tlk->sm->get_rate(tlk->sm, i),
												  &tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
					}
					else {
						tlk->m->p_t(tlk->m,
										bl * tlk->sm->get_rate(tlk->sm, i),
										&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
					}
				}
				else{
					tlk->m->p_t(tlk->m,
									bl * tlk->sm->get_rate(tlk->sm, i),
									&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
				}
#else
				tlk->m->p_t(tlk->m,
								bl * tlk->sm->get_rate(tlk->sm, i),
								&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
#endif
			}
			
			updated = true;
		}
	}
	
	if( !Node_isleaf(n) ){
		bool update_child1 = _calculate_partials( tlk, Node_left(n) );
		bool update_child2 = _calculate_partials( tlk, Node_right(n) );
		
		if( update_child1 || update_child2 ){
			int indx_child1 = Node_id(Node_left(n));
			int indx_child2 = Node_id(Node_right(n));

			if (!tlk->use_upper && tlk->partials[1] != NULL) {
				tlk->current_partials_indexes[Node_id(n)] = 1 - tlk->current_partials_indexes[Node_id(n)];
				tlk->current_partials_indexes[Node_id(n)+Tree_node_count(tlk->tree)] = tlk->current_partials_indexes[Node_id(n)];
			}

			tlk->update_partials( tlk, Node_id(n), indx_child1, indx_child1, indx_child2, indx_child2);
			
			// a bit hackish here, for local optimization of clades
			if( Node_id(n) == tlk->node_id ){
				int nodeID = Node_id(n);
				if(tlk->sm->integrate){
					tlk->integrate_partials(tlk, tlk->partials[tlk->current_partials_indexes[nodeID]][nodeID], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
				}
				else{
					memcpy(tlk->root_partials, tlk->partials[tlk->current_partials_indexes[nodeID]][nodeID], sizeof(double)*tlk->sp->count*tlk->sp->nstate);
				}
				
				tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
				
			}
			updated = true;
		}
	}
	return updated;
}

// update every node
void SingleTreeLikelihood_update_all_nodes( SingleTreeLikelihood *tlk ){
	for (int index = 0; index < Tree_node_count(tlk->tree); index++) {
		tlk->update_nodes[index] = true;
	}
	tlk->update = true;
	tlk->update_upper = true;
	tlk->node_upper = NULL; // does not hurt to set it even if we don't use upper likelihoods
}

// update 1 node
void SingleTreeLikelihood_update_one_node( SingleTreeLikelihood *tlk, const Node *node ){
	tlk->update_nodes[Node_id(node)] = true;
	tlk->update = true;
	tlk->update_upper = true;
}

// update node and its children
//void SingleTreeLikelihood_update_three_nodes( SingleTreeLikelihood *tlk, const int index ){
void SingleTreeLikelihood_update_three_nodes( SingleTreeLikelihood *tlk, const Node *node ){
	
	//Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
	
	tlk->update_nodes[Node_id(node)] = true;
	
	if( !Node_isleaf(node) ){
		tlk->update_nodes[Node_id(node->right)] = true;
		tlk->update_nodes[Node_id(node->left)]  = true;
	}
	tlk->update = true;
	tlk->update_upper = true;
	
#ifdef USE_UPPER_PARTIALS
	//tlk->node_upper = NULL;
#endif
}


/**
 * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
 * Yang (2000) J. Mol. Evol. 51: 423-432
 * <p/>
 * This function looks over the partial likelihoods for each state at each pattern
 * and finds the largest. If this is less than the scalingThreshold (currently set
 * to 1E-40) then it rescales the partials for that pattern by dividing by this number
 * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
 * This is called for every internal node after the partials are calculated so provides
 * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
 * but this sounded like a headache to organize (and he doesn't use the threshold idea
 * which improves the performance quite a bit).
 *
 * @param nodeIndex
 */

void SingleTreeLikelihood_scalePartials( SingleTreeLikelihood *tlk, int nodeIndex, int childIndex1, int childIndex2 ) {
	int u = 0;
	int k,j;
	double* partials = tlk->partials[tlk->current_partials_indexes[nodeIndex]][nodeIndex];
	double* scaling_factors = tlk->scaling_factors[tlk->current_partials_indexes[nodeIndex]][nodeIndex];
	double* scaling_factors1 = childIndex1 >=0 && tlk->partials[0][childIndex1] != NULL ? tlk->scaling_factors[tlk->current_partials_indexes[childIndex1]][childIndex1] : NULL;
	double* scaling_factors2 = childIndex2 >=0 && tlk->partials[0][childIndex2] != NULL ? tlk->scaling_factors[tlk->current_partials_indexes[childIndex2]][childIndex2] : NULL;

	for ( int i = 0; i < tlk->sp->count; i++ ) {
		
		double scaleFactor = 0.0;
		int v = u;
		for ( k = 0; k < tlk->cat_count; k++ ) {
			for ( j = 0; j < tlk->m->nstate; j++ ) {
				if ( partials[v] > scaleFactor ) {
					scaleFactor = partials[v];
				}
				v++;
			}
			v += ( tlk->sp->count - 1 ) * tlk->m->nstate;
		}

		if ( scaleFactor < tlk->scaling_threshold ) {
			
			v = u;
			for ( k = 0; k < tlk->cat_count; k++ ) {
				for ( j = 0; j < tlk->m->nstate; j++ ) {
					partials[v] /= scaleFactor;
					v++;
				}
				v += (tlk->sp->count - 1) * tlk->m->nstate;
			}
			scaling_factors[i] = log(scaleFactor);
		}
		else {
			scaling_factors[i] = 0.0;
		}

		if(scaling_factors1 != NULL){
			scaling_factors[i] += scaling_factors1[i];
		}
		if(scaling_factors2 != NULL){
			scaling_factors[i] += scaling_factors2[i];
		}
		u += tlk->m->nstate;
	}
}

double getLogScalingFactor( const SingleTreeLikelihood *tlk, int pattern ) {
	double log_scale_factor = 0.0;
	if ( tlk->scale ) {
		size_t rootId = Tree_root(tlk->tree)->id;
		return tlk->scaling_factors[tlk->current_partials_indexes[rootId]][rootId][pattern];
	}
	return log_scale_factor;
}

double getLogScalingFactorAtNode( const SingleTreeLikelihood *tlk, int pattern, size_t index ) {
	double log_scale_factor = 0.0;
	if ( tlk->scale && tlk->scaling_factors[tlk->current_partials_indexes[index]][index] != NULL) {
		log_scale_factor = tlk->scaling_factors[tlk->current_partials_indexes[index]][index][pattern];
	}
	return log_scale_factor;
}

// TODO: add calculate_branch_likelihood for codon and proteins and X even
void SingleTreeLikelihood_enable_SSE( SingleTreeLikelihood *tlk, bool value ){
#ifdef SSE3_ENABLED
	tlk->use_SIMD = value;
	if ( tlk->use_SIMD ) {
		//		if( !(tlk->sm->nstate >= 60 || tlk->sm->nstate == 4 ||tlk->sm->nstate == 20 ) ){
		//			tlk->use_SIMD = false;
		//		}
		// Does not work for matrix with odd dimension.
		// If we want to load doubles at indexes 1 and 2 from matrices it will crash
		// same issue with partials, it will try to store at odd indexes
		if( tlk->m->nstate & 1 ){
			tlk->use_SIMD = false;
		}
	}
	
	if ( tlk->use_SIMD ){
		if (tlk->m->nstate == 4 ) {
			tlk->update_partials      = update_partials_4_SSE;
			tlk->update_partials_flexible = update_partials_flexible_4_SSE;
			tlk->integrate_partials   = integrate_partials_4_SSE;
			tlk->node_log_likelihoods = node_log_likelihoods_4_SSE;
			tlk->calculate_per_cat_partials = calculate_branch_partials_4_SSE;
		}
		else if ( tlk->m->nstate == 20 ) {
			tlk->update_partials      = update_partials_20_SSE;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			tlk->calculate_per_cat_partials = calculate_branch_partials_20_SSE;
		}
		else if ( tlk->m->nstate >= 60 ) {
			//			tlk->update_partials      = update_partials_codon_SSE;
			//			tlk->integrate_partials   = integrate_partials_codon_SSE;
			//			tlk->node_log_likelihoods = node_log_likelihoods_codon_SSE;
			
			tlk->update_partials      = update_partials_codon;
			tlk->integrate_partials   = integrate_partials_codon;
			tlk->node_log_likelihoods = node_log_likelihoods_codon;
			
			//            tlk->update_partials      = update_partials_general_transpose;
			//            tlk->integrate_partials   = integrate_partials_general;
			//            tlk->node_log_likelihoods = node_log_likelihoods_general;;
		}
		else {
			tlk->update_partials      = update_partials_general_even_SSE;
			tlk->integrate_partials   = integrate_partials_general_even_SSE;
			tlk->node_log_likelihoods = node_log_likelihoods_general_even_SSE;
			tlk->calculate_per_cat_partials = calculate_branch_partials_even_SSE;
			//TODO: what about the odd case
			// no upper
		}
	}
	else {
		if ( tlk->m->nstate == 4 ) {
			tlk->update_partials      = update_partials_4;
			tlk->update_partials_flexible = update_partials_flexible_4;
			tlk->integrate_partials   = integrate_partials_4;
			tlk->node_log_likelihoods = node_log_likelihoods_4;
			tlk->calculate_per_cat_partials = calculate_branch_partials_4;
		}
		else if ( tlk->m->nstate == 20 ) {
			tlk->update_partials      = update_partials_general;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			tlk->calculate_per_cat_partials = calculate_branch_partials_20;
		}
		else if( tlk->m->nstate >= 60 ){
			tlk->update_partials      = update_partials_codon;
			tlk->integrate_partials   = integrate_partials_codon;
			tlk->node_log_likelihoods = node_log_likelihoods_codon;
		}
		else {
			tlk->update_partials      = update_partials_general;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			tlk->calculate_per_cat_partials = calculate_branch_partials;
		}
	}
	SingleTreeLikelihood_update_all_nodes(tlk);
#endif
}

bool SingleTreeLikelihood_SSE( SingleTreeLikelihood *tlk ){
#ifdef SSE3_ENABLED
	return tlk->use_SIMD;
#else
	return false;
#endif
	
}

const double* get_root_frequencies(SingleTreeLikelihood* tlk){
	return tlk->m->get_frequencies(tlk->m);
}

const double* get_root_frequencies_fixed(SingleTreeLikelihood* tlk){
	assert(tlk->root_frequencies);
	return tlk->root_frequencies;
}


#pragma mark -
#pragma mark Lower Likelihood

double SingleTreeLikelihood_calculate_at_node( SingleTreeLikelihood *tlk, const Node *node ){
    if( tlk->update ) {
		tlk->calculate(tlk);
	}
	
	double *partials   = dvector( tlk->sp->count*tlk->m->nstate );
	double *pattern_lk = dvector(tlk->sp->count);
	
	tlk->integrate_partials(tlk, tlk->partials[tlk->current_matrices_indexes[Node_id(node)]][Node_id(node)], tlk->sm->get_proportions(tlk->sm), partials );
	
	tlk->node_log_likelihoods( tlk, partials, tlk->get_root_frequencies(tlk), pattern_lk);
	
	double lnl = 0;
	for ( int i = 0; i < tlk->sp->count; i++) {
		lnl += pattern_lk[i] * tlk->sp->weights[i];
	}
	free(pattern_lk);
	free(partials);
	
	return lnl;
}


// Should be used with CARE
// compute the likelihood of a subtree, highjacking pattern_lk and root_partials in _calculate_partials
static double _calculate2( SingleTreeLikelihood *tlk ){
	if( !tlk->update ) return tlk->lk;
	
	_calculate_partials( tlk, Tree_node(tlk->tree, tlk->node_id) );
	
	int nodeID = Node_id(Tree_root(tlk->tree));
    tlk->integrate_partials(tlk, tlk->partials[tlk->current_matrices_indexes[nodeID]][nodeID], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
    
    tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
    
	tlk->lk = 0;
	int i = 0;
	
	for ( ; i < tlk->sp->count; i++) {
		tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
	}
	
	if ( tlk->lk == -INFINITY ) {
		fprintf(stderr, "_calculate: rescaling\n");
		SingleTreeLikelihood_use_rescaling(tlk, true );
		
		SingleTreeLikelihood_update_all_nodes( tlk );
		tlk->lk = 0;
		_calculate_partials( tlk, Tree_node(tlk->tree, tlk->node_id) );
		
		int nodeID = Node_id(Tree_root(tlk->tree));
		tlk->integrate_partials(tlk, tlk->partials[tlk->current_matrices_indexes[nodeID]][nodeID], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
        
        tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
		
		for ( i = 0; i < tlk->sp->count; i++) {
			tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
		}
		
	}
	
	for ( i = 0; i < Tree_node_count(tlk->tree); i++) tlk->update_nodes[i] = false;
	tlk->update = false;
	tlk->update_upper = true;
	return tlk->lk;
}

#pragma mark -

bool _update_all_partials_aux(SingleTreeLikelihood *tlk, Node* node, bool* flags){
	bool f = tlk->update_nodes[Node_id(node)];
	if(!Node_isleaf(node)){
		f |= _update_all_partials_aux(tlk, Node_left(node), flags);
		f |= _update_all_partials_aux(tlk, Node_right(node), flags);
	}
	flags[Node_id(node)] = f;
	return f;
}
		

static void _update_upper_partials_efficient(SingleTreeLikelihood *tlk, Node* node, bool* flags, bool p){
	
	if(!Node_isroot(node)){
		Node* parent = Node_parent(node);
		Node* sibling = Node_sibling(node);
		int idNode = Node_id(node);
		int idSibling = Node_id(sibling);
		bool need_update = false;
		
		tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
		
		if(!Node_isroot(parent)){
			int idMatrix = Node_id(parent);
			
			//printf("%s %d %d %d\n", node->name, tlk->upper_partial_indexes[idNode], tlk->upper_partial_indexes[Node_id(parent)], idMatrix);
			need_update = p || flags[idSibling] || flags[Node_id(parent)];
			
			if(need_update){
				tlk->update_partials(tlk,
									 tlk->upper_partial_indexes[idNode],
									 tlk->upper_partial_indexes[Node_id(parent)],
									 idMatrix,
									 idSibling,
									 idSibling);
			}
		}
		// We dont need to calculate upper partials for the children of the root as it is using the lower partials of its sibling
		else{
			tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], idSibling, idSibling, -1, -1);
			need_update = p;
		}
		
		if(!Node_isleaf(node)){
			_update_upper_partials_efficient(tlk, Node_left(node), flags, need_update);
			_update_upper_partials_efficient(tlk, Node_right(node), flags, need_update);
		}
	}
	else{
		tlk->upper_partial_indexes[Node_id(node)] = -1;
		
		if(!Node_isleaf(node)){
			_update_upper_partials_efficient(tlk, Node_left(node), flags, flags[Node_id(Node_right(node))]);
			_update_upper_partials_efficient(tlk, Node_right(node), flags, flags[Node_id(Node_left(node))]);
		}
	}
}
		
// try to update only what's needed. NOT TESTED
void update_all_partials(SingleTreeLikelihood *tlk, Node* node){
	bool* flags = bvector(Tree_node_count(tlk->tree));
	_update_all_partials_aux(tlk, Tree_root(tlk->tree), flags);
	_calculate_partials(tlk, Tree_root(tlk->tree));
	_update_upper_partials_efficient(tlk, Tree_root(tlk->tree), flags, false);
	free(flags);
}
// #define INCLUDE_ROOT_FREQS 1

void inplace_mul_partials_frequencies(size_t dim, const double* frequencies, size_t size, double* partials){
	if(dim == 4){
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		double f[4] __attribute__ ((aligned (16)));
		memcpy(f, frequencies, sizeof(double)*4);
		__m128d f1 = _mm_load_pd(f);
		__m128d f2 = _mm_load_pd(f + 2);
		__m128d* p = (__m128d*)partials;
		for(size_t i = 0; i < size/4; i++){
			*p = _mm_mul_pd(*p, f1);
			p++;
			*p = _mm_mul_pd(*p, f2);
			p++;
		}
#else
		for(size_t i = 0; i < size; i+=4){
			partials[i] *= frequencies[0];
			partials[i+1] *= frequencies[1];
			partials[i+2] *= frequencies[2];
			partials[i+3] *= frequencies[3];
		}
#endif
	}
	else{
		for(size_t i = 0; i < size; ){
			for(size_t j = 0; j < dim; j++, i++){
				partials[i] *= frequencies[j];
			}
		}
	}
}

// calculate upper partials as used by derivatives
void update_upper_partials(SingleTreeLikelihood *tlk, Node* node, bool include_root_freqs){
	int idNode = Node_id(node);
	tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
	
	if(!Node_isroot(node)){
		Node* parent = Node_parent(node);
		Node* sibling = Node_sibling(node);
		int idSibling = Node_id(sibling);
		
		tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
		
		if(!Node_isroot(parent)){
			int idMatrix = Node_id(parent);
			// u_n = (P_p u_p) o (P_s p_s)
			tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], tlk->upper_partial_indexes[Node_id(parent)], idMatrix, idSibling, idSibling);
		}
		else{
			// u_n = P_s u_s
			tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], idSibling, idSibling, -1, -1);
			if(include_root_freqs){
				const double* freqs = tlk->get_root_frequencies(tlk);
				double* partials = tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[idNode]]][tlk->upper_partial_indexes[idNode]];
				size_t dim = tlk->cat_count*tlk->pattern_count*tlk->m->nstate;
				inplace_mul_partials_frequencies(tlk->m->nstate, freqs, dim, partials);
			}
		}
	}
	
	if(!Node_isleaf(node)){
		update_upper_partials(tlk, Node_left(node), include_root_freqs);
		update_upper_partials(tlk, Node_right(node), include_root_freqs);
	}
}

// only update uppers when required : update_nodes[i] == true
void update_upper_partials2(SingleTreeLikelihood *tlk, Node* node){
	int idNode = Node_id(node);
	tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
	
	if(!Node_isroot(node)){
		Node* parent = Node_parent(node);
		Node* sibling = Node_sibling(node);
		int idSibling = Node_id(sibling);
		
		if(tlk->update_nodes[idNode]){
			if(!Node_isroot(parent)){
				int idMatrix = Node_id(parent);
				
				tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], tlk->upper_partial_indexes[Node_id(parent)], idMatrix, idSibling, idSibling);
			}
			else{
				tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], idSibling, idSibling, -1, -1);
			}
		}
	}
	
	if(!Node_isleaf(node)){
		update_upper_partials2(tlk, Node_left(node));
		update_upper_partials2(tlk, Node_right(node));
	}
}


// Calculate derivative of likelihood with respect to branch length at node
// using the new upper likelihood calculation
// use a temporary matrix (the root's)
void calculate_dldt_uppper( SingleTreeLikelihood *tlk, Node *node, double* pattern_dlikelihoods ){
	
	double bl = Node_distance(node);
	const int nodeId = Node_id(node);
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tlk->current_matrices_indexes[tempMatrixId]][tempMatrixId];
	
	if(Tree_root(tlk->tree)->right == node){
		printf("calculate_dldt_uppper: Do not call on root\n");
		exit(1);
	}
	bool node_transpose = tlk->partials[0][Node_id(node)] == NULL;
	for (int i = 0; i < tlk->sm->cat_count; i++) {
		
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		if( tlk->use_SIMD ){
			if( node_transpose ){
				tlk->m->dp_dt_transpose(tlk->m,
											bl * tlk->sm->get_rate(tlk->sm, i),
											&mat[i*tlk->matrix_size]);
			}
			else {
				tlk->m->dp_dt(tlk->m,
								  bl * tlk->sm->get_rate(tlk->sm, i),
								  &mat[i*tlk->matrix_size]);
			}
		}
		else{
			tlk->m->dp_dt(tlk->m,
							  bl * tlk->sm->get_rate(tlk->sm, i),
							  &mat[i*tlk->matrix_size]);
		}
#else
		tlk->m->dp_dt(tlk->m,
						  bl * tlk->sm->get_rate(tlk->sm, i),
						  &mat[i*tlk->matrix_size]);
#endif
		for (int k = 0; k < tlk->matrix_size; k++) {
			mat[i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
		}
	}
	
	
	const int nstate   = tlk->m->nstate;
	const int patternCount = tlk->sp->count;
	double* root_partials = tlk->root_partials  + patternCount*nstate;
	const double* freqs = tlk->get_root_frequencies(tlk);
	
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
	tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
	
	int v = 0;
	
	for ( int k = 0; k < patternCount; k++ ) {
		
		pattern_dlikelihoods[k] = 0;
		for ( int i = 0; i < nstate; i++ ) {
			pattern_dlikelihoods[k] += freqs[i] * root_partials[v];
			v++;
		}
		
		if ( tlk->scale ) {
			//printf("scaling\n");
			pattern_dlikelihoods[k] += getLogScalingFactor( tlk, k);
		}
	}
}

// Assumes that pattern_dlnl pattern_lnl and  are calculated correctly (i.e. _calculate_uppper is not called)
// This function does not change anything regarding partials, matrices, root_partials and pattern_lk
// TODO: Make thread safe
double d2lnldt2_uppper( SingleTreeLikelihood *tlk, Node *node, const double* pattern_likelihoods, const double* pattern_dlikelihoods){
	
	const double bl = Node_distance(node);
	const int nodeId = Node_id(node);
	const int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tlk->current_matrices_indexes[tempMatrixId]][tempMatrixId];
	bool node_transpose = tlk->partials[0][Node_id(node)] == NULL;
	
	for (int i = 0; i < tlk->sm->cat_count; i++) {
		const double rate = tlk->sm->get_rate(tlk->sm, i);
		
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		if( tlk->use_SIMD ){
			if( node_transpose ){
				tlk->m->d2p_d2t_transpose(tlk->m,
											  bl * rate,
											  &mat[i*tlk->matrix_size]);
			}
			else {
				tlk->m->d2p_d2t(tlk->m,
									bl * rate,
									&mat[i*tlk->matrix_size]);
			}
		}
		else{
			tlk->m->d2p_d2t(tlk->m,
								bl * rate,
								&mat[i*tlk->matrix_size]);
		}
#else
		tlk->m->d2p_d2t(tlk->m,
							bl * rate,
							&mat[i*tlk->matrix_size]);
#endif
		for (int k = 0; k < tlk->matrix_size; k++) {
			mat[i*tlk->matrix_size+k] *= rate * rate;
		}
	}
	const int nstate   = tlk->m->nstate;
	const int patternCount = tlk->sp->count;
	double* root_partials = tlk->root_partials  + (patternCount*nstate)*2;
	double* pattern_d2lnl = tlk->pattern_lk + patternCount*3;
	const double* freqs = tlk->get_root_frequencies(tlk);
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
	tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
	
	//tlk->node_log_likelihoods( tlk, root_partials, tlk->m->_freqs, pattern_d2lnl);
	int v = 0;
	double d2lnl = 0;
	
	for ( int k = 0; k < patternCount; k++ ) {
		
		pattern_d2lnl[k] = 0;
		for ( int i = 0; i < nstate; i++ ) {
			pattern_d2lnl[k] += freqs[i] * root_partials[v];
			v++;
		}
		
		if ( tlk->scale ) {
			//printf("scaling\n");
			pattern_d2lnl[k] += getLogScalingFactor( tlk, k);
		}
		//((tlk->pattern_lk[i] * lks[i] - (dfi[Node_id(n)][i]*dfi[Node_id(n)][i])) / (lks[i]*lks[i])) * tlk->sp->weights[i];
		d2lnl += ((pattern_d2lnl[k]*pattern_likelihoods[k] - pattern_dlikelihoods[k]*pattern_dlikelihoods[k]) / (pattern_likelihoods[k]*pattern_likelihoods[k])) * tlk->sp->weights[k];
	}
	
	return d2lnl;
}

double calculate_dlnl_dQ( SingleTreeLikelihood *tlk, int index, const double* pattern_likelihoods ){
	
	Node **nodes = Tree_nodes(tlk->tree);
	
	const int rootId = Node_id(Tree_root(tlk->tree));
	const int rightId = Node_id(Node_right(Tree_root(tlk->tree)));
	
	const int nstate   = tlk->m->nstate;
	const int patternCount = tlk->sp->count;
	bool time_mode = Tree_is_time_mode(tlk->tree);
	
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[0][tempMatrixId];
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	int mat_len = tlk->matrix_size*tlk->sm->cat_count;
    double* mat_transpose =  aligned16_malloc( mat_len * sizeof(double) );
#endif
	double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->m->nstate;
	double* root_partials2 = tlk->root_partials  + tlk->sp->count*tlk->m->nstate*2;
	double* pattern_dlnl = tlk->pattern_lk + tlk->sp->count*2; // pattern_likelihoods points to tlk->pattern_lk + tlk->sp->count
	memset(pattern_dlnl, 0, sizeof(double)*tlk->sp->count);
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[rootId]];
    
    double dlnldQ = 0;
    tlk->m->dQ_need_update = true;
	
	size_t rateCount = tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K - 1;
	if(tlk->m->rates_simplex != NULL && !tlk->m->grad_wrt_reparam){
		rateCount++;
	}
    
	const double* freqs = tlk->get_root_frequencies(tlk);

    // Frequencies
    if(index >= rateCount){
		size_t freqIndex = index - rateCount;
		if(freqs != tlk->root_frequencies){
			double* dphi = dvector(tlk->m->nstate);
			if(tlk->m->grad_wrt_reparam){
				tlk->m->simplex->gradient(tlk->m->simplex, freqIndex, dphi);
			}
			else{
				dphi[freqIndex] = 1.0;
			}
			
			size_t v = 0;
			if(tlk->scale){
				for ( size_t k = 0; k < patternCount; k++ ) {
					double likelihood = 0;
					double dpartial = 0;
					for (size_t i = 0; i < nstate; i++, v++) {
						dpartial += dphi[i]*tlk->root_partials[v];
						likelihood += freqs[i]*tlk->root_partials[v];
					}
					pattern_dlnl[k] = dpartial/likelihood;
				}
			}
			else{
				for ( size_t k = 0; k < patternCount; k++ ) {
					pattern_dlnl[k] = 0;
					for (size_t i = 0; i < nstate; i++, v++) {
						pattern_dlnl[k] += dphi[i]*tlk->root_partials[v];
					}
				}
			}
			free(dphi);
		}
		
        for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
            Node *n = nodes[j];
            const int nodeId = Node_id(n);
            
            if ( nodeId == rootId || (nodeId == rightId && !time_mode)) {
                continue;
            }
            
			double bl;
			if (time_mode) {
				bl = Node_time_elapsed(n) * tlk->bm->get(tlk->bm, n);
			}
			else{
				bl = Node_distance(n);
			}
			
			for (int i = 0; i < tlk->sm->cat_count; i++) {
				tlk->m->dPdp(tlk->m, index, &mat[i*tlk->matrix_size], bl* tlk->sm->get_rate(tlk->sm, i));
			}
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			if( tlk->use_SIMD ){
				bool node_transpose = tlk->partials[0][Node_id(n)] == NULL;
				if(node_transpose){
					for(int l = 0; l < tlk->sm->cat_count; l++){
						for(int ii = 0; ii < nstate; ii++){
							for(int jj = 0; jj < nstate; jj++){
								mat_transpose[l * tlk->matrix_size+ii*nstate+jj] = mat[l * tlk->matrix_size+jj*nstate+ii];
							}
						}
					}
					memcpy(mat, mat_transpose, sizeof(double)*tlk->matrix_size*tlk->sm->cat_count);
				}
			}
#endif
			
			tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
			tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
			
			size_t v = 0;
			if(tlk->include_root_freqs){
				for ( int k = 0; k < patternCount; k++ ) {
					for (int jj = 0; jj < nstate; jj++) {
						pattern_dlnl[k] += root_partials[v];
						v++;
					}
				}
			}
			else{
				if(tlk->scale){
					tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, nodeId );
					tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials2 );
					for ( size_t k = 0; k < patternCount; k++ ) {
						double dlikelihood = 0;
						double likelihood = 0;
						for ( size_t i = 0; i < nstate; i++, v++ ) {
							dlikelihood += freqs[i] * root_partials[v];
							likelihood += freqs[i] * root_partials2[v];
						}
						pattern_dlnl[k] += dlikelihood/likelihood;
					}
				}
				else{
					for ( size_t k = 0; k < patternCount; k++ ) {
						for ( size_t i = 0; i < nstate; i++ ) {
							pattern_dlnl[k] += freqs[i] * root_partials[v];
							v++;
						}
					}
				}
			}
        }
    }
    // Rates
    else{
        for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
            Node *n = nodes[j];
            const int nodeId = Node_id(n);
            
            if ( nodeId == rootId || (nodeId == rightId && !time_mode)) {
                continue;
            }
            
			double bl;
			if (time_mode) {
				bl = Node_time_elapsed(n) * tlk->bm->get(tlk->bm, n);
			}
			else{
				bl = Node_distance(n);
			}
            
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                tlk->m->dPdp(tlk->m, index, &mat[i*tlk->matrix_size], bl* tlk->sm->get_rate(tlk->sm, i));
            }
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD ){
				bool node_transpose = tlk->partials[0][Node_id(n)] == NULL;
                if(node_transpose){
                    for(int l = 0; l < tlk->sm->cat_count; l++){
                        for(int ii = 0; ii < nstate; ii++){
                            for(int jj = 0; jj < nstate; jj++){
                                mat_transpose[l * tlk->matrix_size+ii*nstate+jj] = mat[l * tlk->matrix_size+jj*nstate+ii];
                            }
                        }
                    }
                    memcpy(mat, mat_transpose, sizeof(double)*tlk->matrix_size*tlk->sm->cat_count);
                }
            }
#endif

            tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
			tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
            int v = 0;
			int i = 0;
			const double* freqs = tlk->get_root_frequencies(tlk);

			if(tlk->include_root_freqs){
				if(tlk->scale){
					tlk->calculate_per_cat_partials(tlk, spare_partials, nodeId, tlk->upper_partial_indexes[nodeId], nodeId );
					tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials2 );
					for ( size_t k = 0; k < patternCount; k++ ) {
						double dlikelihood = 0;
						double likelihood = 0;
						for ( i = 0; i < nstate; i++, v++ ) {
							dlikelihood += root_partials[v];
							likelihood += root_partials2[v];
						}
						pattern_dlnl[k] += dlikelihood/likelihood;
						// double logP = log(likelihood) + getLogScalingFactorAtNode(tlk, k, nodeId) + getLogScalingFactorAtNode(tlk, k, tlk->upper_partial_indexes[nodeId] );
						// printf("%d] logP %f %f (%f)\n", nodeId, log(likelihood), logP, tlk->lk);
					}
				}
				else{
					for ( size_t k = 0; k < patternCount; k++ ) {
						for ( i = 0; i < nstate; i++, v++ ) {
							pattern_dlnl[k] += root_partials[v];
						}
					}
				}
			}
			else{
				if(tlk->scale){
					tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, nodeId );
					tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials2 );
					for ( size_t k = 0; k < patternCount; k++ ) {
						double dlikelihood = 0;
						double likelihood = 0;
						for ( i = 0; i < nstate; i++, v++ ) {
							dlikelihood += freqs[i] * root_partials[v];
							likelihood += freqs[i] * root_partials2[v];
						}
						pattern_dlnl[k] += dlikelihood/likelihood;
					}
				}
				else{
					for ( size_t k = 0; k < patternCount; k++ ) {
						for ( i = 0; i < nstate; i++ ) {
							pattern_dlnl[k] += freqs[i] * root_partials[v];
							v++;
						}
					}
				}
			}
        }
    }
	if(tlk->scale){
		for ( size_t k = 0; k < patternCount; k++ ) {
			dlnldQ += pattern_dlnl[k] * tlk->sp->weights[k];
		}
	}
	else{
		for ( size_t k = 0; k < patternCount; k++ ) {
			dlnldQ += pattern_dlnl[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
		}
	}
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	free(mat_transpose);
#endif
	return dlnldQ;
}

#pragma mark -
#pragma mark Upper Likelihood


// should be used when we optimize distances
// should notbe called at the root
// The last node in post order always need an update and the partials at the root too
double _calculate_uppper( SingleTreeLikelihood *tlk, Node *node ){
	// change of node
	if ( tlk->node_upper != node ) {
		if ( tlk->node_upper == NULL || ( Node_isleaf(node) && Node_sibling(node) != tlk->node_upper &&  Node_right(node->parent) == node ) || ( !Node_isleaf(node) && Node_right(node) != tlk->node_upper) ) {
			//  make sure lower partials are calculated
			_calculate_simple(tlk);
			update_upper_partials(tlk, Tree_root(tlk->tree), false);
		}
		else {
			Node *n = tlk->node_upper;
			
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			bool node_transpose = tlk->partials[0][Node_id(n)] == NULL;
			if( tlk->use_SIMD && node_transpose ){
				for (int i = 0; i < tlk->sm->cat_count; i++) {
					tlk->m->p_t_transpose(tlk->m,
											  Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
											  &tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
				}
			}
			else{
				for (int i = 0; i < tlk->sm->cat_count; i++) {
					tlk->m->p_t(tlk->m,
									Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
									&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
				}
			}
#else
			for (int i = 0; i < tlk->sm->cat_count; i++) {
				tlk->m->p_t(tlk->m,
								Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
								&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
			}
#endif
			
			// update lower partials of the parent
			tlk->update_partials(tlk, Node_id(Node_parent(n)), Node_id(n), Node_id(n), Node_id(Node_sibling(n)), Node_id(Node_sibling(n)) );
			
			// update upper partials of its sibling and its sibling's descendents
			// although they technically don't need to be updated if node is a right node since its sibling&descendents won't be reused in postorder

			update_upper_partials(tlk, Node_sibling(n), false);
			
		}
	}
	//		SingleTreeLikelihood_update_all_nodes(tlk);
	//		tlk->calculate(tlk);
	//		calculate_upper(tlk, Tree_root(tlk->tree));
	
	const int nodeId = Node_id(node);
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tlk->current_matrices_indexes[tempMatrixId]][tempMatrixId];
	
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	bool node_transpose = tlk->partials[0][Node_id(node)] == NULL;
	if( tlk->use_SIMD && node_transpose ){
		for (int i = 0; i < tlk->sm->cat_count; i++) {
			tlk->m->p_t_transpose(tlk->m,
									  Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
									  &mat[i*tlk->matrix_size]);
		}
	}
	else{
		for (int i = 0; i < tlk->sm->cat_count; i++) {
			tlk->m->p_t(tlk->m,
							Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
							&mat[i*tlk->matrix_size]);
		}
	}
#else
	for (int i = 0; i < tlk->sm->cat_count; i++) {
		tlk->m->p_t(tlk->m,
						Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
						&mat[i*tlk->matrix_size]);
	}
#endif
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
	if(tlk->sm->integrate){
		tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
	}
	else{
		memcpy(tlk->root_partials, spare_partials, sizeof(double)*tlk->sp->count*tlk->sp->nstate);
	}
	tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
	
	double lk = 0;
	for ( int i = 0; i < tlk->sp->count; i++) {
		lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
	}
	
	// set for update later
	//SingleTreeLikelihood_update_one_node(tlk, node);
	return lk;
}


// should be used for debugging
double calculate_log_likelihood_from_preorder(SingleTreeLikelihood* tlk, size_t index_pre, size_t index_post, size_t index_matrix, bool root_frequencies_included){
	double* spare_partials = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
	double* root_partials = (double*)aligned16_malloc( tlk->root_partials_size * sizeof(double) );
	tlk->calculate_per_cat_partials(tlk, spare_partials, index_pre, index_post, index_matrix );
	tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );

	if(!root_frequencies_included){
		// return tlk->node_log_likelihoods( tlk, root_partials, tlk->get_root_frequencies(tlk), pattern_lk);
	}
	double lk = 0;
	size_t v = 0;
	for ( size_t k = 0; k < tlk->pattern_count; k++ ) {
		double temp = 0;
		for ( size_t i = 0; i < tlk->m->nstate; i++ ) {
			temp += root_partials[v++];
		}
		lk += log(temp)* tlk->sp->weights[k];
	}
	free(spare_partials);
	free(root_partials);
	return lk;
}

#pragma mark - gradients

void gradient_cat_branch_lengths_aux(SingleTreeLikelihood* tlk, double* partials, int nodeId, double* cat_dlikelihoodsNode, const double* pattern_likelihoods){
	const int nstate   = tlk->m->nstate;
	const int patternCount = tlk->sp->count;
	size_t catCount = tlk->sm->cat_count;
	const double* freqs = tlk->get_root_frequencies(tlk);
	double* partialsPtr = partials;
	if(tlk->scale){
		double* spare_partials2 = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
		tlk->calculate_per_cat_partials(tlk, spare_partials2, tlk->upper_partial_indexes[nodeId], nodeId, nodeId );
		double* partialsPtr2 = spare_partials2;
		for ( size_t j = 0; j < catCount; j++ ){
			cat_dlikelihoodsNode[j] = 0;
			for ( size_t k = 0; k < patternCount; k++ ) {
				double temp = freqs[0] * *partialsPtr++;
				double temp2 = freqs[0] * *partialsPtr2++;
				for ( size_t i = 1; i < nstate; i++ ) {
					temp += freqs[i] * *partialsPtr++;
					temp2 += freqs[i] * *partialsPtr2++;
				}
				cat_dlikelihoodsNode[j] += temp/temp2 * tlk->sp->weights[k];
			}
		}
		free(spare_partials2);
	}
	else{
		for ( size_t j = 0; j < catCount; j++ ){
			cat_dlikelihoodsNode[j] = 0;
			for ( size_t k = 0; k < patternCount; k++ ) {
				double temp = freqs[0] * *partialsPtr++;
				for ( size_t i = 1; i < nstate; i++ ) {
					temp += freqs[i] * *partialsPtr++;
				}
				cat_dlikelihoodsNode[j] += temp/pattern_likelihoods[k] * tlk->sp->weights[k];
			}
		}
	}
}

void gradient_cat_branch_lengths_aux_freq_included(SingleTreeLikelihood* tlk, double* partials, int nodeId, double* cat_dlikelihoodsNode, const double* pattern_likelihoods){
	const int nstate   = tlk->m->nstate;
	const int patternCount = tlk->sp->count;
	size_t catCount = tlk->sm->cat_count;
	const double* freqs = tlk->get_root_frequencies(tlk);
	double* partialsPtr = partials;
	if(tlk->scale){
		double* spare_partials2 = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
		tlk->calculate_per_cat_partials(tlk, spare_partials2, tlk->upper_partial_indexes[nodeId], nodeId, nodeId );
		double* partialsPtr2 = spare_partials2;
		for ( size_t j = 0; j < catCount; j++ ){
			cat_dlikelihoodsNode[j] = 0;
			for ( size_t k = 0; k < patternCount; k++ ) {
				double temp = *partialsPtr++;
				double temp2 = *partialsPtr2++;
				for ( size_t i = 1; i < nstate; i++ ) {
					temp += *partialsPtr++;
					temp2 += *partialsPtr2++;
				}
				cat_dlikelihoodsNode[j] += temp/temp2 * tlk->sp->weights[k];
			}
		}
		free(spare_partials2);
	}
	else{
		for ( size_t j = 0; j < catCount; j++ ){
			cat_dlikelihoodsNode[j] = 0;
			for ( size_t k = 0; k < patternCount; k++ ) {
				double temp = *partialsPtr++;
				for ( size_t i = 1; i < nstate; i++ ) {
					temp += *partialsPtr++;
				}
				cat_dlikelihoodsNode[j] += temp/pattern_likelihoods[k] * tlk->sp->weights[k];
			}
		}
	}
}

// calculate derivative of log likelihood wrt to each branch length for each category
// dimension of branch_grandient: nodeCount x category_count
void gradient_cat_branch_lengths( SingleTreeLikelihood *tlk, double* branch_grandient, const double* pattern_likelihoods ){
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tlk->current_matrices_indexes[tempMatrixId]][tempMatrixId];
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	size_t nodeCount = Tree_node_count(tlk->tree);
	const int nstate   = tlk->m->nstate;
	const int patternCount = tlk->sp->count;
	const int catCount = tlk->sm->cat_count;
	const double* freqs = tlk->get_root_frequencies(tlk);

	for(size_t i = 0; i < nodeCount; i++){
		Node* node = Tree_node(tlk->tree, i);
		if(Node_isroot( node)) continue;
		double bl = Node_distance(node);
		const int nodeId = Node_id(node);
		if (tlk->bm != NULL) {
			bl = Node_time_elapsed(node) * tlk->bm->get(tlk->bm, node);
		}

		for (size_t c = 0; c < catCount; c++) {

#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			if( tlk->use_SIMD ){
				bool node_transpose = tlk->partials[0][Node_id(node)] == NULL;
				if( node_transpose ){
					tlk->m->dp_dt_transpose(tlk->m,
												bl * tlk->sm->get_rate(tlk->sm, c),
												&mat[c*tlk->matrix_size]);
				}
				else {
					tlk->m->dp_dt(tlk->m,
									  bl * tlk->sm->get_rate(tlk->sm, c),
									  &mat[c*tlk->matrix_size]);
				}
			}
			else{
				tlk->m->dp_dt(tlk->m,
								  bl * tlk->sm->get_rate(tlk->sm, c),
								  &mat[c*tlk->matrix_size]);
			}
#else
			tlk->m->dp_dt(tlk->m,
							  bl * tlk->sm->get_rate(tlk->sm, c),
							  &mat[c*tlk->matrix_size]);
#endif
		}

		tlk->calculate_per_cat_partials(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );

		double* partialsPtr = spare_partials;
		double* cat_dlikelihoodsNode = branch_grandient + nodeId*tlk->sm->cat_count;

		if(tlk->include_root_freqs){
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			if(tlk->m->nstate == 4){
				__m128d* p = (__m128d*)spare_partials;
				__m128d* p2 = (__m128d*)(spare_partials+2);
				double temp2[2] __attribute__ ((aligned (16)));
				if(tlk->scale){
					double* spare_partials2 = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
					tlk->calculate_per_cat_partials(tlk, spare_partials2, tlk->upper_partial_indexes[nodeId], nodeId, nodeId );
					__m128d* p1 = (__m128d*)spare_partials2;
					__m128d* p12 = (__m128d*)(spare_partials2+2);
					for ( size_t j = 0; j < catCount; j++ ){
						cat_dlikelihoodsNode[j] = 0;
						for ( size_t k = 0; k < patternCount; k++ ) {
							_mm_store_pd(temp2, _mm_add_pd(*p, *p2));
							p += 2;
							p2 += 2;
							double temp = temp2[0] + temp2[1];
							_mm_store_pd(temp2, _mm_add_pd(*p1, *p12));
							p1 += 2;
							p12 += 2;
							cat_dlikelihoodsNode[j] += temp/(temp2[0] + temp2[1]) * tlk->sp->weights[k];
						}
					}
					free(spare_partials2);
				}
				else{
					for ( size_t j = 0; j < catCount; j++ ){
						cat_dlikelihoodsNode[j] = 0;
						for ( size_t k = 0; k < patternCount; k++ ) {
							_mm_store_pd(temp2, _mm_add_pd(*p, *p2));
							p += 2;
							p2 += 2;
							double temp = temp2[0] + temp2[1];
							cat_dlikelihoodsNode[j] += temp/pattern_likelihoods[k] * tlk->sp->weights[k];
						}
					}
				}
			}
			else{
				gradient_cat_branch_lengths_aux_freq_included(tlk, spare_partials, nodeId, cat_dlikelihoodsNode, pattern_likelihoods);
			}
#else
			gradient_cat_branch_lengths_aux_freq_included(tlk, spare_partials, nodeId, cat_dlikelihoodsNode, pattern_likelihoods);
#endif
		}
		else{
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			if(tlk->m->nstate == 4){
				__m128d* p = (__m128d*)spare_partials;
				__m128d f1 = _mm_load_pd(freqs);
				__m128d f2 = _mm_load_pd(freqs+2);
				double temp[2] __attribute__ ((aligned (16)));
				__m128d temp2;
				if(tlk->scale){
					double* spare_partials2 = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
					tlk->calculate_per_cat_partials(tlk, spare_partials2, tlk->upper_partial_indexes[nodeId], nodeId, nodeId );
					__m128d* p2 = (__m128d*)spare_partials2;
					double temp3[2] __attribute__ ((aligned (16)));
					for ( size_t j = 0; j < catCount; j++ ){
						cat_dlikelihoodsNode[j] = 0;
						for ( size_t k = 0; k < patternCount; k++ ) {
							temp2 = _mm_mul_pd(f1, *p);
							p++;
							_mm_store_pd(temp, _mm_add_pd(temp2, _mm_mul_pd(f2, *p)));
							p++;
							temp2 = _mm_mul_pd(f1, *p2);
							p2++;
							_mm_store_pd(temp3, _mm_add_pd(temp2, _mm_mul_pd(f2, *p2)));
							p2++;
							cat_dlikelihoodsNode[j] += (temp[0]+temp[1])/(temp3[0]+temp3[1]) * tlk->sp->weights[k];
						}
					}	
					free(spare_partials2);
				}
				else{
					for ( size_t j = 0; j < catCount; j++ ){
						cat_dlikelihoodsNode[j] = 0;
						for ( size_t k = 0; k < patternCount; k++ ) {
							__m128d temp2 = _mm_mul_pd(f1, *p);
							p++;
							_mm_store_pd(temp, _mm_add_pd(temp2, _mm_mul_pd(f2, *p)));
							p++;
							cat_dlikelihoodsNode[j] += (temp[0]+temp[1])/pattern_likelihoods[k] * tlk->sp->weights[k];
						}
					}
				}
			}
			else{
				gradient_cat_branch_lengths_aux(tlk, spare_partials, nodeId, cat_dlikelihoodsNode, pattern_likelihoods);
			}
#else
			gradient_cat_branch_lengths_aux(tlk, spare_partials, nodeId, cat_dlikelihoodsNode, pattern_likelihoods);
#endif
		}
	}
}

void gradient_pinv_sitemodel(SingleTreeLikelihood* tlk, const double* branch_gradient, const double* branch_lengths, double* gradient){
	double* pattern_likelihoods = tlk->pattern_lk + tlk->sp->count;
	size_t nodeCount = Tree_node_count(tlk->tree);
	size_t stateCount = tlk->m->nstate;
	size_t patternCount = tlk->sp->count;
	size_t catCount = tlk->sm->cat_count;
	double pinv_grad = 0;
	assert(catCount == 2);
	for(size_t i = 0; i < nodeCount; i++){
		pinv_grad += branch_gradient[i*2 + 1] * branch_lengths[i];
	}
	
	double discrete_grad[2];
	discrete_grad[1] = pinv_grad;
	pinv_grad = 0;
	
	double* root_partials = tlk->partials[tlk->current_partials_indexes[Tree_root(tlk->tree)->id]][Tree_root(tlk->tree)->id];
	double* partials0 = root_partials;
	double* partials1 = root_partials + stateCount*patternCount;
	const double* freqs = tlk->get_root_frequencies(tlk);

	for ( size_t k = 0; k < patternCount; k++ ) {
		double L = 0;
		for ( size_t i = 0; i < stateCount; i++ ) {
			L += freqs[i] * (*partials0++ - *partials1++);
		}
		pinv_grad += L/pattern_likelihoods[k] * tlk->sp->weights[k];
	}

	discrete_grad[0] = pinv_grad;
	
	gradient[0] = tlk->sm->derivative(tlk->sm, discrete_grad, Parameters_at(tlk->sm->proportions->parameters, 0));
}

void gradient_pinv_W_sitemodel(SingleTreeLikelihood* tlk, const double* branch_gradient, const double* branch_lengths, double* gradient){
	assert(branch_gradient[Tree_root(tlk->tree)->id*tlk->sm->cat_count] == 0.0);
	double* pattern_likelihoods = tlk->pattern_lk + tlk->sp->count;
	size_t nodeCount = Tree_node_count(tlk->tree);
	size_t stateCount = tlk->m->nstate;
	size_t patternCount = tlk->sp->count;
	size_t catCount = tlk->sm->cat_count;
	size_t variableCount = catCount - 1;
	double pinv_grad = 0;
	double* discrete_grad = dvector(catCount);
	
	// j == 0 invariant site
	double* root_partials = tlk->partials[tlk->current_partials_indexes[Tree_root(tlk->tree)->id]][Tree_root(tlk->tree)->id];
	double* partials0 = root_partials;
	const double* freqs = tlk->get_root_frequencies(tlk);

	for ( size_t k = 0; k < patternCount; k++ ) {
		double Lk = 0;
		for ( size_t i = 0; i < stateCount; i++ ) {
			double temp = 0;
			for(size_t j = 1; j < catCount; j++){
				temp += root_partials[stateCount*patternCount*j + k*stateCount + i];
			}
			Lk += freqs[i] * (*partials0++ - temp/variableCount);
		}
		pinv_grad += Lk/pattern_likelihoods[k] * tlk->sp->weights[k];
	}
	discrete_grad[0] = pinv_grad;
	
	// j > 0
	for(size_t i = 0; i < nodeCount; i++){
		double branch_length = branch_lengths[i];
		for(size_t j = 1; j < catCount; j++){
			discrete_grad[j] += branch_gradient[i*catCount + j] * branch_length;
		}
	}
	gradient[0] = tlk->sm->derivative(tlk->sm, discrete_grad, Parameters_at(tlk->sm->proportions->parameters, 0));
	free(discrete_grad);
}

void gradient_shape_W_sitemodel(SingleTreeLikelihood* tlk, const double* branch_gradient, const double* branch_lengths, double* gradient){
	assert(branch_gradient[Tree_root(tlk->tree)->id*tlk->sm->cat_count] == 0.0);
	size_t catCount = tlk->sm->cat_count;
	size_t nodeCount = Tree_node_count(tlk->tree);
	double* discrete_grad = dvector(catCount);
	// Category 0 does not depend on shape
	for(size_t i = 0; i < nodeCount; i++){
		double branch_length = branch_lengths[i];
		size_t j = (tlk->sm->proportions == NULL ? 0 : 1);
		for( ; j < catCount; j++){
			discrete_grad[j] += branch_gradient[i*catCount + j] * branch_length;
		}
	}
	gradient[0] = tlk->sm->derivative(tlk->sm, discrete_grad, Parameters_at(tlk->sm->rates, 0));
	free(discrete_grad);
}

void gradient_discrete_sitemodel(SingleTreeLikelihood* tlk, const double* branch_gradient, const double* branch_lengths, double* gradient){
	size_t offset = 0;
	
	// derivative wrt shape of Weibull
	if(Parameters_count(tlk->sm->rates) == 1){
		gradient_shape_W_sitemodel(tlk, branch_gradient, branch_lengths, gradient);
		offset++;
	}
	
	// derivative wrt pinv
	if(tlk->sm->proportions != NULL){
		if(Parameters_count(tlk->sm->rates) == 0){
			gradient_pinv_sitemodel(tlk, branch_gradient, branch_lengths, gradient+offset);
		}
		else{
			gradient_pinv_W_sitemodel(tlk, branch_gradient, branch_lengths, gradient+offset);
		}
	}
}

void gradient_clock(SingleTreeLikelihood* tlk, const double* branch_gradient, double* gradient){
	size_t nodeCount = Tree_node_count(tlk->tree);
	if (Parameters_count(tlk->bm->rates) == 1) {
		gradient[0] = 0;
		for(size_t i = 0; i < nodeCount; i++){
			Node* node = Tree_node(tlk->tree, i);
			if(!Node_isroot(node)){
				gradient[0] += branch_gradient[node->id] * Node_time_elapsed(node);
			}
		}
	}
	else{
		size_t tipCount = Tree_tip_count(tlk->tree);
		for(size_t i = 0; i < nodeCount; i++){
			Node* node = Tree_node(tlk->tree, i);
			if(!Node_isroot(node)){
				size_t index = Node_isleaf(node) ? node->class_id : node->class_id + tipCount;
				gradient[index] = branch_gradient[node->id] * Node_time_elapsed(node);
			}
		}
	}
}

void gradient_PMatrix(SingleTreeLikelihood* tlk, const double* pattern_likelihoods, double* gradient){
	size_t parameter_count = tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K;
	if(tlk->m->simplex != NULL){
		parameter_count += tlk->m->simplex->K;
		if(tlk->m->grad_wrt_reparam) parameter_count--;
	}
	if(tlk->m->rates_simplex != NULL && tlk->m->grad_wrt_reparam){
		parameter_count--;
	}
	for(size_t i = 0; i < parameter_count; i++){
		gradient[i] = calculate_dlnl_dQ(tlk, i, pattern_likelihoods);
	}
}

void gradient_PMatrix_rates(SingleTreeLikelihood* tlk, const double* pattern_likelihoods, double* gradient){
	size_t parameter_count = tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K;
	if(tlk->m->rates_simplex != NULL && tlk->m->grad_wrt_reparam){
		parameter_count--;
	}
	for(size_t i = 0; i < parameter_count; i++){
		gradient[i] = calculate_dlnl_dQ(tlk, i, pattern_likelihoods);
	}
}

void gradient_PMatrix_frequencies(SingleTreeLikelihood* tlk, const double* pattern_likelihoods, double* gradient){
	size_t start = tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K;
	size_t parameter_count = tlk->m->simplex->K;
	if(tlk->m->grad_wrt_reparam) parameter_count--;

	if(tlk->m->rates_simplex != NULL && tlk->m->grad_wrt_reparam){
		start--;
	}
	for(size_t i = 0; i < parameter_count; i++){
		gradient[i] = calculate_dlnl_dQ(tlk, i+start, pattern_likelihoods);
	}
}

void gradient_branch_length_from_cat(SingleTreeLikelihood* tlk, const double* cat_branch_gradient, double* gradient){
	size_t nodeCount = Tree_node_count(tlk->tree);
	size_t catCount = tlk->sm->cat_count;
	Node** nodes = Tree_nodes(tlk->tree);
	double* weights = tlk->sm->get_proportions(tlk->sm);
	double* rates = tlk->sm->cat_rates;

	for(size_t i = 0; i < nodeCount; i++){
		size_t nodeId = nodes[i]->id;
		for (size_t j = 0; j < catCount; j++) {
			gradient[nodeId] += cat_branch_gradient[nodeId*catCount + j] * weights[j] * rates[j];
		}
	}
}

void gradient_branch_length_from_cat_inplace(SingleTreeLikelihood* tlk, double* cat_branch_gradient){
	size_t nodeCount = Tree_node_count(tlk->tree);
	size_t catCount = tlk->sm->cat_count;
	Node** nodes = Tree_nodes(tlk->tree);
	double* weights = tlk->sm->get_proportions(tlk->sm);
	double* rates = tlk->sm->cat_rates;

	for(size_t i = 0; i < nodeCount; i++){
		size_t nodeId = nodes[i]->id;
		cat_branch_gradient[nodeId] = cat_branch_gradient[nodeId*catCount] * weights[0] * rates[0];
		for (size_t j = 1; j < catCount; j++) {
			cat_branch_gradient[nodeId] += cat_branch_gradient[nodeId*catCount + j] * weights[j] * rates[j];
		}
	}
}

void gradient_heights(SingleTreeLikelihood* tlk, const double* branch_gradient, double* gradient){
	size_t nodeCount = Tree_node_count(tlk->tree);
	Node** nodes = Tree_get_nodes(tlk->tree, PREORDER);
	for(size_t i = 1; i < nodeCount; i++){
		Node* node = nodes[i];
		double nodeGradient = branch_gradient[node->id] * tlk->bm->get(tlk->bm, node);
		if(!Node_isleaf(node)){
			gradient[node->class_id] = -nodeGradient;
		}
		gradient[node->parent->class_id] += nodeGradient;
	}
}


// \partial log(L)/\partial r_i &= \partial log(L)/\partial h_i \partial h_i\partial r_i
// &= \sum_j \partial log(L)/\partial b_i  \partial b_i/\partial h_i \partial h_i\partial r_i
void gradient_ratios(SingleTreeLikelihood* tlk, const double* branch_gradient, double* gradient){
	double* height_gradient = dvector(Tree_tip_count(tlk->tree) - 1);
	gradient_heights(tlk, branch_gradient, height_gradient);
	//jvp and jacobian
	Tree_node_transform_jvp(tlk->tree, height_gradient, gradient);
	if(tlk->include_jacobian){
		Tree_node_transform_jacobian_gradient(tlk->tree, gradient);
	}
//	node_transform_jvp(tlk->tree, height_gradient, gradient, true);
	free(height_gradient);
}

// derivatives are with respect to the constrained values
void central_finite_differences_simplex(Model* model, Simplex* simplex, double epsilon, double* gradient){
	const double* const_freqs = simplex->get_values(simplex);
	double* freqs = clone_dvector(const_freqs, simplex->K);
	for(size_t i = 0; i < simplex->K; i++){
		double v = freqs[i];

		freqs[i] = v + epsilon;
		simplex->set_values(simplex, freqs);
		double pp = model->logP(model);
		
		freqs[i] = v - epsilon;
		simplex->set_values(simplex, freqs);
		double mm = model->logP(model);
		
		freqs[i] = v;
		simplex->set_values(simplex, freqs);
		
		gradient[i] = (pp - mm)/(2.0*epsilon);
	}
	free(freqs);
}

void central_finite_differences_parameters(Model* model, Parameters* parameters, double epsilon, double* gradient){
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		gradient[i] = Model_first_derivative(model, Parameters_at(parameters, i), epsilon);
	}
}

// compute gradient wrt every parameter
// derivative are concatenated in this order:
// (branch length or ratios), site model, clock model, susbtitution model
void TreeLikelihood_calculate_gradient( Model *model, double* grads ){
	SingleTreeLikelihood* tlk = model->obj;
	double* pattern_likelihoods = tlk->pattern_lk + tlk->sp->count;
	for (int i = 0; i < tlk->sp->count; i++) {
		pattern_likelihoods[i] = exp(tlk->pattern_lk[i]);
	}
	
	bool prepare_tree = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_TREE_MODEL;
	bool prepare_site_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SITE_MODEL;
	bool prepare_branch_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_BRANCH_MODEL;
	bool prepare_substitution_model = tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_UNCONSTRAINED || tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL;
	bool prepare_substitution_model_rates = prepare_substitution_model || (tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_RATES);
	bool prepare_substitution_model_frequencies = prepare_substitution_model || (tlk->prepared_gradient & TREELIKELIHOOD_FLAG_SUBSTITUTION_MODEL_FREQUENCIES);
	
	size_t nodeCount = Tree_node_count(tlk->tree);
	size_t catCount = tlk->sm->cat_count;
	double* cat_branch_gradient = dvector(nodeCount*catCount);
	double* branch_lengths = dvector(nodeCount);
	// calculate per category gradient and gradient of branch lengths
	Node** nodes = Tree_nodes(tlk->tree);
	
	bool time_mode = Tree_is_time_mode(tlk->tree);
	if (time_mode) {
		for(size_t i = 0; i < nodeCount; i++){
			branch_lengths[nodes[i]->id] = Node_time_elapsed(nodes[i]) * tlk->bm->get(tlk->bm, nodes[i]);
		}
	}
	else{
		for(size_t i = 0; i < nodeCount; i++){
			branch_lengths[nodes[i]->id] = Node_distance(nodes[i]);
		}
	}
	if(tlk->sm->mu != NULL){
		double mu = Parameter_value(tlk->sm->mu);
		for(size_t i = 0; i < nodeCount; i++){
			branch_lengths[i] *= mu;
		}
	}
	size_t offset = 0;

	// calculate per category gradient and gradient of branch lengths
	if (prepare_branch_model || prepare_tree || prepare_site_model) {
		gradient_cat_branch_lengths(tlk, cat_branch_gradient, pattern_likelihoods);
	}
	if (!time_mode) {
		size_t rightNodeID = Tree_root(tlk->tree)->right->id;
		for (size_t j = 0; j < catCount; j++) {
			cat_branch_gradient[rightNodeID*catCount + j] = 0;
		}
		branch_lengths[rightNodeID] = 0;
	}
	
	double grad_sitemodel[2];
	if(catCount > 1){
		if (prepare_site_model) {
			gradient_discrete_sitemodel(tlk, cat_branch_gradient, branch_lengths, grad_sitemodel);
		}

		if (prepare_branch_model || prepare_tree) {
			gradient_branch_length_from_cat_inplace(tlk, cat_branch_gradient);
		}
	}
	
	if (prepare_tree && time_mode) {
		gradient_ratios(tlk, cat_branch_gradient, grads);
		offset += Tree_tip_count(tlk->tree) - 1;
	}
	else if (prepare_tree){
		memcpy(grads, cat_branch_gradient, sizeof(double)*nodeCount);
		offset += nodeCount;
	}
	
	if (prepare_site_model) {
		size_t local_offset = 0;
		if(Parameters_count(tlk->sm->rates) == 1){
			grads[offset++] = grad_sitemodel[local_offset++];
		}
		
		// derivative wrt pinv
		if(tlk->sm->proportions != NULL){
			grads[offset++] = grad_sitemodel[local_offset];
		}
		
		// derivative wrt mu
		if(tlk->sm->mu != NULL){
			double grad = 0;
			if (time_mode) {
				for(size_t i = 0; i < nodeCount; i++){
					grad += cat_branch_gradient[nodes[i]->id] * Node_time_elapsed(nodes[i]) * tlk->bm->get(tlk->bm, nodes[i]);
				}
			}
			else{
				for(size_t i = 0; i < nodeCount; i++){
					grad += cat_branch_gradient[nodes[i]->id] * Node_distance(nodes[i]);
				}
			}
			grads[offset++] = grad;
		}
	}
	

	if (prepare_branch_model && time_mode) {
		gradient_clock(tlk, cat_branch_gradient, grads + offset);
		offset += Parameters_count(tlk->bm->rates);
	}

	if(prepare_substitution_model){
		gradient_PMatrix(tlk, pattern_likelihoods, grads+offset);
	}
	else{
		Model** models = (Model**)model->data;
		double epsilon = models[1]->epsilon;

		if(prepare_substitution_model_rates){
			//TODO: frequency parameters
			if(tlk->m->dPdp == NULL || tlk->m->modeltype == NONREVERSIBLE){
				Parameters* params = tlk->m->rates_simplex == NULL ? tlk->m->rates : tlk->m->rates_simplex->parameters;
				for(size_t i = 0; i < Parameters_count(params); i++){
					grads[offset+i] = Model_first_derivative(model, Parameters_at(params, i), 0.000001);
				}
			}
			else{
				if(epsilon > 0.0){
					if(!tlk->m->grad_wrt_reparam && tlk->m->rates_simplex != NULL){
						central_finite_differences_simplex(model, tlk->m->rates_simplex, epsilon, grads+offset);
					}
					else{
						Parameters* params = tlk->m->rates_simplex == NULL ? tlk->m->rates : tlk->m->rates_simplex->parameters;
						central_finite_differences_parameters(model, params, epsilon, grads+offset);
					}
				}
				else{
					gradient_PMatrix_rates(tlk, pattern_likelihoods, grads+offset);
				}
			}
			offset += tlk->m->rates_simplex == NULL ? Parameters_count(tlk->m->rates) : tlk->m->rates_simplex->K;
			if(tlk->m->grad_wrt_reparam) offset--;
		}
		if(prepare_substitution_model_frequencies){
			if(epsilon > 0.0){
				if(!tlk->m->grad_wrt_reparam){
					central_finite_differences_simplex(model, tlk->m->simplex, epsilon, grads+offset);
				}
				else{
					Parameters* params = tlk->m->simplex->parameters;
					central_finite_differences_parameters(model, params, epsilon, grads+offset);
				}
			}
			else{
				gradient_PMatrix_frequencies(tlk, pattern_likelihoods, grads+offset);
			}
		}
	}

	free(branch_lengths);
	free(cat_branch_gradient);
}
