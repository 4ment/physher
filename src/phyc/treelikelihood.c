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
#include <unistd.h> // for sleep

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

#ifdef USE_SSE
#include <xmmintrin.h>
#endif

#ifdef AVX_ENABLED
#include <immintrin.h>
#endif

#define OPTIMIZATION_PRECISION 0.01


// MARK:  Private function declaration


static bool _calculate_partials( SingleTreeLikelihood *tlk, Node *n  );
static double _calculate( SingleTreeLikelihood *tlk );

double _calculate_uppper( SingleTreeLikelihood *tlk, Node *node );

void calculate_dldp_uppper_reparam( SingleTreeLikelihood *tlk, Node *node, double* pattern_dlikelihoods );
double calculate_dl_dr(SingleTreeLikelihood* tlk, double* pattern_likelihoods);
double calculate_dl_dr_discrete(SingleTreeLikelihood* tlk, Node* node, double* pattern_likelihoods);

#pragma mark -
#pragma mark TreeLikelihoodModel

void _treelikelihood_handle_change( Model *self, Model *model, int index ){
	SingleTreeLikelihood *tlk = (SingleTreeLikelihood*)self->obj;

	//printf("%s %d\n", model->name, index);
	if ( model->type == MODEL_TREE ) {
//		printf("node index %d\n", index);
		//SingleTreeLikelihood_update_one_node(tlk, index);
		tlk->update_nodes[index] = true;
		tlk->update = true;
		tlk->update_upper = true;
	}
	else if ( model->type == MODEL_BRANCHMODEL ) {
		if(index == -1){
			SingleTreeLikelihood_update_all_nodes(tlk);
		}
		else{
			tlk->update_nodes[index] = true;
			tlk->update = true;
			tlk->update_upper = true;
		}
	}
	else if ( model->type == MODEL_SITEMODEL ) {
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
	memcpy(tlk->stored_matrices_indexes, tlk->current_matrices_indexes, sizeof(unsigned)*Tree_node_count(tlk->tree)*2);
	memcpy(tlk->stored_partials_indexes, tlk->current_partials_indexes, sizeof(unsigned)*Tree_node_count(tlk->tree)*2);
	Model** models = (Model**)self->data;
	models[0]->store(models[0]); // tree
	models[1]->store(models[1]); // sitemodel
	if(models[2] != NULL){
		models[2]->store(models[2]); // branchmodel
	}
	self->storedLogP = self->lp;
	tlk->stored_lk = tlk->lk;
}

static void _singleTreeLikelihood_restore(Model* self){
	Model** models = (Model**)self->data;
	models[0]->restore(models[0]); // tree
	models[1]->restore(models[1]); // sitemodel
	if(models[2] != NULL){
		models[2]->restore(models[2]); // branchmodel
	}
	self->lp = self->storedLogP;
}

double _singleTreeLikelihood_logP(Model *self){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	self->lp = tlk->calculate(tlk);
	return self->lp;
}

double _singleTreeLikelihood_full_logP(Model *self){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	SingleTreeLikelihood_update_all_nodes(tlk);
	self->lp = tlk->calculate(tlk);
	return self->lp;
}

double calculate_log_jacobian(Node* node, double* logP){
	if(!Node_isleaf(node)){
		double l = calculate_log_jacobian(node->left, logP);
		double r = calculate_log_jacobian(node->right, logP);
		double max = dmax(l, r);
		if(!Node_isroot(node)){
			*logP += log(Node_height(Node_parent(node)) - max);
		}
		return max;
	}
	else{
		return Node_height(node);
	}
}

double _singleTreeLikelihood_reparametrized_logP(Model *self){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	self->lp = tlk->calculate(tlk);
	if(Tree_homochronous(tlk->tree)){
		for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
			Node* node = Tree_node(tlk->tree, i);
			if(!Node_isroot(node) && !Node_isleaf(node)){
				self->lp += log(Node_height(Node_parent(node)));
			}
		}
	}
	else{
		calculate_log_jacobian(Tree_root(tlk->tree), &self->lp);
	}
	return self->lp;
}

double _singleTreeLikelihood_reparametrized_full_logP(Model *self){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
		SingleTreeLikelihood_update_all_nodes(tlk);
	self->lp = tlk->calculate(tlk);
	if(Tree_homochronous(tlk->tree)){
		for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
			Node* node = Tree_node(tlk->tree, i);
			if(!Node_isroot(node) && !Node_isleaf(node)){
				self->lp += log(Node_height(Node_parent(node)));
			}
		}
	}
	else{
		calculate_log_jacobian(Tree_root(tlk->tree), &self->lp);
	}
	return self->lp;
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

double _singleTreeLikelihood_dlogP(Model *self, const Parameter* p){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	
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
		//SingleTreeLikelihood_update_all_nodes(tlk);
		
		double logP = tlk->calculate(tlk); // make sure it is updated
		if (isnan(logP) || isinf(logP)) {
			return logP;
		}
		update_upper_partials(tlk, Tree_root(tlk->tree));
		
		for (int i = 0; i < tlk->sp->count; i++) {
			pattern_likelihoods[i] = exp(tlk->pattern_lk[i]);
		}
		tlk->update_upper = false;
	}
	
	double dlogP = 0;
	
	if (p->model == MODEL_TREE) {
		double* pattern_dlikelihoods = tlk->pattern_lk + 2*tlk->sp->count;
		Node* node = NULL;
        double adjustement = 0;
		// reparameterized heights
		if (p->id < 0) {
			int realIndex = -p->id-1;
			node = Tree_node(tlk->tree, realIndex);
			calculate_dldp_uppper_reparam(tlk, node, pattern_dlikelihoods);
            double* descendant = dvector(Tree_node_count(tlk->tree));
            unsigned *map = get_reparam_map(tlk->tree);
            Parameters* reparams = get_reparams(tlk->tree);
            double* lowers = Tree_lowers(tlk->tree);
            _calculate_dlog_jacobian(node, node, &adjustement, descendant, map, reparams, lowers);
            free(lowers);
            free(descendant);
		}
		// distances
		else{
			node = Tree_node(tlk->tree, p->id);
			calculate_dldt_uppper(tlk, node, pattern_dlikelihoods);
		}
		dlogP = dlnldt_uppper(tlk, pattern_likelihoods, pattern_dlikelihoods) + adjustement;
//		printf("dlogP: %f finite diff: %f (%d %s) %f == %d -- %e\n", dlogP, Model_first_derivative(self, p, 0.00001), node->id, p->name, p->value, node->depth, Node_height(node));
	}
	else if(p->model == MODEL_SUBSTITUTION){
		// Substitution model parameter
		int idx = 0;
		for (; idx < Parameters_count(tlk->sm->m->rates); idx++) {
			if(Parameters_at(tlk->sm->m->rates, idx) == p) break;
		}
		if(idx >= Parameters_count(tlk->sm->m->rates)){
			for (int ii = 0; ii < Parameters_count(tlk->sm->m->simplex->parameters); idx++, ii++) {
				if(Parameters_at(tlk->sm->m->simplex->parameters, ii) == p) break;
			}
		}
		
		if(idx < Parameters_count(tlk->sm->m->rates)+Parameters_count(tlk->sm->m->simplex->parameters)){
			dlogP = calculate_dlnl_dQ( tlk, idx, pattern_likelihoods );
//			printf("dlogP: %f finite diff: %f (%d %s)\n", dlogP, Model_first_derivative(self, p, 0.0001), idx, p->name);
		}
	}
	else if(p->model == MODEL_SITEMODEL){
		// Rate heterogeneity (Weibull only)
		if (Parameters_count(tlk->sm->rates) != 0 && p == Parameters_at(tlk->sm->rates, 0)) {
			dlogP = calculate_dlnl_dWeibull( tlk, pattern_likelihoods );
//			printf("dlogP: %f finite diff: %f (%d %s)\n", dlogP2, Model_first_derivative(self, p, 0.001), idx, p->name);
		}
	}
	else if(p->model == MODEL_BRANCHMODEL){
        if (tlk->bm->name == CLOCK_STRICT) {
            dlogP = calculate_dl_dr(tlk, pattern_likelihoods);
        }
        else if(tlk->bm->name == CLOCK_DISCRETE){
            Node* node = Tree_node(tlk->tree, p->id);
            dlogP = calculate_dl_dr_discrete(tlk, node, pattern_likelihoods);
        }
        else{
            fprintf(stderr, "derivative for clock not yet implemented\n");
            exit(2);
        }
//		printf("dlogP: %f finite diff: %f (%s)\n", dlogP, Model_first_derivative(self, p, 0.00001), p->name);
	}
	
	if(isnan(dlogP)){
		SingleTreeLikelihood_update_all_nodes(tlk);
		//		Tree_print_parameters(tlk->tree);
		//		printf("%f\n", tlk->calculate(tlk));
		//		printf("%f\n", Model_first_derivative(self, p, 0.001));
		//		printf("%s %f\n", p->name, p->value);
		//		exit(10);
	}
	
	// derivative wrt root height of the log determinant of the jacobian
//	double count = Tree_node_count(tlk->tree) - Tree_tip_count(tlk->tree) - 1;
//	dlogP += (double)count/Node_height(Tree_root(tlk->tree));
	
	return dlogP;
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
		update_upper_partials(tlk, Tree_root(tlk->tree));
		
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
		update_upper_partials(tlk, Tree_root(tlk->tree));
		
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
    
    for (int i = 0; i < tlk->sm->cat_count; i++) {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
        if( tlk->use_SIMD ){
            if( Node_isleaf(node1) ){
                tlk->sm->m->dp_dt_transpose(tlk->sm->m,
                                            bl1 * tlk->sm->get_rate(tlk->sm, i),
                                            &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
            }
            else {
                tlk->sm->m->dp_dt(tlk->sm->m,
                                  bl1 * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
            }
            if( Node_isleaf(node2) ){
                tlk->sm->m->dp_dt_transpose(tlk->sm->m,
                                            bl2 * tlk->sm->get_rate(tlk->sm, i),
                                            &tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size]);
            }
            else {
                tlk->sm->m->dp_dt(tlk->sm->m,
                                  bl2 * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size]);
            }
        }
        else{
            tlk->sm->m->dp_dt(tlk->sm->m,
                              bl1 * tlk->sm->get_rate(tlk->sm, i),
                              &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
            
            tlk->sm->m->dp_dt(tlk->sm->m,
                              bl2 * tlk->sm->get_rate(tlk->sm, i),
                              &tlk->matrices[tlk->current_matrices_indexes[Node_id(node2)]][Node_id(node2)][i*tlk->matrix_size]);
        }
#else
        tlk->sm->m->dp_dt(tlk->sm->m,
                          bl1 * tlk->sm->get_rate(tlk->sm, i),
                          &tlk->matrices[tlk->current_matrices_indexes[Node_id(node1)]][Node_id(node1)][i*tlk->matrix_size]);
        
        tlk->sm->m->dp_dt(tlk->sm->m,
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
    
    const int nstate   = tlk->sm->nstate;
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
	if(self->ref_count == 1){
		//printf("Free treelikelihood model %s\n", self->name);
		SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
		int count = (tlk->bm==NULL?2:3);
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
	Model* msm = list[1];
	Model* mbm = list[2];
	Model *mtreeclone = NULL;
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
		clonetlk = clone_SingleTreeLikelihood_with(tlk, (Tree*)mtreeclone->obj, (SiteModel*)msmclone->obj, tlk->sp, (BranchModel*)mbmclone->obj);
	}
	else{
		clonetlk = clone_SingleTreeLikelihood_with(tlk, (Tree*)mtreeclone->obj, (SiteModel*)msmclone->obj, tlk->sp, NULL);
	}
	Model* clone = new_TreeLikelihoodModel(self->name, clonetlk, mtreeclone, msmclone, mbmclone);
	Hashtable_add(hash, clone->name, clone);
	mtreeclone->free(mtreeclone);
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
Model * new_TreeLikelihoodModel( const char* name, SingleTreeLikelihood *tlk,  Model *tree, Model *sm, Model *bm ){
	Model *model = new_Model(MODEL_TREELIKELIHOOD,name, tlk);

	tree->listeners->add( tree->listeners, model );
	if(bm != NULL)bm->listeners->add( bm->listeners, model );
	sm->listeners->add( sm->listeners, model );
	model->handle_restore = _treelikelihood_handle_restore;
	
	model->logP = _singleTreeLikelihood_logP;
	model->full_logP = _singleTreeLikelihood_full_logP;
	model->dlogP = _singleTreeLikelihood_dlogP;
	model->d2logP = _singleTreeLikelihood_d2logP;
	model->ddlogP = _singleTreeLikelihood_ddlogP;
	model->update = _treelikelihood_handle_change;
	model->free = _treeLikelihood_model_free;
	model->clone = _treeLikelihood_model_clone;
	model->data = (Model**)malloc(sizeof(Model*)*3);
	model->store = _singleTreeLikelihood_store;
	model->restore = _singleTreeLikelihood_restore;
	Model** list = (Model**)model->data;
	list[0] = tree;
	list[1] = sm;
	list[2] = bm;
	tree->ref_count++;
	sm->ref_count++;
	if(bm != NULL) bm->ref_count++;
	return model;
}

Model * new_TreeLikelihoodModel_from_json(json_node*node, Hashtable*hash){
	char* allowed[] = {
		"branchmodel",
		"reparameterized",
		"root_frequencies",
		"sitemodel",
		"sitepattern",
		"sse",
		"tree"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* patterns_node = get_json_node(node, "sitepattern");
	json_node* tree_node = get_json_node(node, "tree");
	json_node* sm_node = get_json_node(node, "sitemodel");
	json_node* bm_node = get_json_node(node, "branchmodel");
	bool reparameterized = get_json_node_value_bool(node, "reparameterized", false);
	
	Model* mtree = NULL;
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
	
	SingleTreeLikelihood* tlk = new_SingleTreeLikelihood((Tree*)mtree->obj, (SiteModel*)msm->obj, patterns, bm);
	char* id = get_json_node_value_string(node, "id");
	Model* model = new_TreeLikelihoodModel(id, tlk, mtree, msm, mbm);
	mtree->free(mtree);
	msm->free(msm);
	if(mbm != NULL) mbm->free(mbm);
	
	bool root_freqs_unknown = get_json_node_value_bool(node, "root_frequencies", false);
	if(root_freqs_unknown){
		tlk->get_root_frequencies = get_root_frequencies_fixed;
		tlk->root_frequencies = dvector(tlk->sm->nstate);
		for (int i = 0; i < tlk->sm->nstate; i++) {
			tlk->root_frequencies[i] = 1;
		}
	}
	
	bool use_sse = get_json_node_value_bool(node, "sse", true);
	if (!use_sse) {
		SingleTreeLikelihood_enable_SSE(tlk, false);
	}
	
	if(reparameterized){
		model->logP = _singleTreeLikelihood_reparametrized_logP;
		model->full_logP = _singleTreeLikelihood_reparametrized_full_logP;
	}
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
            if( Node_isleaf( nodes[i] ) ){
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
            if( Node_isleaf( nodes[i] ) ){
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

SingleTreeLikelihood * new_SingleTreeLikelihood( Tree *tree, SiteModel *sm, SitePattern *sp, BranchModel *bm ){
	SingleTreeLikelihood *tlk = (SingleTreeLikelihood *)malloc( sizeof(SingleTreeLikelihood));
	assert(tlk);
	
	tlk->tree = tree;
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
	
	Node **nodes = Tree_get_nodes( tree, POSTORDER );
	
	tlk->partials_size = sp->count * sm->nstate * tlk->cat_count;
	
	
	tlk->matrix_dim = Tree_node_count(tree);
	tlk->matrix_size = sm->nstate * sm->nstate;
	
	tlk->upper_partial_indexes = ivector(Tree_node_count(tree));
	
	tlk->root_partials_size = sp->count*sm->nstate*3;
	tlk->pattern_lk_size = sp->count*4;
	
	tlk->partials_dim = Tree_node_count(tree)*2; // allocate *2 for upper likelihoods
	
    // odd number of state
//    if( sm->nstate & 1 ){
//        tlk->partials_size += sp->count * sm->cat_count;
//        tlk->matrix_size   += sm->nstate;
//    }
	
    allocate_storage(tlk, 0);
    allocate_storage(tlk, 1); // used by mcmc with store/restore
	
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
	tlk->calculate_branch_likelihood = calculate_branch_partials;
	tlk->update_partials_flexible = NULL;
	
	if ( sm->nstate == 4 ) {
		tlk->update_partials      = update_partials_4;
		tlk->update_partials_flexible = update_partials_flexible_4;
		tlk->integrate_partials   = integrate_partials_4;
		tlk->node_log_likelihoods = node_log_likelihoods_4;
		tlk->calculate_branch_likelihood = calculate_branch_partials_4;
	}
	else if( sm->nstate == 20 ){
		tlk->update_partials      = update_partials_general;
		tlk->integrate_partials   = integrate_partials_general;
		tlk->node_log_likelihoods = node_log_likelihoods_general;
		tlk->calculate_branch_likelihood = calculate_branch_partials_20;
	}
	else if( sm->nstate >= 60 ){
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
	if ( sm->nstate == 4 ) {
		tlk->update_partials      = update_partials_4_SSE;
		tlk->update_partials_flexible = update_partials_flexible_4_SSE;
		tlk->integrate_partials   = integrate_partials_4_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_4_SSE;
		tlk->calculate_branch_likelihood = calculate_branch_partials_4_SSE;
        tlk->use_SIMD = true;
	}
	else if ( sm->nstate == 20 ) {
		tlk->update_partials      = update_partials_20_SSE;
		tlk->integrate_partials   = integrate_partials_general;
		tlk->node_log_likelihoods = node_log_likelihoods_general;
		tlk->calculate_branch_likelihood = calculate_branch_partials_20_SSE;
        tlk->use_SIMD = true;
	}
    // Does not work for matrix with odd dimension.
    // If we want to load doubles at indexes 1 and 2 from matrices it will crash
    // same issue with partials, it will try to store at odd indexes
	else if( sm->nstate >= 60 && !(sm->nstate & 1) ){
		tlk->update_partials      = update_partials_codon_SSE;
		tlk->integrate_partials   = integrate_partials_codon_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_codon_SSE;
        tlk->use_SIMD = true;
	}
    else if(!(sm->nstate & 1)){
        tlk->update_partials      = update_partials_general_even_SSE;
        tlk->integrate_partials   = integrate_partials_general_even_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_general_even_SSE;
		tlk->calculate_branch_likelihood = calculate_branch_partials_even_SSE;
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
	
	free(tlk);
}

void free_SingleTreeLikelihood( SingleTreeLikelihood *tlk ){
	if ( tlk->bm != NULL ) tlk->bm->free( tlk->bm, false );
	if ( tlk->sm != NULL ) free_SiteModel( tlk->sm );
	if ( tlk->sp != NULL ) free_SitePattern( tlk->sp );
	free_SingleTreeLikelihood_internals(tlk);
}

SingleTreeLikelihood * clone_SingleTreeLikelihood( SingleTreeLikelihood *tlk ){
	SitePattern *sp = clone_SitePattern(tlk->sp);
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
	return clone_SingleTreeLikelihood_with( tlk, tree, sm, sp, bm  );
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


SingleTreeLikelihood * clone_SingleTreeLikelihood_with( SingleTreeLikelihood *tlk, Tree *tree, SiteModel *sm, SitePattern *sp, BranchModel *bm){
	SingleTreeLikelihood *newtlk = (SingleTreeLikelihood *)malloc( sizeof(SingleTreeLikelihood));
	assert(newtlk);
	
	newtlk->tree = tree;
	newtlk->sm   = sm;
	newtlk->sp   = sp;
    newtlk->bm = bm;
	
	newtlk->cat_count = tlk->cat_count;
	newtlk->pattern_count = tlk->pattern_count;
	
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
	newtlk->calculate_branch_likelihood = tlk->calculate_branch_likelihood;
	
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
		if( !Node_isleaf(nodes[i]) ){
			memcpy(newtlk->partials[0][Node_id( nodes[i] )], tlk->partials[0][Node_id( nodes[i] )], tlk->partials_size * sizeof(double));
		}
	}
	for ( ; i < tlk->partials_dim; i++ ) {
		memcpy(newtlk->partials[0][i], tlk->partials[0][i], tlk->partials_size * sizeof(double));
	}

	for ( int i = 0; i < tlk->matrix_dim; i++ ) {
		memcpy(newtlk->matrices[0][i], tlk->matrices[0][i], tlk->matrix_size*tlk->sm->cat_count * sizeof(double));
	}
    
//    if(tlk->partials[1] != NULL){
        allocate_storage(newtlk, 1);
        memcpy(newtlk->stored_matrices_indexes, tlk->stored_matrices_indexes, 2*nodeCount*sizeof(unsigned));
        memcpy(newtlk->stored_partials_indexes, tlk->stored_partials_indexes, 2*nodeCount*sizeof(unsigned));
        i = 0;
        for ( ; i < nodeCount; i++ ) {
            if( !Node_isleaf(nodes[i]) ){
                memcpy(newtlk->partials[1][Node_id( nodes[i] )], tlk->partials[1][Node_id( nodes[i] )], tlk->partials_size * sizeof(double));
            }
        }
        for ( ; i < tlk->partials_dim; i++ ) {
            memcpy(newtlk->partials[1][i], tlk->partials[1][i], tlk->partials_size * sizeof(double));
        }

        for ( int i = 0; i < tlk->matrix_dim; i++ ) {
            memcpy(newtlk->matrices[1][i], tlk->matrices[1][i], tlk->matrix_size*tlk->sm->cat_count * sizeof(double));
        }
//    }
	
	newtlk->root_partials = (double*)malloc( newtlk->root_partials_size*sizeof(double) );
    assert(newtlk->root_partials);
	memcpy(newtlk->root_partials, tlk->root_partials, newtlk->root_partials_size * sizeof(double));
    
	newtlk->sp->ref_count++;
	
	newtlk->root_frequencies = NULL;
	if(tlk->root_frequencies != NULL){
		newtlk->root_frequencies = clone_dvector(tlk->root_frequencies, tlk->sm->nstate);
	}
	newtlk->get_root_frequencies = tlk->get_root_frequencies;

	newtlk->reparametrized = tlk->reparametrized;
	
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
        tlk->scaling_factors[1] = NULL;
        if(tlk->partials[1] != NULL){
            tlk->scaling_factors[1] = dmatrix(tlk->partials_dim, tlk->sp->count );
        }
	}
}

int SingleTreeLikelihood_df_count( const SingleTreeLikelihood *stlk ){
	int df = 0;
	Node **nodes = Tree_get_nodes(stlk->tree, POSTORDER);
	if( stlk->bm == NULL ){
        // root has no distance
        // if the root is not rooted there should be one df less
		for (int i = 0; i < Tree_node_count(stlk->tree)-1; i++) {
			if ( Parameter_estimate( nodes[i]->distance ) ) {
				df++;
			}
		}
	}
	else {
		df = Parameters_count( stlk->bm->rates );
		for (int i = 0; i < Tree_node_count(stlk->tree); i++) {
			if ( Parameter_estimate( nodes[i]->height ) ) {
				df++;
			}
		}
	}

	// Frequencies
	if(Parameters_estimate(stlk->sm->m->simplex->parameters, 0)){
		df += stlk->sm->nstate-1;
	}
	// Rate bias
	if(Parameters_count(stlk->sm->m->rates) > 0){
		for (int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
			df += Parameters_estimate(stlk->sm->m->rates, i);
		}
	}
	
	// SiteModel
	if(Parameters_count(stlk->sm->rates) > 0){
		for (int i = 0; i < Parameters_count(stlk->sm->rates); i++) {
			df += Parameters_estimate(stlk->sm->rates, i);
		}
	}
	return df;
}


// Rearrange the partials
void SingleTreeLikelihood_rearrange_partials( SingleTreeLikelihood *tlk ){
    
    // assign a partial vector to internal nodes without partials, take it from a tip
    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
	for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
		if ( !Node_isleaf(nodes[i]) && tlk->partials[0][i] == NULL ) {
			for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
				if( Node_isleaf(nodes[j]) && tlk->partials[0][j] != NULL ){
					tlk->partials[i] = tlk->partials[j];
					tlk->partials[j] = NULL;
					break;
				}
			}
		}
	}
    SingleTreeLikelihood_update_all_nodes(tlk);
}

// Swap nodes in Tree and rearrange partials and matrices in SingleTreeLikelihood
// Update the postorder array of nodes in Tree
// Update the taxa/pattern mapping
void SingleTreeLikelihood_rearrange( SingleTreeLikelihood *tlk, Node *node1, Node *node2 ){
    
    Tree_rearrange( tlk->tree, node1, node2 );
    
    SingleTreeLikelihood_rearrange_partials(tlk);
}

void SingleTreeLikelihood_copy_partials( SingleTreeLikelihood *src, SingleTreeLikelihood *dst ){
    Node **nodes_src = Tree_get_nodes(src->tree, POSTORDER);
    Node **nodes_dst = Tree_get_nodes(dst->tree, POSTORDER);
    int nNodes = Tree_node_count(src->tree);
    
    for ( int i = 0; i < nNodes; i++ ) {
        if( !Node_isleaf(nodes_dst[i]) ){
            memcpy(dst->partials[nodes_dst[i]->postorder_idx], src->partials[nodes_src[i]->postorder_idx], src->partials_size*sizeof(double));
        }
        if( !Node_isroot(nodes_src[i]) ){
            memcpy(dst->matrices[nodes_dst[i]->postorder_idx], src->matrices[nodes_src[i]->postorder_idx], src->matrix_size*src->sm->cat_count*sizeof(double));
        }
    }
    
    memcpy(dst->root_partials, src->root_partials, src->sp->count*src->sm->nstate * sizeof(double));
    dst->lk = src->lk;
    memcpy(dst->update_nodes,dst->update_nodes, nNodes*sizeof(bool) );
    dst->update = src->update;
    error("should not be here3\n");
}

// Should not be used if the root is fixed
void SingleTreeLikelihood_add_height( SingleTreeLikelihood *tlk, Node *node, double value ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( Node_isroot(node) || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
		SingleTreeLikelihood_add_height(tlk, node->left, value);
		SingleTreeLikelihood_add_height(tlk, node->right, value);
		Node_set_height(node, Node_height(node)+value);
		SingleTreeLikelihood_update_three_nodes(tlk, node);
	}
}


// Scale nodes above calibrations if any
// if the root is fixed we still try to scale the other nodes
void SingleTreeLikelihood_scaler( SingleTreeLikelihood *tlk, Node *node, const double scaler ){
	if ( node == NULL || Node_isleaf(node) ) {
		return;
	}
	
	Constraint *cnstr = node->height->cnstr;
	if ( Node_isroot(node) || !(Constraint_lower_fixed(cnstr) || Constraint_upper_fixed(cnstr)) ) {
		SingleTreeLikelihood_scaler(tlk, node->left, scaler);
		SingleTreeLikelihood_scaler(tlk, node->right, scaler);
		
        if ( Node_isroot(node) &&  Constraint_lower_fixed(cnstr) && Constraint_upper_fixed(cnstr)){
            return;
        }
        
        if ( scaler > 1) {
            Node_set_height(node, Node_height( node )*scaler);
            SingleTreeLikelihood_update_three_nodes(tlk, node);
        }
        else if( scaler < 1 ){
            double height = Node_height( node );
            double max_height_son = dmax( Node_height( Node_left(node)), Node_height( Node_right(node) ) );
            double bl = (height  - max_height_son)*scaler + max_height_son;
            Node_set_height(node, bl);
            SingleTreeLikelihood_update_three_nodes(tlk, node);
        }
        
	}
}

void SingleTreeLikelihood_scale_root( SingleTreeLikelihood *tlk, const double scaler ){
	Node *node = Tree_root(tlk->tree);
		
    if ( scaler > 1) {
        Node_set_height(node, Node_height( node )*scaler);
        SingleTreeLikelihood_update_all_nodes(tlk);
    }
    else if( scaler < 1 ){
        double height = Node_height( node );
        double max_height_son = dmax( Node_height( Node_left(node)), Node_height( Node_right(node) ) );
        double bl = (height  - max_height_son)*scaler + max_height_son;
        Node_set_height(node, bl);
        SingleTreeLikelihood_update_all_nodes(tlk);
    }
}

//TODO: Only works for codon and nucleotide without SSE or AVX
void SingleTreeLikelihood_set_nthreads( SingleTreeLikelihood *tlk, int nthreads ){
    #if defined (_OPENMP) && !((SSE3_ENABLED) || (AVX_ENABLED))
    //#if defined _OPENMP

     if( tlk->sm->nstate >= 60 || tlk->sm->nstate == 4 ){
        if( nthreads > 1 ){
            tlk->nthreads = nthreads;
            if( tlk->sm->nstate == 4 ){
                tlk->update_partials = update_partials_4_openmp;
            }
            else {
                tlk->update_partials = update_partials_codon_openmp;
            }            
        }
        else {
            tlk->nthreads = 1;
            if( tlk->sm->nstate == 4 ){
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
	update_upper_partials(tlk, Tree_root(tlk->tree));
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
	
	for (int i = 0; i < tlk->sm->cat_count; i++) {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		if( tlk->use_SIMD ){
			if( Node_isleaf(n) ){
				tlk->sm->m->p_t_transpose(tlk->sm->m,
										  bl * tlk->sm->get_rate(tlk->sm, i),
										  &tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
			}
			else {
				tlk->sm->m->p_t(tlk->sm->m,
								bl * tlk->sm->get_rate(tlk->sm, i),
								&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
			}
		}
		else{
			tlk->sm->m->p_t(tlk->sm->m,
							bl * tlk->sm->get_rate(tlk->sm, i),
							&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
		}
#else
		tlk->sm->m->p_t(tlk->sm->m,
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

			for (int i = 0; i < tlk->sm->cat_count; i++) {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
				if( tlk->use_SIMD ){
					if( Node_isleaf(n) ){
						tlk->sm->m->p_t_transpose(tlk->sm->m,
												  bl * tlk->sm->get_rate(tlk->sm, i),
												  &tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
					}
					else {
						tlk->sm->m->p_t(tlk->sm->m,
										bl * tlk->sm->get_rate(tlk->sm, i),
										&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
					}
				}
				else{
					tlk->sm->m->p_t(tlk->sm->m,
									bl * tlk->sm->get_rate(tlk->sm, i),
									&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
				}
#else
				tlk->sm->m->p_t(tlk->sm->m,
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

void SingleTreeLikelihood_scalePartials( SingleTreeLikelihood *tlk, int nodeIndex ) {
	int u = 0;
	int k,j;
	double* partials = tlk->partials[tlk->current_partials_indexes[nodeIndex]][nodeIndex];
	double* scaling_factors = tlk->scaling_factors[tlk->current_partials_indexes[nodeIndex]][nodeIndex];
	for ( int i = 0; i < tlk->sp->count; i++ ) {
		
		double scaleFactor = 0.0;
		int v = u;
		for ( k = 0; k < tlk->cat_count; k++ ) {
			for ( j = 0; j < tlk->sm->nstate; j++ ) {
				if ( partials[v] > scaleFactor ) {
					scaleFactor = partials[v];
				}
				v++;
			}
			v += ( tlk->sp->count - 1 ) * tlk->sm->nstate;
		}
		
		if ( scaleFactor < tlk->scaling_threshold ) {
			
			v = u;
			for ( k = 0; k < tlk->cat_count; k++ ) {
				for ( j = 0; j < tlk->sm->nstate; j++ ) {
					partials[v] /= scaleFactor;
					v++;
				}
				v += (tlk->sp->count - 1) * tlk->sm->nstate;
			}
			scaling_factors[i] = log(scaleFactor);
		}
		else {
			scaling_factors[i] = 0.0;
		}
		u += tlk->sm->nstate;
	}
}

double getLogScalingFactor( const SingleTreeLikelihood *tlk, int pattern ) {
	double log_scale_factor = 0.0;
	if ( tlk->scale ) {
		for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
			log_scale_factor += tlk->scaling_factors[tlk->current_partials_indexes[i]][i][pattern];
		}
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
		if( tlk->sm->nstate & 1 ){
			tlk->use_SIMD = false;
		}
	}
	
	if ( tlk->use_SIMD ){
		if (tlk->sm->nstate == 4 ) {
			tlk->update_partials      = update_partials_4_SSE;
			tlk->update_partials_flexible = update_partials_flexible_4_SSE;
			tlk->integrate_partials   = integrate_partials_4_SSE;
			tlk->node_log_likelihoods = node_log_likelihoods_4_SSE;
			tlk->calculate_branch_likelihood = calculate_branch_partials_4_SSE;
		}
		else if ( tlk->sm->nstate == 20 ) {
			tlk->update_partials      = update_partials_20_SSE;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			tlk->calculate_branch_likelihood = calculate_branch_partials_20_SSE;
		}
		else if ( tlk->sm->nstate >= 60 ) {
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
			tlk->calculate_branch_likelihood = calculate_branch_partials_even_SSE;
			//TODO: what about the odd case
			// no upper
		}
	}
	else {
		if ( tlk->sm->nstate == 4 ) {
			tlk->update_partials      = update_partials_4;
			tlk->update_partials_flexible = update_partials_flexible_4;
			tlk->integrate_partials   = integrate_partials_4;
			tlk->node_log_likelihoods = node_log_likelihoods_4;
			tlk->calculate_branch_likelihood = calculate_branch_partials_4;
		}
		else if ( tlk->sm->nstate == 20 ) {
			tlk->update_partials      = update_partials_general;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			tlk->calculate_branch_likelihood = calculate_branch_partials_20;
		}
		else if( tlk->sm->nstate >= 60 ){
			tlk->update_partials      = update_partials_codon;
			tlk->integrate_partials   = integrate_partials_codon;
			tlk->node_log_likelihoods = node_log_likelihoods_codon;
		}
		else {
			tlk->update_partials      = update_partials_general;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			tlk->calculate_branch_likelihood = calculate_branch_partials;
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
	return tlk->sm->m->get_frequencies(tlk->sm->m);
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
	
	double *partials   = dvector( tlk->sp->count*tlk->sm->nstate );
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

// calculate upper partials as used by derivatives
void update_upper_partials(SingleTreeLikelihood *tlk, Node* node){
	int idNode = Node_id(node);
	tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
	
	if(!Node_isroot(node)){
		Node* parent = Node_parent(node);
		Node* sibling = Node_sibling(node);
		int idSibling = Node_id(sibling);
		
		tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
		
		if(!Node_isroot(parent)){
			int idMatrix = Node_id(parent);
	
			tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], tlk->upper_partial_indexes[Node_id(parent)], idMatrix, idSibling, idSibling);
		}
		else{
			tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], idSibling, idSibling, -1, -1);
		}
	}
	
	if(!Node_isleaf(node)){
		update_upper_partials(tlk, Node_left(node));
		update_upper_partials(tlk, Node_right(node));
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

void setup_upper_indexes(SingleTreeLikelihood *tlk, Node* node){
	
	if(!Node_isroot(node)){
		Node* parent = Node_parent(node);
		Node* sibling = Node_right(parent);
		int idNode = Node_id(node);
		
		if (sibling == node) {
			sibling = Node_left(parent);
		}
		
		if(!Node_isroot(parent)){
			Node* grandParent = Node_parent(parent);
			
			int idMatrix = Node_id(parent);
			
			// The sons of the right node of the root are going to use the lower partials of the left node of the root
			if(Node_isroot(grandParent) && Node_right(grandParent) == parent){
				idMatrix = Node_id(Node_left(grandParent));
			}
			
			tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
		}
		// We dont need to calculate upper partials for the children of the root as it is using the lower partials of its sibling
		// Left node of the root
		else if(Node_left(parent) == node){
			tlk->upper_partial_indexes[idNode] = Node_id(Node_right(parent));
		}
		else{
			tlk->upper_partial_indexes[idNode] = Node_id(Node_left(parent));
		}
	}
	else{
		tlk->upper_partial_indexes[Node_id(node)] = -1;
	}
	
	if(!Node_isleaf(node)){
		setup_upper_indexes(tlk, Node_left(node));
		setup_upper_indexes(tlk, Node_right(node));
	}
}

// Derivative of the log likelihood function
// Assumes that root_partials are calculated correctly (i.e. _calculate_uppper is not called)
// This function does not change anything regarding partials, matrices, root_partials and pattern_lk
// TODO: Make thread safe
double dlnldt_uppper( SingleTreeLikelihood *tlk, const double* pattern_likelihoods, const double* pattern_dlikelihoods ){
	const int patternCount = tlk->sp->count;
	double dlnl = 0;
	
	for ( int k = 0; k < patternCount; k++ ) {
		dlnl += pattern_dlikelihoods[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
	}
	
	return dlnl;
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
	
	for (int i = 0; i < tlk->sm->cat_count; i++) {
		
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		if( tlk->use_SIMD ){
			if( Node_isleaf(node) ){
				tlk->sm->m->dp_dt_transpose(tlk->sm->m,
											bl * tlk->sm->get_rate(tlk->sm, i),
											&mat[i*tlk->matrix_size]);
			}
			else {
				tlk->sm->m->dp_dt(tlk->sm->m,
								  bl * tlk->sm->get_rate(tlk->sm, i),
								  &mat[i*tlk->matrix_size]);
			}
		}
		else{
			tlk->sm->m->dp_dt(tlk->sm->m,
							  bl * tlk->sm->get_rate(tlk->sm, i),
							  &mat[i*tlk->matrix_size]);
		}
#else
		tlk->sm->m->dp_dt(tlk->sm->m,
						  bl * tlk->sm->get_rate(tlk->sm, i),
						  &mat[i*tlk->matrix_size]);
#endif
		for (int k = 0; k < tlk->matrix_size; k++) {
			mat[i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
		}
	}
	
	
	const int nstate   = tlk->sm->nstate;
	const int patternCount = tlk->sp->count;
	double* root_partials = tlk->root_partials  + patternCount*nstate;
	const double* freqs = tlk->get_root_frequencies(tlk);
	
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
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
	
	for (int i = 0; i < tlk->sm->cat_count; i++) {
		const double rate = tlk->sm->get_rate(tlk->sm, i);
		
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		if( tlk->use_SIMD ){
			if( Node_isleaf(node) ){
				tlk->sm->m->d2p_d2t_transpose(tlk->sm->m,
											  bl * rate,
											  &mat[i*tlk->matrix_size]);
			}
			else {
				tlk->sm->m->d2p_d2t(tlk->sm->m,
									bl * rate,
									&mat[i*tlk->matrix_size]);
			}
		}
		else{
			tlk->sm->m->d2p_d2t(tlk->sm->m,
								bl * rate,
								&mat[i*tlk->matrix_size]);
		}
#else
		tlk->sm->m->d2p_d2t(tlk->sm->m,
							bl * rate,
							&mat[i*tlk->matrix_size]);
#endif
		for (int k = 0; k < tlk->matrix_size; k++) {
			mat[i*tlk->matrix_size+k] *= rate * rate;
		}
	}
	const int nstate   = tlk->sm->nstate;
	const int patternCount = tlk->sp->count;
	double* root_partials = tlk->root_partials  + (patternCount*nstate)*2;
	double* pattern_d2lnl = tlk->pattern_lk + patternCount*3;
	const double* freqs = tlk->get_root_frequencies(tlk);
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
	tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
	
	//tlk->node_log_likelihoods( tlk, root_partials, tlk->sm->m->_freqs, pattern_d2lnl);
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


void calculate_dpartials_dfreqs(SingleTreeLikelihood *tlk, int index, double* dpartials, const double* pattern_likelihoods){
    Node **nodes = Tree_nodes(tlk->tree);
    
    const int rootId = Node_id(Tree_root(tlk->tree));
    const int rightId = Node_id(Node_right(Tree_root(tlk->tree)));
    
    const int nstate   = tlk->sm->nstate;
    const int patternCount = tlk->sp->count;
    
    int tempMatrixId = Node_id(Tree_root(tlk->tree));
    double* mat = tlk->matrices[tlk->current_matrices_indexes[tempMatrixId]][tempMatrixId];
    double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->sm->nstate;
    memset(dpartials, 0, sizeof(double)*patternCount);
    
    const int freqIndex = index - Parameters_count(tlk->sm->m->rates);
	const double* freqs = tlk->get_root_frequencies(tlk);
    double* transpose = dvector(nstate*nstate*tlk->sm->cat_count);

    int v = 0;
    for ( int k = 0; k < patternCount; k++ ) {
        dpartials[k] = 0;
        for (int jj = 0; jj < nstate-1; jj++) {
            if(jj == freqIndex){
                dpartials[k] += (freqs[nstate-1] - freqs[jj]*freqs[nstate-1])*tlk->root_partials[v];
            }
            else{
                dpartials[k] -= freqs[jj]*freqs[nstate-1]*tlk->root_partials[v];
            }
            v++;
        }
        dpartials[k] -= freqs[nstate-1]*freqs[nstate-1]*tlk->root_partials[v];
        v++;
    }
    
    for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
        Node *n = nodes[j];
        const int nodeId = Node_id(n);
        
        if ( nodeId == rootId || nodeId == rightId ) {
            continue;
        }
        
        const double bl = Node_distance(n);
        
        for (int i = 0; i < tlk->sm->cat_count; i++) {
    
            if(tlk->sm->m->dPdp != NULL){
                tlk->sm->m->dPdp(tlk->sm->m, index, &mat[i*tlk->matrix_size], bl* tlk->sm->get_rate(tlk->sm, i));
            }
            else{
                // Use finite differences
                fprintf(stderr, "Not yet implemented: calculate_dpartials_dQ\n");
                exit(1);
            }
        }
        
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
        if( tlk->use_SIMD ){
            if(Node_isleaf(n)){
                for(int l = 0; l < tlk->sm->cat_count; l++){
                    for(int ii = 0; ii < nstate; ii++){
                        for(int jj = 0; jj < nstate; jj++){
                            transpose[l * tlk->matrix_size+ii*4+jj] = mat[l * tlk->matrix_size+jj*4+ii];
                        }
                    }
                }
                memcpy(mat, transpose, sizeof(double)*tlk->matrix_size*tlk->sm->cat_count);
            }
        }
#endif
		double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[rootId]];
        tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
        tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
		
        v = 0;
        for ( int k = 0; k < patternCount; k++ ) {
            for (int jj = 0; jj < nstate; jj++) {
                dpartials[k] += freqs[jj]*root_partials[v];
                v++;
            }
        }
    }
    free(transpose);
}

double calculate_dlnl_dQ( SingleTreeLikelihood *tlk, int index, const double* pattern_likelihoods ){
	
	Node **nodes = Tree_nodes(tlk->tree);
	
	const int rootId = Node_id(Tree_root(tlk->tree));
	const int rightId = Node_id(Node_right(Tree_root(tlk->tree)));
	
	const int nstate   = tlk->sm->nstate;
	const int patternCount = tlk->sp->count;
	bool time_mode = Tree_is_time_mode(tlk->tree);
	
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[0][tempMatrixId];
	double* mat_transpose = tlk->matrices[1][tempMatrixId];
	double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->sm->nstate;
	double* pattern_dlnl = tlk->pattern_lk + tlk->sp->count*2; // pattern_likelihoods points to tlk->pattern_lk + tlk->sp->count
	memset(pattern_dlnl, 0, sizeof(double)*tlk->sp->count);
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[rootId]];
    
    double dlnldQ = 0;
    tlk->sm->m->dQ_need_update = true;
    
    const int freqIndex = index - Parameters_count(tlk->sm->m->rates);
	const double* freqs = tlk->get_root_frequencies(tlk);
	const double phi1 = Parameters_value(tlk->sm->m->simplex->parameters, 0);
	const double phi2 = Parameters_value(tlk->sm->m->simplex->parameters, 1);
	const double phi3 = Parameters_value(tlk->sm->m->simplex->parameters, 2);
	
    // Frequencies
    if(index >= Parameters_count(tlk->sm->m->rates)){
		double dphi1[4] = {1, -phi2, (phi2-1)*phi3, -(phi2 - 1)*phi3 + phi2 - 1};
		double dphi2[4] = {0, -phi1+1, (phi1-1.0)*phi3, -(phi1 - 1)*phi3 + phi1 - 1};
		double dphi3[4] = {0, 0, (phi1 - 1)*phi2 - phi1 + 1, -(phi1 - 1)*phi2 + phi1 - 1};
		double* dphi;
		if (freqIndex == 0) {
			dphi = dphi1;
		}
		else if (freqIndex == 1) {
			dphi = dphi2;
		}
		else{
			dphi = dphi3;
		}
		
        //calculate_dpartials_dfreqs(tlk, index, dpartials, pattern_likelihoods);
        int v = 0;
        for ( int k = 0; k < patternCount; k++ ) {
            pattern_dlnl[k] = 0;
            for (int jj = 0; jj < nstate; jj++) {
                    pattern_dlnl[k] += dphi[jj]*tlk->root_partials[v];
				
                v++;
            }
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
				tlk->sm->m->dPdp(tlk->sm->m, index, &mat[i*tlk->matrix_size], bl* tlk->sm->get_rate(tlk->sm, i));
			}
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			if( tlk->use_SIMD ){
				if(Node_isleaf(n)){
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
			
			tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
			tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
			
            v = 0;
            for ( int k = 0; k < patternCount; k++ ) {
                for (int jj = 0; jj < nstate; jj++) {
					pattern_dlnl[k] += freqs[jj]*root_partials[v];
                    v++;
                }
            }
        }
        for ( int k = 0; k < patternCount; k++ ) {
            dlnldQ += pattern_dlnl[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
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
                tlk->sm->m->dPdp(tlk->sm->m, index, &mat[i*tlk->matrix_size], bl* tlk->sm->get_rate(tlk->sm, i));
            }
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD ){
                if(Node_isleaf(n)){
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

            tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
			tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
            int v = 0;
			int i = 0;
			const double* freqs = tlk->get_root_frequencies(tlk);
			
            for ( int k = 0; k < patternCount; k++ ) {
                for ( i = 0; i < nstate; i++ ) {
                    
                    pattern_dlnl[k] += freqs[i] * root_partials[v];
                    v++;
                }
                
                if ( tlk->scale ) {
                    //printf("scaling\n");
                    tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
                }
            }
        }
        
        for ( int k = 0; k < patternCount; k++ ) {
            dlnldQ += pattern_dlnl[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
        }
    }
//	SingleTreeLikelihood_update_all_nodes(tlk);
	return dlnldQ;
}

double calculate_dlnl_dWeibull( SingleTreeLikelihood *tlk, const double* pattern_likelihoods ){
	
	Node **nodes = Tree_nodes(tlk->tree);
	
	const int rootId = Node_id(Tree_root(tlk->tree));
	const int rightId = Node_id(Node_right(Tree_root(tlk->tree)));
	
	const int nstate   = tlk->sm->nstate;
	const int patternCount = tlk->sp->count;
	bool time_mode = Tree_is_time_mode(tlk->tree);
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[0][tempMatrixId];
	double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->sm->nstate;
	double* pattern_dlnl = tlk->pattern_lk + tlk->sp->count*2; // pattern_likelihoods points to tlk->pattern_lk + tlk->sp->count
	memset(pattern_dlnl, 0, sizeof(double)*tlk->sp->count);
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[rootId]];
	
	double dlnldW = 0;
	
	const double* freqs = tlk->get_root_frequencies(tlk);
	double* proportions = tlk->sm->get_proportions(tlk->sm);
	double cat = tlk->sm->cat_count;
	double shape = Parameters_value(tlk->sm->rates, 0);
	
	for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
		Node *n = nodes[j];
		const int nodeId = Node_id(n);
		
		if ( nodeId == rootId || (nodeId == rightId && !time_mode) ) {
			continue;
		}
		
		double bl;
		if (time_mode) {
			bl = Node_time_elapsed(n) * tlk->bm->get(tlk->bm, n);
		}
		else{
			bl = Node_distance(n);
		}
		
		double derivsum = 0;
		double sum = 0;
		
		for (int i = 0; i < tlk->sm->cat_count; i++) {
			double prob = (2.0*i + 1.0)/(2.0*cat);
			derivsum += -pow(-log(1.0 - prob), 1.0/shape) * log(-log(1.0 - prob))/shape/shape*proportions[i];
			sum += pow(-log(1.0 - prob), 1.0/shape)*proportions[i];
		}
		
		for (int i = 0; i < tlk->sm->cat_count; i++) {
			
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			if( tlk->use_SIMD ){
				if( Node_isleaf(n) ){
					tlk->sm->m->dp_dt_transpose(tlk->sm->m,
												bl * tlk->sm->get_rate(tlk->sm, i),
												&mat[i*tlk->matrix_size]);
				}
				else {
					tlk->sm->m->dp_dt(tlk->sm->m,
									  bl * tlk->sm->get_rate(tlk->sm, i),
									  &mat[i*tlk->matrix_size]);
				}
			}
			else{
				tlk->sm->m->dp_dt(tlk->sm->m,
								  bl * tlk->sm->get_rate(tlk->sm, i),
								  &mat[i*tlk->matrix_size]);
			}
#else
			tlk->sm->m->dp_dt(tlk->sm->m,
							  bl * tlk->sm->get_rate(tlk->sm, i),
							  &mat[i*tlk->matrix_size]);
#endif
			
			double prob = (2.0*i + 1.0)/(2.0*cat);
			double deriv = -pow(-log(1.0 - prob), 1.0/shape) * log(-log(1.0 - prob))/shape/shape;
			deriv = deriv/sum - (pow(-log(1.0 - prob), 1.0/shape)*derivsum)/sum/sum;
			
			for (int k = 0; k < tlk->matrix_size; k++) {
				mat[i*tlk->matrix_size+k] *= bl*deriv;
			}
		}

		tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
		tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
		int v = 0;
		int i = 0;
		
		for ( int k = 0; k < patternCount; k++ ) {
			for ( i = 0; i < nstate; i++ ) {
				
				pattern_dlnl[k] += freqs[i] * root_partials[v];
				v++;
			}
			
			if ( tlk->scale ) {
				//printf("scaling\n");
				tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
			}
		}
	}
	
	for ( int k = 0; k < patternCount; k++ ) {
		dlnldW += pattern_dlnl[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
	}
	
	return dlnldW;
}

#pragma mark -
#pragma mark derivative heights

// calculate first step to get the derivative of P matrix wrt branch length
// store Q.P in mat
void get_Q_times_P(SingleTreeLikelihood* tlk, const Node* node, double* mat){
	int cat_count = tlk->sm->cat_count;
	double bl = Node_time_elapsed(node) * tlk->bm->get(tlk->bm, node);
	
	for (int i = 0; i < cat_count; i++) {
		
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		if( tlk->use_SIMD ){
			if( Node_isleaf(node) ){
				tlk->sm->m->dp_dt_transpose(tlk->sm->m,
											bl * tlk->sm->get_rate(tlk->sm, i),
											&mat[i*tlk->matrix_size]);
			}
			else {
				tlk->sm->m->dp_dt(tlk->sm->m,
								  bl * tlk->sm->get_rate(tlk->sm, i),
								  &mat[i*tlk->matrix_size]);
			}
		}
		else{
			tlk->sm->m->dp_dt(tlk->sm->m,
							  bl * tlk->sm->get_rate(tlk->sm, i),
							  &mat[i*tlk->matrix_size]);
		}
#else
		tlk->sm->m->dp_dt(tlk->sm->m,
						  bl * tlk->sm->get_rate(tlk->sm, i),
						  &mat[i*tlk->matrix_size]);
#endif
	}
}

// Derivative of P matrix wrt to proportion p
// every node branch should depend on p
// Noderef is the node associated with the proportion

// if it is called it means the derivative is not 0
// Three cases:
// 1) derivative wrt h
// 2) p2p1h - p3p2p1h derivative wrt p3 => index == node_map_index
// 3) p2p1h - p3p2p1h derivative wrt a p but not p3
void derivative_P_wrt_p(SingleTreeLikelihood* tlk, Node* node, const Node* noderef){
	
	int node_id = Node_id(node);
	double* mat = tlk->matrices[1-tlk->current_matrices_indexes[node_id]][node_id];
	int cat_count = tlk->sm->cat_count;
	get_Q_times_P(tlk, node, mat);
	
	unsigned *map = get_reparam_map(tlk->tree);
	Parameters* reparams = get_reparams(tlk->tree);
	
	// derivative of branch length wrt p
	double d = 1;
	if (Node_isroot(noderef)) {
		Node* n = Node_parent(node);
		while(!Node_isroot(n)){
			d *= Parameters_value(reparams, map[Node_id(n)]);
			n = n->parent;
		}
		if(!Node_isleaf(node)){
			d = d - d*Parameters_value(reparams, map[Node_id(node)]);
		}
	}
	else if (noderef == node){
		d = -Node_height(Node_parent(node));
	}
	else{
		double prop = Parameters_value(reparams, map[Node_id(noderef)]);
		if(Node_isleaf(node)){
			d = Node_height(Node_parent(node))/prop;
		}
		else{
			d = Node_time_elapsed(node)/prop;
		}
	}
	
	d *= tlk->bm->get(tlk->bm, node);
	
	for (int i = 0; i < cat_count; i++) {
		double v = d * tlk->sm->get_rate(tlk->sm, i);
		for (int k = 0; k < tlk->matrix_size; k++) {
			mat[i*tlk->matrix_size+k] *= v;
		}
	}
}

void derivative_P_wrt_p_heterochronous(SingleTreeLikelihood* tlk, Node* node, const Node* noderef){
	
	int node_id = Node_id(node);
	double* mat = tlk->matrices[1-tlk->current_matrices_indexes[node_id]][node_id];
	int cat_count = tlk->sm->cat_count;
	get_Q_times_P(tlk, node, mat);
	
	unsigned *map = get_reparam_map(tlk->tree);
	Parameters* reparams = get_reparams(tlk->tree);
	
	// derivative of branch length wrt p
	double d = 1;
	if (Node_isroot(noderef)) {
		Node* n = Node_parent(node);
		while(!Node_isroot(n)){
			d *= Parameters_value(reparams, map[Node_id(n)]);
			n = n->parent;
		}
		if(!Node_isleaf(node)){
			d = d - d*Parameters_value(reparams, map[Node_id(node)]);
		}
		
	}
	else if (noderef == node){
		double* lowers = dvector(Tree_node_count(tlk->tree));
		Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
		for ( int i = 0; i < Tree_node_count(tlk->tree); i++) {
			Node* n = nodes[i];
			if( Node_isleaf(n) ){
				lowers[Node_id(n)] = Node_height(n);
			}
			else{
				lowers[Node_id(n)] = dmax(lowers[Node_left(n)->id], lowers[Node_right(n)->id]);
			}
		}
		d =  -(Node_height(Node_parent(node)) - lowers[Node_id(node)]);
		free(lowers);
	}
	else{
		double* lowers = dvector(Tree_node_count(tlk->tree));
		Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
		for ( int i = 0; i < Tree_node_count(tlk->tree); i++) {
			Node* n = nodes[i];
			if( Node_isleaf(n) ){
				lowers[Node_id(n)] = Node_height(n);
			}
			else{
				lowers[Node_id(n)] = dmax(lowers[Node_left(n)->id], lowers[Node_right(n)->id]);
			}
		}
		
		Node* n = Node_parent(node);
		while(n != noderef){
			d *= Parameters_value(reparams, map[Node_id(n)]);
			n = n->parent;
		}
		d *= Node_height(n->parent) - lowers[Node_id(n)];
		if(!Node_isleaf(node)){
			d = d - d*Parameters_value(reparams, map[Node_id(node)]);
		}
		free(lowers);
	}
	
	d *= tlk->bm->get(tlk->bm, node);
	
    for (int i = 0; i < cat_count; i++) {
        double v = d * tlk->sm->get_rate(tlk->sm, i);
        for (int k = 0; k < tlk->matrix_size; k++) {
            mat[i*tlk->matrix_size+k] *= v;
        }
    }
}

// called on node that is reparameterized with p
void calculate_dldp_uppper_aux(SingleTreeLikelihood *tlk, Node *node, const Node* noderef){
	if(!Node_isroot(node)){
		if (Tree_homochronous(tlk->tree)) {
			derivative_P_wrt_p(tlk, node, noderef);
		}
		else{
			derivative_P_wrt_p_heterochronous(tlk, node, noderef);
		}

//#ifndef NDEBUG
//		if(Tree_homochronous(tlk->tree) && !Node_isleaf(node)){
//			assert(Node_height(node->parent)*Parameters_value(get_reparams(tlk->tree), get_reparam_map(tlk->tree)[node->id]) == Node_height(node) );
//		}
//		else if(!Tree_homochronous(tlk->tree) && !Node_isleaf(node)){
//			double* lowers = dvector(Tree_node_count(tlk->tree));
//			Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
//			for ( int i = 0; i < Tree_node_count(tlk->tree); i++) {
//				Node* node = nodes[i];
//				if( Node_isleaf(node) ){
//					lowers[Node_id(node)] = Node_height(node);
//				}
//				else{
//					lowers[Node_id(node)] = dmax(lowers[Node_left(node)->id], lowers[Node_right(node)->id]);
//				}
//			}
//            if(lowers[Node_id(node)] + (Node_height(node->parent) - lowers[Node_id(node)])*Parameters_value(get_reparams(tlk->tree), get_reparam_map(tlk->tree)[node->id]) != Node_height(node)){
//                fprintf(stderr, "%f %f\n",  lowers[Node_id(node)] + (Node_height(node->parent) - lowers[Node_id(node)])*Parameters_value(get_reparams(tlk->tree), get_reparam_map(tlk->tree)[node->id]), Node_height(node));
//            }
//
//			assert(lowers[Node_id(node)] + (Node_height(node->parent) - lowers[Node_id(node)])*Parameters_value(get_reparams(tlk->tree), get_reparam_map(tlk->tree)[node->id]) == Node_height(node) );
//			free(lowers);
//		}
//#endif
	}
	
	if(!Node_isleaf(node)){
		calculate_dldp_uppper_aux(tlk, node->left, noderef);
		calculate_dldp_uppper_aux(tlk, node->right, noderef);
		
		int node_id = Node_id(node);
		int indx_child1 = Node_id(Node_left(node));
		int indx_child2 = Node_id(Node_right(node));
		bool isleaf_child1 = Node_isleaf(Node_left(node));
		bool isleaf_child2 = Node_isleaf(Node_right(node));
		
		// P matrices and partials of 2 children
		double* mat_1 = tlk->matrices[tlk->current_matrices_indexes[indx_child1]][indx_child1];
		double* mat_2 = tlk->matrices[tlk->current_matrices_indexes[indx_child2]][indx_child2];
		double* partials_1 = tlk->partials[tlk->current_partials_indexes[indx_child1]][indx_child1];
		double* partials_2 = tlk->partials[tlk->current_partials_indexes[indx_child2]][indx_child2];
		
		// derivative of P matrices and partials of 2 children
		double* mat_deriv_1 = tlk->matrices[1-tlk->current_matrices_indexes[indx_child1]][indx_child1];
		double* mat_deriv_2 = tlk->matrices[1-tlk->current_matrices_indexes[indx_child2]][indx_child2];
		double* partials_deriv_1 = tlk->partials[1-tlk->current_partials_indexes[indx_child1]][indx_child1];
		double* partials_deriv_2 = tlk->partials[1-tlk->current_partials_indexes[indx_child2]][indx_child2];
		
		double* partials_deriv = tlk->partials[1-tlk->current_partials_indexes[node_id]][node_id]; // store derivatives for node
		
		int dim = tlk->sp->count*tlk->sm->nstate*tlk->sm->cat_count;
		
		double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
		double* spare_partials2 = tlk->partials[1][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
		
		// child 1
		tlk->update_partials_flexible( tlk, partials_deriv,
									  indx_child1,
									  partials_1,  //partial child 1
									  mat_deriv_1, //derivative P child 1
									  -1,
									  NULL,
									  NULL);

		// not a leaf
		if(!isleaf_child1){
			tlk->update_partials_flexible( tlk, spare_partials,
										  indx_child1,
										  partials_deriv_1, // derivative partial child 1
										  mat_1,            // P child 1
										  -1,
										  NULL,
										  NULL);
			
			for (int k = 0; k < dim; k++) {
				partials_deriv[k] += spare_partials[k];
			}
		}

		tlk->update_partials_flexible( tlk, spare_partials,
									  indx_child2,
									  partials_2, // partial child 2
									  mat_2,      // P child 2
									  -1,
									  NULL,
									  NULL);

		for (int k = 0; k < dim; k++) {
			partials_deriv[k] *= spare_partials[k];
		}
		
		// child 2
		tlk->update_partials_flexible( tlk, spare_partials2,
									  indx_child2,
									  partials_2,  //partial child 2
									  mat_deriv_2, //derivative P child 1
									  -1,
									  NULL,
									  NULL);
		
		// not a leaf
		if(!isleaf_child2){
			tlk->update_partials_flexible( tlk, spare_partials,
										  indx_child2,
										  partials_deriv_2, // derivative partial child 2
										  mat_2,            // P child 2
										  -1,
										  NULL,
										  NULL);
			
			for (int k = 0; k < dim; k++) {
				spare_partials2[k] += spare_partials[k];
			}
		}
		
		tlk->update_partials_flexible( tlk, spare_partials,
									  indx_child1,
									  partials_1, // partial child 1
									  mat_1,      // P child 1
									  -1,
									  NULL,
									  NULL);

		for (int k = 0; k < dim; k++) {
			partials_deriv[k] += spare_partials2[k]*spare_partials[k];
		}
	}
}

void calculate_dldp_uppper_reparam( SingleTreeLikelihood *tlk, Node *node, double* pattern_dlikelihoods ){
	const int nstate   = tlk->sm->nstate;
	const int patternCount = tlk->sp->count;
	double* root_partials = tlk->root_partials  + patternCount*nstate;
	
	// this function return \frac{\partial d_{node}}{\partial d_p}
	calculate_dldp_uppper_aux(tlk, node, node);
	
	if(Node_isroot(node)){
		tlk->integrate_partials(tlk, tlk->partials[1-tlk->current_partials_indexes[Tree_root(tlk->tree)->id]][Tree_root(tlk->tree)->id], tlk->sm->get_proportions(tlk->sm), root_partials );
	}
	else{
		const int nodeId = Node_id(node);
		double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
		double* spare_partials2 = tlk->partials[1][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
		memset(spare_partials, 0, sizeof(double)*patternCount*nstate*tlk->sm->cat_count);
		memset(spare_partials2, 0, sizeof(double)*patternCount*nstate*tlk->sm->cat_count);
		
		// calculate P . \frac{\partial d_{node}}{\partial d_p}
		tlk->update_partials_flexible( tlk, spare_partials,
									  nodeId,
									  tlk->partials[1-tlk->current_partials_indexes[nodeId]][nodeId],
									  tlk->matrices[tlk->current_matrices_indexes[nodeId]][nodeId],
									  -1,
									  NULL,
									  NULL);
		
		// calculate \frac{\partial P_{node}}{\partial d_p} . d_{node}
		tlk->update_partials_flexible( tlk, spare_partials2,
									  nodeId,
									  tlk->partials[tlk->current_partials_indexes[nodeId]][nodeId],
									  tlk->matrices[1-tlk->current_matrices_indexes[nodeId]][nodeId],
									  -1,
									  NULL,
									  NULL);
		
		// now calcuate \frac{\partial \delta_{node}}{\partial d_p} # upper_{node}
		double* upper_partials = tlk->partials[tlk->current_partials_indexes[tlk->upper_partial_indexes[nodeId]]][tlk->upper_partial_indexes[nodeId]];
		
		for (int k = 0; k < patternCount*nstate*tlk->sm->cat_count; k++) {
			spare_partials[k] = (spare_partials[k] + spare_partials2[k])*upper_partials[k];
		}
		
		tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
	}
	
	const double* freqs = tlk->get_root_frequencies(tlk);
	int v = 0;
	for ( int k = 0; k < patternCount; k++ ) {
		
		pattern_dlikelihoods[k] = 0;
		for ( int i = 0; i < nstate; i++ ) {
			pattern_dlikelihoods[k] += freqs[i] * root_partials[v];
			v++;
		}
		
		if ( tlk->scale ) {
			printf("scaling\n");
			pattern_dlikelihoods[k] += getLogScalingFactor( tlk, k);
		}
	}
}
#pragma mark -
#pragma mark Derivative clock

double calculate_dl_dr(SingleTreeLikelihood* tlk, double* pattern_likelihoods){
	
	Node **nodes = Tree_nodes(tlk->tree);
	
	const int rootId = Node_id(Tree_root(tlk->tree));
	
	const int nstate   = tlk->sm->nstate;
	const int patternCount = tlk->sp->count;
	
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[0][tempMatrixId];
	double* mat_transpose = tlk->matrices[1][tempMatrixId];
	double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->sm->nstate;
	double* pattern_dlnl = tlk->pattern_lk + tlk->sp->count*2;
	memset(pattern_dlnl, 0, sizeof(double)*tlk->sp->count);
	
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	
	double dlnldR = 0;
	
	const double* freqs = tlk->get_root_frequencies(tlk);
	double* proportions = tlk->sm->get_proportions(tlk->sm);
	double cat = tlk->sm->cat_count;
	
	
	for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
		Node *n = nodes[j];
		const int nodeId = Node_id(n);
		
		if ( nodeId == rootId ) {
			continue;
		}
		
		// calculate derivative of probability matrix wrt substitution rate
		get_Q_times_P(tlk, n, mat);
		for (int i = 0; i < cat; i++) {
			double v = Node_time_elapsed(n) * tlk->sm->get_rate(tlk->sm, i);
			for (int k = 0; k < tlk->matrix_size; k++) {
				mat[i*tlk->matrix_size+k] *= v;
			}
		}
		
		tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
		tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
		int v = 0;
		int i = 0;
		
		for ( int k = 0; k < patternCount; k++ ) {
			for ( i = 0; i < nstate; i++ ) {
				
				pattern_dlnl[k] += freqs[i] * root_partials[v];
				v++;
			}
			
			if ( tlk->scale ) {
				//printf("scaling\n");
				tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
			}
		}
	}
	
	for ( int k = 0; k < patternCount; k++ ) {
		dlnldR += pattern_dlnl[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
	}
	
	return dlnldR;
}

// one rate per branch
double calculate_dl_dr_discrete(SingleTreeLikelihood* tlk, Node* node, double* pattern_likelihoods){
    const int nstate   = tlk->sm->nstate;
    const int patternCount = tlk->sp->count;
    
    int tempMatrixId = Node_id(Tree_root(tlk->tree));
    double* mat = tlk->matrices[0][tempMatrixId];
    double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->sm->nstate;
    double* pattern_dlnl = tlk->pattern_lk + tlk->sp->count*2;
    memset(pattern_dlnl, 0, sizeof(double)*tlk->sp->count);
    
    double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
    
    double dlnldR = 0;
    
    const double* freqs = tlk->get_root_frequencies(tlk);
    double cat = tlk->sm->cat_count;
    const int nodeId = Node_id(node);
    
    // calculate derivative of probability matrix wrt substitution rate
    get_Q_times_P(tlk, node, mat);
    for (int i = 0; i < cat; i++) {
        double v = Node_time_elapsed(node) * tlk->sm->get_rate(tlk->sm, i);
        for (int k = 0; k < tlk->matrix_size; k++) {
            mat[i*tlk->matrix_size+k] *= v;
        }
    }
    
    tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
    tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), root_partials );
    int v = 0;
    int i = 0;
    
    for ( int k = 0; k < patternCount; k++ ) {
        for ( i = 0; i < nstate; i++ ) {
            pattern_dlnl[k] += freqs[i] * root_partials[v];
            v++;
        }
        
        if ( tlk->scale ) {
            //printf("scaling\n");
            tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
        }
    }
    
    for ( int k = 0; k < patternCount; k++ ) {
        dlnldR += pattern_dlnl[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
    }
    
    return dlnldR;
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
			update_upper_partials(tlk, Tree_root(tlk->tree));
		}
		else {
			Node *n = tlk->node_upper;
			
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
			if( tlk->use_SIMD && Node_isleaf(n) ){
				for (int i = 0; i < tlk->sm->cat_count; i++) {
					tlk->sm->m->p_t_transpose(tlk->sm->m,
											  Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
											  &tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
				}
			}
			else{
				for (int i = 0; i < tlk->sm->cat_count; i++) {
					tlk->sm->m->p_t(tlk->sm->m,
									Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
									&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
				}
			}
#else
			for (int i = 0; i < tlk->sm->cat_count; i++) {
				tlk->sm->m->p_t(tlk->sm->m,
								Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
								&tlk->matrices[tlk->current_matrices_indexes[Node_id(n)]][Node_id(n)][i*tlk->matrix_size]);
			}
#endif
			
			// update lower partials of the parent
			tlk->update_partials(tlk, Node_id(Node_parent(n)), Node_id(n), Node_id(n), Node_id(Node_sibling(n)), Node_id(Node_sibling(n)) );
			
			// update upper partials of its sibling and its sibling's descendents
			// although they technically don't need to be updated if node is a right node since its sibling&descendents won't be reused in postorder

			update_upper_partials(tlk, Node_sibling(n));
			
		}
	}
	//		SingleTreeLikelihood_update_all_nodes(tlk);
	//		tlk->calculate(tlk);
	//		calculate_upper(tlk, Tree_root(tlk->tree));
	
	const int nodeId = Node_id(node);
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tlk->current_matrices_indexes[tempMatrixId]][tempMatrixId];
	
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	if( tlk->use_SIMD && Node_isleaf(node) ){
		for (int i = 0; i < tlk->sm->cat_count; i++) {
			tlk->sm->m->p_t_transpose(tlk->sm->m,
									  Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
									  &mat[i*tlk->matrix_size]);
		}
	}
	else{
		for (int i = 0; i < tlk->sm->cat_count; i++) {
			tlk->sm->m->p_t(tlk->sm->m,
							Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
							&mat[i*tlk->matrix_size]);
		}
	}
#else
	for (int i = 0; i < tlk->sm->cat_count; i++) {
		tlk->sm->m->p_t(tlk->sm->m,
						Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
						&mat[i*tlk->matrix_size]);
	}
#endif
	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
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
