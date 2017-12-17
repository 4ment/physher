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

#include "pooledtreelikelihood.h"

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
double _calculate_uppper_height( SingleTreeLikelihood *tlk, Node *node );
double _calculate_uppper_height_sse( SingleTreeLikelihood *tlk, Node *node );

static bool _calculate_partials_noexp_integrate( SingleTreeLikelihood *tlk, Node *n  );



#pragma mark -
#pragma mark TreeLikelihoodModel

void _treelikelihood_handle_change( Model *self, Model *model, int index ){
	SingleTreeLikelihood *tlk = (SingleTreeLikelihood*)self->obj;
	//printf("%s %d\n", model->name, index);
	if ( strcmp(model->type, "tree") == 0 ) {
//		printf("node index %d\n", index);
		//SingleTreeLikelihood_update_one_node(tlk, index);
		tlk->update_nodes[index] = true;
		tlk->update = true;
	}
	else if ( strcmp(model->type, "branchmodel") == 0 ) {
		if(index == -1){
			SingleTreeLikelihood_update_all_nodes(tlk);
		}
		else{
			tlk->update_nodes[index] = true;
			tlk->update = true;
		}
	}
	else if ( strcmp(model->type, "sitemodel") == 0 ) {
		SingleTreeLikelihood_update_all_nodes(tlk);
	}
	else {
		fprintf(stderr, "%s of type %s\n", model->name, model->type);
		error("Unknown change in SingleLikelihood\n");
	}
}

double _singleTreeLikelihood_logP(Model *self){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	return tlk->calculate(tlk);
}

double _singleTreeLikelihood_dlogP(Model *self, const Parameter* p){
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)self->obj;
	
	
	Node* node = NULL;
	int i = 0;
	for (; i < Tree_node_count(tlk->tree); i++) {
		node = Tree_node(tlk->tree, i);
		if (strcmp(node->distance->name, Parameter_name(p)) == 0) {
			break;
		}
	}
	
	double* pattern_likelihoods = tlk->pattern_lk + tlk->sp->count;
	
	if(tlk->update){
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
					SingleTreeLikelihood_update_all_nodes(tlk);

		double logP = tlk->calculate(tlk); // make sure it is updated
		if (isnan(logP) || isinf(logP)) {
			return logP;
		}
		calculate_upper(tlk, Tree_root(tlk->tree));
		
		for (int i = 0; i < tlk->sp->count; i++) {
			pattern_likelihoods[i] = exp(tlk->pattern_lk[i]);
		}
	}
	if(Tree_node_count(tlk->tree) == i){
		return Model_first_derivative(self, p, 0.001);
	}
	double* pattern_dlikelihoods = tlk->pattern_lk + 2*tlk->sp->count;
	calculate_dldt_uppper(tlk, node, pattern_dlikelihoods);
	
	double dlogP = dlnldt_uppper(tlk, node, pattern_likelihoods, pattern_dlikelihoods);
	if(isnan(dlogP)){
		SingleTreeLikelihood_update_all_nodes(tlk);
//		Tree_print_parameters(tlk->tree);
//		printf("%f\n", tlk->calculate(tlk));
//		printf("%f\n", Model_first_derivative(self, p, 0.001));
//		printf("%s %f\n", p->name, p->value);
//		exit(10);
	}
	
	return dlogP;
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
	return clone;
}

static void _treeLikelihood_model_get_free_parameters(Model* model, Parameters* parameters){
	Model** list = (Model**)model->data;
	Model* mtree = list[0];
	Model* msm = list[1];
	Model* mbm = list[2];
	
	mtree->get_free_parameters(mtree, parameters);
	msm->get_free_parameters(msm, parameters);
	if(mbm != NULL){
		mbm->get_free_parameters(mbm, parameters);
	}
}


// TreeLikelihood listen to the TreeModel, SiteModel, BranchModel
Model * new_TreeLikelihoodModel( const char* name, SingleTreeLikelihood *tlk,  Model *tree, Model *sm, Model *bm ){
	Model *model = new_Model("treelikelihood",name, tlk);

	tree->listeners->add( tree->listeners, model );
	if(bm != NULL)bm->listeners->add( bm->listeners, model );
	sm->listeners->add( sm->listeners, model );

	model->logP = _singleTreeLikelihood_logP;
	model->dlogP = _singleTreeLikelihood_dlogP;
	model->update = _treelikelihood_handle_change;
	model->free = _treeLikelihood_model_free;
	model->clone = _treeLikelihood_model_clone;
	model->get_free_parameters = _treeLikelihood_model_get_free_parameters;
	model->data = (Model**)malloc(sizeof(Model*)*3);
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
	json_node* patterns_node = get_json_node(node, "sitepattern");
	json_node* tree_node = get_json_node(node, "tree");
	json_node* sm_node = get_json_node(node, "sitemodel");
	json_node* bm_node = get_json_node(node, "branchmodel");
	
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
//			sm = new_SiteModel_from_json(sm_node, hash);
		}
		bm = mbm->obj;
	}
	
	SingleTreeLikelihood* tlk = new_SingleTreeLikelihood((Tree*)mtree->obj, (SiteModel*)msm->obj, patterns, bm);
	char* id = get_json_node_value_string(node, "id");
	Model* model = new_TreeLikelihoodModel(id, tlk, mtree, msm, mbm);
	mtree->free(mtree);
	msm->free(msm);
	if(mbm != NULL) mbm->free(mbm);
	
	return model;
}

#pragma mark -
// MARK: SingleTreeLikelihood


SingleTreeLikelihood * new_SingleTreeLikelihood( Tree *tree, SiteModel *sm, SitePattern *sp, BranchModel *bm ){
	SingleTreeLikelihood *tlk = (SingleTreeLikelihood *)malloc( sizeof(SingleTreeLikelihood));
	assert(tlk);
	
	tlk->tree = tree;
	tlk->sm = sm;
	tlk->sp = sp;
	tlk->bm = bm;
	
	tlk->pattern_count = tlk->sp->count;
	tlk->cat_count = tlk->sm->cat_count;
	
	Node **nodes = Tree_get_nodes( tree, POSTORDER );
	
	tlk->partials_size = sp->count * sm->nstate * sm->cat_count;
	
	
	tlk->matrix_dim = Tree_node_count(tree);
	tlk->matrix_size = sm->nstate * sm->nstate;
	
	tlk->upper_partial_indexes = ivector(Tree_node_count(tree));
	
	tlk->root_partials_size = sp->count*sm->nstate*3;
	tlk->pattern_lk_size = sp->count*3;
	
	tlk->partials_dim = Tree_node_count(tree)*2; // allocate *2 for upper likelihoods
	
    // odd number of state
//    if( sm->nstate & 1 ){
//        tlk->partials_size += sp->count * sm->cat_count;
//        tlk->matrix_size   += sm->nstate;
//    }
    
	tlk->partials = (double**)malloc( tlk->partials_dim*sizeof(double*) );
    assert(tlk->partials);
	int i = 0;
	for ( ; i < Tree_node_count(tree); i++ ) {
		if( Node_isleaf( nodes[i] ) ){
			tlk->partials[Node_id( nodes[i] )] = NULL;
		}
		else {
			tlk->partials[Node_id( nodes[i] )] = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
            assert(tlk->partials[Node_id( nodes[i] )]);
			memset(tlk->partials[Node_id( nodes[i] )], 0.0, tlk->partials_size * sizeof(double));
		}
	}
	for ( ; i < tlk->partials_dim; i++ ) {
		tlk->partials[i] = (double*)aligned16_malloc( tlk->partials_size * sizeof(double) );
		assert(tlk->partials[i]);
		memset(tlk->partials[i], 0.0, tlk->partials_size * sizeof(double));
	}
	
	tlk->matrices = (double**)malloc( tlk->matrix_dim*sizeof(double*) );
    assert(tlk->matrices);
    int mat_len = tlk->matrix_size*sm->cat_count;
	for ( int i = 0; i < tlk->matrix_dim; i++ ) {
		tlk->matrices[i] = aligned16_malloc( mat_len * sizeof(double) );
        assert(tlk->matrices[i]);
		memset(tlk->matrices[i], 0.0, mat_len * sizeof(double));
	}
	
	tlk->root_partials = aligned16_malloc( tlk->root_partials_size * sizeof(double) );
    assert(tlk->root_partials);
	
	tlk->update_nodes = bvector(Tree_node_count(tree));
	for (int i = 0; i < Tree_node_count(tree); i++){
		tlk->update_nodes[i] = true;
	}
	tlk->update = true;
	
	tlk->pattern_lk = dvector(tlk->pattern_lk_size);
	tlk->lk = 0.;
	
	tlk->calculate = _calculate;
	tlk->calculate_upper = _calculate_uppper;
	
	tlk->update_partials      = update_partials_general;
	tlk->integrate_partials   = integrate_partials_general;
	tlk->node_log_likelihoods = node_log_likelihoods_general;
	tlk->calculate_branch_likelihood = calculate_branch_likelihood;
	tlk->update_partials_upper = update_partials_upper_general;

	if ( sm->nstate == 4 ) {
		tlk->update_partials      = update_partials_4;
		tlk->integrate_partials   = integrate_partials_4;
		tlk->node_log_likelihoods = node_log_likelihoods_4;
		tlk->calculate_branch_likelihood = calculate_branch_likelihood_4;
		tlk->update_partials_upper = update_partials_upper_4;
		tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_4;
	}
//	else if( sm->nstate == 20 ){
//		tlk->update_partials      = update_partials_20;
//		tlk->integrate_partials   = integrate_partials_20;
//		tlk->node_log_likelihoods = node_log_likelihoods_20;
//	}
	else if( sm->nstate >= 60 ){
		tlk->update_partials      = update_partials_codon;
		tlk->integrate_partials   = integrate_partials_codon;
		tlk->node_log_likelihoods = node_log_likelihoods_codon;
		tlk->update_partials_upper = update_partials_upper_codon;
		tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_codon;
	}
    
	tlk->mapping = ivector(Tree_node_count(tree));
	
	// map node names to sequence names
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
		if( Node_isleaf( nodes[i] ) ){
			tlk->mapping[Node_id(nodes[i])] = get_sequence_index(tlk->sp, nodes[i]->name);
		}
		else tlk->mapping[Node_id( nodes[i] )] = -1;
	}
	
	tlk->scale = false;
	tlk->scaling_factors = NULL;
	tlk->scaling_threshold = 1.E-40;
	
	OptConfig_init(tlk);
    
    
	tlk->hessian = NULL;
    tlk->approx = TREELIKELIHOOD_APPROXIMATION_NONE;
    tlk->hessian_length = 0;
    tlk->lnl_bl = INFINITY;
	
	tlk->node_id = -1;
    
    // Upper likelihood calculation
    // Variables are instanciated when we need them using SingleTreeLikelihood_use_upper
    tlk->use_upper = false;
    tlk->update_partials_upper = NULL;
	tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_general;
    tlk->node_upper = NULL;
    
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
	tlk->use_SIMD = false;
#endif
	
#ifdef SSE3_ENABLED
	if ( sm->nstate == 4 ) {
		tlk->update_partials      = update_partials_4_SSE;
		tlk->integrate_partials   = integrate_partials_4_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_4_SSE;
		tlk->calculate_branch_likelihood = calculate_branch_likelihood_4_SSE;
		tlk->update_partials_upper = update_partials_upper_sse_4;
		tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_sse_4;
        tlk->use_SIMD = true;
	}
	else if ( sm->nstate == 20 ) {
		tlk->update_partials      = update_partials_20_SSE;
		tlk->integrate_partials   = integrate_partials_general;
		tlk->node_log_likelihoods = node_log_likelihoods_general;
		tlk->update_partials_upper = update_partials_upper_sse_20;
		tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_sse_20;
        tlk->use_SIMD = true;
	}
    // Does not work for matrix with odd dimension.
    // If we want to load doubles at indexes 1 and 2 from matrices it will crash
    // same issue with partials, it will try to store at odd indexes
	else if( sm->nstate >= 60 && !sm->nstate & 1 ){
		tlk->update_partials      = update_partials_codon_SSE;
		tlk->integrate_partials   = integrate_partials_codon_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_codon_SSE;
		tlk->update_partials_upper = update_partials_upper_sse_codon;
		tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_codon;
        tlk->use_SIMD = true;
	}
    else if(!(sm->nstate & 1)){
        tlk->update_partials      = update_partials_general_even_SSE;
        tlk->integrate_partials   = integrate_partials_general_even_SSE;
		tlk->node_log_likelihoods = node_log_likelihoods_general_even_SSE;
		tlk->update_partials_upper = update_partials_upper_sse_codon;
		tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_codon;
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
	return tlk;
}

void free_SingleTreeLikelihood_internals( SingleTreeLikelihood *tlk ){
	if(tlk->scaling_factors != NULL ) free_dmatrix( tlk->scaling_factors, Tree_node_count(tlk->tree));
	
	for ( int i = 0; i < tlk->partials_dim; i++ ) {
		if( tlk->partials[i] != NULL ){
			free(tlk->partials[i]);
		}
	}
	free(tlk->partials);
	
	free_dmatrix(tlk->matrices, tlk->matrix_dim);
	
	free(tlk->mapping);
	free(tlk->update_nodes);
	free(tlk->pattern_lk);
	free(tlk->root_partials);
	
	if ( tlk->hessian != NULL ) {
		free(tlk->hessian);
	}
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

void free_SingleTreeLikelihood_share( SingleTreeLikelihood *tlk, bool shared_sitepattern, bool shared_sitemodel ){
	
	if(tlk->scaling_factors != NULL ) free_dmatrix( tlk->scaling_factors, tlk->partials_dim);
	
	for ( int i = 0; i < tlk->partials_dim; i++ ) {
		if( tlk->partials[i] != NULL ){
			free(tlk->partials[i]);
		}
	}
	free(tlk->partials);
	
	free_dmatrix(tlk->matrices, tlk->matrix_dim);
	
	free(tlk->mapping);
	free(tlk->update_nodes);
	free(tlk->pattern_lk);
	free(tlk->root_partials);
    
	if ( tlk->bm != NULL )     tlk->bm->free( tlk->bm, false );
	if ( tlk->tree != NULL )   free_Tree( tlk->tree );
	if ( !shared_sitemodel && tlk->sm != NULL)   free_SiteModel( tlk->sm );
	if ( !shared_sitepattern && tlk->sp != NULL) free_SitePattern( tlk->sp );
	
	if ( tlk->hessian != NULL ) {
		free(tlk->hessian);
	}
	
	free(tlk->upper_partial_indexes);
    if(tlk->root_frequencies != NULL)free(tlk->root_frequencies);
	free(tlk);
}

void free_SingleTreeLikelihood_share2( SingleTreeLikelihood *tlk, bool shared_tree, bool shared_sitemodel, bool shared_sitepattern, bool shared_branchmodel ){
	
	if(tlk->scaling_factors != NULL ) free_dmatrix( tlk->scaling_factors, tlk->partials_dim);
	
	for ( int i = 0; i < tlk->partials_dim; i++ ) {
		if( tlk->partials[i] != NULL ){
			free(tlk->partials[i]);
		}
	}
	free(tlk->partials);

	free_dmatrix(tlk->matrices, tlk->matrix_dim);
	
	free(tlk->mapping);
	free(tlk->update_nodes);
	free(tlk->pattern_lk);
	free(tlk->root_partials);
	
	free(tlk->upper_partial_indexes);

	if ( tlk->bm != NULL && !shared_branchmodel ){
		tlk->bm->free( tlk->bm, shared_tree );
	}
	else {
		if( !shared_tree && tlk->tree != NULL) free_Tree( tlk->tree );
	}

	if ( !shared_sitemodel && tlk->sm != NULL) free_SiteModel( tlk->sm );
	if ( !shared_sitepattern && tlk->sp != NULL) free_SitePattern( tlk->sp );
	
	if ( tlk->hessian != NULL ) {
		free(tlk->hessian);
	}
    if(tlk->root_frequencies != NULL)free(tlk->root_frequencies);
	free(tlk);
}

SingleTreeLikelihood * clone_SingleTreeLikelihood( SingleTreeLikelihood *tlk ){
	SitePattern *sp = clone_SitePattern(tlk->sp);
	SiteModel *sm   = clone_SiteModel(tlk->sm);
	Tree* tree = clone_Tree(tlk->tree);
	BranchModel* bm = NULL;
	if(tlk->bm != NULL){
		bm = clone_BranchModel(tlk->bm, tree);
	}
	return clone_SingleTreeLikelihood_with( tlk, tree, sm, sp, bm  );
}

SingleTreeLikelihood * clone_SingleTreeLikelihood_share( SingleTreeLikelihood *tlk, bool share_sitepattern, bool share_sitemodel ){
	SitePattern *sp = NULL;
	SiteModel *sm = NULL;
	if( share_sitepattern ) sp = tlk->sp;
	else sp = clone_SitePattern(tlk->sp);
	
	if( share_sitemodel ) sm = tlk->sm;
	else sm = clone_SiteModel(tlk->sm);
	Tree* tree = clone_Tree(tlk->tree);
	BranchModel* bm = NULL;
	if(tlk->bm != NULL){
		bm = clone_BranchModel(tlk->bm, tree);
	}
	return clone_SingleTreeLikelihood_with( tlk, tree, sm, sp, bm );
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
	
	newtlk->id = tlk->id;
	
	newtlk->mapping = clone_ivector(tlk->mapping, Tree_node_count(tlk->tree));
		
	newtlk->partials_size = tlk->partials_size;
	newtlk->pattern_lk_size = tlk->pattern_lk_size;
	newtlk->root_partials_size = tlk->root_partials_size;
	
	newtlk->scaling_factors = NULL;
	if ( tlk->scaling_factors != NULL ){
		newtlk->scaling_factors = clone_dmatrix( tlk->scaling_factors, tlk->partials_dim, tlk->sp->count );
	}
	
	newtlk->scale = tlk->scale;
	newtlk->scaling_threshold = tlk->scaling_threshold;
	
	newtlk->matrix_size = tlk->matrix_size;
	
	newtlk->update_nodes = clone_bvector( tlk->update_nodes, Tree_node_count(newtlk->tree));
	newtlk->update = tlk->update;
	
	
	newtlk->pattern_lk = clone_dvector( tlk->pattern_lk, tlk->pattern_lk_size );
	newtlk->lk = tlk->lk;
	
	newtlk->calculate            = tlk->calculate;
	newtlk->update_partials      = tlk->update_partials;
	newtlk->integrate_partials   = tlk->integrate_partials;
	newtlk->node_log_likelihoods = tlk->node_log_likelihoods;
	newtlk->calculate_branch_likelihood = tlk->calculate_branch_likelihood;
	
	OptConfig_copy( &tlk->opt, &newtlk->opt );
    
    
	newtlk->hessian = NULL;
    newtlk->hessian_length = 0;
    newtlk->lnl_bl = INFINITY;
    newtlk->approx = TREELIKELIHOOD_APPROXIMATION_NONE;
	
	if ( tlk->hessian != NULL ) {
		newtlk->hessian = clone_dvector(tlk->hessian, tlk->hessian_length );
        newtlk->hessian_length = tlk->hessian_length;
        newtlk->lnl_bl = tlk->lnl_bl;
        newtlk->approx = tlk->approx;
	}
	
	newtlk->upper_partial_indexes = clone_ivector(tlk->upper_partial_indexes, Tree_node_count(tlk->tree));
	
	newtlk->node_id = tlk->node_id;
    
    newtlk->nthreads = tlk->nthreads;
    
	newtlk->use_upper = false;
    newtlk->calculate_upper = NULL;
    newtlk->update_partials_upper = NULL;
    newtlk->node_log_likelihoods_upper = NULL;
    newtlk->node_upper = NULL;
	
	newtlk->matrix_dim = tlk->matrix_dim;
	newtlk->partials_dim = tlk->partials_dim;
	
	Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
	
	newtlk->use_upper = tlk->use_upper;
	newtlk->calculate_upper = tlk->calculate_upper;
	newtlk->update_partials_upper = tlk->update_partials_upper;
	newtlk->node_log_likelihoods_upper = tlk->node_log_likelihoods_upper;
	
	newtlk->node_upper = NULL;
	
    
    // partials and matrices need to be aligned when compiled with GCC >= 4.7
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
    newtlk->use_SIMD = tlk->use_SIMD;
#endif
	
	newtlk->partials = (double**)malloc( tlk->partials_dim*sizeof(double*) );
    assert(newtlk->partials);
	int i = 0;
	for ( ; i < Tree_node_count(tlk->tree); i++ ) {
		if( Node_isleaf(nodes[i]) ){
			newtlk->partials[Node_id( nodes[i] )] = NULL;
		}
		else {
			newtlk->partials[Node_id( nodes[i] )] = aligned16_malloc( tlk->partials_size * sizeof(double) );
            assert(newtlk->partials[Node_id( nodes[i] )]);
			memcpy(newtlk->partials[Node_id( nodes[i] )], tlk->partials[Node_id( nodes[i] )], tlk->partials_size * sizeof(double));
		}
	}
	for ( ; i < tlk->partials_dim; i++ ) {
		newtlk->partials[i] = aligned16_malloc( tlk->partials_size * sizeof(double) );
		assert(newtlk->partials[i]);
		memcpy(newtlk->partials[i], tlk->partials[i], tlk->partials_size * sizeof(double));
	}
	
	newtlk->matrices = (double**)malloc( tlk->matrix_dim*sizeof(double*) );
    assert(newtlk->matrices);
	for ( int i = 0; i < tlk->matrix_dim; i++ ) {
		newtlk->matrices[i] = aligned16_malloc( tlk->matrix_size*tlk->sm->cat_count * sizeof(double) );
        assert(newtlk->matrices[i]);
		memcpy(newtlk->matrices[i], tlk->matrices[i], tlk->matrix_size*tlk->sm->cat_count * sizeof(double));
	}
	
	newtlk->root_partials = (double*)malloc( newtlk->root_partials_size*sizeof(double) );
    assert(newtlk->root_partials);
	memcpy(newtlk->root_partials, tlk->root_partials, newtlk->root_partials_size * sizeof(double));
    
	newtlk->sp->ref_count++;
	
	if(tlk->root_frequencies != NULL){
		newtlk->root_frequencies = clone_dvector(tlk->root_frequencies, tlk->sm->nstate);
	}
	newtlk->get_root_frequencies = tlk->get_root_frequencies;
	return newtlk;
}

void SingleTreeLikelihood_set_BranchModel( SingleTreeLikelihood *stlk, BranchModel *bm, bool removeTree ){
	SingleTreeLikelihood_remove_BranchModel(stlk, removeTree);
	
	stlk->bm = bm;
	if ( removeTree ){
		stlk->tree = bm->tree;
		Node **nodes = Tree_get_nodes(bm->tree, POSTORDER);
		for ( int i = 0; i < Tree_node_count(bm->tree); i++ ) {
			if( Node_isleaf(nodes[i]) ){
                stlk->mapping[Node_id(nodes[i])] = get_sequence_index(stlk->sp, nodes[i]->name);
			}
			else stlk->mapping[Node_id( nodes[i] )] = -1;
		}
	}
	SingleTreeLikelihood_update_all_nodes( stlk );

    stlk->node_upper = NULL;
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
    if( stlk->use_SIMD ){
        stlk->calculate_upper = _calculate_uppper_height_sse;
    }
    else {
        stlk->calculate_upper = _calculate_uppper_height;
    }
#else
    stlk->calculate_upper = _calculate_uppper_height;
#endif
}

void SingleTreeLikelihood_remove_BranchModel( SingleTreeLikelihood *stlk, bool removeTree ){
	if ( stlk->bm != NULL ) {
		stlk->bm->free(stlk->bm, removeTree);
		stlk->bm = NULL;
	}
    stlk->calculate_upper = _calculate_uppper;
}


void SingleTreeLikelihood_set_Tree( SingleTreeLikelihood *stlk, Tree *tree ){
	SingleTreeLikelihood_remove_Tree(stlk);
	
	stlk->tree = tree;
	stlk->bm->tree = tree;
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if( nodes[i]->left == NULL ){
            stlk->mapping[Node_id(nodes[i])] = get_sequence_index(stlk->sp, nodes[i]->name);
		}
		else stlk->mapping[Node_id( nodes[i] )] = -1;
	}
	SingleTreeLikelihood_update_all_nodes( stlk );
	
}

void SingleTreeLikelihood_remove_Tree( SingleTreeLikelihood *stlk ){
	if ( stlk->tree != NULL ) {
		stlk->bm->tree = NULL;
		free_Tree(stlk->tree);
		stlk->tree = NULL;
	}
}


void SingleTreeLikelihood_use_rescaling( SingleTreeLikelihood *tlk, bool use ){
	tlk->scale = use;
	if ( tlk->scaling_factors == NULL && use ) {
		tlk->scaling_factors = dmatrix(tlk->partials_dim, tlk->sp->count );
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
        if ( !Node_isleaf(nodes[i]) && tlk->partials[i] == NULL ) {
            for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
                if( Node_isleaf(nodes[j]) && tlk->partials[j] != NULL ){
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

// partials:    number of patterns * number of categpries * number of states
// frequencies: number of states
// posteriors:  [number of patterns x number of categories]
double ** SingleTreeLikelihood_posterior_sites( SingleTreeLikelihood *tlk ){
	
	if ( tlk->update ) {
		_calculate(tlk);
	}
	const double *partials = tlk->partials[ Tree_node_count(tlk->tree)-1];
	const double *frequencies = tlk->get_root_frequencies(tlk);
	const double *root_lk   = tlk->pattern_lk; // P[D_i]
	const double *rate_prop = tlk->sm->get_proportions(tlk->sm); // P[R_j]
	
	double **posteriors = dmatrix(tlk->sp->count, tlk->sm->cat_count);
	int v = 0;
	int i = 0;
	
	// P[ D_i| r_i=r_j ]
	for ( int j = 0; j < tlk->sm->cat_count; j++) {
		for ( i = 0; i < tlk->sp->count; i++) {
			// sum pi_x P[ D_i| r_i=r_j, D_i^r=x ]
			for ( int x = 0; x < tlk->sm->nstate; x++ ) {
				posteriors[i][j] += frequencies[x] * partials[v];
				v++;
			}
			// P[ r_i=r_j| D_i ] = P[ D_i| r_i=r_j ] * P[ r_i=r_j ] / P[D_i]
			posteriors[i][j] = posteriors[i][j] *rate_prop[j]/exp(root_lk[i]);
		}
	}
	
	return posteriors;
}


double _calculate( SingleTreeLikelihood *tlk ){
	if (tlk->use_upper) {
		tlk->use_upper = false;
		// In brent the initial point is calculated without changing the parameter
		if(tlk->update){
			for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
				if (tlk->update_nodes[i]) {
					double lk = tlk->calculate_upper(tlk, Tree_node(tlk->tree, i));
					tlk->node_upper = Tree_node(tlk->tree, i);
					tlk->use_upper = true;
//					printf("+ %d %f\n", i, lk);
					tlk->update_nodes[i] = false;
					return lk;
				}
			}
		}
		
		tlk->calculate(tlk);
		calculate_upper(tlk, Tree_root(tlk->tree));
		tlk->use_upper = true;
//		printf("== %f\n", tlk->lk);
		return tlk->lk;
	}
	
	if( !tlk->update ){
		return tlk->lk;
	}
	
	//tlk->update_partials = update_partials_noexp_integrate_4;
	//_calculate_partials_noexp_integrate( tlk, Tree_root(tlk->tree) );
	_calculate_partials( tlk, Tree_root(tlk->tree) );
	
	tlk->integrate_partials(tlk, tlk->partials[Node_id(Tree_root(tlk->tree))], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
	
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
		
		tlk->integrate_partials(tlk, tlk->partials[Node_id(Tree_root(tlk->tree))], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
		
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

bool _calculate_partials( SingleTreeLikelihood *tlk, Node *n  ){
	bool updated = false;
	
	if( tlk->update_nodes[ Node_id(n) ] && !Node_isroot(n) ){
		// a bit hackish here, for local optimization of clades
		if ( tlk->node_id == -1 || Node_id(n) != tlk->node_id ) {
			
			double bl = 0;
			if( tlk->bm == NULL ){
				bl = Node_distance(n);
			}
			else{
				bl = tlk->bm->get(tlk->bm, n) * Node_time_elapsed(n);
				
				if(bl < 0 ){
#ifdef TIMETEST
					fprintf(stderr,"%f right: %d diff %E\n", Node_t(n), (Node_right(Node_parent(n))==n), (Node_height(Node_parent(n))-Node_height(n)));
					
#endif
					fprintf(stderr, "calculate_partials: %s branch length = %E rate = %f height = %f - parent height [%s]= %f (%f)\n", n->name, bl, tlk->bm->get(tlk->bm, n), Node_height(n), n->parent->name, Node_height(Node_parent(n)), Node_distance(n));
					exit(1);
				}
			}
			
			for (int i = 0; i < tlk->sm->cat_count; i++) {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
				if( tlk->use_SIMD ){
					if( Node_isleaf(n) ){
						tlk->sm->m->p_t_transpose(tlk->sm->m,
												  bl * tlk->sm->get_rate(tlk->sm, i),
												  &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
					}
					else {
						tlk->sm->m->p_t(tlk->sm->m,
										bl * tlk->sm->get_rate(tlk->sm, i),
										&tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
					}
				}
				else{
					tlk->sm->m->p_t(tlk->sm->m,
									bl * tlk->sm->get_rate(tlk->sm, i),
									&tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
				}
#else
				tlk->sm->m->p_t(tlk->sm->m,
								bl * tlk->sm->get_rate(tlk->sm, i),
								&tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
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
			
			tlk->update_partials( tlk, Node_id(n), indx_child1, indx_child1, indx_child2, indx_child2);
			
			// a bit hackish here, for local optimization of clades
			if( Node_id(n) == tlk->node_id ){
				tlk->integrate_partials(tlk, tlk->partials[Node_id(n)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
				
				tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
				
			}
			updated = true;
		}
	}
	return updated;
}

bool _calculate_partials_noexp_integrate( SingleTreeLikelihood *tlk, Node *n  ){
	bool updated = false;
	
	if( tlk->update_nodes[ Node_id(n) ] && !Node_isroot(n) ){
		// a bit hackish here, for local optimization of clades
		if ( tlk->node_id == -1 || Node_id(n) != tlk->node_id ) {
			
			double bl = 0;
			if( tlk->bm == NULL ){
				bl = Node_distance(n);
			}
			else{
				bl = tlk->bm->get(tlk->bm, n) * (Node_height(Node_parent(n)) - Node_height(n) );
				
				if(bl < 0 ){
					fprintf(stderr, "calculate_partials: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", n->name, bl, tlk->bm->get(tlk->bm, n), Node_height(n), n->parent->name, Node_height(Node_parent(n)), Node_distance(n));
					exit(1);
				}
			}
			
			if( tlk->sm->m->need_update ){
				tlk->sm->m->update_Q(tlk->sm->m);
			}
			
			// precompute exp(lambda_i*bl). dim = cat_count*nstate
			// it would be better to precompute u_{x,i}exp(lambda_i*bl) instead but dimension would be  dim = cat_count*nstate *nstate
			for ( int c = 0; c < tlk->sm->cat_count; c++ ) {
				for ( int i = 0; i < tlk->sm->nstate ; i++ ) {
					tlk->matrices[Node_id(n)][c*tlk->sm->nstate+i] = exp(tlk->sm->m->eigendcmp->eval[i] * bl * tlk->sm->get_rate(tlk->sm, c) );
				}
			}
			
			updated = true;
		}
	}
	
	if( !Node_isleaf(n) ){
		bool update_child1 = _calculate_partials_noexp_integrate( tlk, Node_left(n) );
		bool update_child2 = _calculate_partials_noexp_integrate( tlk, Node_right(n) );
		
		if( update_child1 || update_child2 ){
			int indx_child1 = Node_id(Node_left(n));
			int indx_child2 = Node_id(Node_right(n));
			
			tlk->update_partials( tlk, Node_id(n), indx_child1, indx_child1, indx_child2, indx_child2);
			
			// a bit hackish here, for local optimization of clades
			if( Node_isroot(n) || Node_id(n) == tlk->node_id ){
				tlk->integrate_partials(tlk, tlk->partials[Node_id(n)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
				
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
	
	tlk->node_upper = NULL; // does not hurt to set it even if we don't use upper likelihoods
}

// update 1 node
void SingleTreeLikelihood_update_one_node( SingleTreeLikelihood *tlk, const Node *node ){
	tlk->update_nodes[Node_id(node)] = true;
	tlk->update = true;
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
	
	for ( int i = 0; i < tlk->sp->count; i++ ) {
		
		double scaleFactor = 0.0;
		int v = u;
		for ( k = 0; k < tlk->sm->cat_count; k++ ) {
			for ( j = 0; j < tlk->sm->nstate; j++ ) {
				if ( tlk->partials[nodeIndex][v] > scaleFactor ) {
					scaleFactor = tlk->partials[nodeIndex][v];
				}
				v++;
			}
			v += ( tlk->sp->count - 1 ) * tlk->sm->nstate;
		}
		
		if ( scaleFactor < tlk->scaling_threshold ) {
			
			v = u;
			for ( k = 0; k < tlk->sm->cat_count; k++ ) {
				for ( j = 0; j < tlk->sm->nstate; j++ ) {
					tlk->partials[nodeIndex][v] /= scaleFactor;
					v++;
				}
				v += (tlk->sp->count - 1) * tlk->sm->nstate;
			}
			tlk->scaling_factors[nodeIndex][i] = log(scaleFactor);
		}
		else {
			tlk->scaling_factors[nodeIndex][i] = 0.0;
		}
		u += tlk->sm->nstate;
	}
}

double getLogScalingFactor( const SingleTreeLikelihood *tlk, int pattern ) {
	double log_scale_factor = 0.0;
	if ( tlk->scale ) {
		for (int i = 0; i < Tree_node_count(tlk->tree); i++) {
			log_scale_factor += tlk->scaling_factors[i][pattern];
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
			tlk->integrate_partials   = integrate_partials_4_SSE;
			tlk->node_log_likelihoods = node_log_likelihoods_4_SSE;
			tlk->calculate_branch_likelihood = calculate_branch_likelihood_4_SSE;
			
			tlk->update_partials_upper      = update_partials_upper_sse_4;
			tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_sse_4;
		}
		else if ( tlk->sm->nstate == 20 ) {
			tlk->update_partials      = update_partials_20_SSE;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			
			tlk->update_partials_upper      = update_partials_upper_sse_20;
			tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_sse_20;
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
			//            tlk->node_log_likelihoods = node_log_likelihoods_general;
			
			tlk->update_partials_upper      = update_partials_upper_sse_codon;
			tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_codon;
		}
		else {
			tlk->update_partials      = update_partials_general_even_SSE;
			tlk->integrate_partials   = integrate_partials_general_even_SSE;
			tlk->node_log_likelihoods = node_log_likelihoods_general_even_SSE;
			tlk->calculate_branch_likelihood = calculate_branch_likelihood;
			//TODO: what about the odd case
			// no upper
		}
	}
	else {
		if ( tlk->sm->nstate == 4 ) {
			tlk->update_partials      = update_partials_4;
			tlk->integrate_partials   = integrate_partials_4;
			tlk->node_log_likelihoods = node_log_likelihoods_4;
			tlk->calculate_branch_likelihood = calculate_branch_likelihood_4;
			
			tlk->update_partials_upper      = update_partials_upper_4;
			tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_4;
		}
		else if ( tlk->sm->nstate == 20 ) {
			tlk->update_partials      = update_partials_general;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			
			tlk->update_partials_upper      = update_partials_upper_20;
			tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_20;
		}
		else if( tlk->sm->nstate >= 60 ){
			tlk->update_partials      = update_partials_codon;
			tlk->integrate_partials   = integrate_partials_codon;
			tlk->node_log_likelihoods = node_log_likelihoods_codon;
			
			tlk->update_partials_upper      = update_partials_upper_codon;
			tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_codon;
		}
		else {
			tlk->update_partials      = update_partials_general;
			tlk->integrate_partials   = integrate_partials_general;
			tlk->node_log_likelihoods = node_log_likelihoods_general;
			tlk->calculate_branch_likelihood = calculate_branch_likelihood;
			
			// no upper??
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
	
	
	tlk->integrate_partials(tlk, tlk->partials[Node_id(node)], tlk->sm->get_proportions(tlk->sm), partials );
	
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
    
    tlk->integrate_partials(tlk, tlk->partials[Node_id(Tree_root(tlk->tree))], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
    
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
		
        tlk->integrate_partials(tlk, tlk->partials[Node_id(Tree_root(tlk->tree))], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
        
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

void TreeLikelihood_set_calculate( SingleTreeLikelihood *tlk, Node *node ){
	if ( node == NULL ) {
		tlk->calculate = _calculate;
        tlk->node_id = -1;
	}
	else {
		tlk->calculate = _calculate2;
        tlk->node_id = Node_id(node);
	}
	
	SingleTreeLikelihood_update_all_nodes(tlk);
}


#pragma mark -
#pragma mark Analytical derivative

// update partials only, does not update matrices
static bool _calculate_partials_noP( SingleTreeLikelihood *tlk, Node *n  ){
	bool updated = false;
	
	if( tlk->update_nodes[ Node_id(n) ] && !Node_isroot(n) ){
			updated = true;
	}
	
	if( !Node_isleaf(n) ){
		bool update_child1 = _calculate_partials_noP( tlk, Node_left(n) );
		bool update_child2 = _calculate_partials_noP( tlk, Node_right(n) );
		
		if( update_child1 || update_child2 ){
			int indx_child1 = Node_id(Node_left(n));
			int indx_child2 = Node_id(Node_right(n));
            
			tlk->update_partials( tlk, Node_id(n), indx_child1, indx_child1, indx_child2, indx_child2);
			updated = true;
		}
	}
	return updated;
}

// Calculate the first derivative with respect to each branch length parameter of the likleihood function
// df [#nodes][#sites]
// fist dimension is indexed by node postorder_idx
void calculate_all_df_dt( SingleTreeLikelihood *tlk, double **df ){
    
    Node **nodes = Tree_nodes(tlk->tree);
    
    double *pmats = dvector(tlk->sm->cat_count * tlk->matrix_size);
    
    for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
        Node *n = nodes[j];
        
        if ( Node_isroot(n) || (Node_isroot(Node_parent(n)) && Node_right(Node_parent(n)) == n )) {
            continue;
        }
        double bl = Node_distance(n);
        
        memcpy(pmats, tlk->matrices[Node_id(n)], tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        
        
        
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD ){
                if( Node_isleaf(n) ){
                    tlk->sm->m->dp_dt_transpose(tlk->sm->m,
                                                bl * tlk->sm->get_rate(tlk->sm, i),
                                                &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
                else {
                    tlk->sm->m->dp_dt(tlk->sm->m,
                                      bl * tlk->sm->get_rate(tlk->sm, i),
                                      &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
            }
            else{
                tlk->sm->m->dp_dt(tlk->sm->m,
                                  bl * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
            }
#else
            tlk->sm->m->dp_dt(tlk->sm->m,
                              bl * tlk->sm->get_rate(tlk->sm, i),
                              &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
#endif
            for (int k = 0; k < tlk->matrix_size; k++) {
                tlk->matrices[Node_id(n)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
            }
        }
       
        
        SingleTreeLikelihood_update_one_node(tlk, n);
        _calculate_partials_noP( tlk, Tree_root(tlk->tree) );
        tlk->integrate_partials(tlk, tlk->partials[Node_id(Tree_root(tlk->tree))], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
        int v = 0;
        int i = 0;
        
        const int nstate   = tlk->sm->nstate;
        const int patternCount = tlk->sp->count;
		const double* freqs = tlk->get_root_frequencies(tlk);
        
        for ( int k = 0; k < patternCount; k++ ) {
            
            tlk->pattern_lk[k] = 0;
            for ( i = 0; i < nstate; i++ ) {
                
                tlk->pattern_lk[k] += freqs[i] * tlk->partials[Node_id(Tree_root(tlk->tree))][v];
                v++;
            }
            //tlk->pattern_lk[k] = log(tlk->pattern_lk[k]);
            
            if ( tlk->scale ) {
                //printf("scaling\n");
                tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
            }
        }
        

        memcpy(df[Node_id(n)], tlk->pattern_lk, tlk->sp->count*sizeof(double));

        memcpy( tlk->matrices[Node_id(n)], pmats, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        SingleTreeLikelihood_update_one_node(tlk, n);
        
    }
    free(pmats);
}

// Calculate the first derivative with respect to each branch length parameter of the loglikelihood function
void calculate_all_dlnl_dt( SingleTreeLikelihood *tlk, double *dlnl ){
    
    Node **nodes = Tree_nodes(tlk->tree);
    
    double *pmats = dvector(tlk->sm->cat_count * tlk->matrix_size);
    double *lnls = dvector(tlk->sp->count);
    tlk->calculate(tlk);
    for ( int i = 0; i < tlk->sp->count; i++) {
        lnls[i] = exp(tlk->pattern_lk[i]);
	}
    
    
    for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
        Node *n = nodes[j];
        //if( !Parameter_estimate(n->distance) ) continue;

        if ( Node_isroot(n) || (Node_isroot(Node_parent(n)) && Node_right(Node_parent(n)) == n )) {
            continue;
        }
        
        double bl = Node_distance(n);
        
        memcpy(pmats, tlk->matrices[Node_id(n)], tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        //print_P(pmats, 4);
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD ){
                if( Node_isleaf(n) ){
                    tlk->sm->m->dp_dt_transpose(tlk->sm->m,
                                                bl * tlk->sm->get_rate(tlk->sm, i),
                                                &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
                else {
                    tlk->sm->m->dp_dt(tlk->sm->m,
                                      bl * tlk->sm->get_rate(tlk->sm, i),
                                      &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
            }
            else{
                tlk->sm->m->dp_dt(tlk->sm->m,
                                  bl * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
            }
#else
            tlk->sm->m->dp_dt(tlk->sm->m,
                              bl * tlk->sm->get_rate(tlk->sm, i),
                              &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
#endif
            for (int k = 0; k < tlk->matrix_size; k++) {
                tlk->matrices[Node_id(n)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
            }
        }
        
         //print_P(tlk->matrices[Node_id(n)], 4);
        SingleTreeLikelihood_update_one_node(tlk, n);
        _calculate_partials_noP( tlk, Tree_root(tlk->tree) );
		
		const int rootId = Node_id(Tree_root(tlk->tree));
        tlk->integrate_partials(tlk, tlk->partials[rootId], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
		//tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->sm->m->_freqs, tlk->pattern_lk);
		int v = 0;
		int i = 0;
		
		const double* freqs = tlk->get_root_frequencies(tlk);
		const int nstate   = tlk->sm->nstate;
		const int patternCount = tlk->sp->count;
		dlnl[Node_id(n)] = 0;
		
		for ( int k = 0; k < patternCount; k++ ) {
			
			tlk->pattern_lk[k] = 0;
			for ( i = 0; i < nstate; i++ ) {
				//TODO: this is changing the pattern_lk
				tlk->pattern_lk[k] += freqs[i] * tlk->root_partials[v];
				v++;
			}
			//tlk->pattern_lk[k] = log(tlk->pattern_lk[k]);
			
			if ( tlk->scale ) {
				//printf("scaling\n");
				tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
			}
			dlnl[Node_id(n)] += tlk->pattern_lk[k]/lnls[k] * tlk->sp->weights[k];
		}

        memcpy( tlk->matrices[Node_id(n)], pmats, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        SingleTreeLikelihood_update_one_node(tlk, n);
        
    }
    free(lnls);
    free(pmats);
}

// Calculate the second derivative with respect to each branch length parameter
// Highly inefficient
void calculate_all_d2lnl_d2t( SingleTreeLikelihood *tlk, double *d2lnl ){
    
    double **dfi = dmatrix(Tree_node_count(tlk->tree), tlk->sp->count);
    calculate_all_df_dt(tlk, dfi);
    
    Node **nodes = Tree_nodes(tlk->tree);
    double *pmats = dvector(tlk->sm->cat_count * tlk->matrix_size);
    double *lks = dvector(tlk->sp->count);
    
    SingleTreeLikelihood_update_all_nodes(tlk);
    tlk->calculate(tlk);
    for ( int i = 0; i < tlk->sp->count; i++) {
		lks[i] = exp(tlk->pattern_lk[i]);
	}
    
    for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
        Node *n = nodes[j];
        
        if ( Node_isroot(n) || (Node_isroot(Node_parent(n)) && Node_right(Node_parent(n)) == n )) {
            continue;
        }
        
        double bl = Node_distance(n);
        
        memcpy(pmats, tlk->matrices[Node_id(n)], tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD ){
                if( Node_isleaf(n) ){
                    tlk->sm->m->d2p_d2t_transpose(tlk->sm->m,
                                                bl * tlk->sm->get_rate(tlk->sm, i),
                                                &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
                else {
                    tlk->sm->m->d2p_d2t(tlk->sm->m,
                                      bl * tlk->sm->get_rate(tlk->sm, i),
                                      &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
            }
            else{
                tlk->sm->m->d2p_d2t(tlk->sm->m,
                                  bl * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
            }
#else
            tlk->sm->m->d2p_d2t(tlk->sm->m,
                              bl * tlk->sm->get_rate(tlk->sm, i),
                              &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
#endif
            
            for (int k = 0; k < tlk->matrix_size; k++) {
                tlk->matrices[Node_id(n)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i) * tlk->sm->get_rate(tlk->sm, i);
            }
            
        }
        
        SingleTreeLikelihood_update_one_node(tlk, n);
        _calculate_partials_noP( tlk, Tree_root(tlk->tree) );
        
        tlk->integrate_partials(tlk, tlk->partials[Node_id(Tree_root(tlk->tree))], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
        int v = 0;
        int i = 0;
        
        const int nstate   = tlk->sm->nstate;
		const int patternCount = tlk->sp->count;
		const double* freqs = tlk->get_root_frequencies(tlk);
		
        for ( int k = 0; k < patternCount; k++ ) {
            
            tlk->pattern_lk[k] = 0;
            for ( i = 0; i < nstate; i++ ) {
                
                tlk->pattern_lk[k] += freqs[i] * tlk->partials[Node_id(Tree_root(tlk->tree))][v];
                v++;
            }
            //tlk->pattern_lk[k] = log(tlk->pattern_lk[k]);
            
            if ( tlk->scale ) {
                //printf("scaling\n");
                tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
            }
        }
        
        for ( int i = 0; i < tlk->sp->count; i++) {
            d2lnl[Node_id(n)] += ((tlk->pattern_lk[i] * lks[i] - (dfi[Node_id(n)][i]*dfi[Node_id(n)][i])) / (lks[i]*lks[i])) * tlk->sp->weights[i];
        }
        
        memcpy( tlk->matrices[Node_id(n)], pmats, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        SingleTreeLikelihood_update_one_node(tlk, n);
        
    }
    free(lks);
    free(pmats);
    free_dmatrix(dfi, Tree_node_count(tlk->tree));
}

// Calculate the Hessian fo branches only
// length(d2lnl) == nodeCount*nodeCount
void calculate_hessian_branches( SingleTreeLikelihood *tlk, double *d2lnl ){
    
    int nodeCount = Tree_node_count(tlk->tree);
    double **dfi = dmatrix(nodeCount, tlk->sp->count);
    calculate_all_df_dt(tlk, dfi);
    Node **nodes = Tree_nodes(tlk->tree);
    double *pmats1 = dvector(tlk->sm->cat_count * tlk->matrix_size);
    double *pmats2 = dvector(tlk->sm->cat_count * tlk->matrix_size);
    double *lks = dvector(tlk->sp->count);
    
    SingleTreeLikelihood_update_all_nodes(tlk);
    tlk->calculate(tlk);
    for ( int i = 0; i < tlk->sp->count; i++) {
        lks[i] = exp(tlk->pattern_lk[i]);
    }
    
    Node *root = Tree_root(tlk->tree);
    Node *unwantedNode = Node_right(root);
    
    for ( int j = 0; j < nodeCount; j++ ) {
        Node *n1 = nodes[j];
        
        if ( n1 == root || n1 == unwantedNode ) {
            continue;
        }
        
        double bl1 = Node_distance(n1);
        
        memcpy(pmats1, tlk->matrices[Node_id(n1)], tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        
        // Diagonal
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD ){
                if( Node_isleaf(n1) ){
                    tlk->sm->m->d2p_d2t_transpose(tlk->sm->m,
                                                  bl1 * tlk->sm->get_rate(tlk->sm, i),
                                                  &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
                }
                else {
                    tlk->sm->m->d2p_d2t(tlk->sm->m,
                                        bl1 * tlk->sm->get_rate(tlk->sm, i),
                                        &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
                }
            }
            else{
                tlk->sm->m->d2p_d2t(tlk->sm->m,
                                    bl1 * tlk->sm->get_rate(tlk->sm, i),
                                    &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
            }
#else
            tlk->sm->m->d2p_d2t(tlk->sm->m,
                                bl1 * tlk->sm->get_rate(tlk->sm, i),
                                &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
#endif
            
            for (int k = 0; k < tlk->matrix_size; k++) {
                tlk->matrices[Node_id(n1)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i) * tlk->sm->get_rate(tlk->sm, i);
            }
            
        }
        Node* root = Tree_root(tlk->tree);
        SingleTreeLikelihood_update_one_node(tlk, n1);
        _calculate_partials( tlk, root );
        
        tlk->integrate_partials(tlk, tlk->partials[Node_id(root)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
        int v = 0;
        int i = 0;
        
        const int nstate   = tlk->sm->nstate;
		const int patternCount = tlk->sp->count;
		const double* freqs = tlk->get_root_frequencies(tlk);
		
        for ( int k = 0; k < patternCount; k++ ) {
            
            tlk->pattern_lk[k] = 0;
            for ( i = 0; i < nstate; i++ ) {
                
                tlk->pattern_lk[k] += freqs[i] * tlk->partials[Node_id(root)][v];
                v++;
            }
            //tlk->pattern_lk[k] = log(tlk->pattern_lk[k]);
            
            if ( tlk->scale ) {
                //printf("scaling\n");
                tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
            }
        }
        
        for ( int i = 0; i < tlk->sp->count; i++) {
            d2lnl[nodeCount*Node_id(n1)+Node_id(n1)] += ((tlk->pattern_lk[i] * lks[i] - (dfi[Node_id(n1)][i]*dfi[Node_id(n1)][i])) / (lks[i]*lks[i])) * tlk->sp->weights[i];
        }
        
        memcpy( tlk->matrices[Node_id(n1)], pmats1, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
        SingleTreeLikelihood_update_one_node(tlk, n1);
        
        
        // Off diagonal
        for ( int k = j+1; k < nodeCount; k++ ) {
            Node *n2 = nodes[k];
            
            if ( n2 == root || n2 == unwantedNode ) {
                continue;
            }
            
            double bl2 = Node_distance(n2);
            
            memcpy(pmats2, tlk->matrices[Node_id(n2)], tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
            
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                Node *n2 = nodes[i];
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
                if( tlk->use_SIMD ){
                    if( Node_isleaf(n1) ){
                        tlk->sm->m->dp_dt_transpose(tlk->sm->m,
                                                    bl1 * tlk->sm->get_rate(tlk->sm, i),
                                                    &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
                    }
                    else {
                        tlk->sm->m->dp_dt(tlk->sm->m,
                                          bl1 * tlk->sm->get_rate(tlk->sm, i),
                                          &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
                    }
                    if( Node_isleaf(n2) ){
                        tlk->sm->m->dp_dt_transpose(tlk->sm->m,
                                                    bl2 * tlk->sm->get_rate(tlk->sm, i),
                                                    &tlk->matrices[Node_id(n2)][i*tlk->matrix_size]);
                    }
                    else {
                        tlk->sm->m->dp_dt(tlk->sm->m,
                                          bl2 * tlk->sm->get_rate(tlk->sm, i),
                                          &tlk->matrices[Node_id(n2)][i*tlk->matrix_size]);
                    }
                }
                else{
                    tlk->sm->m->dp_dt(tlk->sm->m,
                                      bl1 * tlk->sm->get_rate(tlk->sm, i),
                                      &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
                    
                    tlk->sm->m->dp_dt(tlk->sm->m,
                                      bl2 * tlk->sm->get_rate(tlk->sm, i),
                                      &tlk->matrices[Node_id(n2)][i*tlk->matrix_size]);
                }
#else
                tlk->sm->m->dp_dt(tlk->sm->m,
                                  bl1 * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[Node_id(n1)][i*tlk->matrix_size]);
                
                tlk->sm->m->dp_dt(tlk->sm->m,
                                  bl2 * tlk->sm->get_rate(tlk->sm, i),
                                  &tlk->matrices[Node_id(n2)][i*tlk->matrix_size]);
#endif
                
                for (int k = 0; k < tlk->matrix_size; k++) {
                    tlk->matrices[Node_id(n1)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
                    tlk->matrices[Node_id(n2)][i*tlk->matrix_size+k] *= tlk->sm->get_rate(tlk->sm, i);
                }
                
            }
            
            SingleTreeLikelihood_update_one_node(tlk, n1);
            SingleTreeLikelihood_update_one_node(tlk, n2);
            
            Node* root = Tree_root(tlk->tree);
            
            _calculate_partials( tlk, root );
            
            tlk->integrate_partials(tlk, tlk->partials[Node_id(root)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
            int v = 0;
            int i = 0;
            
            const int nstate   = tlk->sm->nstate;
			const int patternCount = tlk->sp->count;
			const double* freqs = tlk->get_root_frequencies(tlk);
			
            for ( int k = 0; k < patternCount; k++ ) {
                
                tlk->pattern_lk[k] = 0;
                for ( i = 0; i < nstate; i++ ) {
                    
                    tlk->pattern_lk[k] += freqs[i] * tlk->partials[Node_id(root)][v];
                    v++;
                }
                //tlk->pattern_lk[k] = log(tlk->pattern_lk[k]);
                
                if ( tlk->scale ) {
                    //printf("scaling\n");
                    tlk->pattern_lk[k] += getLogScalingFactor( tlk, k);
                }
            }
            
            for ( int i = 0; i < tlk->sp->count; i++) {
                d2lnl[nodeCount*Node_id(n1)+Node_id(n2)] += ((tlk->pattern_lk[i] * lks[i] - (dfi[Node_id(n1)][i]*dfi[Node_id(n2)][i])) / (lks[i]*lks[i])) * tlk->sp->weights[i];
            }
            
            d2lnl[nodeCount*Node_id(n2)+Node_id(n1)] = d2lnl[nodeCount*Node_id(n1)+Node_id(n2)];
            
            memcpy( tlk->matrices[Node_id(n1)], pmats1, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
            SingleTreeLikelihood_update_one_node(tlk, n1);
            
            memcpy( tlk->matrices[Node_id(n2)], pmats2, tlk->sm->cat_count *tlk->matrix_size*sizeof(double));
            SingleTreeLikelihood_update_one_node(tlk, n2);
        }
    }
    free(lks);
    free(pmats1);
    free(pmats2);
    free_dmatrix(dfi, Tree_node_count(tlk->tree));
}

#pragma mark -

// calculate upper partials as used by derivatives
void calculate_upper(SingleTreeLikelihood *tlk, Node* node){
	
	if(!Node_isroot(node)){
		Node* parent = Node_parent(node);
		Node* sibling = Node_right(parent);
		int idNode = Node_id(node);
		
		if (sibling == node) {
			sibling = Node_left(parent);
		}
		int idSibling = Node_id(sibling);
		
		if(!Node_isroot(parent)){
			Node* grandParent = Node_parent(parent);
			
			int idMatrix = Node_id(parent);
			
			// The sons of the right node of the root are going to use the lower partials of the left node of the root
			if(Node_isroot(grandParent) && Node_right(grandParent) == parent){
				idMatrix = Node_id(Node_left(grandParent));
			}
			
			tlk->upper_partial_indexes[idNode] = idNode + Tree_node_count(tlk->tree);
			
			//printf("%s %d %d %d\n", node->name, tlk->upper_partial_indexes[idNode], tlk->upper_partial_indexes[Node_id(parent)], idMatrix);
			tlk->update_partials(tlk, tlk->upper_partial_indexes[idNode], tlk->upper_partial_indexes[Node_id(parent)], idMatrix, idSibling, idSibling);
		}
		// We dont need to calculate upper partials for the children of the root as it is using the lower partials of its sibling
		else{
			tlk->upper_partial_indexes[idNode] = idSibling;
		}
	}
	else{
		tlk->upper_partial_indexes[Node_id(node)] = -1;
	}
	
	if(!Node_isleaf(node)){
		calculate_upper(tlk, Node_left(node));
		calculate_upper(tlk, Node_right(node));
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
double dlnldt_uppper( SingleTreeLikelihood *tlk, Node *node, const double* pattern_likelihoods, const double* pattern_dlikelihoods ){
	const int patternCount = tlk->sp->count;
	double dlnl = 0;
	
	for ( int k = 0; k < patternCount; k++ ) {
		dlnl += pattern_dlikelihoods[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
	}
	
	return dlnl;
}

// Derivative of the likelihood function
// Assumes that root_partials are calculated correctly (i.e. _calculate_uppper is not called)
// This function does not change anything regarding partials, matrices, root_partials and pattern_lk
// TODO: Make thread safe
double dldt_uppper2( SingleTreeLikelihood *tlk, Node *node ){
	
	double bl = Node_distance(node);
	const int nodeId = Node_id(node);
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tempMatrixId];
	
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
	
	double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->sm->nstate;
	double* pattern_dlnl = tlk->pattern_lk + tlk->sp->count;
	
	tlk->calculate_branch_likelihood(tlk, root_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
	
	//tlk->node_log_likelihoods( tlk, root_partials, tlk->sm->m->_freqs, pattern_dlnl);
	int v = 0;
	
	const int nstate   = tlk->sm->nstate;
	const int patternCount = tlk->sp->count;
	const double* freqs = tlk->get_root_frequencies(tlk);
	double dl = 0;
	
	for ( int k = 0; k < patternCount; k++ ) {
		
		pattern_dlnl[k] = 0;
		for ( int i = 0; i < nstate; i++ ) {
			pattern_dlnl[k] += freqs[i] * root_partials[v];
			v++;
		}
		
		if ( tlk->scale ) {
			//printf("scaling\n");
			pattern_dlnl[k] += getLogScalingFactor( tlk, k);
		}
		dl += pattern_dlnl[k] * tlk->sp->weights[k];
	}
	
	return dl;
}

// Calculate derivative of likelihood with respect to branch length at node
// using the new upper likelihood calculation
// use a temporary matrix (the root's)
void calculate_dldt_uppper( SingleTreeLikelihood *tlk, Node *node, double* pattern_dlikelihoods ){
	
	double bl = Node_distance(node);
	const int nodeId = Node_id(node);
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tempMatrixId];
	
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
	
	tlk->calculate_branch_likelihood(tlk, root_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
	
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
	double* mat = tlk->matrices[tempMatrixId];
	
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
	double* pattern_d2lnl = tlk->pattern_lk + patternCount*2;
	const double* freqs = tlk->get_root_frequencies(tlk);
	tlk->calculate_branch_likelihood(tlk, root_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
	
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
    double* mat = tlk->matrices[tempMatrixId];
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
        tlk->calculate_branch_likelihood(tlk, root_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
        
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
	
	int tempMatrixId = Node_id(Tree_root(tlk->tree));
	double* mat = tlk->matrices[tempMatrixId];
	double* root_partials = tlk->root_partials  + tlk->sp->count*tlk->sm->nstate;
	double* pattern_dlnl = tlk->pattern_lk + tlk->sp->count;
	memset(pattern_dlnl, 0, sizeof(double)*tlk->sp->count);
    
    double dlnldQ = 0;
    tlk->sm->m->dQ_need_update = true;
    
    const int freqIndex = index - Parameters_count(tlk->sm->m->rates);
	const double* freqs = tlk->get_root_frequencies(tlk);
    double* transpose = dvector(nstate*nstate*tlk->sm->cat_count);
    
    // Frequencies
    if(index >= Parameters_count(tlk->sm->m->rates)){
        double *dpartials = dvector(patternCount);
        //calculate_dpartials_dfreqs(tlk, index, dpartials, pattern_likelihoods);
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
            
            tlk->calculate_branch_likelihood(tlk, root_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
            
            v = 0;
            for ( int k = 0; k < patternCount; k++ ) {
                for (int jj = 0; jj < nstate; jj++) {
                    dpartials[k] += freqs[jj]*root_partials[v];
                    v++;
                }
            }
        }
        for ( int k = 0; k < patternCount; k++ ) {
            dlnldQ += dpartials[k]/pattern_likelihoods[k] * tlk->sp->weights[k];
        }
        free(dpartials);
    }
    // Rates
    else{
        for ( int j = 0; j < Tree_node_count(tlk->tree); j++ ) {
            Node *n = nodes[j];
            const int nodeId = Node_id(n);
            
            if ( nodeId == rootId || nodeId == rightId ) {
                continue;
            }
            
            const double bl = Node_distance(n);
            
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                tlk->sm->m->dPdp(tlk->sm->m, index, &mat[i*tlk->matrix_size], bl* tlk->sm->get_rate(tlk->sm, i));
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
                
            
            tlk->calculate_branch_likelihood(tlk, root_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
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
    free(transpose);
	SingleTreeLikelihood_update_all_nodes(tlk);
	return dlnldQ;
}


#pragma mark -
#pragma mark Upper Likelihood

void SingleTreeLikelihood_set_upper_function( SingleTreeLikelihood *tlk, calculate_upper_t function ){
    if( function == NULL ){
        if( tlk->bm == NULL ){
            tlk->calculate_upper = _calculate_uppper;
        }
        else {
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if ( tlk->use_SIMD ){
                tlk->calculate_upper = _calculate_uppper_height_sse;
            }
            else {
                tlk->calculate_upper = _calculate_uppper_height;
            }
#else
            tlk->calculate_upper = _calculate_uppper_height;
#endif
        }
    }
    else {
        tlk->calculate_upper = function;
    }
}

// Should NOT be used before setting or removing a BranchModel
void SingleTreeLikelihood_use_upper( SingleTreeLikelihood *tlk, bool use_upper ){
    if( use_upper ){
		tlk->calculate_upper = _calculate_uppper;
        
        tlk->node_upper = NULL;
        
        SingleTreeLikelihood_set_upper_function(tlk, NULL);
        
        
		if ( tlk->sm->nstate == 4 ) {
			tlk->update_partials_upper      = update_partials_upper_4;
			tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_4;
			tlk->calculate_branch_likelihood = calculate_branch_likelihood_4;
			
#ifdef SSE3_ENABLED
            if ( tlk->use_SIMD ){
                tlk->update_partials_upper      = update_partials_upper_sse_4;
                tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_sse_4;
				tlk->calculate_branch_likelihood = calculate_branch_likelihood_4_SSE;
            }
#endif
        }
        else if( tlk->sm->nstate == 20 ){
#ifdef SSE3_ENABLED
            if ( tlk->use_SIMD ){
                tlk->update_partials_upper      = update_partials_upper_sse_20;
                tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_sse_20;
            }
            else {
                tlk->update_partials_upper      = update_partials_upper_20;
                tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_20;
            }
#else
            tlk->update_partials_upper      = update_partials_upper_20;
            tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_20;
#endif
        }
        else if( tlk->sm->nstate >= 60 ){
            tlk->update_partials_upper      = update_partials_upper_codon;
            tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_codon;
        }
        else {
            tlk->update_partials_upper      = update_partials_upper_general;
            tlk->node_log_likelihoods_upper = node_log_likelihoods_upper_general;
			tlk->calculate_branch_likelihood = calculate_branch_likelihood;
        }
    }

    tlk->use_upper = use_upper;
}


// Update all the upper partials
static void _preorder_traverse( SingleTreeLikelihood *tlk, Node *node ){
    if( !Node_isroot(node) ){
        tlk->update_partials_upper(tlk, node);
    }
    
    if( !Node_isleaf(node) ){
        _preorder_traverse(tlk, Node_left(node));
        _preorder_traverse(tlk, Node_right(node));
    }
}


// should be used when we optimize heights
// should notbe called on tips
double _calculate_uppper_height( SingleTreeLikelihood *tlk, Node *node ){
    
    //if( !tlk->update ) return tlk->lk;
    
    double bl;
    if( !Node_isroot(node) ){
        // Calculate probability matrix of parent
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            bl = tlk->bm->get(tlk->bm, node) * (Node_height(Node_parent(node)) - Node_height(node) );
            if(bl < 0 )
                fprintf(stderr, "calculate_partials_upper1: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->name, bl, tlk->bm->get(tlk->bm, node), Node_height(node), node->parent->name, Node_height(Node_parent(node)), Node_distance(node));
            
            tlk->sm->m->p_t(tlk->sm->m,
                            bl * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
        }
    }
    
    for (int i = 0; i < tlk->sm->cat_count; i++) {
        // Calculate probability matrix of left child
        bl = tlk->bm->get(tlk->bm, node->left) * (Node_height(node) - Node_height(node->left) );
        if(bl < 0 )
            fprintf(stderr, "calculate_partials_upper2: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->left->name, bl, tlk->bm->get(tlk->bm, node->left), Node_height(node->left), node->name, Node_height(node), Node_distance(node->left));
        
        tlk->sm->m->p_t(tlk->sm->m,
                        bl * tlk->sm->get_rate(tlk->sm, i),
                        &tlk->matrices[Node_id(node->left)][i*tlk->matrix_size]);
        
        // Calculate probability matrix of right child
        bl = tlk->bm->get(tlk->bm, node->right) * (Node_height(node) - Node_height(node->right) );
        if(bl < 0 )
            fprintf(stderr, "calculate_partials_upper3: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->right->name, bl, tlk->bm->get(tlk->bm, node->right), Node_height(node->right), node->name, Node_height(node), Node_distance(node->right));
        
        tlk->sm->m->p_t(tlk->sm->m,
                        bl * tlk->sm->get_rate(tlk->sm, i),
                        &tlk->matrices[Node_id(node->right)][i*tlk->matrix_size]);
    }
    
    
    // We have changed of node
    if ( tlk->node_upper != node ) {
        // Calculate all the upper partial likelihoods
        if ( tlk->node_upper == NULL ) {
            _preorder_traverse(tlk, Tree_root(tlk->tree));
        }
        // Update probability matrix upper and partial likelihoods of the previous node
        // not going to work if there is fixed node
        // should do this step recursively from the index_upper up to postorder_idx
        else {
            Node *n = tlk->node_upper;
            
            
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                bl = tlk->bm->get(tlk->bm, n) * (Node_height(Node_parent(n)) - Node_height(n) );
                if(bl < 0 )
                    fprintf(stderr, "calculate_partials_upper1: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->name, bl, tlk->bm->get(tlk->bm, node), Node_height(node), node->parent->name, Node_height(Node_parent(node)), Node_distance(node));
                
                tlk->sm->m->p_t(tlk->sm->m,
                                bl * tlk->sm->get_rate(tlk->sm, i),
                                &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                
                bl = tlk->bm->get(tlk->bm, n->left) * (Node_height(n) - Node_height(n->left) );
                if(bl < 0 )
                    fprintf(stderr, "calculate_partials_upper2: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->left->name, bl, tlk->bm->get(tlk->bm, node->left), Node_height(node->left), node->name, Node_height(node), Node_distance(node->left));
                
                tlk->sm->m->p_t(tlk->sm->m,
                                bl * tlk->sm->get_rate(tlk->sm, i),
                                &tlk->matrices[Node_id(n->left)][i*tlk->matrix_size]);
                
                bl = tlk->bm->get(tlk->bm, n->right) * (Node_height(n) - Node_height(n->right) );
                if(bl < 0 )
                    fprintf(stderr, "calculate_partials_upper3: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->right->name, bl, tlk->bm->get(tlk->bm, node->right), Node_height(node->right), node->name, Node_height(node), Node_distance(node->right));
                
                tlk->sm->m->p_t(tlk->sm->m,
                                bl * tlk->sm->get_rate(tlk->sm, i),
                                &tlk->matrices[Node_id(n->right)][i*tlk->matrix_size]);
            }
            
            
            // update its lower partials
			int indx_child1 = Node_id(Node_left(n));
			int indx_child2 = Node_id(Node_right(n));
			tlk->update_partials( tlk, Node_id(n), indx_child1, indx_child1, indx_child2, indx_child2);
			
            // update lower partials of its parent
            if( !Node_isroot(n) ){
                tlk->update_partials(tlk, Node_id(Node_parent(n)), Node_id(n), Node_id(n), Node_id(Node_sibling(n)), Node_id(Node_sibling(n)) );
            }
            
            // update upper partials of its sibling and its sibling's descendents
            _preorder_traverse(tlk, Node_sibling(n));
        }
        tlk->update_nodes[Node_id(node)]        = false;
        tlk->update_nodes[Node_id(node->left)]  = false;
        tlk->update_nodes[Node_id(node->right)] = false;
    }
    
    tlk->node_upper = node;
	
	tlk->update_partials(tlk, Node_id(node), Node_id(Node_left(node)), Node_id(Node_left(node)), Node_id(Node_right(node)), Node_id(Node_right(node)) );
    // we don't need to update the upper partials of the children at this stage as they are not needed
    
    
    tlk->lk = 0;
	int i = 0;
    
    if ( Node_isroot(node) ) {
        tlk->integrate_partials(tlk, tlk->partials[Node_id(node)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
        tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
        for ( ; i < tlk->sp->count; i++) {
            tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
        }
    }
    else {
        tlk->node_log_likelihoods_upper( tlk, node );
		
		double* pattern_lk = tlk->pattern_lk+ tlk->sp->count;
        for ( ; i < tlk->sp->count; i++) {
            tlk->lk += log(pattern_lk[i]) * tlk->sp->weights[i];
        }
    }
    
    
    //TODO
	if ( tlk->lk == -INFINITY ) {
		fprintf(stderr, "_calculate: rescaling\n");
        //		SingleTreeLikelihood_use_rescaling(tlk, true );
        //
        //		SingleTreeLikelihood_update_all_nodes( tlk );
        //		tlk->lk = 0;
        //		_calculate_partials( tlk, Tree_root(tlk->tree) );
        //
        //
        //		for ( i = 0; i < tlk->sp->count; i++) {
        //			tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
        //		}
		
	}
	
	for ( i = 0; i < Tree_node_count(tlk->tree); i++) tlk->update_nodes[i] = false;
	
    //printf("%s %f\n", Node_name(node), tlk->lk);
	return tlk->lk;
}

double _calculate_uppper_height_sse( SingleTreeLikelihood *tlk, Node *node ){
    
    //if( !tlk->update ) return tlk->lk;
    
    double bl;
    if( !Node_isroot(node) ){
        // Calculate probability matrix of nodesim
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            bl = tlk->bm->get(tlk->bm, node) * Node_time_elapsed(node);
            if(bl < 0 )
                fprintf(stderr, "calculate_partials_upper1: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->name, bl, tlk->bm->get(tlk->bm, node), Node_height(node), node->parent->name, Node_height(Node_parent(node)), Node_distance(node));
            
            tlk->sm->m->p_t(tlk->sm->m,
                            bl * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
        }
    }
    
    // Calculate probability matrix of left child
    if( Node_isleaf(node->left) ){
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            bl = tlk->bm->get(tlk->bm, node->left) * Node_time_elapsed(node->left);
            if(bl < 0 )
                fprintf(stderr, "calculate_partials_upper2: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->left->name, bl, tlk->bm->get(tlk->bm, node->left), Node_height(node->left), node->name, Node_height(node), Node_distance(node->left));
            
            tlk->sm->m->p_t_transpose(tlk->sm->m,
                            bl * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(node->left)][i*tlk->matrix_size]);
        }
    }
    else {
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            bl = tlk->bm->get(tlk->bm, node->left) * Node_time_elapsed(node->left);
            if(bl < 0 )
                fprintf(stderr, "calculate_partials_upper2: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->left->name, bl, tlk->bm->get(tlk->bm, node->left), Node_height(node->left), node->name, Node_height(node), Node_distance(node->left));
            
            tlk->sm->m->p_t(tlk->sm->m,
                            bl * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(node->left)][i*tlk->matrix_size]);
        }
    }
    
    
    // Calculate probability matrix of right child
    if( Node_isleaf(node->right) ){
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            bl = tlk->bm->get(tlk->bm, node->right) * Node_time_elapsed(node->right);
            if(bl < 0 )
                fprintf(stderr, "calculate_partials_upper3: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->right->name, bl, tlk->bm->get(tlk->bm, node->right), Node_height(node->right), node->name, Node_height(node), Node_distance(node->right));
            
            tlk->sm->m->p_t_transpose(tlk->sm->m,
                            bl * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(node->right)][i*tlk->matrix_size]);
        }
    }
    else {
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            bl = tlk->bm->get(tlk->bm, node->right) * Node_time_elapsed(node->right);
            if(bl < 0 )
                fprintf(stderr, "calculate_partials_upper3: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->right->name, bl, tlk->bm->get(tlk->bm, node->right), Node_height(node->right), node->name, Node_height(node), Node_distance(node->right));
            
            tlk->sm->m->p_t(tlk->sm->m,
                            bl * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(node->right)][i*tlk->matrix_size]);
        }
    }
    
    
    
    
    // We have changed of node
    if ( tlk->node_upper != node ) {
        // Calculate all the upper partial likelihoods
        if ( tlk->node_upper == NULL ) {
            tlk->calculate(tlk);
            _preorder_traverse(tlk, Tree_root(tlk->tree));
        }
        // Update probability matrix upper and partial likelihoods of the previous node
        // not going to work if there is fixed node
        // should do this step recursively from the index_upper up to postorder_idx
        else {
            Node *n = tlk->node_upper;
            
            for ( int i = 0; i < tlk->sm->cat_count; i++) {
                bl = tlk->bm->get(tlk->bm, n) * (Node_height(Node_parent(n)) - Node_height(n) );
                if(bl < 0 )
                    fprintf(stderr, "calculate_partials_upper1: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->name, bl, tlk->bm->get(tlk->bm, node), Node_height(node), node->parent->name, Node_height(Node_parent(node)), Node_distance(node));
                
                tlk->sm->m->p_t(tlk->sm->m,
                                bl * tlk->sm->get_rate(tlk->sm, i),
                                &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
            }
            
            if( Node_isleaf(n->left)){
                for (int i = 0; i < tlk->sm->cat_count; i++) {
                    
                    bl = tlk->bm->get(tlk->bm, n->left) * (Node_height(n) - Node_height(n->left) );
                    if(bl < 0 )
                        fprintf(stderr, "calculate_partials_upper2: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->left->name, bl, tlk->bm->get(tlk->bm, node->left), Node_height(node->left), node->name, Node_height(node), Node_distance(node->left));
                    
                    tlk->sm->m->p_t_transpose(tlk->sm->m,
                                    bl * tlk->sm->get_rate(tlk->sm, i),
                                    &tlk->matrices[Node_id(n->left)][i*tlk->matrix_size]);
                }
                
            }
            else {
                for (int i = 0; i < tlk->sm->cat_count; i++) {
                    
                    bl = tlk->bm->get(tlk->bm, n->left) * (Node_height(n) - Node_height(n->left) );
                    if(bl < 0 )
                        fprintf(stderr, "calculate_partials_upper2: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->left->name, bl, tlk->bm->get(tlk->bm, node->left), Node_height(node->left), node->name, Node_height(node), Node_distance(node->left));
                    
                    tlk->sm->m->p_t(tlk->sm->m,
                                    bl * tlk->sm->get_rate(tlk->sm, i),
                                    &tlk->matrices[Node_id(n->left)][i*tlk->matrix_size]);
                }
            }
            
            
            if( Node_isleaf(n->right)){
                
                for (int i = 0; i < tlk->sm->cat_count; i++) {
                    bl = tlk->bm->get(tlk->bm, n->right) * (Node_height(n) - Node_height(n->right) );
                    if(bl < 0 )
                        fprintf(stderr, "calculate_partials_upper3: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->right->name, bl, tlk->bm->get(tlk->bm, node->right), Node_height(node->right), node->name, Node_height(node), Node_distance(node->right));
                    
                    tlk->sm->m->p_t_transpose(tlk->sm->m,
                                    bl * tlk->sm->get_rate(tlk->sm, i),
                                    &tlk->matrices[Node_id(n->right)][i*tlk->matrix_size]);
                }
            }
            else {
                for (int i = 0; i < tlk->sm->cat_count; i++) {                    
                    bl = tlk->bm->get(tlk->bm, n->right) * (Node_height(n) - Node_height(n->right) );
                    if(bl < 0 )
                        fprintf(stderr, "calculate_partials_upper3: %s branch length = %f rate = %f height = %f - parent height [%s]= %f (%f)\n", node->right->name, bl, tlk->bm->get(tlk->bm, node->right), Node_height(node->right), node->name, Node_height(node), Node_distance(node->right));
                    
                    tlk->sm->m->p_t(tlk->sm->m,
                                    bl * tlk->sm->get_rate(tlk->sm, i),
                                    &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
            }
            
            
            
            
			// update its lower partials
			tlk->update_partials(tlk, Node_id(n), Node_id(Node_left(n)), Node_id(Node_left(n)), Node_id(Node_right(n)), Node_id(Node_right(n)) );
			
            // update lower partials of its parent
            if( !Node_isroot(n) ){
                tlk->update_partials(tlk, Node_id(Node_parent(n)), Node_id(n), Node_id(n), Node_id(Node_sibling(n)), Node_id(Node_sibling(n)) );
            }
            
            // update upper partials of its sibling and its sibling's descendents
            _preorder_traverse(tlk, Node_sibling(n));
        }
        tlk->update_nodes[Node_id(node)]        = false;
        tlk->update_nodes[Node_id(node->left)]  = false;
        tlk->update_nodes[Node_id(node->right)] = false;
    }
    
    tlk->node_upper = node;
    
    tlk->update_partials(tlk, Node_id(node), Node_id(Node_left(node)), Node_id(Node_left(node)), Node_id(Node_right(node)), Node_id(Node_right(node)) );
    // we don't need to update the upper partials of the children at this stage as they are not needed
    
    
    tlk->lk = 0;
	int i = 0;
    
    if ( Node_isroot(node) ) {
        tlk->integrate_partials(tlk, tlk->partials[Node_id(node)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
        tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
        for ( ; i < tlk->sp->count; i++) {
            tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
        }
    }
    else {
        tlk->node_log_likelihoods_upper( tlk, node );
		
		double* pattern_lk = tlk->pattern_lk+ tlk->sp->count;
        for ( ; i < tlk->sp->count; i++) {
            tlk->lk += log(pattern_lk[i]) * tlk->sp->weights[i];
        }
    }
    
    
    //TODO
	if ( tlk->lk == -INFINITY ) {
		fprintf(stderr, "_calculate: rescaling\n");
        //		SingleTreeLikelihood_use_rescaling(tlk, true );
        //
        //		SingleTreeLikelihood_update_all_nodes( tlk );
        //		tlk->lk = 0;
        //		_calculate_partials( tlk, Tree_root(tlk->tree) );
        //
        //
        //		for ( i = 0; i < tlk->sp->count; i++) {
        //			tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
        //		}
		
	}
	
	for ( i = 0; i < Tree_node_count(tlk->tree); i++) tlk->update_nodes[i] = false;
	
    //printf("%s %f\n", Node_name(node), tlk->lk);
	return tlk->lk;
}



// should be used when we optimize distances
// should notbe called at the root
// The last node in post order always need an update and the partials at the root too
double _calculate_uppper( SingleTreeLikelihood *tlk, Node *node ){
	
    //printf("_calculate_uppper %s\n", node->name);
    //if( !tlk->update ) return tlk->lk;
	
	if(true){
		
		/*if(tlk->update_upper){
			
			//  make sure lower partials are calculated
			// a call to this function set update_upper to true if there was an update
			tlk->calculate(tlk);
			
			//calculate_upper(tlk, Tree_root(tlk->tree));
			tlk->update_upper = false;
		}*/
		
		// change of node
		if ( tlk->node_upper != node ) {
			if ( tlk->node_upper == NULL || ( Node_isleaf(node) && Node_sibling(node) != tlk->node_upper &&  Node_right(node->parent) == node ) || ( !Node_isleaf(node) && Node_right(node) != tlk->node_upper) ) {
				//  make sure lower partials are calculated
				tlk->calculate(tlk);
				calculate_upper(tlk, Tree_root(tlk->tree));
			}
			else {
				Node *n = tlk->node_upper;
				
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
				if( tlk->use_SIMD && Node_isleaf(n) ){
					for (int i = 0; i < tlk->sm->cat_count; i++) {
						tlk->sm->m->p_t_transpose(tlk->sm->m,
												  Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
												  &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
					}
				}
				else{
					for (int i = 0; i < tlk->sm->cat_count; i++) {
						tlk->sm->m->p_t(tlk->sm->m,
										Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
										&tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
					}
				}
#else
				for (int i = 0; i < tlk->sm->cat_count; i++) {
					tlk->sm->m->p_t(tlk->sm->m,
									Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
									&tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
				}
#endif
				
				// update lower partials of the parent
				tlk->update_partials(tlk, Node_id(Node_parent(n)), Node_id(n), Node_id(n), Node_id(Node_sibling(n)), Node_id(Node_sibling(n)) );
				
				// update upper partials of its sibling and its sibling's descendents
				calculate_upper(tlk, Node_sibling(n));
				
			}
		}
//		SingleTreeLikelihood_update_all_nodes(tlk);
//		tlk->calculate(tlk);
//		calculate_upper(tlk, Tree_root(tlk->tree));
		
		const int nodeId = Node_id(node);
        int tempMatrixId = Node_id(Tree_root(tlk->tree));
        double* mat = tlk->matrices[tempMatrixId];
		
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
		tlk->calculate_branch_likelihood(tlk, tlk->root_partials, tlk->upper_partial_indexes[nodeId], nodeId, tempMatrixId );
		
		tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
		
		double lk = 0;
		for ( int i = 0; i < tlk->sp->count; i++) {
			lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
		}
		
		// set for update later
		//SingleTreeLikelihood_update_one_node(tlk, node);
		return lk;
	}
	
    if ( tlk->node_upper != node ) {
        if ( tlk->node_upper == NULL || ( Node_isleaf(node) && Node_sibling(node) != tlk->node_upper &&  Node_right(node->parent) == node ) || ( !Node_isleaf(node) && Node_right(node) != tlk->node_upper) ) {
            //  make sure lower partials are calculated
            tlk->calculate(tlk);
            _preorder_traverse(tlk, Tree_root(tlk->tree));
        }
        else {
//            if( tlk->node_upper > node->postorder_idx ){
//                printf("_calculate_uppper %s %s\n", Node_name(node), Node_name(Tree_get_node(tlk->tree, POSTORDER, tlk->index_upper)));
//                exit(1);
//            }
//            if( Node_isleaf(node) ){
//                if( Node_parent(node)->right == node && Node_sibling(node)->postorder_idx != tlk->index_upper ) error("_calculate_uppper leaf and not sibling\n");
//            }
//            else{
//                if( Node_right(node)->postorder_idx != tlk->index_upper ) error("_calculate_uppper\n");
//            }
            
            Node *n = tlk->node_upper;
            
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD && Node_isleaf(n) ){
                for (int i = 0; i < tlk->sm->cat_count; i++) {
                    tlk->sm->m->p_t_transpose(tlk->sm->m,
                                              Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
                                              &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
            }
            else{
                for (int i = 0; i < tlk->sm->cat_count; i++) {
                    tlk->sm->m->p_t(tlk->sm->m,
                                    Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
                                    &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                }
            }
#else
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                tlk->sm->m->p_t(tlk->sm->m,
                                Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
                                &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
            }
#endif
            
            // update lower partials of the parent
            tlk->update_partials(tlk, Node_id(Node_parent(n)), Node_id(n), Node_id(n), Node_id(Node_sibling(n)), Node_id(Node_sibling(n)) );
            
            // update upper partials of its sibling and its sibling's descendents
            _preorder_traverse(tlk, Node_sibling(n));
            
        }
    }
    
    
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
    if( tlk->use_SIMD && Node_isleaf(node) ){
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            tlk->sm->m->p_t_transpose(tlk->sm->m,
                                      Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
                                      &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
        }
    }
    else{
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            tlk->sm->m->p_t(tlk->sm->m,
                            Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
        }
    }
#else
    for (int i = 0; i < tlk->sm->cat_count; i++) {
        tlk->sm->m->p_t(tlk->sm->m,
                        Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
                        &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
    }
#endif
    
    tlk->node_upper = node;
    
    
    tlk->node_log_likelihoods_upper( tlk, node );
	
	tlk->lk = 0;
	int i = 0;
	
	double* pattern_lk = tlk->pattern_lk+ tlk->sp->count;
	for ( ; i < tlk->sp->count; i++) {
		tlk->lk += log(pattern_lk[i]) * tlk->sp->weights[i];
	}
    
    //TODO
	if ( tlk->lk == -INFINITY ) {
		fprintf(stderr, "_calculate_upper: rescaling\n");
        //		SingleTreeLikelihood_use_rescaling(tlk, true );
        //
        //		SingleTreeLikelihood_update_all_nodes( tlk );
        //		tlk->lk = 0;
        //		_calculate_partials( tlk, Tree_root(tlk->tree) );
        //
        //
        //		for ( i = 0; i < tlk->sp->count; i++) {
        //			tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
        //		}
		
	}
    
	return tlk->lk;
}

// Can be called on a tip but certainly not on the root
// The likelihood is calculated at node's parent
double calculate_uppper_2nodes( SingleTreeLikelihood *tlk, Node *node ){
	
    //if( !tlk->update ) return tlk->lk;
	
    if( Tree_topology_changed(tlk->tree) ){
        printf("reorder calculate upper\n");
        SingleTreeLikelihood_rearrange_partials(tlk);// this function makes the tree update its internal strutures (e.g. postorder...)
    }
	
	
    if ( tlk->node_upper != node ) {
        if ( tlk->node_upper == NULL || Node_right(node) != tlk->node_upper ) {
			
            //  make sure lower partials are calculated
            tlk->calculate(tlk);
            _preorder_traverse(tlk, Tree_root(tlk->tree));
        }
        else {
//            if( Node_isleaf(node) ){
//                if( Node_parent(node)->right == node && Node_sibling(node)->postorder_idx != tlk->index_upper ) error("_calculate_uppper leaf and not sibling\n");
//            }
//            else{
//                if( Node_right(node)->postorder_idx != tlk->index_upper ) error("_calculate_uppper\n");
//            }
            
            Node *n = tlk->node_upper;
            
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
            if( tlk->use_SIMD ){
                if( Node_isleaf(n) ){
                    for (int i = 0; i < tlk->sm->cat_count; i++) {
                        tlk->sm->m->p_t_transpose(tlk->sm->m,
                                                  Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
                                                  &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                    }
                }
                else{
                    for (int i = 0; i < tlk->sm->cat_count; i++) {
                        tlk->sm->m->p_t(tlk->sm->m,
                                        Node_distance(n) * tlk->sm->get_rate(tlk->sm, i),
                                        &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                    }
                }
                
                if( Node_isleaf( Node_left(n) ) ){
                    for (int i = 0; i < tlk->sm->cat_count; i++) {
                        tlk->sm->m->p_t_transpose(tlk->sm->m,
                                                  Node_distance(Node_left(n)) * tlk->sm->get_rate(tlk->sm, i),
                                                  &tlk->matrices[Node_id(Node_left(n))][i*tlk->matrix_size]);
                    }
                }
                else{
                    for (int i = 0; i < tlk->sm->cat_count; i++) {
                        tlk->sm->m->p_t(tlk->sm->m,
                                        Node_distance(Node_left(n)) * tlk->sm->get_rate(tlk->sm, i),
                                        &tlk->matrices[Node_id(Node_left(n))][i*tlk->matrix_size]);
                    }
                }
                
                if( Node_isleaf( Node_right(n) ) ){
                    for (int i = 0; i < tlk->sm->cat_count; i++) {
                        tlk->sm->m->p_t_transpose(tlk->sm->m,
                                                  Node_distance(Node_right(n)) * tlk->sm->get_rate(tlk->sm, i),
                                                  &tlk->matrices[Node_id(Node_right(n))][i*tlk->matrix_size]);
                    }
                }
                else{
                    for (int i = 0; i < tlk->sm->cat_count; i++) {
                        tlk->sm->m->p_t(tlk->sm->m,
                                        Node_distance(Node_right(n)) * tlk->sm->get_rate(tlk->sm, i),
                                        &tlk->matrices[Node_id(Node_right(n))][i*tlk->matrix_size]);
                    }
                }
            }
            else {
                for (int i = 0; i < tlk->sm->cat_count; i++) {
                    double r = tlk->sm->get_rate(tlk->sm, i);
                    
                    tlk->sm->m->p_t(tlk->sm->m,
                                    Node_distance(n) * r,
                                    &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                    
                    tlk->sm->m->p_t(tlk->sm->m,
                                    Node_distance(Node_left(n)) * r,
                                    &tlk->matrices[Node_id(Node_left(n))][i*tlk->matrix_size]);
                    
                    tlk->sm->m->p_t(tlk->sm->m,
                                    Node_distance(Node_right(n)) * r,
                                    &tlk->matrices[Node_id(Node_right(n))][i*tlk->matrix_size]);
                }
            }
#else
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                double r = tlk->sm->get_rate(tlk->sm, i);
                
                tlk->sm->m->p_t(tlk->sm->m,
                                Node_distance(n) * r,
                                &tlk->matrices[Node_id(n)][i*tlk->matrix_size]);
                
                tlk->sm->m->p_t(tlk->sm->m,
                                Node_distance(n->left) * r,
                                &tlk->matrices[Node_id(n->left)][i*tlk->matrix_size]);
                
                tlk->sm->m->p_t(tlk->sm->m,
                                Node_distance(n->right) * r,
                                &tlk->matrices[Node_id(n->right)][i*tlk->matrix_size]);
            }
#endif
            
            // update lower partials of the parent
            tlk->update_partials(tlk, Node_id(Node_parent(n)), Node_id(n), Node_id(n), Node_id(Node_sibling(n)), Node_id(Node_sibling(n)) );
            
            // update upper partials of its sibling and its sibling's descendents
            _preorder_traverse(tlk, Node_sibling(n));
            
        }
    }
    
    Node *parent  = Node_parent(node);
    
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
    if( tlk->use_SIMD ){
        if( Node_isleaf(node) ){
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                tlk->sm->m->p_t_transpose(tlk->sm->m,
                                          Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
                                          &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
            }
        }
        else{
            for (int i = 0; i < tlk->sm->cat_count; i++) {
                tlk->sm->m->p_t(tlk->sm->m,
                                Node_distance(node) * tlk->sm->get_rate(tlk->sm, i),
                                &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
            }
        }
        
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            tlk->sm->m->p_t(tlk->sm->m,
                            Node_distance(parent) * tlk->sm->get_rate(tlk->sm, i),
                            &tlk->matrices[Node_id(parent)][i*tlk->matrix_size]);
        }
        
        
    }
    else {
        Node *parent  = Node_parent(node);
        double r;
        for (int i = 0; i < tlk->sm->cat_count; i++) {
            r = tlk->sm->get_rate(tlk->sm, i);
            
            tlk->sm->m->p_t(tlk->sm->m,
                            Node_distance(node) * r,
                            &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
            
            tlk->sm->m->p_t(tlk->sm->m,
                            Node_distance(parent) * r,
                            &tlk->matrices[Node_id(parent)][i*tlk->matrix_size]);
        }
    }
#else
    
    
    double r;
    for (int i = 0; i < tlk->sm->cat_count; i++) {
        r = tlk->sm->get_rate(tlk->sm, i);
        
        tlk->sm->m->p_t(tlk->sm->m,
                        Node_distance(node) * r,
                        &tlk->matrices[Node_id(node)][i*tlk->matrix_size]);
        
        tlk->sm->m->p_t(tlk->sm->m,
                        Node_distance(parent) * r,
                        &tlk->matrices[Node_id(parent)][i*tlk->matrix_size]);
    }
#endif
    
    tlk->node_upper = node;
    
    tlk->update_partials(tlk, Node_id(parent), Node_id(node), Node_id(node), Node_id(Node_sibling(node)), Node_id(Node_sibling(node)) );
    
    
    tlk->node_log_likelihoods_upper( tlk, parent );
	
	tlk->lk = 0;
	int i = 0;
	double* pattern_lk = tlk->pattern_lk+ tlk->sp->count;
	for ( ; i < tlk->sp->count; i++) {
		tlk->lk += log(pattern_lk[i]) * tlk->sp->weights[i];
	}
    
    //TODO
	if ( tlk->lk == -INFINITY ) {
		fprintf(stderr, "_calculate: rescaling\n");
        //		SingleTreeLikelihood_use_rescaling(tlk, true );
        //
        //		SingleTreeLikelihood_update_all_nodes( tlk );
        //		tlk->lk = 0;
        //		_calculate_partials( tlk, Tree_root(tlk->tree) );
        //
        //
        //		for ( i = 0; i < tlk->sp->count; i++) {
        //			tlk->lk += tlk->pattern_lk[i] * tlk->sp->weights[i];
        //		}
		
	}
    
    for ( i = 0; i < Tree_node_count(tlk->tree); i++) tlk->update_nodes[i] = false;
    
	return tlk->lk;
}



#pragma mark -
#pragma mark Second order Taylor series approximation

void SingleTreeLikelihood_fisher_information( SingleTreeLikelihood *tlk, double *hessian ) {
    int dim = Tree_node_count(tlk->tree)-2;
    SingleTreeLikelihood_Hessian(tlk, hessian, NULL);
    
    for ( int i = 0; i < dim*dim; i++) {
		hessian[i] = -hessian[i];
	}
}

void SingleTreeLikelihood_covariance( SingleTreeLikelihood *tlk, double *hessian ) {
    int dim = Tree_node_count(tlk->tree)-2;
    SingleTreeLikelihood_fisher_information(tlk, hessian);
    
    inverse2(hessian, dim);
}

/*
 * Approximation of the Hessian using the second-order central difference approximation
 * d^2f(x)/dx_idx_j
 */
bool SingleTreeLikelihood_Hessian( SingleTreeLikelihood *tlk, double *hessian, double *gradient ) {
	
    //PooledTreeLikelihood *pool = new_PooledTreeLikelihood(tlk, 10, true, false );
    
	Tree *tree = tlk->tree;
	int dim = Tree_node_count(tree)-2;// ignore the root and the right child because it is 0 anyway
	double eps = 0.00001;
	
	double *e = dvector(dim);
	
	int i = 0;
	Node **nodes = Tree_get_nodes(tree, POSTORDER);
	for ( i = 0; i < dim; i++) {
		e[i] = eps * Node_distance( nodes[i] );
	}
    
	double pp, pm, mp, mm;
	
	double lk = tlk->lnl_bl = tlk->calculate(tlk);
	
	int j = 0;
    double d1,d2;
	int k = 0;
    
	for ( i = 0; i < dim; i++ ){
        
        d1 = Node_distance( nodes[i] );
        // +
        Node_set_distance( nodes[i], d1 + e[i] );
        SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
        pp = tlk->calculate(tlk);
        
        // -
        Node_set_distance( nodes[i], d1 - e[i] );
        SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
        mm = tlk->calculate(tlk);
        
        
        Node_set_distance( nodes[i], d1 );
        SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
        
        hessian[i*dim+i] = (pp + mm -2.*lk)/(e[i]*e[i]);
        
        if ( gradient != NULL ) {
            gradient[i*dim+i] = (pp - mm)/(4.*e[i]);
        }
        
//        if(hessian[i*dim+i] >= 0 ) {
//            printf("\nPositive second derivative: %s %f\n", Node_name(nodes[i]), d1 );
//            free(e);
//            return false;
//        }
        
		for ( j = i+1; j < dim; j++ ){
            
            d2 = Node_distance( nodes[j] );
            
            // + +
            Node_set_distance( nodes[i], d1 + e[i] );
            Node_set_distance( nodes[j], d2 + e[j] );
            
            SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
            SingleTreeLikelihood_update_one_node(tlk, nodes[j]);
            
            pp = tlk->calculate(tlk);
            
            // - -
            Node_set_distance( nodes[i], d1 - e[i] );
            Node_set_distance( nodes[j], d2 - e[j] );
            
            SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
            SingleTreeLikelihood_update_one_node(tlk, nodes[j]);
            
            mm = tlk->calculate(tlk);
            
            // + -
            Node_set_distance( nodes[i], d1 + e[i] );
            SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
            pm = tlk->calculate(tlk);
            
            // - +
            Node_set_distance( nodes[i], d1 - e[i] );
            Node_set_distance( nodes[j], d2 + e[j] );
            
            SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
            SingleTreeLikelihood_update_one_node(tlk, nodes[j]);
            mp = tlk->calculate(tlk);
            
            hessian[j*dim+i] = hessian[i*dim+j] = (pp + mm - pm - mp) / (4.0*e[i]*e[j]);
            
            //fprintf(stderr, "%f %f %f %f %.10f %.50f\n", pp, mm, pm, mp, 4*(e[i]*e[j]), hessian[i*dim+j]);
            
            Node_set_distance( nodes[i], d1 );
            Node_set_distance( nodes[j], d2 );
            
            SingleTreeLikelihood_update_one_node(tlk, nodes[i]);
            SingleTreeLikelihood_update_one_node(tlk, nodes[j]);
            k++;
            fprintf(stderr, "%d/%d\r",k,(dim*(dim-1)/2));
            
//            if(hessian[i*dim+j] >= 0 ) {
//                printf("\nPositive second derivative: %s %f %s %f\n", Node_name(nodes[i]), d1, Node_name(nodes[j]), d2 );
//                free(e);
//                return false;
//            }
        }
		
		
	}
	fprintf(stderr, "\n");
	free(e);
    return true;
}

/*
 * Approximation of the a diagonal Hessian using the second-order central difference approximation or forward difference if the parameter is at the boundry
 * d^2f(x)/dx_idx_j
 */
bool SingleTreeLikelihood_Hessian_diag( SingleTreeLikelihood *tlk, double *hessian, int *len ) {
	
	Tree *tree = tlk->tree;
	int dim = Tree_node_count(tree);// ignore the root and the right child because it is 0 anyway
	//double eps = 0.00001;
    double eps = 0.001;
	
	Node **nodes = Tree_nodes(tree);
    
	double p, m;
	
	double lnl = tlk->lnl_bl = tlk->calculate(tlk);
	
    double d,e;
	
    
	for ( int i = 0; i < dim; i++ ){
        Node *n = nodes[i];
        
        if ( Node_isroot(n) || (Node_isroot(Node_parent(n)) && Node_right(Node_parent(n)) == n )) {
            continue;
        }
        
        d = Node_distance( n );
        
        //if( d > 1e-7 ){
            
            e = eps * d;
            
            // +
            Node_set_distance( n, d + e );
            SingleTreeLikelihood_update_one_node(tlk, n);
            p = tlk->calculate(tlk);
            
            // -
            Node_set_distance( n, d - e );
            SingleTreeLikelihood_update_one_node(tlk, n);
            m = tlk->calculate(tlk);
            
            
            Node_set_distance( n, d );
            SingleTreeLikelihood_update_one_node(tlk, n);
            
            hessian[i] = (p + m -2*lnl)/(e*e);
            
            
            if(hessian[i] >= 0 ) {
                printf("Non-negative second derivative: %s bl: %e d2: %e\n", Node_name(n), d, hessian[i] );
                //return false;
            }
        
//        }
//        else {
            //Node_set_sticky(nodes[i], true);
            //            e = 0.00001*d;//1.0/tlk->sp->nsites;
            //
            //            // +
            //            Node_set_distance( nodes[i], d + e );
            //            SingleTreeLikelihood_update_one_node(tlk, i);
            //            p = tlk->calculate(tlk);
            //
            //            // ++
            //            Node_set_distance( nodes[i], d + 2*e );
            //            SingleTreeLikelihood_update_one_node(tlk, i);
            //            double pp = tlk->calculate(tlk);
            //
            //
            //            Node_set_distance( nodes[i], d );
            //            SingleTreeLikelihood_update_one_node(tlk, i);
            //
            //            hessian[i] = (pp + -2*p + lnl)/(e*e);
            //            if(hessian[i] >= 0 ) {
            //                printf("Positive second derivative: %s bl: %e d2: %e [%e %e %e] e: %e\n", Node_name(nodes[i]), d, hessian[i],pp, p, lnl,e );
            //                //return false;
            //            }
            
        //}
        
	}
    
    return true;
}


#pragma mark -
// MARK: Optimizable

void Optimizable_init( Optimizable *opt, bool optimize, opt_algorithm method, const int max_iteration, const double tolx, const double tolfx ){
	(*opt).optimize      = optimize;
	(*opt).method        = method;
	(*opt).max_iteration = max_iteration;
	(*opt).tolx          = tolx;
	(*opt).tolfx         = tolfx;
}

void Optimizable_copy( const Optimizable *src, Optimizable *dst ){
	(*dst).optimize      = (*src).optimize;
	(*dst).method        = (*src).method;
	(*dst).max_iteration = (*src).max_iteration;
	(*dst).tolx          = (*src).tolx;
	(*dst).tolfx         = (*src).tolfx;
}

void OptConfig_init( SingleTreeLikelihood *tlk ){
	Optimizable_init(&tlk->opt.freqs, false, OPT_BRENT, 100, 0.001, OPTIMIZATION_PRECISION);
	Optimizable_init(&tlk->opt.relative_rates, false, OPT_BRENT, 100, 0.001, OPTIMIZATION_PRECISION);
	Optimizable_init(&tlk->opt.bl, false, OPT_BRENT, 100, 0.001, OPTIMIZATION_PRECISION);
	Optimizable_init(&tlk->opt.rates, false, OPT_BRENT, 100, 0.0001, OPTIMIZATION_PRECISION);
	Optimizable_init(&tlk->opt.heights, false, OPT_BRENT, 100, 0.0001, OPTIMIZATION_PRECISION);
    tlk->opt.topology_optimize = 0;
    tlk->opt.topology_alogrithm = TREE_SEARCH_NNI;
    tlk->opt.topology_threads = 1;
	tlk->opt.max_rounds = 10000;
	tlk->opt.precision = 0.001;
	tlk->opt.verbosity = 1;
    tlk->opt.interruptible = false;
}

void OptConfig_copy( const OptConfig *src, OptConfig *dst ){
	Optimizable_copy(&src->freqs, &dst->freqs);
	Optimizable_copy(&src->relative_rates, &dst->relative_rates);
	Optimizable_copy(&src->bl, &dst->bl);
	Optimizable_copy(&src->rates, &dst->rates);
	Optimizable_copy(&src->heights, &dst->heights);
	
    dst->topology_optimize  = src->topology_optimize;
    dst->topology_alogrithm = src->topology_alogrithm;
    dst->topology_threads   = src->topology_threads;
    
	dst->max_rounds = src->max_rounds;
	dst->precision  = src->precision;
	dst->verbosity  = src->verbosity;
	dst->interruptible  = src->interruptible;
}





/* Work out the Hessian for the likelihood function. Only branch lengths are considered as variable.
 This function is very much inspired from Jeff Thorne's 'hessian' function in his program 'estbranches'. */

//double * Hessian2( TreeLikelihood *tlk ) {
//	
//	Tree *tree = tlk->list[0]->tree;
//	
//	double zero_zero;
//	int n_ok_edges;
//	int i,j;
//
//	double lk;
//	double lnL1,lnL2;
//	
//	int dim = Tree_node_count(tree)-1;// 2*tree->n_otu-3;
//	double eps = 0.01;
//	
//	double *hessian = dvector(dim*dim);
//	
//	double *ori_bl      = dvector(dim);
//	double *plus_plus   = dvector(dim*dim);
//	double *minus_minus = dvector(dim*dim);
//	double *plus_minus  = dvector(dim*dim);
//	double *plus_zero   = dvector(dim);
//	double *minus_zero  = dvector(dim);
//	double *inc         = dvector(dim);
//	double *buff        = dvector(dim*dim);
//
//	int *ok_edges    = ivector(dim);
//	bool *is_ok      = bvector(dim);
//	
//	double lnL = lnL1 = lnL2 = UNLIKELY;
//	
//	//Lk(tree);
//	
//	Node **nodes = Tree_get_nodes(tree, POSTORDER);
//	
//	for ( i = 0; i < dim; i++ ) ori_bl[i] = Node_distance( nodes[i] );
//	
//	
//	n_ok_edges = 0;
//	for ( i = 0; i < dim; i++ ) {
//		if( ori_bl[i]*(1.-eps) > BL_MIN) {	  
//			inc[i] = eps * ori_bl[i];
//			ok_edges[n_ok_edges] = i;
//			n_ok_edges++;
//			is_ok[i] = true;
//		}
//		else {
//			inc[i] = -1.0;
//			is_ok[i] = false;
//		}
//    }
//	
//	/* zero zero */  
//	zero_zero = tlk->calculate(tlk);
//	
//	/* plus zero */  
//	for ( i = 0; i < dim; i++ ) {
//		if ( is_ok[i] ) {
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) + inc[i] );
//			TreeLikelihood_update_one_node(tlk, i);
//			plus_zero[i] = tlk->calculate(tlk); //Lk_At_Given_Edge(tree->t_edges[i],tree);
//			Node_set_distance( nodes[i], ori_bl[i] );
//			TreeLikelihood_update_one_node(tlk, i);
//		}
//    }
//	
//	
//	/* minus zero */  
//	for ( i = 0; i < dim; i++ ) {
//		if ( is_ok[i] ) {
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) - inc[i] );
//			TreeLikelihood_update_one_node(tlk, i);
//			minus_zero[i] = tlk->calculate(tlk); //Lk_At_Given_Edge(tree->t_edges[i],tree);
//			Node_set_distance( nodes[i], ori_bl[i] );
//			TreeLikelihood_update_one_node(tlk, i);
//		}
//    }
//	
//	//for ( i = 0; i < dim; i++ ) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
//	
//	/* plus plus  */  
//	for ( i = 0; i < dim; i++ ) {
//		if( is_ok[i] ) {
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) + inc[i] );
//			TreeLikelihood_update_one_node(tlk, i); //Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
//			
//			for ( int j = 0; j < 3; j++ ){
//				if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
//					Recurr_Hessian(tree->t_edges[i]->left, tree->t_edges[i]->left->v[j], 1, inc, plus_plus+i*dim, is_ok, tree);
//			}
//			
//			for ( int j = 0; j < 3; j++ ){
//				if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
//					Recurr_Hessian(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],1,inc,plus_plus+i*dim,is_ok,tree);
//			}
//			
//			Node_set_distance( nodes[i], ori_bl[i] );
//			Lk(tree);
//		}
//    }
//	
//	/* plus minus */  
//	for ( i = 0; i < dim; i++ ) {
//		if( is_ok[i] ) {
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) + inc[i] );
//			Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
//			
//			for ( int j = 0; j < 3; j++ ){
//				if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
//					Recurr_Hessian(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
//			}
//			
//			for ( int j = 0; j < 3; j++ ){
//				if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
//					Recurr_Hessian(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],-1,inc,plus_minus+i*dim,is_ok,tree);
//			}
//			
//			Node_set_distance( nodes[i], ori_bl[i] );
//			Lk(tree);
//		}
//    }
//	
//	/* minus minus */  
//	for ( i = 0; i < dim; i++ ) {
//		if(is_ok[i]) {
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) - inc[i] );
//			
//			if(tree->t_edges[i]->l < BL_MIN) {
//				PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//				Exit("\n");
//			}
//			
//			Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
//			
//			for ( int j = 0; j < 3; j++ ){
//				if((!tree->t_edges[i]->left->tax) && (tree->t_edges[i]->left->v[j] != tree->t_edges[i]->rght))
//					Recurr_Hessian(tree->t_edges[i]->left,tree->t_edges[i]->left->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
//			}
//			
//			for ( int j = 0; j < 3; j++ ){
//				if((!tree->t_edges[i]->rght->tax) && (tree->t_edges[i]->rght->v[j] != tree->t_edges[i]->left))
//					Recurr_Hessian(tree->t_edges[i]->rght,tree->t_edges[i]->rght->v[j],-1,inc,minus_minus+i*dim,is_ok,tree);
//			}
//			
//			Node_set_distance( nodes[i], ori_bl[i] );
//			Lk(tree);
//		}
//    }
//	
//	
//	for ( i = 0; i < dim; i++ ) {
//		if( is_ok[i] ) {
//			hessian[i*dim+i] = (plus_zero[i]-2*zero_zero+minus_zero[i])/(pow(inc[i],2));
//			
//			for( j = i+1; j < dim; j++ ) {
//				if ( is_ok[j] ) {
//					hessian[i*dim+j] = (plus_plus[i*dim+j]-plus_minus[i*dim+j]-plus_minus[j*dim+i]+minus_minus[i*dim+j]) / (4*inc[i]*inc[j]);
//					hessian[j*dim+i] = hessian[i*dim+j];
//				}
//			}
//		}
//    }
//	
//	for ( i = 0; i < n_ok_edges; i++ ) {
//		for ( j = 0; j < n_ok_edges; j++ ) {
//			buff[i*n_ok_edges+j] = -1.0*hessian[ok_edges[i]*dim+ok_edges[j]];
//		}
//    }
//	
//	/*if( !Matinv(buff,n_ok_edges,n_ok_edges) ) {
//		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//		Exit("\n");      
//    }*/
//	
//	for ( i = 0; i < n_ok_edges; i++ ) {
//		for ( j = 0; j < n_ok_edges; j++ ) {
//			hessian[ ok_edges[i]*dim+ok_edges[j] ] = buff[i*n_ok_edges+j];
//		}
//    }
//	
//	/* Approximate variance for very short branches */
//	for ( i = 0; i < dim; i++ ) {
//		if( inc[i] < 0.0 ) {
//			lnL  = tree->c_lnL;
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) + eps );
//			lnL1 = Lk_At_Given_Edge(tree->t_edges[i],tree);
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) + eps );
//			lnL2 = Lk_At_Given_Edge(tree->t_edges[i],tree);
//			
//			hessian[i*dim+i] = (lnL2 - 2*lnL1 + lnL) / pow(eps,2);
//			hessian[i*dim+i] = -1.0 / hessian[i*dim+i];
//		}
//	}
//	
//	for ( i = 0; i < dim; i++ ) {
//		if ( hessian[i*dim+i] < MIN_VAR_BL ) {
//			lnL  = tree->c_lnL;
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) + eps );
//			lnL1 = Lk_At_Given_Edge(tree->t_edges[i],tree);
//			Node_set_distance( nodes[i], Node_distance( nodes[i] ) + eps );
//			lnL2 = Lk_At_Given_Edge(tree->t_edges[i],tree);
//			
//			hessian[i*dim+i] = (lnL2 - 2*lnL1 + lnL) / pow(eps,2);
//			hessian[i*dim+i] = -1.0 / hessian[i*dim+i];
//		}
//	}
//	
//	for ( i = 0; i < dim; i++ ) {
//		if ( hessian[i*dim+i] < 0.0 ) {
//			PhyML_Printf("\n. l=%G var=%G", Node_distance( nodes[i] ), hessian[i*dim+i] );
//			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//			Exit("\n");
//		}
//	}
//	
//	for ( i = 0; i < dim; i++ ) {
//		if( hessian[i*dim+i] < MIN_VAR_BL ) {
//			PhyML_Printf("\n. l=%G var=%G", Node_distance( nodes[i] ), hessian[i*dim+i]);
//			hessian[i*dim+i] = MIN_VAR_BL;
//			PhyML_Printf("\n. Numerical precision issues may alter this analysis...");
//		}
//	}
//	
//	
//	/*if( !Matinv(hessian,dim,dim) ){
//		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//		Exit("\n");
//    }*/
//	
//	for ( i = 0; i < dim*dim; i++ ) hessian[i] *= -1.0;
//	
//	for ( i = 0; i < dim; i++ ) {
//		for ( j = 0; j < dim; j++ ) {
//			if(fabs(hessian[i*dim+j]-hessian[j*dim+i]) > 1.E-3) {
//				PhyML_Printf("\n. Hessian not symmetrical.");
//				PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//				Exit("\n");
//			}
//			hessian[i*dim+j] = (hessian[i*dim+j] + hessian[j*dim+i]) / 2.; 
//			hessian[j*dim+i] = hessian[i*dim+j];  
//		}
//    }
//	
//		
//	free(ori_bl);
//	free(plus_plus);
//	free(minus_minus);
//	free(plus_zero);
//	free(minus_zero);
//	free(plus_minus);
//	free(inc);
//	free(buff);
//	free(ok_edges);
//	free(is_ok);
//	
//	return hessian;
//	
//}

/*char * SingleTreeLikelihood_stringify( SingleTreeLikelihood *stlk ){
	StringBuffer *buffer = new_StringBuffer(1000);
	
	buffer = SingleTreeLikelihood_bufferize( buffer, stlk);
	
	char *final = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return final;
 }
 
 
 StringBuffer * SingleTreeLikelihood_bufferize( StringBuffer *buffer, SingleTreeLikelihood *stlk ){
	buffer = StringBuffer_append_format(buffer, "(SingleTreeLikelihood:\n(id:\"stlk%d\")\n", stlk->id );
	buffer = Tree_Bufferize( buffer, stlk->tree );
	buffer = StringBuffer_append_char(buffer, '\n');
	
	buffer = SiteModel_bufferize( buffer, stlk->sm );
	buffer = StringBuffer_append_char(buffer, '\n');
	
	buffer = BranchModel_bufferize( buffer, stlk->bm );
	buffer = StringBuffer_append_char(buffer, '\n');
	
	buffer = SitePattern_bufferize( buffer, stlk->sp );
	buffer = StringBuffer_append_char(buffer, '\n');
	
	buffer = StringBuffer_append_char(buffer, ')');
	
	return buffer;
 }
 
 
 void * SingleTreeLikelihood_SML_to_object( ObjectStore *store, SMLNode node ){
	fprintf(stderr, "SingleTreeLikelihood_SML_to_object\n");
	
	
	SMLNode sm_n = SML_get_element( node, "SiteModel");
	SiteModel    *sm = SiteModel_SML_to_object( store, sm_n);
	
	SitePattern *sp = NULL;
	fprintf(stderr, "SingleTreeLikelihood_SML_to_object:SitePattern\n");
	SMLNode sp_n = SML_get_element( node, "SitePattern");
	if( sp_n != NULL ){
 sp = SitePattern_SML_to_object(store, sp_n);
	}
 else {
 error("Need a sitepattern\n");
 }
	
	BranchModel *bm = NULL;
	SMLNode bm_n = SML_get_element( node, "BranchModel");
	if( bm_n != NULL ){
 bm = BranchModel_SML_to_object( store, bm_n);
	}
	
	
	SMLNode tree_n = SML_get_element( node, "Tree");
	Tree *tree = Tree_SML_to_object(store, tree_n);
	
	SingleTreeLikelihood *stlk = new_SingleTreeLikelihood(tree, sm, sp, bm);
	
	return stlk;
 }*/
