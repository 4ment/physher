//
//  test_clone.c
//  physher
//  Created by Mathieu Fourment on 13/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>
#include <strings.h>

#include "minunit.h"

#include "phyc/treelikelihood.h"
#include "phyc/compoundmodel.h"
#include "phyc/parameters.h"
#include "phyc/mjson.h"
#include "phyc/hashtable.h"
#include "phyc/filereader.h"
#include "phyc/vb.h"
#include "phyc/random.h"
#include "phyc/optimizer.h"

char* test_jc69_clone(){
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	// for cloning
	Hashtable* hash2 = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	char* content = load_file("jc69.json");
	json_node* json = create_json_tree(content);
	free(content);
	
	json_node* child = json->children[0];
	Model* model = new_TreeLikelihoodModel_from_json(child, hash);
	double logP = model->logP(model);
	Model* clone = model->clone(model, hash2);
	model->free(model);
	double logP2 = clone->logP(clone);
	mu_assert(logP == logP2, "logP not matching");
	clone->free(clone);
	free_Hashtable(hash);
	free_Hashtable(hash2);
	json_free_tree(json);
	
	return NULL;
}

char* test_f81_clone(){
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	// for cloning
	Hashtable* hash2 = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	char* content = load_file("f81.json");
	json_node* json = create_json_tree(content);
	free(content);
	
	json_node* child = json->children[0];
	Model* model = new_TreeLikelihoodModel_from_json(child, hash);
	double logP = model->logP(model);
	Model* clone = model->clone(model, hash2);
	model->free(model);
	double logP2 = clone->logP(clone);
	mu_assert(logP == logP2, "logP not matching");
	clone->free(clone);
	free_Hashtable(hash);
	free_Hashtable(hash2);
	json_free_tree(json);
	
	return NULL;
}

char* test_hky_clone(){
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	// for cloning
	Hashtable* hash2 = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	char* content = load_file("hky.json");
	json_node* json = create_json_tree(content);
	free(content);
	
	json_node* child = json->children[0];
	Model* model = new_TreeLikelihoodModel_from_json(child, hash);
	double logP = model->logP(model);
	Model* clone = model->clone(model, hash2);
	model->free(model);
	double logP2 = clone->logP(clone);
	mu_assert(logP == logP2, "logP not matching");
	clone->free(clone);
	free_Hashtable(hash);
	free_Hashtable(hash2);
	json_free_tree(json);
	
	return NULL;
}

char* test_gtr_clone(){
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	Hashtable* hash2 = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	char* content = load_file("gtr-c2-priors.json");
	json_node* json = create_json_tree(content);
	free(content);
	
	json_node* child = json->children[0];
	Model* model = new_CompoundModel_from_json(child, hash);
	char* id = get_json_node_value_string(child, "id");
	Hashtable_add(hash, id, model);
	double logP = model->logP(model);
	
	Model* clone = model->clone(model, hash2);
	
	json_node* run_node = get_json_node(json, "physher");
	child = run_node->children[0];
	Optimizer* opt = new_Optimizer_from_json(child, hash);
	double temp;
	opt_optimize(opt, NULL, &temp);
	free_Optimizer(opt);
	model->free(model);
	
	double logP2 = clone->logP(clone);
//	fprintf(stderr, "logP %f %f\n", logP, logP2);
	mu_assert(logP == logP2, "logP not matching");
	clone->free(clone);
	free_Hashtable(hash);
	free_Hashtable(hash2);
	json_free_tree(json);
	
	return NULL;
}

char* test_variational_clone(){
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	init_genrand(1);
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(r, 1);
	char* rand_key = "RANDOM_GENERATOR!@";
	Hashtable_add(hash, rand_key, r);
	
	// for cloning
	Hashtable* hash2 = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	char* content = load_file("jc69-vb.json");
	json_node* json = create_json_tree(content);
	free(content);
	
	json_node* child = json->children[0];
	Model* posterior = new_CompoundModel_from_json(child, hash);
	Hashtable_add(hash, get_json_node_value_string(child, "id"), posterior);
	double logP = posterior->logP(posterior);
	
	json_node* child1 = json->children[1];
	Model* var = new_Variational_from_json(child1, hash);
	posterior->free(posterior);
	
	
	
	
	Model* clone = var->clone(var, hash2);
	
	variational_t* clone_var = clone->obj;
	Model* clone_posterior = clone_var->posterior;
	double logP2 = clone_posterior->logP(clone_posterior);
	mu_assert(logP == logP2, "logP not matching");
	
	init_genrand(1);
	gsl_rng_set(r, 1);
	double logQ = var->logP(var);
	var->free(var);
	
	init_genrand(1);
	gsl_rng_set(r, 1);
	double logQ2 = clone->logP(clone);
	fprintf(stderr, "%f %f\n", logQ, logQ2);
// 	mu_assert(logQ == logQ2, "logQ not matching");
	
	clone->free(clone);
	free_Hashtable(hash);
	free_Hashtable(hash2);
	json_free_tree(json);
	gsl_rng_free(r);
	
	return NULL;
}

char* test_variational_blocks_clone(){
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	init_genrand(1);
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(r, 1);
	char* rand_key = "RANDOM_GENERATOR!@";
	Hashtable_add(hash, rand_key, r);
	
	// for cloning
	Hashtable* hash2 = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	char* content = load_file("jc69-vb-blocks.json");
	json_node* json = create_json_tree(content);
	free(content);
	
	json_node* child = json->children[0];
	Model* posterior = new_CompoundModel_from_json(child, hash);
	Hashtable_add(hash, get_json_node_value_string(child, "id"), posterior);
	double logP = posterior->logP(posterior);
	
	json_node* child1 = json->children[1];
	Model* var = new_Variational_from_json(child1, hash);
	posterior->free(posterior);
	
	
	
	
	Model* clone = var->clone(var, hash2);
	
	variational_t* clone_var = clone->obj;
	Model* clone_posterior = clone_var->posterior;
	double logP2 = clone_posterior->logP(clone_posterior);
	mu_assert(logP == logP2, "logP not matching");
	
	init_genrand(1);
	gsl_rng_set(r, 1);
	double logQ = var->logP(var);
	var->free(var);
	
	init_genrand(1);
	gsl_rng_set(r, 1);
	double logQ2 = clone->logP(clone);
	fprintf(stderr, "%f %f\n", logQ, logQ2);
// 	mu_assert(logQ == logQ2, "logQ not matching");
	
	clone->free(clone);
	free_Hashtable(hash);
	free_Hashtable(hash2);
	json_free_tree(json);
	gsl_rng_free(r);
	
	return NULL;
}

char *all_tests(){
	mu_suite_start();
	mu_run_test(test_jc69_clone);
	mu_run_test(test_f81_clone);
	mu_run_test(test_hky_clone);
	mu_run_test(test_gtr_clone);
	mu_run_test(test_variational_clone);
	mu_run_test(test_variational_blocks_clone);
 
     return NULL;
 }
  
RUN_TESTS(all_tests);
