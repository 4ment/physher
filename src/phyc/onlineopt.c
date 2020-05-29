//
//  onlineopt.c
//  physher
//
//  Created by mathieu on 28/5/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "onlineopt.h"

#include "treelikelihood.h"
#include "compoundmodel.h"

double _calculate_insertion_at_node(OnlineOptimizer* online, Node* focal, int newTaxonIndex, int newNodeIndex){
	SingleTreeLikelihood* tlk = online->treelikelihood->obj;
	tlk->new_taxon_bl = online->new_taxon_bl;
	tlk->new_node_id = newNodeIndex;
	tlk->new_taxon_id = newTaxonIndex;
	tlk->node_upper = focal;
	return tlk->calculate(tlk); // Model should be not updated at this stage
}

void exhaustive_online_unrooted(OnlineOptimizer* online){
	SingleTreeLikelihood* tlk = online->treelikelihood->obj;
	Model** models = online->treelikelihood->data;
	Tree* tree = tlk->tree;
	
	for (size_t index = 0; index < online->indexes_count; index++) {
		int tipCount = Tree_tip_count(tree);
		int nodeCount = Tree_node_count(tree);
		char* taxonName = tlk->sp->names[online->indexes[index]];
		//update everything but should it be the whole model like compound
		online->treelikelihood->logP(online->treelikelihood);
		tlk->mapping[tipCount] = get_sequence_index(tlk->sp, taxonName);
		Node** nodes = Tree_nodes(tree);// topology should update and so is this array
		Node* bestFocal = NULL;
		double bestScore = -INFINITY;
		for (size_t i = 0; i < nodeCount; i++) {
			if(!Node_isroot(nodes[i])){
				double score = _calculate_insertion_at_node(online, nodes[i], tipCount, nodeCount);
				if (score > bestScore) {
					bestScore = score;
					bestFocal = nodes[i];
				}
			}
		}
		TreeModel_insert_taxon(models[0], bestFocal, taxonName);
		tlk->new_taxon_id = tlk->new_node_id = -1;
		tlk->use_upper = false;
		tlk->node_upper = NULL;
	}
}

OnlineOptimizer * new_OnlineOptimizer( Model *model, Model *mlikelihood ){
    OnlineOptimizer *opt = (OnlineOptimizer*)malloc(sizeof(OnlineOptimizer));
	opt->model = model;
	opt->treelikelihood = mlikelihood;
	opt->indexes = NULL;
	opt->indexes_count = 0;
	opt->optimize = false;
	opt->new_taxon_bl = 0.1;
    return opt;
}

void _free_OnlineOptimizer( OnlineOptimizer *opt ){
    if( opt->indexes != NULL )free(opt->indexes);
	opt->treelikelihood->free(opt->treelikelihood);
	opt->model->free(opt->model);
    free(opt);
}


double _OnlineOptimizer_optimize(OnlineOptimizer* online){
	exhaustive_online_unrooted(online);
	
	double logP = online->treelikelihood->logP(online->treelikelihood);
	return logP;
}

OnlineOptimizer* new_OnlineOptimizer_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"branchlength",
		"criterion",// exhaustive...
		"model",
		"treelikelihood",
		"verbosity"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* criterion = get_json_node_value_string(node, "criterion");
	int verbosity = get_json_node_value_int(node, "verbosity", 1);
	double default_branch_length = get_json_node_value_double(node, "branchlength", 0.1);
	json_node* tlk_node = get_json_node(node, "treelikelihood"); // treelikelihood with the tree
	json_node* model_node = get_json_node(node, "model"); // model to optimize
	
	if(tlk_node == NULL){
		fprintf(stderr, "treelikleihood must be provided for optimizing in online inference");
		exit(2);
	}
	Model* likelihood = NULL;
	if (tlk_node->node_type == MJSON_OBJECT) {
		likelihood = new_TreeLikelihoodModel_from_json(tlk_node, hash);
		char* id = get_json_node_value_string(tlk_node, "id");
		Hashtable_add(hash, id, likelihood);
	}
	else if(tlk_node->node_type == MJSON_STRING){
		char* ref = (char*)tlk_node->value;
		likelihood = Hashtable_get(hash, ref+1);
		likelihood->ref_count++;
	}
	else{
		exit(10);
	}
	
	
	Model* model = NULL; // full model (could be treelikelihood)
	if (model_node->node_type == MJSON_OBJECT) {
		json_node* type_node = get_json_node(model_node, "type");
		char* id = get_json_node_value_string(model_node, "id");
		
		if (strcasecmp((char*)type_node->value, "compound") == 0) {
			model = new_CompoundModel_from_json(model_node, hash);
		}
		else if(strcasecmp((char*)type_node->value, "treelikelihood") == 0){
			model = new_TreeLikelihoodModel_from_json(model_node, hash);
		}
		else{
			exit(10);
		}
		Hashtable_add(hash, id, model);
	}
	else if(model_node->node_type == MJSON_STRING){
		char* ref = (char*)model_node->value;
		model = Hashtable_get(hash, ref+1);
		model->ref_count++;
	}
	else{
		exit(10);
	}
	
	OnlineOptimizer *opt = new_OnlineOptimizer(model, likelihood);
	opt->indexes = NULL;
	opt->indexes_count = 0;
	opt->new_taxon_bl = default_branch_length;
	opt->optimize = _OnlineOptimizer_optimize;
	opt->free = _free_OnlineOptimizer;
	
	return opt;
}
