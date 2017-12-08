//
//  opvb.c
//  physher
//
//  Created by Mathieu Fourment on 7/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "opvb.h"

#include "tree.h"
#include "vb.h"
#include "matrix.h"

void operator_vb_store(Operator* op){
	Parameters* ps = op->x;
	for(int i = 0; i < Parameters_count(ps); i++){
		Parameter_store(Parameters_at(ps, i));
	}
}

void operator_vb_restore(Operator* op){
	Tree* tree = op->models[1]->obj;
	Node* node1 = Tree_node(tree, op->indexes[0]);
	Node* node2 = Tree_node(tree, op->indexes[1]);
	Node* node3 = Tree_node(tree, op->indexes[2]);
	Parameter_restore(node1->distance);
	Parameter_restore(node2->distance);
	Parameter_restore(node3->distance);
}

void operator_vb_restore_1(Operator* op){
	Tree* tree = op->models[1]->obj;
	Node* node = Tree_node(tree, op->indexes[0]);
	Parameter_restore(node->distance);
}

bool operator_vb(Operator* op, double* logHR){
	variational_t* var = op->models[0]->obj;
	size_t dim = Parameters_count(var->parameters);
	Tree* tree = op->models[1]->obj;
	Node* root = Tree_root(tree);
	Node* left_root = Tree_root(tree)->left;
	Node* right_root = Tree_root(tree)->right;
	int index;
	Node* node;
	double* values = dvector(dim);
	Parameters* ps = new_Parameters(3);
	if(Node_isleaf(right_root)){
		do {
			index = random_int(Tree_node_count(tree)-1);
			node = Tree_node(tree, index);
		}while(node == root || node == left_root || Node_isleaf(node));
	}
	else{
		do {
			index = random_int(Tree_node_count(tree)-1);
			node = Tree_node(tree, index);
		}while(node == root || node == right_root || Node_isleaf(node));
	}
	Parameters_add(ps, node->distance);
	Parameters_add(ps, Node_left(node)->distance);
	Parameters_add(ps, Node_right(node)->distance);
	
	op->indexes[0] = Node_id(node);
	op->indexes[1] = Node_id(Node_left(node));
	op->indexes[2] = Node_id(Node_right(node));
	
	var->sample_some(var, ps, values);
	
	Parameter_set_value(node->distance, values[0]);
	Parameter_set_value(Node_left(node)->distance, values[1]);
	Parameter_set_value(Node_right(node)->distance, values[2]);
	
//		Parameters_print(var->parameters);
	free_Parameters(ps);
	free(values);
	*logHR = 0;
	return true;
}

bool operator_vb_1(Operator* op, double* logHR){
	variational_t* var = op->models[0]->obj;
	Tree* tree = op->models[1]->obj;
	Node* root = Tree_root(tree);
	Node* right_root = Tree_root(tree)->right;
	int index;
	Node* node;
	double value;
	Parameters* ps = new_Parameters(1);
	do {
		index = random_int(Tree_node_count(tree)-1);
		node = Tree_node(tree, index);
	}while(node == root || node == right_root);
		
	Parameters_add(ps, node->distance);
	
	op->indexes[0] = Node_id(node);
	
	var->sample_some(var, ps, &value);
	
	Parameter_set_value(node->distance, value);
	
	//	Parameters_print(var->parameters);
	free_Parameters(ps);
	*logHR = 0;
	return true;
}

Operator* new_VariationalOperator_from_json(json_node* node, Hashtable* hash){
	Operator* op = malloc(sizeof(Operator));
	const char* id_string = get_json_node_value_string(node, "id");
	const char* algorithm_string = get_json_node_value_string(node, "algorithm");
	const char* x_string = get_json_node_value_string(node, "x");
	op->x = NULL;
	op->weight = get_json_node_value_double(node, "weight", 1);
	op->name = String_clone(id_string);
	op->parameters = NULL;

	char* ref = get_json_node_value_string(node, "x");
	char* var_ref = get_json_node_value_string(node, "var");
	op->model_count = 1;
	if (ref != NULL && var_ref != NULL) {
		op->model_count++;
	}
	op->models = malloc(op->model_count*sizeof(Model*));
	// variational model
	op->models[0] = Hashtable_get(hash, var_ref+1);
	op->models[0]->ref_count++;
	// tree model or parameters
	if(ref != NULL){
		op->models[1] = Hashtable_get(hash, ref+1);
		op->models[1]->ref_count++;
		Tree* tree = op->models[1]->obj;
		op->x = new_Parameters(Tree_node_count(tree));
		for (int i = 0; i < Tree_node_count(tree); i++) {
			Parameters_add(op->x, Tree_node(tree, i)->distance);
		}
	}
	op->propose = operator_vb_1;
	op->store = operator_vb_store;
	op->restore = operator_vb_restore_1;
	op->optimize = NULL;
	op->parameters = dvector(1);
	op->parameters[0] = 1000;
	op->indexes = ivector(3);
	
	op->rejected_count = 0;
	op->accepted_count = 0;
	return op;
}
