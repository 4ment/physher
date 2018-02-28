//
//  operator.c
//  physher
//
//  Created by Mathieu Fourment on 4/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "operator.h"

#include <string.h>
#include <strings.h>

#include "simplex.h"
#include "random.h"
#include "gamma.h"
#include "dirichlet.h"
#include "tree.h"
#include "matrix.h"
#include "treesearch.h"
#include "opvb.h"
#include "ophmc.h"

double tune(int accepted, int count, double target, double tuner, double min, double max, bool inverse){
	double delta = 1./sqrt(count);
	//delta = fmin(0.01, delta);
	double logTuner = log(tuner);
	double accRatio = (double)accepted/count;
	if (inverse) {
		delta = -delta;
	}
	if(accRatio > target){
		logTuner += delta;
	}
	else {
		logTuner -= delta;
	}
	double t = exp(logTuner);
	
	if (t <= min || t >= max)
		t = tuner;
	return t;
}

double getDeltaP(int count, double logAlpha, double target){
	return ((1.0 / count) * (exp(fmin(logAlpha, 0.)) - target));
}

double optimizeScaleFactor(double scaleFactor, int count, double logAlpha, double target ){
	double delta = getDeltaP(count, logAlpha, target);
	delta += log(1./scaleFactor-1.0);
	return 1./(exp(delta)+1);
}

void operator_store(Operator* op){
	Parameters* ps = op->x;
	if(op->models!= NULL){
		ps = ((Simplex*)op->models[0]->obj)->parameters;
	}
	for(int i = 0; i < Parameters_count(ps); i++){
		Parameter_store(Parameters_at(ps, i));
	}
}

void operator_snni_store(Operator* op){
	
}

void operator_restore(Operator* op){
	Parameters* ps = op->x;
	if(op->models != NULL){
		ps = ((Simplex*)op->models[0]->obj)->parameters;
	}
	for(int i = 0; i < Parameters_count(ps); i++){
		Parameter_restore(Parameters_at(ps, i));
	}
}

void operator_snni_restore(Operator* op){
	Tree* tree = op->models[0]->obj;
	Node* node1 = Tree_node(tree, op->indexes[0]);
	Node* node2 = Tree_node(tree, op->indexes[1]);
	NNI_move(tree, node1, node2);
	node1->distance->listeners->fire(node1->distance->listeners, NULL, Node_id(node1));
	node2->distance->listeners->fire(node2->distance->listeners, NULL, Node_id(node2));
}


bool operator_scaler(Operator* op, double* logHR){
	op->index = random_int(Parameters_count(op->x)-1);
	Parameter* p = Parameters_at(op->x, op->index);
	double scaler_scaleFactor = op->parameters[0];
	double v = Parameter_value(p);
	double vv = v;
	double s = (scaler_scaleFactor + (random_double() * ((1.0 / scaler_scaleFactor) - scaler_scaleFactor)));
	vv *= s;
	*logHR = -log(s);
	if ( vv > Parameter_upper(p) || vv < Parameter_lower(p ) ) {
		op->rejected_count++;
		return false;
	}
	Parameter_set_value(p, vv);
	return true;
}

bool operator_slider(Operator* op, double* logHR){
	op->index = random_int(Parameters_count(op->x)-1);
	Parameter* p = Parameters_at(op->x, op->index);
	double v = Parameter_value(p);
	double vv = v;
	double w = (random_double() - 0.5)*op->parameters[0];//slider_delta;
	vv += w;
	//	printf("%s %f %f ", Parameter_name(p), w, Parameter_value(p));
	bool ok = false;
	do{
		if ( vv > Parameter_upper(p) ) {
			vv = 2.0*Parameter_upper(p) - vv;
		}
		else if ( vv < Parameter_lower(p ) ) {
			vv = 2.0*Parameter_lower(p) - vv;
		}
		else{
			ok = true;
		}
	}while(!ok);
	
	*logHR = 0;
	
	Parameter_set_value(p, vv);
	return true;
}

bool operator_dirichlet(Operator* op, double* logHR){
	Simplex* simplex = op->models[0]->obj;
	double alpha = op->parameters[0];
	double sum = 0;
	double* scaledOld = dvector(simplex->K);
	double* newScaled = dvector(simplex->K);
	double* newValues = dvector(simplex->K);
	
	for (int i = 0; i < simplex->K; i++) {
		scaledOld[i] = alpha*simplex->values[i];
	}
	
	rdirichlet(newValues, simplex->K, scaledOld);
	
	for (int i = 0; i < simplex->K; i++) {
		newScaled[i] = alpha*newValues[i];
	}
	
	for (int i = 0; i < simplex->K; i++) {
		if (isnan(newValues[i])) {
			free(scaledOld);
			free(newValues);
			free(newScaled);
			op->rejected_count++;
			return false;
		}
	}
	simplex->set_values(simplex, newValues);
	
	double f = ddirchletln(newValues, simplex->K, scaledOld);
	double b = ddirchletln(simplex->values, simplex->K, newScaled);
	*logHR = b-f;
	free(scaledOld);
	free(newValues);
	free(newScaled);
	return true;
}

// stochastic NNI
bool operator_sNNI(Operator* op, double* logHR){
	Tree* tree = op->models[0]->obj;
	Node* root = Tree_root(tree);
	Node* left_root = Tree_root(tree)->left;
	Node* right_root = Tree_root(tree)->right;
	int index;
	Node* node;
	// right is a leaf so left is not but we do not do nni around left
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
	Node* sibling = Node_sibling(node);
	
	if( Node_isroot(Node_parent(node)) ){
		sibling = Node_left(sibling);
	}
	
	Node* node_swap = node->left;
	int indexSwap = random_int(1);
	
	if(indexSwap == 1){
		node_swap = node->right;
	}
	
	NNI_move(tree, sibling, node_swap);
	
	op->indexes[0] = Node_id(sibling);
	op->indexes[1] = Node_id(node_swap);
	sibling->distance->listeners->fire(sibling->distance->listeners, NULL, Node_id(sibling));
	node_swap->distance->listeners->fire(sibling->distance->listeners, NULL, Node_id(node_swap));
	*logHR = 0;
	return true;
}

void operator_scaler_optimize(Operator* op, double alpha){
	op->parameters[0] = optimizeScaleFactor(op->parameters[0], op->accepted_count+op->rejected_count+1, alpha, 0.24);
}

void operator_slider_optimize(Operator* op, double alpha){
	op->parameters[0] = tune(op->accepted_count, op->accepted_count+op->rejected_count, 0.24, op->parameters[0], 0.0001, 20, false);
}

void operator_dirichlet_optimize(Operator* op, double alpha){
	op->parameters[0] = tune(op->accepted_count, op->accepted_count+op->rejected_count, 0.24, op->parameters[0], 0.1, 10000, true);
}

Operator* new_Operator_from_json(json_node* node, Hashtable* hash){
	
	const char* algorithm_string = get_json_node_value_string(node, "algorithm");
	if (strcasecmp(algorithm_string, "vb") == 0) {
		return new_VariationalOperator_from_json(node, hash);
	}
	else if (strcasecmp(algorithm_string, "hmc") == 0) {
		return new_HMCOperator_from_json(node, hash);
	}
	
	Operator* op = malloc(sizeof(Operator));
	const char* id_string = get_json_node_value_string(node, "id");
	const char* x_string = get_json_node_value_string(node, "x");
	op->x = NULL;
	op->weight = get_json_node_value_double(node, "weight", 1);
	op->name = String_clone(id_string);
	op->store = operator_store;
	op->restore = operator_restore;
	op->parameters = NULL;
	op->indexes = NULL;
	op->models = NULL;
	
	if (strcasecmp(algorithm_string, "scaler") == 0) {
		op->x = new_Parameters(1);
		get_parameters_references2(node, hash, op->x, "x");
		op->propose = operator_scaler;
		op->optimize = operator_scaler_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 0.9;
		
	}
	else if (strcasecmp(algorithm_string, "slider") == 0) {
		op->x = new_Parameters(1);
		get_parameters_references2(node, hash, op->x, "x");
		op->propose = operator_slider;
		op->optimize = operator_slider_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 0.001;
	}
	else if (strcasecmp(algorithm_string, "dirichlet") == 0) {
		char* ref = get_json_node_value_string(node, "x");
		op->model_count = 1;
		op->models = malloc(op->model_count*sizeof(Model*));
		op->models[0] = Hashtable_get(hash, ref+1);
		op->models[0]->ref_count++;
		op->propose = operator_dirichlet;
		op->optimize = operator_dirichlet_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 1000;
	}
	else if (strcasecmp(algorithm_string, "nni") == 0) {
		char* ref = get_json_node_value_string(node, "x");
		op->model_count = 1;
		op->models = malloc(op->model_count*sizeof(Model*));
		op->models[0] = Hashtable_get(hash, ref+1);
		op->models[0]->ref_count++;
		op->propose = operator_sNNI;
		op->store = operator_snni_store;
		op->restore = operator_snni_restore;
		op->optimize = NULL;
		op->parameters = dvector(1);
		op->indexes = ivector(2);
		op->parameters[0] = 1000;
	}
	
	op->rejected_count = 0;
	op->accepted_count = 0;
	op->failure_count = 0;
	return op;
}

void free_Operator(Operator* op){
	if(op->parameters != NULL) free(op->parameters);
	if(op->models != NULL){
		for (int i = 0; i < op->model_count; i++) {
			op->models[i]->free(op->models[i]);
		}
		free(op->models);
	}
	if(op->indexes != NULL) free(op->indexes);
	free_Parameters(op->x);
	free(op->name);
	free(op);
}
