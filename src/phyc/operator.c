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
#include "gamma.h"
#include "dirichlet.h"
#include "tree.h"
#include "matrix.h"
#include "treesearch.h"
// #include "opvb.h"
#include "ophmc.h"
#include "discreteoperator.h"
#include "demographicmodels.h"

#include <gsl/gsl_randist.h>

double tune(long accepted, long count, double target, double tuner, double min, double max, bool inverse){
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

double getDeltaP(long count, double logAlpha, double target){
	return ((1.0 / count) * (exp(fmin(logAlpha, 0.)) - target));
}

double optimizeScaleFactor(double scaleFactor, long count, double logAlpha, double target ){
	double delta = getDeltaP(count, logAlpha, target);
	delta += log(1./scaleFactor-1.0);
	return 1./(exp(delta)+1.0);
}

bool operator_uniform_height(Operator* op, double* logHR){
	Tree* tree = op->models[0]->obj;
	Node* internal = Tree_root(tree);
	int nodeCount = Tree_node_count(tree);
	while (Node_isroot(internal) || Node_isleaf(internal)) {
		size_t idx = gsl_rng_uniform_int(op->rng, nodeCount);
		internal = Tree_node(tree, idx);
	}
	double lower = fmax(Node_height(Node_left(internal)), Node_height(Node_right(internal)));
	double upper = Node_height(Node_parent(internal));
	double newValue = (gsl_rng_uniform(op->rng) * (upper - lower)) + lower;
	Node_set_height(internal, newValue);
	*logHR = 0;
	return true;
}

bool operator_uniform(Operator* op, double* logHR){
	size_t index = gsl_rng_uniform_int(op->rng, Parameters_count(op->x));
	Parameter* p = Parameters_at(op->x, index);
	double lower = Parameter_lower(p);
	double upper = Parameter_upper(p);
	double newValue = (gsl_rng_uniform(op->rng) * (upper - lower)) + lower;
	Parameter_set_value(p, newValue);
	*logHR = 0;
	return true;
}

bool operator_interval_scaler(Operator* op, double* logHR){
	double scaler_scaleFactor = op->parameters[0];
	double s = (scaler_scaleFactor + (gsl_rng_uniform(op->rng) * ((1.0 / scaler_scaleFactor) - scaler_scaleFactor)));
	
	Coalescent* coal = op->models[0]->obj;
	Node* internal = Tree_root(coal->tree);
	int nodeCount = Tree_node_count(coal->tree);
	while (Node_isleaf(internal)) {
		size_t idx = gsl_rng_uniform_int(op->rng, nodeCount);
		internal = Tree_node(coal->tree, idx);
	}
	int endPoint = 0;
	for( ; endPoint < coal->n; endPoint++){
		if(coal->nodes[endPoint]->index == Node_id(internal)) break;
	}
	double oldInterval = Node_height(internal);
	if (endPoint > 0) {
		oldInterval -= Node_height(Tree_node(coal->tree, coal->nodes[endPoint-1]->index));
	}
	double intervalIncrement = s*oldInterval - oldInterval;
	
	for(int i = endPoint; i < coal->n; i++){
		Node* node = Tree_node(coal->tree, coal->nodes[i]->index);
		Node_set_height(node, Node_height(node)+intervalIncrement);
	}

	*logHR = -log(s);
	return true;
}

bool operator_scaler(Operator* op, double* logHR){
	if(op->x != NULL){
		size_t index = 0;
		if(Parameters_count(op->x) != 1){
			index = gsl_rng_uniform_int(op->rng, Parameters_count(op->x));
		}
		Parameter* p = Parameters_at(op->x, index);
		double scaler_scaleFactor = op->parameters[0];
		double v = Parameter_value(p);
		double vv = v;
		double s = (scaler_scaleFactor + (gsl_rng_uniform(op->rng)  * ((1.0 / scaler_scaleFactor) - scaler_scaleFactor)));
		vv *= s;
		*logHR = -log(s);
		
		if ( vv > Parameter_upper(p) || vv < Parameter_lower(p ) ) {
			op->rejected_count++;
			return false;
		}
		Parameter_set_value(p, vv);
	}
	else{
		Tree* tree = op->models[0]->obj;
		double scaler_scaleFactor = op->parameters[0];
		double s = (scaler_scaleFactor + (gsl_rng_uniform(op->rng) * ((1.0 / scaler_scaleFactor) - scaler_scaleFactor)));
		
		if (op->all) {
			Node** nodes = Tree_nodes(tree);
			int nodeCount = Tree_node_count(tree);
			for(int i = 0; i < nodeCount; i++){
				if(!Node_isleaf(nodes[i])){
					double newValue = s*Node_height(nodes[i]);
					Node_set_height(nodes[i], newValue);
				}
			}
			int internalNodes = nodeCount - Tree_tip_count(tree);
			*logHR = log(s)*(internalNodes - 2);
		}
		else{
			Node* root = Tree_root(tree);
			double newValue = s*Node_height(root);
			double lower = fmax(Node_height(Node_left(root)), Node_height(Node_right(root)));
			if(newValue < lower){
				op->rejected_count++;
				return false;
			}
			Node_set_height(root, newValue);
			*logHR = -log(s);
		}
	}
	return true;
}

bool operator_slider(Operator* op, double* logHR){
	size_t index = gsl_rng_uniform_int(op->rng, Parameters_count(op->x));
	Parameter* p = Parameters_at(op->x, index);
	double v = Parameter_value(p);
	double vv = v;
	double w = (gsl_rng_uniform(op->rng) - 0.5)*op->parameters[0];//slider_delta;
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


bool operator_random_walk_unif(Operator* op, double* logHR){
	size_t index = gsl_rng_uniform_int(op->rng, Parameters_count(op->x));
	Parameter* p = Parameters_at(op->x, index);
	double v = Parameter_value(p);
	double vv = v;
	double w = (gsl_rng_uniform(op->rng) * 2.0 - 1.0)*op->parameters[0];
	vv += w;
	
	if ( vv > Parameter_upper(p) || vv < Parameter_lower(p) ) {
		return false;
	}
	
	*logHR = 0;
	
	Parameter_set_value(p, vv);
	return true;
}

bool operator_simplex_exchange(Operator* op, double* logHR){
	Simplex* simplex = op->models[0]->obj;
	double lambda = op->parameters[0];
	long idx1 = gsl_rng_uniform_int(op->rng, simplex->K);
	long idx2 = idx1;
	while(idx1 == idx2){
		idx2 = gsl_rng_uniform_int(op->rng, simplex->K);
	}
	
	const double* oldValues = simplex->get_values(simplex);
	double* newValues = clone_dvector(oldValues, simplex->K);
	
	double w = gsl_rng_uniform(op->rng)*lambda;
	
	newValues[idx1] += w;
	newValues[idx2] -= w;
	
	for(int i = 0; i < simplex->K; i++){
		if (newValues[i] < 0.0001) {
			free(newValues);
			op->rejected_count++;;
			op->failure_count++;
			return false;
		}
	}
	
	for (int i = 0; i < simplex->K; i++) {
		if (isnan(newValues[i])) {
			free(newValues);
			op->rejected_count++;;
			op->failure_count++;
			return false;
		}
	}
	simplex->set_values(simplex, newValues);
	free(newValues);
	return true;
}

bool operator_beta(Operator* op, double* logHR){
	double alpha = op->parameters[0];
	double v = Parameters_value(op->x, 0);
	double newValue = gsl_ran_beta(op->rng, alpha*v+1.0, alpha*(1.0-v)+1.0);
	if (newValue == 1.0 || newValue == 0.0) {
		return false;
	}
	Parameters_set_value(op->x, 0, newValue);
	*logHR = log(gsl_ran_beta_pdf(v, alpha*newValue+1.0, alpha*(1.0-newValue)+1.0)/gsl_ran_beta_pdf(newValue, alpha*v+1.0, alpha*(1.0-v))+1.0);
	return true;
}

bool operator_dirichlet(Operator* op, double* logHR){
	Simplex* simplex = op->models[0]->obj;
	double alpha = op->parameters[0];
	double* scaledOld = dvector(simplex->K);
	double* newScaled = dvector(simplex->K);
	double* newValues = dvector(simplex->K);
	
	double* oldValues = clone_dvector(simplex->get_values(simplex), simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		scaledOld[i] = alpha*oldValues[i];
	}
	
	bool ok = false;
	int failures = 0;
	while (!ok && failures != 50) {
		gsl_ran_dirichlet(op->rng, simplex->K, scaledOld, newValues); //rdirichlet(newValues, simplex->K, scaledOld);
		for(int i = 0; i < simplex->K; i++){
			ok = true;
			if (newValues[i] < 0.00001) {
				ok = false;
				failures++;
				break;
			}
		}
	}
	if(failures == 50){
		free(oldValues);
		free(scaledOld);
		free(newValues);
		free(newScaled);
		op->rejected_count++;
		op->failure_count++;
		return false;
	}
	
	for (int i = 0; i < simplex->K; i++) {
		newScaled[i] = alpha*newValues[i];
	}
	
	for (int i = 0; i < simplex->K; i++) {
		if (isnan(newValues[i])) {
			free(oldValues);
			free(scaledOld);
			free(newValues);
			free(newScaled);
			op->rejected_count++;
			op->failure_count++;
			return false;
		}
	}
	simplex->set_values(simplex, newValues);

	double f = gsl_ran_dirichlet_lnpdf(simplex->K, scaledOld, newValues);// ddirchletln(newValues, simplex->K, scaledOld);
	double b = gsl_ran_dirichlet_lnpdf(simplex->K, newScaled, oldValues); //ddirchletln(oldValues, simplex->K, newScaled);

	*logHR = b-f;
	free(oldValues);
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
	size_t index;
	Node* node;
	// right is a leaf so left is not but we do not do nni around left
	if(Node_isleaf(right_root)){
		do {
			index = gsl_rng_uniform_int(op->rng, Tree_node_count(tree));
			node = Tree_node(tree, index);
		}while(node == root || node == left_root || Node_isleaf(node));
	}
	else{
		do {
			index = gsl_rng_uniform_int(op->rng, Tree_node_count(tree));
			node = Tree_node(tree, index);
		}while(node == root || node == right_root || Node_isleaf(node));
	}
	Node* sibling = Node_sibling(node);
	
	if( Node_isroot(Node_parent(node)) ){
		sibling = Node_left(sibling);
	}
	
	Node* node_swap = node->left;
	size_t indexSwap = gsl_rng_uniform_int(op->rng, 2);
	
	if(indexSwap == 1){
		node_swap = node->right;
	}
	
	NNI_move(tree, sibling, node_swap);
	
	*logHR = 0;
	return true;
}

void operator_scaler_optimize(Operator* op, double alpha){
	long count = op->accepted_count+op->rejected_count - op->tuning_delay;
	if(count >= 0){
		op->parameters[0] = optimizeScaleFactor(op->parameters[0], count+1, alpha, 0.24);
	}
}

void operator_slider_optimize(Operator* op, double alpha){
	long count = op->accepted_count+op->rejected_count - op->tuning_delay;
	if(count >= 0){
		op->parameters[0] = tune(op->accepted_count, op->accepted_count+op->rejected_count, 0.24, op->parameters[0], 0.0001, 20, false);
	}
}

void operator_dirichlet_optimize(Operator* op, double alpha){
	long count = op->accepted_count+op->rejected_count - op->tuning_delay;
	if(count >= 0){
		op->parameters[0] = tune(op->accepted_count, op->accepted_count+op->rejected_count, 0.24, op->parameters[0], 0.01, 10000, true);
	}
}


void operator_beta_optimize(Operator* op, double alpha){
	long count = op->accepted_count+op->rejected_count - op->tuning_delay;
	if(count >= 0){
		op->parameters[0] = tune(op->accepted_count, op->accepted_count+op->rejected_count, 0.24, op->parameters[0], 0.01, 10000, true);
	}
}
void operator_exchange_optimize(Operator* op, double alpha){
	long count = op->accepted_count+op->rejected_count - op->tuning_delay;
	if(count >= 0){
		double delta = getDeltaP(count+1, alpha, 0.24) + log(op->parameters[0]);
		op->parameters[0] = exp(delta);
	}
}

void operator_discrete_exchange_optimize(Operator* op, double alpha){
	long count = op->accepted_count+op->rejected_count - op->tuning_delay;
	if(count >= 0){
		double delta = getDeltaP(count+1, alpha, 0.24);
		if (delta > DBL_MAX || delta < DBL_MIN) {
			delta = 0;
		}
		delta += log(op->parameters[0]);
		op->parameters[0] = exp(delta);
		op->parameters[0] = fmax(0.5000000001, op->parameters[0]);
	}
}

Operator* new_Operator_from_json(json_node* node, Hashtable* hash){
	
	const char* algorithm_string = get_json_node_value_string(node, "algorithm");
	/*if (strcasecmp(algorithm_string, "vb") == 0) {
		return new_VariationalOperator_from_json(node, hash);
	}
	else*/ if (strcasecmp(algorithm_string, "hmc") == 0) {
		return new_HMCOperator_from_json(node, hash);
	}
	
	char* allowed[] = {
		"algorithm",
		"all",
		"coalescent",
		"delay",
		"parameters",
		"tree",
		"weight",
		"x"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* p_node = get_json_node(node, "parameters");
	
	Operator* op = malloc(sizeof(Operator));
	const char* id_string = get_json_node_value_string(node, "id");
	const char* x_string = get_json_node_value_string(node, "x");
	op->x = NULL;
	op->weight = get_json_node_value_double(node, "weight", 1);
	op->name = String_clone(id_string);
	op->parameters = NULL;
	op->models = NULL;
	op->optimize = NULL;
	op->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	//op->tuning_delay = get_json_node_value_size_t(node, "delay", 10000);
	op->all = get_json_node_value_bool(node, "all", false);
	
	if (strcasecmp(algorithm_string, "scaler") == 0) {
		char* ref_tree = get_json_node_value_string(node, "tree");
		char* ref_coalescent = get_json_node_value_string(node, "coalescent");
		if(ref_coalescent != NULL){
			op->model_count = 1;
			op->models = malloc(op->model_count*sizeof(Model*));
			op->models[0] = Hashtable_get(hash, ref_coalescent+1);
			op->models[0]->ref_count++;
			op->propose = operator_interval_scaler;
		}
		else{
			if(ref_tree != NULL){
				op->model_count = 1;
				op->models = malloc(op->model_count*sizeof(Model*));
				op->models[0] = Hashtable_get(hash, ref_tree+1);
				op->models[0]->ref_count++;
			}
			else{
				op->x = new_Parameters(1);
				get_parameters_references2(node, hash, op->x, "x");
			}
			op->propose = operator_scaler;
		}
		op->optimize = operator_scaler_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 0.9;
		if(p_node != NULL){
			if (p_node->node_type == MJSON_PRIMITIVE){
				op->parameters[0] = get_json_node_value_double(node, "parameters", 0.9);
			}
			else if(p_node->node_type == MJSON_ARRAY){
				json_node* child = p_node->children[0];
				op->parameters[0] = get_json_node_value_double(child, "parameters", 0.9);
				
			}
		}
	}
	else if (strcasecmp(algorithm_string, "slider") == 0) {
		op->x = new_Parameters(1);
		get_parameters_references2(node, hash, op->x, "x");
		op->propose = operator_slider;
		op->optimize = operator_slider_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 0.001;
		if(p_node != NULL){
			if (p_node->node_type == MJSON_PRIMITIVE){
				op->parameters[0] = get_json_node_value_double(node, "parameters", 0.001);
			}
			else if(p_node->node_type == MJSON_ARRAY){
				json_node* child = p_node->children[0];
				op->parameters[0] = get_json_node_value_double(child, "parameters", 0.001);
				
			}
		}
	}
	else if (strcasecmp(algorithm_string, "randomwalk") == 0) {
		op->x = new_Parameters(1);
		get_parameters_references2(node, hash, op->x, "x");
		op->propose = operator_random_walk_unif;
		op->optimize = operator_slider_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = 0.75;
		if(p_node != NULL){
			if (p_node->node_type == MJSON_PRIMITIVE){
				op->parameters[0] = get_json_node_value_double(node, "parameters", 0.75);
			}
			else if(p_node->node_type == MJSON_ARRAY){
				json_node* child = p_node->children[0];
				op->parameters[0] = get_json_node_value_double(child, "parameters", 0.75);
				
			}
		}
	}
	else if (strcasecmp(algorithm_string, "beta") == 0) {
		op->x = new_Parameters(1);
		get_parameters_references2(node, hash, op->x, "x");
		op->propose = operator_beta;
		op->optimize = operator_beta_optimize;
		op->parameters = dvector(1);
		op->parameters[0] = get_json_node_value_double(node, "parameters", 10);
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
		op->parameters[0] = 100;
		if(p_node != NULL){
			if (p_node->node_type == MJSON_PRIMITIVE){
				op->parameters[0] = get_json_node_value_double(node, "parameters", 100);
			}
			else if(p_node->node_type == MJSON_ARRAY){
				json_node* child = p_node->children[0];
				op->parameters[0] = get_json_node_value_double(child, "parameters", 100);
				
			}
		}
	}
	else if (strcasecmp(algorithm_string, "exchange") == 0) {
		char* ref = get_json_node_value_string(node, "x");
		op->model_count = 1;
		op->models = malloc(op->model_count*sizeof(Model*));
		op->models[0] = Hashtable_get(hash, ref+1);
		op->models[0]->ref_count++;
		op->parameters = dvector(1);
		
		if(op->models[0]->type == MODEL_DISCRETE_PARAMETER){
			op->propose = operator_discrete_exchange;
			op->optimize = operator_discrete_exchange_optimize;
			op->parameters[0] = 1;
		}
		else{
			op->propose = operator_simplex_exchange;
			op->optimize = operator_exchange_optimize;
			op->parameters[0] = 0.001;
		}
		if(p_node != NULL){
			if (p_node->node_type == MJSON_PRIMITIVE){
				op->parameters[0] = get_json_node_value_double(node, "parameters", 0.01);
			}
			else if(p_node->node_type == MJSON_ARRAY){
				json_node* child = p_node->children[0];
				op->parameters[0] = get_json_node_value_double(child, "parameters", 0.01);
				
			}
		}
	}
	else if (strcasecmp(algorithm_string, "nni") == 0) {
		char* ref = get_json_node_value_string(node, "x");
		op->model_count = 1;
		op->models = malloc(op->model_count*sizeof(Model*));
		op->models[0] = Hashtable_get(hash, ref+1);
		op->models[0]->ref_count++;
		op->propose = operator_sNNI;
		op->parameters = dvector(1);
		op->parameters[0] = 1000;
		if(p_node != NULL){
			if (p_node->node_type == MJSON_PRIMITIVE){
				op->parameters[0] = get_json_node_value_double(node, "parameters", 1000);
			}
			else if(p_node->node_type == MJSON_ARRAY){
				json_node* child = p_node->children[0];
				op->parameters[0] = get_json_node_value_double(child, "parameters", 1000);
				
			}
		}
	}
	else if (strcasecmp(algorithm_string, "uniform") == 0) {
		char* ref_tree = get_json_node_value_string(node, "tree");
		if (ref_tree == NULL) {
			op->x = new_Parameters(1);
			get_parameters_references2(node, hash, op->x, "x");
			op->propose = operator_uniform;
		}
		else{
			op->model_count = 1;
			op->models = malloc(op->model_count*sizeof(Model*));
			op->models[0] = Hashtable_get(hash, ref_tree+1);
			op->models[0]->ref_count++;
			op->propose = operator_uniform_height;
		}
	}
	else if (strcasecmp(algorithm_string, "bitflip") == 0) {
		char* ref = get_json_node_value_string(node, "x");
		op->model_count = 1;
		op->models = malloc(op->model_count*sizeof(Model*));
		op->models[0] = Hashtable_get(hash, ref+1);
		op->models[0]->ref_count++;
		op->propose = operator_discrete_bitflip;
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
	free_Parameters(op->x);
	free(op->name);
	free(op);
}
