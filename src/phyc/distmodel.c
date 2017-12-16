//
//  distmodel.c
//  physher
//
//  Created by Mathieu Fourment on 3/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "distmodel.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "tree.h"

#include "exponential.h"
#include "gamma.h"
#include "dirichlet.h"

static void _free_dist(DistributionModel*dm){
	if(dm->x != NULL) free_Parameters(dm->x);
	if(dm->parameters != NULL) free_Parameters(dm->parameters);
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	if(dm->simplex != NULL) free_Simplex(dm->simplex);
	// freeing data is left to the user
	free(dm);
}

static DistributionModel* _clone_dist(DistributionModel* dm){
	DistributionModel* clone = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(clone);
	clone->x = NULL;
	clone->simplex = NULL;
	
	clone->parameters = new_Parameters(Parameters_count(dm->parameters));
	for (int i = 0; i < Parameters_count(dm->parameters); i++) {
		Parameters_move(clone->parameters, clone_Parameter(Parameters_at(dm->parameters, i)));
	}
	if(dm->x != NULL){
		clone->x = new_Parameters(Parameters_count(dm->x));
	}
	else if(dm->simplex != NULL){
		clone->simplex = clone_Simplex(dm->simplex);
	}
	clone->logP = dm->logP;
	clone->dlogP = dm->dlogP;
	clone->clone = dm->clone;
	clone->free = dm->free;
	clone->tempp = NULL;
	clone->tempx = NULL;
	if(dm->tempp != NULL) clone->tempp = clone_dvector(dm->tempp, Parameters_count(dm->parameters));
	if(dm->tempx != NULL) clone->tempx = clone_dvector(dm->tempx, Parameters_count(dm->x));
	return clone;
}

DistributionModel* clone_DistributionModel_with_parameters(DistributionModel* dm, const Parameters* params, const Parameters* x, Simplex* simplex){
	DistributionModel* clone = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(clone);
	clone->parameters = NULL;
	if(Parameters_count(params) > 0){
		clone->parameters = new_Parameters(Parameters_count(params));
		for (int i = 0; i < Parameters_count(params); i++) {
			Parameters_add(clone->parameters, Parameters_at(params, i));
		}
	}
	clone->x = NULL;
	clone->simplex = NULL;
	if(x != NULL){
		clone->x = new_Parameters(Parameters_count(x));
		for (int i = 0; i < Parameters_count(x); i++) {
			Parameters_add(clone->x, Parameters_at(x, i));
		}
	}
	else if(simplex != NULL){
		clone->simplex = simplex;
	}
	clone->logP = dm->logP;
	clone->dlogP = dm->dlogP;
	clone->clone = dm->clone;
	clone->free = dm->free;
	clone->tempp = NULL;
	clone->tempx = NULL;
	if(dm->tempp != NULL) clone->tempp = clone_dvector(dm->tempp, Parameters_count(dm->parameters));
	if(dm->tempx != NULL) clone->tempx = clone_dvector(dm->tempx, Parameters_count(dm->x));
	return clone;
}

DistributionModel* new_DistributionModel(const Parameters* p, const Parameters* x){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
	dm->parameters = NULL;
	if(p != NULL){
		dm->parameters = new_Parameters(Parameters_count(p));
		Parameters_add_parameters(dm->parameters, p);
	}
	dm->x = new_Parameters(Parameters_count(x));
	Parameters_add_parameters(dm->x, x);
	dm->simplex = NULL;
	dm->logP = NULL;
	dm->dlogP = NULL;
	dm->free = _free_dist;
	dm->clone = _clone_dist;
	dm->data = NULL;
	dm->tempx = NULL;
	dm->tempp = NULL;
	return dm;
}

DistributionModel* new_DistributionModelSimplex(Parameters* p, Simplex* simplex){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
	dm->parameters = p;
	dm->x = NULL;
	dm->simplex = simplex;
	dm->logP = NULL;
	dm->dlogP = NULL;
	dm->free = _free_dist;
	dm->clone = _clone_dist;
	dm->data = NULL;
	dm->tempx = NULL;
	dm->tempp = NULL;
	return dm;
}

double DistributionModel_log_gamma(DistributionModel* dm){
	double alpha = Parameters_value(dm->parameters, 0);
	double beta = Parameters_value(dm->parameters, 1);
	double logP = -gammln(alpha) * Parameters_count(dm->x);
	double beta_alpha = pow(beta, alpha);
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		double x = Parameters_value(dm->x, i);
		logP += log(beta_alpha*pow(x,alpha-1.0)) - beta*x;
//		logP += dloggamma(Parameters_value(dm->x, i), alpha, beta);
	}
	return logP;
}

double DistributionModel_dlog_gamma(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			return (Parameters_value(dm->parameters, 0)-1.0)/Parameter_value(p) - Parameters_value(dm->parameters, 1);
		}
	}
	return 0;
}

double DistributionModel_log_exp(DistributionModel* dm){
	double lambda = Parameters_value(dm->parameters, 0);
	double logP = log(lambda) * Parameters_count(dm->x);
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		//logP += log(dexp(Parameters_value(dm->x, i), Parameters_value(dm->parameters, 0)));
		logP -= lambda * Parameters_value(dm->x, i);
	}
	return logP;
}

double DistributionModel_dlog_exp(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			return -Parameters_value(dm->parameters, 0);
		}
	}
	return 0;
}

// Flat dirichlet
double DistributionModel_log_flat_dirichlet(DistributionModel* dm){
	return log(ddirchlet_flat(dm->simplex->K));
}

double DistributionModel_dlog_flat_dirichlet(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

// Dirichlet
double DistributionModel_log_dirichlet(DistributionModel* dm){
	if(dm->x != NULL){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			dm->tempx[i] = Parameters_value(dm->x, i);
		}
	}
	else if(dm->simplex != NULL){
		return log(ddirchlet(dm->simplex->get_values(dm->simplex), dm->simplex->K, dm->tempp));
	}
	else{
		assert(0);
	}
	return log(ddirchlet(dm->tempx, Parameters_count(dm->x), dm->tempp));
}

//TODO: implement
double DistributionModel_dlog_dirichlet(DistributionModel* dm, const Parameter* p){
	exit(1);
	return 0;
}

DistributionModel* new_IndependantGammaDistributionModel(const double shape, const double rate, const Parameters* x){
	Parameters* ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter("gamma.shape", shape, NULL));
	Parameters_move(ps, new_Parameter("gamma.rate", rate, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->logP = DistributionModel_log_gamma;
	dm->dlogP = DistributionModel_dlog_gamma;
	dm->clone = _clone_dist;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantGammaDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->logP = DistributionModel_log_gamma;
	dm->dlogP = DistributionModel_dlog_gamma;
	dm->clone = _clone_dist;
	return dm;
}

DistributionModel* new_IndependantExpDistributionModel(const double lambda, const Parameters* x){
	Parameters* ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter("exp.lambda", lambda, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->logP = DistributionModel_log_exp;
	dm->dlogP = DistributionModel_dlog_exp;
	dm->clone = _clone_dist;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantExpDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->logP = DistributionModel_log_exp;
	dm->dlogP = DistributionModel_dlog_exp;
	dm->clone = _clone_dist;
	return dm;
}

DistributionModel* new_FlatDirichletDistributionModel(Simplex* simplex){
	DistributionModel* dm = new_DistributionModelSimplex(NULL, simplex);
	dm->logP = DistributionModel_log_flat_dirichlet;
	dm->dlogP = DistributionModel_dlog_flat_dirichlet;
	dm->clone = _clone_dist;
	return dm;
}

DistributionModel* new_DirichletDistributionModel(const double* alpha, Simplex* simplex){
	Parameters* ps = new_Parameters(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		Parameters_move(ps, new_Parameter("dirichlet.", alpha[i], NULL));
	}
	DistributionModel* dm = new_DistributionModelSimplex(ps, simplex);// len(x)==len(alpha)
	dm->logP = DistributionModel_log_dirichlet;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->clone = _clone_dist;
	dm->tempx = dvector(simplex->K);
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		dm->tempp[i] = Parameters_value(dm->parameters, i);
	}
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_DirichletDistributionModel_with_parameters(const Parameters* parameters, Simplex* simplex){
	Parameters* ps = new_Parameters(simplex->K);
	Parameters_add_parameters(ps, parameters);
	DistributionModel* dm = new_DistributionModelSimplex(ps, simplex);// len(x)==len(alpha)
	dm->logP = DistributionModel_log_dirichlet;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->clone = _clone_dist;
	dm->tempx = dvector(simplex->K);
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		dm->tempp[i] = Parameters_value(dm->parameters, i);
	}
	return dm;
}

//MARK: tree prior

double logFactorial(int n){
    double logF = 0;
    for (int i = 2; i < n; i++) {
        logF += log(i);
    }
    return logF;
}


double DistributionModel_log_uniform_tree(DistributionModel* dm){
    Tree* tree = ((Model*)dm->data)->obj;
    int n = Tree_tip_count(tree);
    return -logFactorial(n*2-5) + logFactorial(n-3) + (n-3)*log(2);
}

double DistributionModel_dlog_uniform_tree(DistributionModel* dm, const Parameter* p){
    return 0.0;
}

DistributionModel* new_UniformTreeDistribution(Model* tree){
    DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
    assert(dm);
    dm->parameters = NULL;
    dm->x = NULL;
    dm->simplex = NULL;
    dm->logP = NULL;
    dm->dlogP = NULL;
    dm->free = _free_dist;
    dm->clone = _clone_dist;
    dm->data = tree;
    dm->tempx = NULL;
    dm->tempp = NULL;
    dm->logP = DistributionModel_log_uniform_tree;
    dm->dlogP = DistributionModel_dlog_uniform_tree;
    dm->clone = _clone_dist;
    return dm;
}

static double _dist_model_logP(Model *self){
	DistributionModel* cm = (DistributionModel*)self->obj;
	return cm->logP(cm);
}

static double _dist_model_dlogP(Model *self, const Parameter* p){
	DistributionModel* cm = (DistributionModel*)self->obj;
	for (int i = 0; i < Parameters_count(cm->x); i++) {
		if(Parameters_at(cm->x, i) == p){
			return cm->dlogP(cm, p);
		}
	}
	return 0;
}

static void _dist_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free distribution model %s\n", self->name);
		DistributionModel* cm = (DistributionModel*)self->obj;
		Model* msimplex = (Model*)self->data;
		if(msimplex != NULL){
			msimplex->free(msimplex);
		}
		if(cm->x != NULL) free_Parameters(cm->x);
		if(cm->parameters != NULL) free_Parameters(cm->parameters);
		if(cm->tempx != NULL) free(cm->tempx);
		if(cm->tempp != NULL) free(cm->tempp);
		free(cm);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _dist_model_clone( Model *self, Hashtable* hash ){
	if(Hashtable_exists(hash, self->name)){
		return Hashtable_get(hash, self->name);
	}
	DistributionModel* dm = (DistributionModel*)self->obj;
	
	Model* msimplex = (Model*)self->data;
	Model* msimplexclone = NULL;
	
	Parameters* x = NULL;
	if(dm->x != NULL){
		x = new_Parameters(Parameters_count(dm->x));
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			char* name = Parameters_name(dm->x, i);
			if(Hashtable_exists(hash, name)){
				Parameters_add(x, Hashtable_get(hash, name));
			}
			else{
				Parameter* p = clone_Parameter(Parameters_at(dm->x, i));
				Parameters_move(x, p);
				Hashtable_add(hash, name, p);
			}
		}
	}
	else if(dm->simplex != NULL){
		if (Hashtable_exists(hash, msimplex->name)) {
			msimplexclone = Hashtable_get(hash, msimplex->name);
			msimplexclone->ref_count++; // it is decremented at the end using free
		}
		else{
			msimplexclone = msimplex->clone(msimplex, hash);
			Hashtable_add(hash, msimplexclone->name, msimplexclone);
		}
		
	}
	Parameters* params = NULL;
	
	// Flat dirichlet does not have parameters
	if(dm->parameters != NULL){
		params = new_Parameters(Parameters_count(dm->parameters));
		for (int i = 0; i < Parameters_count(dm->parameters); i++) {
			char* name = Parameters_name(dm->parameters, i);
			if(Hashtable_exists(hash, name)){
				Parameters_add(params, Hashtable_get(hash, name));
			}
			else{
				Parameter* p = clone_Parameter(Parameters_at(dm->parameters, i));
				Parameters_move(params, p);
				Hashtable_add(hash, name, p);
			}
		}
	}
	Simplex* s = NULL;
	if(msimplexclone != NULL){
		s = (Simplex*)msimplexclone->obj;
	}
	DistributionModel* dmclone = clone_DistributionModel_with_parameters(dm, params, x, s);
	
	free_Parameters(params);
	free_Parameters(x);
	Model* clone = new_DistributionModel3(self->name, dmclone, msimplexclone);
	Hashtable_add(hash, clone->name, clone);
	if(msimplexclone != NULL){
		msimplexclone->free(msimplexclone);
	}
	
	return clone;
}

static void _dist_model_get_free_parameters(Model* model, Parameters* parameters){
	DistributionModel* dm = (DistributionModel*)model->obj;
	
	Model* msimplex = (Model*)model->data;
}

Model* new_DistributionModel2(const char* name, DistributionModel* dm){
	Model *model = new_Model("distribution",name, dm);
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	model->get_free_parameters = _dist_model_get_free_parameters;
	return model;
}

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* simplex){
	Model *model = new_Model("distribution",name, dm);
	model->data = simplex;
	if(simplex != NULL) simplex->ref_count++;
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	model->get_free_parameters = _dist_model_get_free_parameters;
	return model;
}

Model* new_TreeDistributionModel(const char* name, DistributionModel* dm, Model* tree){
    Model *model = new_Model("distribution",name, dm);
    model->data = tree;
    tree->ref_count++;
    model->logP = _dist_model_logP;
    model->dlogP = _dist_model_dlogP;
    model->free = _dist_model_free;
    model->clone = _dist_model_clone;
    model->get_free_parameters = _dist_model_get_free_parameters;
    return model;
}

char** get_parameters(json_node* node, Hashtable* hash, Parameters* parameters){
	json_node* parameters_node = get_json_node(node, "parameters");
	char** param_names = NULL;
	
	if(parameters_node->node_type == MJSON_ARRAY){
		for (int i = 0; i < parameters_node->child_count; i++) {
			json_node* child = parameters_node->children[i];
			// it's a ref
			if (child->node_type == MJSON_STRING) {
				char* ref = (char*)child->value;
				Parameter* p = Hashtable_get(hash, ref+1);
				Parameters_add(parameters, p);
			}
			// it's a value
			else if(child->node_type == MJSON_PRIMITIVE){
				double v = atof((char*)child->value);
				Parameters_move(parameters, new_Parameter("anonymous", v, NULL));
			}
			else{
				exit(1);
			}
		}
	}
	else if(parameters_node->node_type == MJSON_OBJECT){
		param_names = malloc(sizeof(char*)*parameters_node->child_count);
		for(int i = 0; i < parameters_node->child_count; i++){
			json_node* p_node = parameters_node->children[i];
			Parameter* p = new_Parameter_from_json(p_node, hash);
			Parameters_move(parameters, p);
			param_names[i] = p_node->key;
		}
	}
	// it's a ref
	else if(parameters_node->node_type == MJSON_STRING){
		char* ref = (char*)parameters_node->value;
		Parameter* p = Hashtable_get(hash, ref+1);
		Parameters_add(parameters, p);
	}
	else{
		exit(1);
	}
	return param_names;
}

void get_x(json_node* node, Hashtable* hash, Parameters* x){
	json_node* x_node = get_json_node(node, "x");
	
	if(x_node->node_type == MJSON_ARRAY){
		for (int i = 0; i < x_node->child_count; i++) {
			json_node* child = x_node->children[i];
			// it's a ref
			if (child->node_type == MJSON_STRING) {
				char* ref = (char*)child->value;
				Parameter* p = Hashtable_get(hash, ref+1);
				Parameters_add(x, p);
			}
			else{
				exit(1);
			}
		}
	}
	// it's a ref
	else if(x_node->node_type == MJSON_STRING){
		char* ref = (char*)x_node->value;
		Parameter* p = Hashtable_get(hash, ref+1);
		Parameters_add(x, p);
	}
	else{
		exit(1);
	}
}

Model* get_simplex(json_node* node, Hashtable* hash){
	json_node* x_node = get_json_node(node, "x");
	char* ref = (char*)x_node->value;
	return Hashtable_get(hash, ref+1);
}

Model* new_DistributionModel_from_json(json_node* node, Hashtable* hash){
	char* d_string = get_json_node_value_string(node, "distribution");
	char* id = get_json_node_value_string(node, "id");
	json_node* tree_node = get_json_node(node, "tree");
	
	Parameters* parameters = NULL;
	Parameters* x = NULL;
	DistributionModel* dm = NULL;
	Model* model = NULL;
	char** param_names = NULL;
	
	if (strcasecmp(d_string, "exponential") == 0) {
		parameters = new_Parameters(1);
		x = new_Parameters(1);
		param_names = get_parameters(node, hash, parameters);
		
		if (tree_node != NULL) {
			char* ref = (char*)tree_node->value;
			Model* mtree = Hashtable_get(hash, ref+1);
			Tree* tree = mtree->obj;
			for (int i = 0; i < Tree_node_count(tree); i++) {
				Node* n = Tree_node(tree, i);
				if (!Node_isroot(n) && !(Node_isroot(n) && Node_right(Node_parent(n)) == n)) {
					Parameters_add(x, n->distance);
				}
			}
		}
		else{
			get_x(node, hash, x);
		}
		dm = new_IndependantExpDistributionModel_with_parameters(parameters, x);
		model = new_DistributionModel2(id, dm);
	}
	else if (strcasecmp(d_string, "dirichlet") == 0) {
		parameters = new_Parameters(1);
		get_parameters(node, hash, parameters);
		Model* simplex = get_simplex(node, hash);
		int i = 0;
		for ( ; i < Parameters_count(parameters); i++) {
			if(Parameters_value(parameters, i) != 1) break;
		}
		if(i == Parameters_count(parameters)){
			dm = new_FlatDirichletDistributionModel((Simplex*)simplex->obj);
		}
		else{
			dm = new_DirichletDistributionModel_with_parameters(parameters, (Simplex*)simplex->obj);
		}
		Model* model = new_DistributionModel3(id, dm, simplex);
		
		free_Parameters(parameters);
		free_Parameters(x);
		return model;
    }
    else if (strcasecmp(d_string, "gamma") == 0) {
        parameters = new_Parameters(1);
        x = new_Parameters(1);
        param_names = get_parameters(node, hash, parameters);
        get_x(node, hash, x);
        dm = new_IndependantGammaDistributionModel_with_parameters(parameters, x);
        model = new_DistributionModel2(id, dm);
    }
    else if (strcasecmp(d_string, "topology") == 0) {
        char* ref = get_json_node_value_string(node, "tree");
        Model* mtree = Hashtable_get(hash, ref+1);
        dm = new_UniformTreeDistribution(mtree);
        model = new_TreeDistributionModel(id, dm, mtree);
    }
	else{
		printf("%s\n", d_string);
		exit(10);
	}
	
	free_Parameters(parameters);
	free_Parameters(x);
	if(param_names != NULL) free(param_names);
	
	return model;
}
