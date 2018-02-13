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

#include "filereader.h"
#include "statistics.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


static double _DistributionModel_dlog_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

static double _DistributionModel_d2log_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

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
	clone->d2logP = dm->d2logP;
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
	clone->d2logP = dm->d2logP;
	clone->clone = dm->clone;
	clone->free = dm->free;
	clone->tempp = NULL;
	clone->tempx = NULL;
	if(dm->tempp != NULL) clone->tempp = clone_dvector(dm->tempp, Parameters_count(dm->parameters));
	if(dm->tempx != NULL) clone->tempx = clone_dvector(dm->tempx, Parameters_count(dm->x));
	clone->lp = dm->lp;
	clone->need_update = dm->need_update;
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
	dm->d2logP = NULL;
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
	dm->d2logP = NULL;
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

double DistributionModel_d2log_gamma(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			return -(Parameters_value(dm->parameters, 0)-1.0)/Parameter_value(p)/Parameter_value(p);
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

double DistributionModel_d2log_exp(DistributionModel* dm, const Parameter* p){
	return 0;
}

// Flat dirichlet
double DistributionModel_log_flat_dirichlet(DistributionModel* dm){
	return log(ddirchlet_flat(dm->simplex->K));
}

double DistributionModel_dlog_flat_dirichlet(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_d2log_flat_dirichlet(DistributionModel* dm, const Parameter* p){
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

//TODO: implement
double DistributionModel_d2log_dirichlet(DistributionModel* dm, const Parameter* p){
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
	dm->d2logP = DistributionModel_d2log_gamma;
	dm->clone = _clone_dist;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantGammaDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->logP = DistributionModel_log_gamma;
	dm->dlogP = DistributionModel_dlog_gamma;
	dm->d2logP = DistributionModel_d2log_gamma;
	dm->clone = _clone_dist;
	return dm;
}

DistributionModel* new_IndependantExpDistributionModel(const double lambda, const Parameters* x){
	Parameters* ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter("exp.lambda", lambda, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->logP = DistributionModel_log_exp;
	dm->dlogP = DistributionModel_dlog_exp;
	dm->d2logP = DistributionModel_d2log_exp;
	dm->clone = _clone_dist;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantExpDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->logP = DistributionModel_log_exp;
	dm->dlogP = DistributionModel_dlog_exp;
	dm->d2logP = DistributionModel_d2log_exp;
	dm->clone = _clone_dist;
	return dm;
}

DistributionModel* new_FlatDirichletDistributionModel(Simplex* simplex){
	DistributionModel* dm = new_DistributionModelSimplex(NULL, simplex);
	dm->logP = DistributionModel_log_flat_dirichlet;
	dm->dlogP = DistributionModel_dlog_flat_dirichlet;
	dm->d2logP = DistributionModel_d2log_flat_dirichlet;
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
	dm->d2logP = DistributionModel_d2log_dirichlet;
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
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->clone = _clone_dist;
	dm->tempx = dvector(simplex->K);
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		dm->tempp[i] = Parameters_value(dm->parameters, i);
	}
	return dm;
}

//MARK: Multivariate Normal

typedef struct gsl_multivariate_normal_wrapper_t{
	gsl_vector* mu;
	gsl_matrix* L;
	gsl_vector* x; // used for sampling or pdf
	gsl_vector * work;
	gsl_rng* rng;
	bool transform;
}gsl_multivariate_normal_wrapper_t;

// if dst is NULL we assign directly the sampled values to dm->x
static void _sample_multivariate_normal(DistributionModel* dm, double* dst){
	gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
	
	gsl_ran_multivariate_gaussian(wrapper->rng, wrapper->mu, wrapper->L, wrapper->x);
	
	if(wrapper->transform){
		if(dst == NULL){
			for (int j = 0; j < Parameters_count(dm->x); j++) {
				Parameters_set_value(dm->x, j, exp(gsl_vector_get(wrapper->x, j)));
			}
		}
		else{
			for (int j = 0; j < Parameters_count(dm->x); j++) {
				dst[j] = exp(gsl_vector_get(wrapper->x, j));
			}
		}
	}
	else{
		if(dst == NULL){
			for (int j = 0; j < Parameters_count(dm->x); j++) {
				Parameters_set_value(dm->x, j, gsl_vector_get(wrapper->x, j));
			}
		}
		else{
			for (int j = 0; j < Parameters_count(dm->x); j++) {
				dst[j] = gsl_vector_get(wrapper->x, j);
			}
		}
	}
}

static double _multivariate_normal_logP(DistributionModel* dm){
	gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
	
	double logJocobian = 0;
	if(wrapper->transform){
		for (int j = 0; j < Parameters_count(dm->x); j++) {
			double y = log(Parameters_value(dm->x, j));
			logJocobian -= y;
			gsl_vector_set(wrapper->x, j,  y);
		}
	}
	else{
		for (int j = 0; j < Parameters_count(dm->x); j++) {
			gsl_vector_set(wrapper->x, j,  Parameters_value(dm->x, j));
		}
	}
	double logP;
	gsl_ran_multivariate_gaussian_log_pdf (wrapper->x,  wrapper->mu, wrapper->L, &logP, wrapper->work);

	return logP + logJocobian;
}

static void _free_dist_gsl_multivariate_normal(DistributionModel*dm){
	if(dm->x != NULL) free_Parameters(dm->x);
	if(dm->parameters != NULL) free_Parameters(dm->parameters);
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	if(dm->simplex != NULL) free_Simplex(dm->simplex);
	gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
	gsl_vector_free(wrapper->mu);
	gsl_vector_free(wrapper->x);
	gsl_vector_free(wrapper->work);
	gsl_matrix_free(wrapper->L);
	free(wrapper);
	free(dm);
}

DistributionModel* new_MultivariateNormalDistributionModel_with_parameters(const double* mu, const double* sigma, const Parameters* x, gsl_rng* rng){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
	dm->parameters = NULL;
	dm->x = new_Parameters(Parameters_count(x));
	Parameters_add_parameters(dm->x, x);
	dm->simplex = NULL;
	dm->free = _free_dist_gsl_multivariate_normal;
	dm->clone = _clone_dist;
	dm->tempx = NULL;
	dm->tempp = NULL;
	dm->logP = _multivariate_normal_logP;
	dm->dlogP = _DistributionModel_dlog_0;
	dm->d2logP = _DistributionModel_d2log_0;
	dm->sample = _sample_multivariate_normal;
	dm->clone = _clone_dist;
	dm->need_update = true;
	
	size_t dim = Parameters_count(x);
	gsl_multivariate_normal_wrapper_t* wrapper = (gsl_multivariate_normal_wrapper_t*)malloc(sizeof(gsl_multivariate_normal_wrapper_t));
	wrapper->mu = gsl_vector_calloc(dim);
	//wrapper->Sigma = gsl_matrix_calloc(dim, dim);
	wrapper->L = gsl_matrix_calloc(dim, dim);
	wrapper->x = gsl_vector_calloc(dim);
	wrapper->work = gsl_vector_calloc(dim);

	for (int i = 0; i < dim; i++) {
		gsl_vector_set(wrapper->mu, i, mu[i]);
		for (int j = 0; j < dim; j++) {
			gsl_matrix_set(wrapper->L, i, j, sigma[i*dim+j]);
		}
	}
	gsl_linalg_cholesky_decomp1(wrapper->L);
	wrapper->rng = rng;
	dm->data = wrapper;
	return dm;
}

//MARK: tree prior

double DistributionModel_log_uniform_tree(DistributionModel* dm){
	if(dm->need_update){
		Tree* tree = ((Model*)dm->data)->obj;
		int n = Tree_tip_count(tree);
		dm->lp = logDoubleFactorial(n*2-5);
		dm->need_update = false;
	}
	return dm->lp;
}

DistributionModel* new_UniformTreeDistribution(Model* tree){
    DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
    assert(dm);
    dm->parameters = NULL;
    dm->x = NULL;
    dm->simplex = NULL;
    dm->free = _free_dist;
    dm->clone = _clone_dist;
    dm->data = tree;
    dm->tempx = NULL;
    dm->tempp = NULL;
	dm->logP = DistributionModel_log_uniform_tree;
	dm->dlogP = _DistributionModel_dlog_0;
	dm->d2logP = _DistributionModel_d2log_0;
    dm->clone = _clone_dist;
	dm->need_update = true;
    return dm;
}

static void _dist_model_store(Model* self){
	self->storedLogP = self->lp;
}

static void _dist_model_restore(Model* self){
	self->lp = self->storedLogP;
}

static double _dist_model_logP(Model *self){
	DistributionModel* cm = (DistributionModel*)self->obj;
	self->lp = cm->logP(cm);
	return self->lp;
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

static double _dist_model_d2logP(Model *self, const Parameter* p){
	DistributionModel* cm = (DistributionModel*)self->obj;
	for (int i = 0; i < Parameters_count(cm->x); i++) {
		if(Parameters_at(cm->x, i) == p){
			return cm->d2logP(cm, p);
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
	clone->store = self->store;
	clone->restore = self->restore;
	clone->storedLogP = self->storedLogP;
	clone->lp = self->lp;
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
	model->d2logP = _dist_model_d2logP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	model->get_free_parameters = _dist_model_get_free_parameters;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
	return model;
}

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* simplex){
	Model *model = new_Model("distribution",name, dm);
	model->data = simplex;
	if(simplex != NULL) simplex->ref_count++;
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->d2logP = _dist_model_d2logP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	model->get_free_parameters = _dist_model_get_free_parameters;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
	return model;
}

Model* new_TreeDistributionModel(const char* name, DistributionModel* dm, Model* tree){
    Model *model = new_Model("distribution",name, dm);
    model->data = tree;
    tree->ref_count++;
    model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->d2logP = _dist_model_d2logP;
    model->free = _dist_model_free;
    model->clone = _dist_model_clone;
    model->get_free_parameters = _dist_model_get_free_parameters;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
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
				if (!Node_isroot(n) && !(Node_isroot(Node_parent(n)) && Node_right(Node_parent(n)) == n)) {
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
	else if(strcasecmp(d_string, "empirical") == 0){
		char* file = get_json_node_value_string(node, "file");
		size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
		
		x = new_Parameters(1);
		get_parameters_references2(node, hash, x, "x");
		size_t paramCount = Parameters_count(x);
		
		int count = 0;
		char *ptr = NULL;
		double *temp = NULL;
		int l;
		Vector** vec = malloc(sizeof(Vector*)*paramCount);
		for(int i = 0; i < paramCount; i++){
			vec[i] = new_Vector(1000);
		}
		
		FileReader *reader = new_FileReader(file, 1000);
		reader->read_line(reader);
		ptr = reader->line;
		l = 0;
		char** header = String_split_char(ptr, '\t', &l);
		bool* in = bvector(l);
		
		for (int j = 0; j < paramCount; j++) {
			for (int i = 0; i < l; i++) {
				if(strcmp(Parameters_name(x, j), header[i]) == 0){
					in[i] = true;
					break;
				}
			}
		}
		for (int i = 0; i < l; i++) {
			free(header[i]);
		}
		free(header);
		
		while ( reader->read_line(reader) ) {
			StringBuffer_trim(reader->buffer);
			
			if ( reader->buffer->length == 0){
				continue;
			}
			if ( count >= burnin){
				ptr = reader->line;
				l = 0;
				temp = String_split_char_double( ptr, '\t', &l );
				int index = 0;
				for (int i = 0; i < l; i++) {
					if(in[i]){
						Vector_push(vec[index], log(temp[i]));
						index++;
					}
				}
				free(temp);
			}
			count++;
		}
		free_FileReader(reader);

		int n = Vector_length(vec[0]);
		double* mu = dvector(paramCount);
		double* sigma = dvector(paramCount*paramCount);
		
		for (int i = 0; i < paramCount; i++) {
			mu[i] = mean(Vector_data(vec[i]), n);
		}
		
		// Calculate covariance matrix
		for (int i = 0; i < paramCount; i++) {
			double* pp = Vector_data(vec[i]);
			sigma[i*paramCount+i] = variance(pp, n, mu[i]);
			for (int j = i+1; j < paramCount; j++) {
				double* pp2 = Vector_data(vec[j]);
				sigma[i*paramCount+j] = sigma[j*paramCount+i] = covariance(pp, pp2, mu[i], mu[j], n);
			}
		}
		
		dm = new_MultivariateNormalDistributionModel_with_parameters(mu, sigma, x, Hashtable_get(hash, "RANDOM_GENERATOR!@"));
		model = new_DistributionModel2(id, dm);
		
		gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
		wrapper->transform = true;
		
		free(mu);
		free(sigma);
		for (int i = 0; i < paramCount; i++) {
			free_Vector(vec[i]);
		}
		free(vec);
		free(in);
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
