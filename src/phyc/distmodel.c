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
		Parameters_move(clone->parameters, clone_Parameter(Parameters_at(dm->parameters, i), true));
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
	clone->parameters = new_Parameters(Parameters_count(params));
	for (int i = 0; i < Parameters_count(params); i++) {
		Parameters_add(clone->parameters, Parameters_at(params, i));
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
		clone->simplex = clone_Simplex(simplex);
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

DistributionModel* new_DistributionModel(Parameters* p, const Parameters* x){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
	dm->parameters = p;
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
	double logP = 0;
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		logP += dloggamma(Parameters_value(dm->x, i), Parameters_value(dm->parameters, 0), Parameters_value(dm->parameters, 1));
	}
	return logP;
}

double DistributionModel_dlog_gamma(DistributionModel* dm, const Parameter* p){
	return (Parameters_value(dm->parameters, 0)-1.0)/Parameter_value(p) - Parameters_value(dm->parameters, 1);
}

double DistributionModel_log_exp(DistributionModel* dm){
	double logP = 0;
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		logP += log(dexp(Parameters_value(dm->x, i), Parameters_value(dm->parameters, 0)));
	}
	return logP;
}

double DistributionModel_dlog_exp(DistributionModel* dm, const Parameter* p){
	return -Parameters_value(dm->parameters, 0);
}

// Flat dirichlet
double DistributionModel_log_flat_dirichlet(DistributionModel* dm){
	return log(ddirchlet_flat(Parameters_count(dm->x)));
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
	return dm;
}

DistributionModel* new_IndependantExpDistributionModel(const double lambda, const Parameters* x){
	Parameters* ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter("exp.lambda", lambda, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
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
		printf("Free distribution model %s\n", self->name);
		DistributionModel* cm = (DistributionModel*)self->obj;
		cm->free(cm);
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
				Parameter* p = clone_Parameter(Parameters_at(dm->x, i), true);
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
				Parameter* p = clone_Parameter(Parameters_at(dm->parameters, i), true);
				Parameters_move(params, p);
				Hashtable_add(hash, name, p);
			}
		}
	}
	
	DistributionModel* dmclone = clone_DistributionModel_with_parameters(dm, params, x, (Simplex*)msimplexclone->obj);
	
	free_Parameters(params);
	free_Parameters(x);
	Model* clone = new_DistributionModel3(self->name, dmclone, msimplexclone);
	Hashtable_add(hash, clone->name, clone);
	msimplexclone->free(msimplexclone);
	
	return clone;
}

Model* new_DistributionModel2(const char* name, DistributionModel* dm){
	Model *model = new_Model(name, dm);
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	return model;
}

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* simplex){
	Model *model = new_Model(name, dm);
	model->data = simplex;
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	return model;
}