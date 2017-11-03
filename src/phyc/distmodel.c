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
	if(dm->x != NULL) free_Parameters_soft(dm->x);
	if(dm->parameters != NULL) free_Parameters(dm->parameters);
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	// freeing data is left to the user
	free(dm);
}

static DistributionModel* _clone_dist(DistributionModel* dm){
	
	DistributionModel* clone = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(clone);
	clone->parameters = new_Parameters(Parameters_count(dm->parameters));
	for (int i = 0; i < Parameters_count(dm->parameters); i++) {
		Parameters_set_value(clone->parameters, i, Parameters_value(dm->parameters, i));
	}
	clone->x = new_Parameters(Parameters_count(dm->x));
	clone->logP = dm->logP;
	clone->dlogP = dm->dlogP;
	clone->clone = dm->clone;
	clone->free = dm->free;
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

double DistributionModel_log_flat_dirichlet(DistributionModel* dm){
	return log(ddirchlet_flat(Parameters_count(dm->x)));
}

double DistributionModel_dlog_flat_dirichlet(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_log_dirichlet(DistributionModel* dm){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		dm->tempx[i] = Parameters_value(dm->x, i);
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

DistributionModel* new_FlatDirichletDistributionModel(const Parameters* x){
	DistributionModel* dm = new_DistributionModel(NULL, x);
	dm->logP = DistributionModel_log_flat_dirichlet;
	dm->dlogP = DistributionModel_dlog_flat_dirichlet;
	dm->clone = _clone_dist;
	return dm;
}

DistributionModel* new_DirichletDistributionModel(const double* alpha, const Parameters* x){
	Parameters* ps = new_Parameters(Parameters_count(x));
	for (int i = 0; i < Parameters_count(x); i++) {
		Parameters_move(ps, new_Parameter("dirichlet.", alpha[i], NULL));
	}
	DistributionModel* dm = new_DistributionModel(ps, x);// len(x)==len(alpha)
	dm->logP = DistributionModel_log_dirichlet;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->clone = _clone_dist;
	dm->tempx = dvector(Parameters_count(x));
	dm->tempp = dvector(Parameters_count(x));
	for (int i = 0; i < Parameters_count(x); i++) {
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

static Model* _dist_model_clone( Model *self ){
	DistributionModel* cm = (DistributionModel*)self->obj;
	DistributionModel* cmclone = cm->clone(cm);
	return new_DistributionModel2(self->name, cmclone);
}

Model* new_DistributionModel2(const char* name, DistributionModel* cm){
	Model *model = new_Model(name, cm);
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	return model;
}
