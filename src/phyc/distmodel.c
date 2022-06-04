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
#include <strings.h>
#include <math.h>

#include "matrix.h"
#include "tree.h"

#include "exponential.h"
#include "gamma.h"

#include "filereader.h"
#include "statistics.h"
#include "parametersio.h"
#include "utilsio.h"

#include "distbeta.h"
#include "distbetaprime.h"
#include "distexp.h"
#include "distgamma.h"
#include "distcauchy.h"
#include "gmrf.h"
#include "distkumaraswamy.h"
#include "distdirichlet.h"
#include "distnormal.h"
#include "distmultinormal.h"
#include "distlognormal.h"
#include "distoneonx.h"
#include "ctmcscale.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>


double DistributionModel_dlog_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_d2log_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_ddlog_0(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	return 0.0;
}

static void _DistributionModel_error_sample(DistributionModel* dm, double* samples){
	fprintf(stderr, "Sample method not implemented\n");
	exit(1);
}

static double _DistributionModel_error_sample_evaluate(DistributionModel* dm){
	fprintf(stderr, "Sample and evaluate method not implemented\n");
	exit(1);
}

// do not free tree and simplex since they are managed by their model
static void _free_partial_distribution(DistributionModel*dm){
	if(dm->x != NULL) free_Parameters(dm->x);
    if(dm->parameters != NULL){
        for(size_t i = 0; i < dm->parameter_count; i++){
            free_Parameters(dm->parameters[i]);
        }
        free(dm->parameters);
    }
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	if(dm->gradient != NULL) free(dm->gradient);
	// freeing data is left to the user
	free(dm);
}

//static void _free_full_distribution(DistributionModel*dm){
//    _free_partial_distribution(dm);
//	if(dm->simplex != NULL) free_Simplex(dm->simplex);
//	if(dm->tree != NULL) free_Tree(dm->tree);
//	// freeing data is left to the user
//	free(dm);
//}

static DistributionModel* DistributionMode_clone(DistributionModel* dm){
    Parameters** parameters = NULL;
    Parameters* x = NULL;
    Simplex* simplex = NULL;
    
    if(dm->parameters != NULL){
        parameters = malloc(sizeof(Parameters*)*dm->parameter_count);
        for(size_t i = 0; i < dm->parameter_count; i++){
            parameters[i] = new_Parameters(Parameters_count(dm->parameters[i]));
            for (int j = 0; j < Parameters_count(dm->parameters[i]); j++) {
                Parameters_move(parameters[i], clone_Parameter(Parameters_at(dm->parameters[i], j)));
            }
        }
    }
	if(dm->x != NULL){
		x = new_Parameters(Parameters_count(dm->x));
	}
	else if(dm->simplex != NULL){
		simplex = clone_Simplex(dm->simplex);
	}
    
	return clone_DistributionModel_with_parameters(dm, parameters, x, simplex);
}

DistributionModel* clone_DistributionModel_with_parameters(DistributionModel* dm, Parameters** params, const Parameters* x, Simplex* simplex){
	DistributionModel* clone = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(clone);
	clone->type = dm->type;
	clone->parameters = NULL;
    clone->parameter_count = dm->parameter_count;
    
	if(dm->parameter_count > 0){
        clone->parameters = malloc(sizeof(Parameters*)*dm->parameter_count);
        for(size_t i = 0; i < dm->parameter_count; i++){
            clone->parameters[i] = new_Parameters(Parameters_count(dm->parameters[i]));
            for (int j = 0; j < Parameters_count(dm->parameters[i]); j++) {
                Parameters_add(clone->parameters[i], Parameters_at(params[i], j));
            }
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
	clone->ddlogP = dm->ddlogP;
	clone->sample = dm->sample;
	clone->sample_evaluate = dm->sample_evaluate;
	clone->logP_with_values = dm->logP_with_values;
	clone->clone = dm->clone;
	clone->free = dm->free;
	clone->tempp = NULL;
	clone->tempx = NULL;
    if(dm->tempp != NULL){
        size_t size = 0;
        for(size_t i = 0; i < dm->parameter_count; i++){
            size += Parameters_count(dm->parameters[i]);
        }
        clone->tempp = clone_dvector(dm->tempp, size);
    }
	if(dm->tempx != NULL) clone->tempx = clone_dvector(dm->tempx, Parameters_count(dm->x));
	clone->lp = dm->lp;
    clone->stored_lp = dm->stored_lp;
	clone->need_update = dm->need_update;
	clone->shift = dm->shift;
	clone->type = dm->type;
	clone->parameterization = dm->parameterization;
	clone->rng = dm->rng;
    clone->data = NULL;
	clone->gradient = NULL;
	clone->gradient_length = dm->gradient_length;
	clone->need_update = dm->need_update;
	clone->need_update_gradient = dm->need_update_gradient;
	clone->prepared_gradient = dm->prepared_gradient;
	if(dm->gradient != NULL){
		clone->gradient = clone_dvector(dm->gradient, dm->gradient_length);
	}
	return clone;
}

DistributionModel* new_DistributionModel(Parameters** p, size_t dim, Parameters* x){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
    dm->x = NULL;
	dm->parameters = NULL;
    dm->parameter_count = 0;
	
    if(p != NULL){
        dm->parameters = malloc(sizeof(Parameters*)*dim);
        dm->parameter_count = dim;
        for(size_t i = 0; i < dim; i++){
            dm->parameters[i] = new_Parameters(Parameters_count(p[i]));
            Parameters_add_parameters(dm->parameters[i], p[i]);
            Parameters_set_name2(dm->parameters[i], Parameters_name2(p[i]));
        }
    }
    dm->x = x;
//    if(x != NULL){
//        dm->x = new_Parameters(Parameters_count(x));
//        Parameters_add_parameters(dm->x, x);
//        Parameters_set_name2(dm->x, Parameters_name2(x));
//    }
	dm->simplex = NULL;
	dm->tree = NULL;
	dm->logP = NULL;
	dm->dlogP = NULL;
	dm->d2logP = NULL;
	dm->ddlogP = NULL;
	dm->sample = _DistributionModel_error_sample;
	dm->sample_evaluate = _DistributionModel_error_sample_evaluate;
	dm->free = _free_partial_distribution;
	dm->clone = DistributionMode_clone;
	dm->data = NULL;
	dm->tempx = NULL;
	dm->tempp = NULL;
	dm->need_update = true;
	
	dm->prepared_gradient = 0;;
	dm->gradient = NULL;
	dm->gradient_length = 0;
	dm->need_update_gradient = true;
	return dm;
}


//MARK: copula

typedef struct gsl_gaussian_copula_wrapper_t{
	gsl_vector* mu;
	gsl_matrix* L;
	gsl_vector* x; // used for sampling or pdf
	gsl_vector * work;
	gsl_matrix* cor;
	gsl_permutation* p;
	double det;
}gsl_gaussian_copula_wrapper_t;

static void _gaussian_copula_gamma_sample(DistributionModel* dm, double* dst){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	gsl_ran_multivariate_gaussian(dm->rng, wrapper->mu, wrapper->L, wrapper->x);
	
	for (int i = 0; i < paramCount; i++) {
		double u = gsl_cdf_gaussian_P(gsl_vector_get(wrapper->x, i), 1.0); // uniform
		double z = gsl_cdf_gamma_Pinv(u, Parameters_value(dm->parameters[0], i), 1.0/Parameters_value(dm->parameters[1], i)); // gamma distributed
		if(dst == NULL){
			Parameters_set_value(dm->x, i, z);
		}
		else{
			dst[i] = z;
		}
	}
}

static void _gaussian_copula_lognormal_sample(DistributionModel* dm, double* dst){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	gsl_ran_multivariate_gaussian(dm->rng, wrapper->mu, wrapper->L, wrapper->x);
	
	for (int i = 0; i < paramCount; i++) {
		double u = gsl_cdf_gaussian_P(gsl_vector_get(wrapper->x, i), 1.0); // uniform
		double z = gsl_cdf_lognormal_Pinv(u, Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i)); // lognormal distributed
		if(dst == NULL){
			Parameters_set_value(dm->x, i, z);
		}
		else{
			dst[i] = z;
		}
	}
}

double _gaussian_copula_gamma_logP_with_values(DistributionModel* dm, const double* values){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	double logP = -log(sqrt(wrapper->det));
	for (int i = 0; i < paramCount; i++) {
		logP += log(gsl_ran_gamma_pdf(values[i], Parameters_value(dm->parameters[0], i), 1.0/Parameters_value(dm->parameters[1], i)));
	}
	return logP;
}

static double _gaussian_copula_lognormal_logP_with_values(DistributionModel* dm, const double* values){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	double logP = -log(sqrt(wrapper->det));
	for (int i = 0; i < paramCount; i++) {
		logP += log(gsl_ran_lognormal_pdf(values[i], Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i)));
	}
	return logP;
}

static double _gaussian_copula_gamma_logP(DistributionModel* dm){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	double logP = -log(sqrt(wrapper->det));
	for (int i = 0; i < paramCount; i++) {
		logP += log(gsl_ran_gamma_pdf(Parameters_value(dm->x, i), Parameters_value(dm->parameters[0], i), 1.0/Parameters_value(dm->parameters[1], i)));
	}
	return logP;
}

static double _gaussian_copula_lognormal_logP(DistributionModel* dm){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	double logP = -log(sqrt(wrapper->det));
	for (int i = 0; i < paramCount; i++) {
		logP += log(gsl_ran_lognormal_pdf(Parameters_value(dm->x, i), Parameters_value(dm->parameters[0], i), Parameters_value(dm->parameters[1], i)));
	}
	return logP;
}

static DistributionModel* _clone_gaussian_copula_gamma(DistributionModel* dm){
	DistributionModel* clone = DistributionMode_clone(dm);
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)malloc(sizeof(gsl_gaussian_copula_wrapper_t));
	gsl_gaussian_copula_wrapper_t* wrapper_orig = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	wrapper->mu = gsl_vector_calloc(paramCount);
	wrapper->L = gsl_matrix_calloc(paramCount, paramCount);
	wrapper->cor = gsl_matrix_calloc(paramCount, paramCount);
	wrapper->x = gsl_vector_calloc(paramCount);
	wrapper->work = gsl_vector_calloc(paramCount);
	wrapper->p = gsl_permutation_alloc(paramCount);
	
	gsl_matrix_memcpy(wrapper->cor, wrapper_orig->cor);
	gsl_matrix_memcpy(wrapper->L, wrapper_orig->L);
	gsl_vector_memcpy(wrapper->mu, wrapper_orig->mu);
	gsl_vector_memcpy(wrapper->x, wrapper_orig->x);
	gsl_vector_memcpy(wrapper->work, wrapper_orig->work);
	gsl_permutation_memcpy(wrapper->p, wrapper_orig->p);
	wrapper->det = wrapper_orig->det;
	clone->data = wrapper;
	return clone;
}

static void _free_gaussian_copula_gamma(DistributionModel*dm){
	if(dm->x != NULL) free_Parameters(dm->x);
    if(dm->parameters != NULL){
        for(size_t i = 0; i < dm->parameter_count; i++){
            free_Parameters(dm->parameters[i]);
        }
        free(dm->parameters);
    }
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	if(dm->simplex != NULL) free_Simplex(dm->simplex);
	gsl_gaussian_copula_wrapper_t* wrapper = dm->data;
	gsl_vector_free(wrapper->mu);
	gsl_vector_free(wrapper->x);
	gsl_vector_free(wrapper->work);
	gsl_matrix_free(wrapper->L);
	gsl_matrix_free(wrapper->cor);
	gsl_permutation_free(wrapper->p);
	free(wrapper);
	free(dm);
}

static double _gaussian_copula_gamma_determinant(DistributionModel* dm){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	
	for (size_t i = 0; i < paramCount; i++) {
		for (size_t j = 0; j < paramCount; j++) {
			gsl_matrix_set(wrapper->cor, i, j, Parameters_value(dm->parameters[1], i*paramCount+j));
		}
	}
	
	int signum;
	gsl_linalg_LU_decomp(wrapper->cor, wrapper->p , &signum);
	return gsl_linalg_LU_det(wrapper->cor, signum);
}

DistributionModel* new_CopulaDistributionModel_with_parameters(Parameters** p, size_t dim, const Parameters* x){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
	dm->parameters = NULL;
    dm->parameter_count = 0;

	if(p != NULL){
        dm->parameters = malloc(sizeof(Parameters*)*dim);
        dm->parameter_count = dim;
        for(size_t i = 0; i < dim; i++){
            dm->parameters[i] = new_Parameters(Parameters_count(p[i]));
            Parameters_add_parameters(dm->parameters[i], p[i]);
        }
	}
	dm->x = new_Parameters(Parameters_count(x));
	Parameters_add_parameters(dm->x, x);
	dm->simplex = NULL;
	dm->tree = NULL;
	dm->free = _free_gaussian_copula_gamma;
	dm->clone = _clone_gaussian_copula_gamma;
	dm->tempx = NULL;
	dm->tempp = NULL;
	dm->logP = _gaussian_copula_gamma_logP;
	dm->logP_with_values = _gaussian_copula_gamma_logP_with_values;
	//TODO: copula methods
	dm->dlogP = NULL;
	dm->d2logP = NULL;
	dm->ddlogP = NULL;
	dm->sample = NULL;
	dm->need_update = true;
	
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)malloc(sizeof(gsl_gaussian_copula_wrapper_t));
	size_t paramCount = Parameters_count(x);
	wrapper->mu = gsl_vector_calloc(paramCount);
	wrapper->L = gsl_matrix_calloc(paramCount, paramCount);
	wrapper->cor = gsl_matrix_calloc(paramCount, paramCount);
	wrapper->x = gsl_vector_calloc(paramCount);
	wrapper->work = gsl_vector_calloc(paramCount);
	wrapper->p = gsl_permutation_alloc(paramCount);
	
	for (int i = 0; i < paramCount; i++) {
		gsl_vector_set(wrapper->mu, i, 0);
		for (int j = 0; j < paramCount; j++) {
			gsl_matrix_set(wrapper->cor, i, j, Parameters_value(dm->parameters[1], i*paramCount+j));
		}
	}
	gsl_matrix_memcpy(wrapper->L, wrapper->cor);
	gsl_linalg_cholesky_decomp1(wrapper->L);
	
	dm->data = wrapper;
	wrapper->det = _gaussian_copula_gamma_determinant(dm);
	return dm;
}

//MARK: tree prior

double DistributionModel_log_uniform_tree(DistributionModel* dm){
	if(dm->need_update){
		int n = Tree_tip_count(dm->tree);
		dm->lp = -logDoubleFactorial(n*2-5);
		dm->need_update = false;
	}
	return dm->lp;
}

DistributionModel* new_UniformTreeDistribution(Tree* tree){
    DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
    assert(dm);
	dm->type = DISTRIBUTION_UNIFORM;
    dm->parameters = NULL;
    dm->x = NULL;
    dm->simplex = NULL;
    dm->free = _free_partial_distribution;
    dm->clone = DistributionMode_clone;
    dm->tree = tree;
    dm->tempx = NULL;
    dm->tempp = NULL;
	dm->logP = DistributionModel_log_uniform_tree;
	dm->dlogP = DistributionModel_dlog_0;
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->need_update = true;
    return dm;
}

//MARK: Model functions

void _dist_model_handle_change( Model *self, Model *model, int index ){
	DistributionModel* dm = self->obj;
	dm->need_update = true;
	dm->need_update_gradient = true;
	self->listeners->fire( self->listeners, self, index );
}

void _dist_model_handle_restore( Model *self, Model *model, int index ){
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _dist_model_store(Model* self){
	self->storedLogP = self->lp;
	DistributionModel* dm = self->obj;
	dm->stored_lp = dm->lp;
    for(size_t i = 0; i < dm->parameter_count; i++){
        Parameters_store(dm->parameters[i]);
    }
	if(dm->simplex == NULL){
		Parameters_store(dm->x);
	}
	else{
		Model* msimplex = self->data;
		msimplex->store(msimplex);
	}
}

static void _dist_model_restore(Model* self){
	self->lp = self->storedLogP;
	DistributionModel* dm = self->obj;
	dm->lp = dm->stored_lp;
	bool changed = false;
	Parameter*p = NULL;
	// restore the parameters of the model
    for(size_t i = 0; i < dm->parameter_count; i++){
        for (size_t j = 0; j < Parameters_count(dm->parameters[i]); j++) {
            p = Parameters_at(dm->parameters[i], j);
            if (Parameter_changed(p)) {
                changed = true;
                Parameter_restore_quietly(p);
            }
        }
    }
	if (changed) {
		p->listeners->fire_restore(p->listeners, NULL, p->id);
	}
	
	// restore the domain
	changed = false;
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		p = Parameters_at(dm->x, i);
		if (Parameter_changed(p)) {
			changed = true;
			Parameter_restore_quietly(p);
		}
	}
	if (changed) {
		p->listeners->fire_restore(p->listeners, NULL, p->id);
	}
	
	if(dm->simplex != NULL){
		Model* msimplex = self->data;
		msimplex->restore(msimplex);
	}
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
	for (int i = 0; i < cm->parameter_count; i++) {
        for (int j = 0; j < Parameters_count(cm->parameters[i]); j++) {
            if(Parameters_at(cm->parameters[i], j) == p){
                return cm->dlogP(cm, p);
            }
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
	for (int i = 0; i < cm->parameter_count; i++) {
        for (int j = 0; j < Parameters_count(cm->parameters[i]); j++) {
            if(Parameters_at(cm->parameters[i], j) == p){
                return cm->d2logP(cm, p);
            }
        }
    }
	return 0;
}

static double _dist_model_ddlogP(Model *self, const Parameter* p1, const Parameter* p2){
	DistributionModel* cm = (DistributionModel*)self->obj;
	bool found1 = false;
	bool found2 = false;
	for (int i = 0; i < Parameters_count(cm->x); i++) {
		if(Parameters_at(cm->x, i) == p1){
			found1 = true;
		}
		else if(Parameters_at(cm->x, i) == p2){
			found2 = true;
		}
	}
	if(found1 && found2) return cm->ddlogP(cm, p1, p2);
	
	found1 = false;
	found2 = false;
    for (int i = 0; i < cm->parameter_count; i++) {
        for (int j = 0; j < Parameters_count(cm->parameters[i]); j++) {
            if(Parameters_at(cm->parameters[i], j) == p1){
                found1 = true;
            }
            else if(Parameters_at(cm->parameters[i], j) == p2){
                found2 = true;
            }
        }
    }
	if(found1 && found2) return cm->ddlogP(cm, p1, p2);
	return 0;
}

static void _dist_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free distribution model %s\n", self->name);
		DistributionModel* cm = (DistributionModel*)self->obj;
		// tree or simplex
		if(self->data != NULL){
			Model* model = self->data;
			model->free(model);
		}
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
	Parameters** params = NULL;
	
	// Flat dirichlet does not have parameters
	if(dm->parameters != NULL){
        params = malloc(sizeof(Parameters*)*dm->parameter_count);
		for (int i = 0; i < dm->parameter_count; i++) {
            params[i] = new_Parameters(Parameters_count(dm->parameters[i]));
            for (int j = 0; j < Parameters_count(dm->parameters[i]); j++) {
                char* name = Parameters_name(dm->parameters[i], j);
                if(Hashtable_exists(hash, name)){
                    Parameters_add(params[i], Hashtable_get(hash, name));
                }
                else{
                    Parameter* p = clone_Parameter(Parameters_at(dm->parameters[i], j));
                    Parameters_move(params[i], p);
                    Hashtable_add(hash, name, p);
                }
            }
		}
	}
	Simplex* s = NULL;
	if(msimplexclone != NULL){
		s = (Simplex*)msimplexclone->obj;
	}
	DistributionModel* dmclone = clone_DistributionModel_with_parameters(dm, params, x, s);
	
    for (int i = 0; i < dm->parameter_count; i++) {
        free_Parameters(params[i]);
    }
    free(params);
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
	clone->sample = self->sample;
	clone->samplable = self->samplable;
    clone->need_update = self->need_update;
	return clone;
}

static void _dist_model_sample(Model* model, double* samples, double* logP){
	DistributionModel* dm = (DistributionModel*)model->obj;
	dm->sample(dm, samples);
	if(logP != NULL){
		*logP = dm->logP_with_values(dm, samples);
	}
}

static double _dist_model_sample_evaluate(Model* model){
	DistributionModel* dm = (DistributionModel*)model->obj;
	model->lp = dm->sample_evaluate(dm);
	return model->lp;
}

Model* new_DistributionModel2(const char* name, DistributionModel* dm){
	Model *model = new_Model(MODEL_DISTRIBUTION,name, dm);
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->d2logP = _dist_model_d2logP;
	model->ddlogP = _dist_model_ddlogP;
	model->free = _dist_model_free;
	model->clone = _dist_model_clone;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
	model->update = _dist_model_handle_change;
	model->handle_restore = _dist_model_handle_restore;
	model->sample = _dist_model_sample;
	model->sample_evaluate = _dist_model_sample_evaluate;
	model->samplable = false;
	
	for ( size_t i = 0; i < dm->parameter_count; i++ ) {
        for ( size_t j = 0; j < Parameters_count(dm->parameters[i]); j++ ) {
            Parameters_at(dm->parameters[i], j)->listeners->add( Parameters_at(dm->parameters[i], j)->listeners, model );
        }
	}
	for ( int i = 0; i < Parameters_count(dm->x); i++ ) {
		Parameters_at(dm->x, i)->listeners->add( Parameters_at(dm->x, i)->listeners, model );
	}
	return model;
}

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* simplex){
	Model *model = new_DistributionModel2(name, dm);
	model->data = simplex;
    if(simplex != NULL){
        simplex->ref_count++;
        simplex->listeners->add(simplex->listeners, model);
    }
	return model;
}

double _gamma_wrapper(void* data, double x){
	double* params = (double*)data;
	return gsl_cdf_gamma_P(x, params[0], 1.0/params[1]);
}

double _lognormal_wrapper(void* data, double x){
	double* params = (double*)data;
	return gsl_cdf_lognormal_P(x, params[0], params[1]);
}

double _exponential_wrapper(void* data, double x){
	double* params = (double*)data;
	return gsl_cdf_exponential_P(x, params[0]);
}

double ks(double* x, size_t length, double(*cdf)(void*, double), void* data){
	qsort(x, length, sizeof(double), qsort_asc_dvector);
	double fn = 0;
	double fo = 0;
	double d = 0;
	for (size_t i = 0; i < length; i++) {
		fn = (double)i/length;
		double ff = cdf(data, x[i]);
		double dt = dmax(fabs(fo-ff), fabs(fn-ff));
		if(dt > d){
			d = dt;
		}
		fo = fn;
	}
	return d;
}


Parameters** distmodel_get_parameters(const char* who, json_node* parameters_node, Hashtable* hash, size_t *dim){
    Parameters** parameters = malloc(sizeof(Parameters*)*parameters_node->child_count);
    *dim = parameters_node->child_count;
    for (int i = 0; i < parameters_node->child_count; i++) {
        parameters[i] = new_Parameters(1);
        json_node* p_node = parameters_node->children[i];
        // it is a reference
        if(p_node->node_type == MJSON_STRING){
            char* ref = p_node->value;
            if(ref[0] == '&'){
                Parameter* p = Hashtable_get(hash, ref+1);
                Parameters_add(parameters[i], p);
            }
            else if (ref[0] == '%'){
                Parameters* ps = Hashtable_get(hash, ref+1);
                Parameters_add_parameters(parameters[i], ps);
            }
            else{
                fprintf(stderr, "Distribution with ID %s cannot access ref: %s", who, ref);
                exit(1);
            }
            Parameters_set_name2(parameters[i], ref+1);
        }
        else if(p_node->node_type == MJSON_OBJECT){
            // it is a multi dimensional parameter
            if(get_json_node_value_int(p_node, "dimension", 0) != 0){
                Parameters* multi_parameter = new_MultiParameter_from_json(p_node, hash);
                Parameters_add_parameters(parameters[i], multi_parameter);
                Parameters_set_name2(parameters[i], Parameters_name2(multi_parameter));

                for (int i = 0; i < Parameters_count(multi_parameter); i++) {
                    Hashtable_add(hash, Parameters_name(multi_parameter, i), Parameters_at(multi_parameter, i));
                }
                free_Parameters(multi_parameter);
                
                Hashtable_add(hash, Parameters_name2(parameters[i]), parameters[i]);
            }
            // single parameter
            else{
                Parameter* parameter = new_Parameter_from_json(p_node, hash);
                Parameters_move(parameters[i], parameter);
                Parameters_set_name2(parameters[i], Parameter_name(parameter));

                Hashtable_add(hash, Parameter_name(parameter), parameter);
            }
        }
        else{
            fprintf(stderr, "%s cannot read parameters\n", who);
            fprintf(stderr, "with id: %s", get_json_node_value_string(p_node, "id"));
            fflush(stderr);
            exit(2);
        }
    }
    return parameters;
}

void distmodel_get_ref(const char* ref, Hashtable* hash, Parameters* parameters){
	if (ref[0] == '&') {
		Parameter* p = Hashtable_get(hash, ref+1);
		Parameters_add(parameters, p);
	}
	else if (ref[0] == '%') {
		// slicing
		if (ref[strlen(ref)-1] == ']') {
			get_parameters_slice(ref+1, parameters, hash);
		}
		else{
			Parameters* ps = Hashtable_get(hash, ref+1);
			Parameters_add_parameters(parameters, ps);
		}
	}
	// simplex
	else if (ref[0] == '$') {
		Model* msimplex = Hashtable_get(hash, ref+1);
		Simplex* simplex = msimplex->obj;
		Parameters_add_parameters(parameters, simplex->parameters);
	}
}

Parameters* distmodel_get_x(const char* who, json_node* node, Hashtable* hash){
//void get_parameters_from_node(json_node* node, Hashtable* hash, Parameters* parameters){
    Parameters* parameters = new_Parameters(1);
    // it's a ref
    if(node->node_type == MJSON_STRING){
        char* ref = (char*)node->value;
		distmodel_get_ref(ref, hash, parameters);
        Parameters_set_name2(parameters, ref+1);
    }
	else if(node->node_type == MJSON_ARRAY){
		for (int i = 0; i < node->child_count; i++) {
			json_node* child = node->children[i];
			char* ref = (char*)child->value;
			distmodel_get_ref(ref, hash, parameters);
		}
	}
    else if(node->node_type == MJSON_OBJECT){
        json_node* p_node_dimension = get_json_node(node, "dimension");
        if(p_node_dimension != NULL ){
            Parameters* multi_parameter = new_MultiParameter_from_json(node, hash);
            Parameters_add_parameters(parameters, multi_parameter);
            Parameters_set_name2(parameters, Parameters_name2(multi_parameter));

            for (int i = 0; i < Parameters_count(multi_parameter); i++) {
                Hashtable_add(hash, Parameters_name(multi_parameter, i), Parameters_at(multi_parameter, i));
            }
            free_Parameters(multi_parameter);
            
            Hashtable_add(hash, Parameters_name2(parameters), parameters);
        }
        else{
            Parameter* p = new_Parameter_from_json(node, hash);
            Parameters_move(parameters, p);
            Parameters_set_name2(parameters, Parameter_name(p));

            Hashtable_add(hash, Parameter_name(p), p);
        }
    }
    else{
        fprintf(stderr, "Do not recognize node type of %s", node->key);
        exit(1);
    }
    return parameters;
}

Model* new_DistributionModel_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"burnin",
		"distribution",
		"file",
		"from",
		"margin",
		"parameters",
		"parameterization",
		"posterior",
		"tree",
		"x"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* d_string = get_json_node_value_string(node, "distribution");
	
	if (strcasecmp(d_string, "beta") == 0) {
        return new_BetaDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "betaprime") == 0) {
        return new_BetaPrimeDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "exponential") == 0) {
        return new_ExponentialDistributionModel_from_json(node, hash);
    }
	else if (strcasecmp(d_string, "gamma") == 0) {
		return new_GammaDistributionModel_from_json(node, hash);
	}
	else if (strcasecmp(d_string, "gmrf") == 0) {
		return new_GMRFModel_from_json(node, hash);
	}
	else if (strcasecmp(d_string, "dirichlet") == 0){
		return new_DirichletDistributionModel_from_json(node, hash);
	}
	else if (strcasecmp(d_string, "cauchy") == 0){
		return new_CauchyDistributionModel_from_json(node, hash);
	}
    else if (strcasecmp(d_string, "kumaraswamy") == 0){
        return new_KumaraswamyDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "lognormal") == 0){
        return new_LogNormalDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "multivariatenormal") == 0){
        return new_MultivariateNormalDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "normal") == 0 || strcasecmp(d_string, "gaussian") == 0){
        return new_NormalDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "halfnormal") == 0){
        return new_HalfNormalDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "oneonx") == 0){
        return new_OneOnXDistributionModel_from_json(node, hash);
    }
    else if (strcasecmp(d_string, "ctmcscale") == 0){
        return new_CTMCScaleModel_from_json(node, hash);
    }
    
	char* id = get_json_node_value_string(node, "id");
	
	Parameters** parameters = NULL;
    size_t parameters_dim = 0;
	Parameters* x = new_Parameters(1);
	DistributionModel* dm = NULL;
	Model* model = NULL;
	
	char* file = get_json_node_value_string(node, "file");
	Vector** samples = NULL;
	if (file != NULL){
		get_parameters_references2(node, hash, x, "x");
		size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
		samples = read_log_for_parameters_t(file, burnin, x);
	}
	
	if (strcasecmp(d_string, "topology") == 0) {
        char* ref = get_json_node_value_string(node, "tree");
        Model* mtree = Hashtable_get(hash, ref+1);
        dm = new_UniformTreeDistribution(mtree->obj);
        model = new_DistributionModel3(id, dm, mtree);
    }
	else if(strcasecmp(d_string, "copula") == 0){
		char* margin_str = get_json_node_value_string(node, "margin");
//		parameters = new_Parameters(1);
		
		if (samples != NULL) {
			int n = Vector_length(samples[0]);
			size_t paramCount = Parameters_count(x);
			double* means = dvector(paramCount);
			double* variances = dvector(paramCount);
			
			for (int i = 0; i < paramCount; i++) {
				const double* vec = Vector_data(samples[i]);
				means[i] = mean(vec, n);
				variances[i] = variance(vec, n, means[i]);
			}
			
//			size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
//			Vector* probs = read_log_column_with_id(file, burnin, "posterior");
//			
//			double** mat = dmatrix(paramCount, 6);
//			double*a = dvector(n);
//			double* params = dvector(2);
//			for (int i = 0; i < paramCount; i++) {
//				memcpy(a, Vector_data(samples[i]), sizeof(double)*n);
//				
//				double mode = 0;
//				double max = Vector_at(probs, 0);
//				for (int j = 1; j < n; j++) {
//					if(Vector_at(probs, j) > max){
//						max = Vector_at(probs, j);
//						mode = Vector_at(samples[i], j);
//					}
//				}
//				
//				params[0] = means[i]*means[i]/variances[i];
//				params[1] = means[i]/variances[i];
//				double mode_g = (params[0]-1.0)/params[1];
//				mat[i][0] = ks(a, n, _gamma_wrapper, params);
//				for (int j = 0; j < n; j++){
////					mat[i][3] += log(_gamma_wrapper(params, (params[0]-1.0)/params[1]));
//					mat[i][3] += log(_gamma_wrapper(params, mode));
//				}
//				mat[i][3] = 4.0 -2.0*mat[i][3];
//				
//				double mean2 = means[i]*means[i];
//				params[0] = log(mean2/sqrt(variances[i] + mean2));
//				params[1] = sqrt(log(1.0 + variances[i]/mean2));
//				double mode_ln = exp(params[0]-params[1]*params[1]);
//				mat[i][1] = ks(a, n, _lognormal_wrapper, params);
//				for (int j = 0; j < n; j++){
////					mat[i][4] += log(_lognormal_wrapper(params, exp(params[0]-params[1]*params[1])));
//					mat[i][4] += log(_lognormal_wrapper(params, mode));
//				}
//				mat[i][4] = 4.0 -2.0*mat[i][4];
//				
//				params[0] = means[i];
//				mat[i][2] = ks(a, n, _exponential_wrapper, params);
//				for (int j = 0; j < n; j++){
//					mat[i][5] += log(_exponential_wrapper(params, mode));
//				}
//				mat[i][5] = 2.0 -2.0*mat[i][5];
//
//				printf("%f %f %f %f | %f %f %f %s %s | %f %f %f %s\n", Parameters_value(x, i), mode, mode_g, mode_ln,
//					   mat[i][0], mat[i][1], mat[i][2], (mat[i][0] > mat[i][2] ? "*":""), (mat[i][0] > mat[i][1] ? "+":""),
//					   mat[i][3], mat[i][4], mat[i][5], (mat[i][3] > mat[i][4] ? "+":""));
//			}
			
			
//			for (int i = 0; i < paramCount; i++) {
//				printf("%f %f | %f %f %f %s %s | %f %f %f %s\n", Parameters_value(x, i), mode, mat[i][0], mat[i][1], mat[i][2], (mat[i][0] > mat[i][2] ? "*":""), (mat[i][0] > mat[i][1] ? "+":""), mat[i][3], mat[i][4], mat[i][5], (mat[i][3] > mat[i][4] ? "+":""));
//			}
			
			json_node* from_node = get_json_node(node, "from");
            parameters = malloc(sizeof(Parameters*)*3);
            parameters[0] = new_Parameters(paramCount);
            parameters[1] = new_Parameters(paramCount);
            parameters[2] = new_Parameters(paramCount*paramCount);
			// calculate margins parameters
			if(from_node == NULL){
				if(strcasecmp(margin_str, "gamma") == 0){
					for (int i = 0; i < paramCount; i++) {
						Parameters_move(parameters[0], new_Parameter("alpha", means[i]*means[i]/variances[i], NULL));
						Parameters_move(parameters[1], new_Parameter("beta", means[i]/variances[i], NULL));
					}
				}
				else if(strcasecmp(margin_str, "lognormal") == 0){
					for (int i = 0; i < paramCount; i++) {
						double mean2 = means[i]*means[i];
						Parameters_move(parameters[0], new_Parameter("mu", log(mean2/sqrt(variances[i] + mean2)), NULL));
						Parameters_move(parameters[1], new_Parameter("s", sqrt(log(1.0 + variances[i]/mean2)), NULL));
					}
				}
			}
			else{
				char* ref = (char*)from_node->value;
				Model* simpleModel = Hashtable_get(hash, ref+1);
				DistributionModel* simpleDM = simpleModel->obj;
				for (int i = 0; i < paramCount; i++) {
                    Parameters_move(parameters[0], clone_Parameter(Parameters_at(simpleDM->parameters[i], i)));
                    Parameters_move(parameters[1], clone_Parameter(Parameters_at(simpleDM->parameters[i], i)));
				}
			}
		
			double* cov = dvector(paramCount*paramCount);
			// Calculate sample covariance matrix
			for (int i = 0; i < paramCount; i++) {
				const double* pp = Vector_data(samples[i]);
				cov[i*paramCount+i] = covariance(pp, pp, means[i], means[i], n);
				for (int j = i+1; j < paramCount; j++) {
					const double* pp2 = Vector_data(samples[j]);
					cov[i*paramCount+j] = cov[j*paramCount+i] = covariance(pp, pp2, means[i], means[j], n);
				}
			}
			
			// sample correlation matrix
			for (int i = 0; i < paramCount; i++) {
				double covii = sqrt(cov[i*paramCount + i]);
				for (int j = 0; j < paramCount; j++) {
					double covjj = sqrt(cov[j*paramCount + j]);
					double covij = cov[i*paramCount + j];
					Parameters_move(parameters[2], new_Parameter("ij", covij/(covii*covjj), NULL));
				}
			}
			free(cov);
			
			dm = new_CopulaDistributionModel_with_parameters(parameters, 3, x);
			
			if(strcasecmp(margin_str, "lognormal") == 0){
				dm->logP = _gaussian_copula_lognormal_logP;
				dm->logP_with_values = _gaussian_copula_lognormal_logP_with_values;
				dm->sample = _gaussian_copula_lognormal_sample;
			}
			
			free(means);
			free(variances);
		}
		else{
			get_parameters_references2(node, hash, x, "x");
            
            size_t paramCount = Parameters_count(x);
			parameters = malloc(sizeof(Parameters*)*3);
            parameters[0] = new_Parameters(paramCount);
            parameters[1] = new_Parameters(paramCount);
            parameters[2] = new_Parameters(paramCount*paramCount);
			
			json_node* from_node = get_json_node(node, "from");
			char* ref = (char*)from_node->value;
			Model* simpleModel = Hashtable_get(hash, ref+1);
			DistributionModel* simpleDM = simpleModel->obj;
			for (int i = 0; i < paramCount*2; i++) {
				Parameters_move(parameters[i], clone_Parameter(Parameters_at(simpleDM->parameters[i], i)));
			}
			
			json_node* posterior_node = get_json_node(node, "posterior");
			char* refm = (char*)posterior_node->value;
			Model* mpost = Hashtable_get(hash, refm+1);
			gsl_matrix* H = gsl_matrix_alloc(paramCount, paramCount);
			gsl_matrix* L = gsl_matrix_alloc(paramCount, paramCount);
			gsl_permutation* perm = gsl_permutation_alloc(paramCount);
			for (int i = 0; i < paramCount; i++) {
				double mapi = Parameters_value(x, i);
				double d2logP = mpost->d2logP(mpost, Parameters_at(x, i));
				gsl_matrix_set(H, i, i, d2logP);
				for (int j = i+1; j < paramCount; j++) {
					double didj = mpost->ddlogP(mpost, Parameters_at(x, i), Parameters_at(x, j));
					gsl_matrix_set(H, i, j, didj);
					gsl_matrix_set(H, j, i, didj);
				}
			}
			int signum;
			gsl_linalg_LU_decomp (H, perm, &signum);
			
			gsl_linalg_LU_invert (H, perm, L);
			for (int i = 0; i < paramCount; i++) {
				double covii = sqrt(-gsl_matrix_get(L, i, i));
				for (int j = 0; j < paramCount; j++) {
					double covjj = sqrt(-gsl_matrix_get(L, j, j));
					double covij = -gsl_matrix_get(L, i, j);
					Parameters_move(parameters[2], new_Parameter("ij", covij/(covii*covjj), NULL));
				}
			}
			gsl_permutation_free(perm);
			gsl_matrix_free(H);
			gsl_matrix_free(L);
			
			dm = new_CopulaDistributionModel_with_parameters(parameters, 3, x);
			
			if(strcasecmp(margin_str, "lognormal") == 0){
				dm->logP = _gaussian_copula_lognormal_logP;
				dm->logP_with_values = _gaussian_copula_lognormal_logP_with_values;
				dm->sample = _gaussian_copula_lognormal_sample;
			}
		}
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->sample_evaluate = _dist_model_sample_evaluate;
		model->samplable = true;
	}
	else if(strcasecmp(d_string, "oneonx") == 0){
		get_parameters_references2(node, hash, x, "x");
		dm = new_OneOnXDistributionModel(x);
		model = new_DistributionModel2(id, dm);
	}
	else{
		printf("Distribution unknown: %s\n", d_string);
		exit(10);
	}
	
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	if (samples != NULL) {
		size_t paramCount = Parameters_count(x);
		for (int i = 0; i < paramCount; i++) {
			free_Vector(samples[i]);
		}
		free(samples);
	}
    if(parameters != NULL){
        for(int i = 0; i < parameters_dim; i++) free_Parameters(parameters[i]);
        free(parameters);
    }
	free_Parameters(x);
	
	return model;
}
