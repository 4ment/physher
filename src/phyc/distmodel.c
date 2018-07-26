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

#include "matrix.h"
#include "tree.h"

#include "exponential.h"
#include "gamma.h"
#include "dirichlet.h"

#include "filereader.h"
#include "statistics.h"
#include "parametersio.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>


static double _DistributionModel_dlog_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

static double _DistributionModel_d2log_0(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

static double _DistributionModel_ddlog_0(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	return 0.0;
}

static void _DistributionModel_error_sample(DistributionModel* dm, double* samples){
	fprintf(stderr, "Sample method not implemented\n");
	exit(1);
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
	clone->ddlogP = dm->ddlogP;
	clone->sample = dm->sample;
	clone->clone = dm->clone;
	clone->free = dm->free;
	clone->logP_with_values= dm->logP_with_values;
	clone->tempp = NULL;
	clone->tempx = NULL;
	if(dm->tempp != NULL) clone->tempp = clone_dvector(dm->tempp, Parameters_count(dm->parameters));
	if(dm->tempx != NULL) clone->tempx = clone_dvector(dm->tempx, Parameters_count(dm->x));
	clone->rng = dm->rng;
	clone->data = NULL;
	clone->lp = dm->lp;
	clone->need_update = dm->need_update;
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
	clone->ddlogP = dm->ddlogP;
	clone->sample = dm->sample;
	clone->logP_with_values = dm->logP_with_values;
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
	dm->ddlogP = NULL;
	dm->sample = _DistributionModel_error_sample;
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
	dm->ddlogP = NULL;
	dm->sample = _DistributionModel_error_sample;
	dm->free = _free_dist;
	dm->clone = _clone_dist;
	dm->data = NULL;
	dm->tempx = NULL;
	dm->tempp = NULL;
	return dm;
}

double DistributionModel_log_gamma(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			double x = Parameters_value(dm->x, i);
//			logP += gammln(alpha) - log(pow(beta, alpha)*pow(x,alpha-1.0)) - beta*x;
			logP += log(gsl_ran_gamma_pdf(x, alpha, 1.0/beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		logP = -gammln(alpha) * Parameters_count(dm->x);
		double beta_alpha = pow(beta, alpha);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i);
			logP += log(beta_alpha*pow(x,alpha-1.0)) - beta*x;
	//		logP += dloggamma(Parameters_value(dm->x, i), alpha, beta);
		}
	}
	return logP;
}

double DistributionModel_log_gamma_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			double x = values[i];
			//			logP += gammln(alpha) - log(pow(beta, alpha)*pow(x,alpha-1.0)) - beta*x;
			logP += log(gsl_ran_gamma_pdf(x, alpha, 1.0/beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		logP = -gammln(alpha) * Parameters_count(dm->x);
		double beta_alpha = pow(beta, alpha);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = values[i];
			logP += log(beta_alpha*pow(x,alpha-1.0)) - beta*x;
			//		logP += dloggamma(Parameters_value(dm->x, i), alpha, beta);
		}
	}
	return logP;
}

double DistributionModel_dlog_gamma(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters, 0);
			double beta = Parameters_value(dm->parameters, 1);
			if(Parameters_count(dm->parameters) > 2){
				alpha = Parameters_value(dm->parameters, i*2);
				beta = Parameters_value(dm->parameters, i*2+1);
			}
			return (alpha-1.0)/Parameter_value(p) - beta;
			
		}
	}
	return 0;
}

double DistributionModel_d2log_gamma(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters, 0);
			if(Parameters_count(dm->parameters) > 2){
				alpha = Parameters_value(dm->parameters, i*2);
			}
			return -(alpha-1.0)/Parameter_value(p)/Parameter_value(p);
		}
	}
	return 0;
}

static void DistributionModel_gamma_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_gamma(dm->rng, Parameters_value(dm->parameters, i*2), 1.0/Parameters_value(dm->parameters, i*2+1));
		}
	}
	else{
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_gamma(dm->rng, Parameters_value(dm->parameters, 0), 1.0/Parameters_value(dm->parameters, 1));
		}
	}
}

double DistributionModel_lognormal_logP(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double mu = Parameters_value(dm->parameters, i*2);
			double sigma = Parameters_value(dm->parameters, i*2+1);
			double x = Parameters_value(dm->x, i);
			logP += log(gsl_ran_lognormal_pdf(x, mu, sigma));
		}
	}
	else{
		double mu = Parameters_value(dm->parameters, 0);
		double sigma = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i);
			logP += log(gsl_ran_lognormal_pdf(x, mu, sigma));
		}
	}
	return logP;
}

double DistributionModel_lognormal_logP_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double mu = Parameters_value(dm->parameters, i*2);
			double sigma = Parameters_value(dm->parameters, i*2+1);
			logP += log(gsl_ran_lognormal_pdf(values[i], mu, sigma));
		}
	}
	else{
		double mu = Parameters_value(dm->parameters, 0);
		double sigma = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(gsl_ran_lognormal_pdf(values[i], mu, sigma));
		}
	}
	return logP;
}

double DistributionModel_lognormal_dlogP(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double mu = Parameters_value(dm->parameters, 0);
			double sigma = Parameters_value(dm->parameters, 1);
			if(Parameters_count(dm->parameters) > 2){
				mu = Parameters_value(dm->parameters, i*2);
				sigma = Parameters_value(dm->parameters, i*2+1);
			}
			error("DistributionModel_lognormal_dlogP not yet implemented");
			return INFINITY;
			
		}
	}
	return 0;
}

double DistributionModel_lognormal_d2logP(DistributionModel* dm, const Parameter* p){
	error("DistributionModel_lognormal_d2logP not yet implemented");
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double alpha = Parameters_value(dm->parameters, 0);
			if(Parameters_count(dm->parameters) > 2){
				alpha = Parameters_value(dm->parameters, i*2);
			}
			return -(alpha-1.0)/Parameter_value(p)/Parameter_value(p);
		}
	}
	return 0;
}

static void DistributionModel_lognormal_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters, i*2), Parameters_value(dm->parameters, i*2+1));
		}
	}
	else{
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters, 0), Parameters_value(dm->parameters, 1));
		}
	}
}

double DistributionModel_log_exp(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = Parameters_value(dm->parameters, i);
			logP += lambda - lambda * Parameters_value(dm->x, i);
		}
	}
	else{
		double lambda = Parameters_value(dm->parameters, 0);
		logP = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			//logP += log(dexp(Parameters_value(dm->x, i), Parameters_value(dm->parameters, 0)));
			logP -= lambda * Parameters_value(dm->x, i);
		}
	}
	return logP;
}

double DistributionModel_log_exp_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double lambda = Parameters_value(dm->parameters, i);
			logP += lambda - lambda * values[i];
		}
	}
	else{
		double lambda = Parameters_value(dm->parameters, 0);
		logP = log(lambda) * Parameters_count(dm->x);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			//logP += log(dexp(Parameters_value(dm->x, i), Parameters_value(dm->parameters, 0)));
			logP -= lambda * values[i];
		}
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

static void DistributionModel_exp_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 1){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = rexp(Parameters_value(dm->parameters, i));
		}
	}
	else{
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = rexp(Parameters_value(dm->parameters, 0));
		}
	}
}

// Flat dirichlet
double DistributionModel_log_flat_dirichlet(DistributionModel* dm){
	return log(ddirchlet_flat(dm->simplex->K));
}

double DistributionModel_log_flat_dirichlet_with_values(DistributionModel* dm, const double* values){
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

double DistributionModel_log_dirichlet_with_values(DistributionModel* dm, const double* values){
	return log(ddirchlet(values, Parameters_count(dm->x), dm->tempp));
}
static void DistributionModel_dirichlet_sample(DistributionModel* dm, double* samples){
	rdirichlet(samples, dm->simplex->K, dm->tempp);
}

//TODO: implement
double DistributionModel_dlog_dirichlet(DistributionModel* dm, const Parameter* p){
	// find corresponding alpha
	// return (alpha-1.0)/Parameter_value(p);
	exit(1);
	return 0;
}

//TODO: implement
double DistributionModel_d2log_dirichlet(DistributionModel* dm, const Parameter* p){
	// find corresponding alpha
	// return -(alpha-1.0)/(Parameter_value(p)*Parameter_value(p));
	exit(1);
	return 0;
}

DistributionModel* new_IndependantGammaDistributionModel(const double shape, const double rate, const Parameters* x){
	Parameters* ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter("gamma.shape", shape, NULL));
	Parameters_move(ps, new_Parameter("gamma.rate", rate, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->logP = DistributionModel_log_gamma;
	dm->logP_with_values = DistributionModel_log_gamma_with_values;
	dm->dlogP = DistributionModel_dlog_gamma;
	dm->d2logP = DistributionModel_d2log_gamma;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_gamma_sample;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantGammaDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->logP = DistributionModel_log_gamma;
	dm->logP_with_values = DistributionModel_log_gamma_with_values;
	dm->dlogP = DistributionModel_dlog_gamma;
	dm->d2logP = DistributionModel_d2log_gamma;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_gamma_sample;
	return dm;
}


DistributionModel* new_IndependantLognormalDistributionModel(const double shape, const double rate, const Parameters* x){
	Parameters* ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter("gamma.shape", shape, NULL));
	Parameters_move(ps, new_Parameter("gamma.rate", rate, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->logP = DistributionModel_lognormal_logP;
	dm->logP_with_values = DistributionModel_lognormal_logP_with_values;
	dm->dlogP = DistributionModel_lognormal_dlogP;
	dm->d2logP = DistributionModel_lognormal_d2logP;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_lognormal_sample;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantLognormalDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->logP = DistributionModel_lognormal_logP;
	dm->logP_with_values = DistributionModel_lognormal_logP_with_values;
	dm->dlogP = DistributionModel_lognormal_dlogP;
	dm->d2logP = DistributionModel_lognormal_d2logP;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_lognormal_sample;
	return dm;
}

DistributionModel* new_IndependantExpDistributionModel(const double lambda, const Parameters* x){
	Parameters* ps = new_Parameters(1);
	Parameters_move(ps, new_Parameter("exp.lambda", lambda, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->logP = DistributionModel_log_exp;
	dm->logP_with_values = DistributionModel_log_exp_with_values;
	dm->dlogP = DistributionModel_dlog_exp;
	dm->d2logP = DistributionModel_d2log_exp;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_exp_sample;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantExpDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->logP = DistributionModel_log_exp;
	dm->logP_with_values = DistributionModel_log_exp_with_values;
	dm->dlogP = DistributionModel_dlog_exp;
	dm->d2logP = DistributionModel_d2log_exp;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_exp_sample;
	return dm;
}

DistributionModel* new_FlatDirichletDistributionModel(Simplex* simplex){
	DistributionModel* dm = new_DistributionModelSimplex(NULL, simplex);
	dm->logP = DistributionModel_log_flat_dirichlet;
	dm->logP_with_values = DistributionModel_log_flat_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_flat_dirichlet;
	dm->d2logP = DistributionModel_d2log_flat_dirichlet;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i <simplex->K; i++) {
		dm->tempp[i] = 1;
	}
	return dm;
}

DistributionModel* new_DirichletDistributionModel(const double* alpha, Simplex* simplex){
	Parameters* ps = new_Parameters(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		Parameters_move(ps, new_Parameter("dirichlet.", alpha[i], NULL));
	}
	DistributionModel* dm = new_DistributionModelSimplex(ps, simplex);// len(x)==len(alpha)
	dm->logP = DistributionModel_log_dirichlet;
	dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_dirichlet_sample;
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
	dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = _DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_dirichlet_sample;
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


static double _DistributionModel_dlog_mvn(DistributionModel* dm, const Parameter* p){
	fprintf(stderr, "%s is not implemented (line %d of file %s)\n", __func__, __LINE__, __FILE__);
	exit(1);
	return 0.0;
}

static double _DistributionModel_d2log_mvn(DistributionModel* dm, const Parameter* p){
	fprintf(stderr, "%s is not implemented (line %d of file %s)\n", __func__, __LINE__, __FILE__);
	exit(1);
	return 0.0;
}

static double _DistributionModel_ddlog_mvn(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	fprintf(stderr, "%s is not implemented (line %d of file %s)\n", __func__, __LINE__, __FILE__);
	exit(1);
	return 0.0;
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

DistributionModel* new_MultivariateNormalDistributionModel_with_parameters(const double* mu, const double* sigma, const Parameters* x){
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
	dm->dlogP = _DistributionModel_dlog_mvn;
	dm->d2logP = _DistributionModel_d2log_mvn;
	dm->ddlogP = _DistributionModel_ddlog_mvn;
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
	dm->data = wrapper;
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
		double z = gsl_cdf_gamma_Pinv(u, Parameters_value(dm->parameters, i*2), 1.0/Parameters_value(dm->parameters, i*2+1)); // gamma distributed
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
		double z = gsl_cdf_lognormal_Pinv(u, Parameters_value(dm->parameters, i*2), Parameters_value(dm->parameters, i*2+1)); // lognormal distributed
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
		logP += log(gsl_ran_gamma_pdf(values[i], Parameters_value(dm->parameters, i*2), 1.0/Parameters_value(dm->parameters, i*2+1)));
	}
	return logP;
}

static double _gaussian_copula_lognormal_logP_with_values(DistributionModel* dm, const double* values){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	double logP = -log(sqrt(wrapper->det));
	for (int i = 0; i < paramCount; i++) {
		logP += log(gsl_ran_lognormal_pdf(values[i], Parameters_value(dm->parameters, i*2), Parameters_value(dm->parameters, i*2+1)));
	}
	return logP;
}

static double _gaussian_copula_gamma_logP(DistributionModel* dm){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	double logP = -log(sqrt(wrapper->det));
	for (int i = 0; i < paramCount; i++) {
		logP += log(gsl_ran_gamma_pdf(Parameters_value(dm->x, i), Parameters_value(dm->parameters, i*2), 1.0/Parameters_value(dm->parameters, i*2+1)));
	}
	return logP;
}

static double _gaussian_copula_lognormal_logP(DistributionModel* dm){
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)dm->data;
	size_t paramCount = Parameters_count(dm->x);
	double logP = -log(sqrt(wrapper->det));
	for (int i = 0; i < paramCount; i++) {
		logP += log(gsl_ran_lognormal_pdf(Parameters_value(dm->x, i), Parameters_value(dm->parameters, i*2), Parameters_value(dm->parameters, i*2+1)));
	}
	return logP;
}

static DistributionModel* _clone_gaussian_copula_gamma(DistributionModel* dm){
	DistributionModel* clone = _clone_dist(dm);
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
	if(dm->parameters != NULL) free_Parameters(dm->parameters);
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
	size_t offset = paramCount*2;
	
	for (size_t i = 0; i < paramCount; i++) {
		for (size_t j = 0; j < paramCount; j++) {
			gsl_matrix_set(wrapper->cor, i, j, Parameters_value(dm->parameters, offset+i*paramCount+j));
		}
	}
	
	int signum;
	gsl_linalg_LU_decomp(wrapper->cor, wrapper->p , &signum);
	return gsl_linalg_LU_det(wrapper->cor, signum);
}

DistributionModel* new_CopulaDistributionModel_with_parameters(const Parameters* p, const Parameters* x){
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
	dm->free = _free_gaussian_copula_gamma;
	dm->clone = _clone_gaussian_copula_gamma;
	dm->tempx = NULL;
	dm->tempp = NULL;
	dm->logP = _gaussian_copula_gamma_logP;
	dm->logP_with_values = _gaussian_copula_gamma_logP_with_values;
	//TODO: copula methods
	dm->dlogP = _DistributionModel_dlog_mvn;
	dm->d2logP = _DistributionModel_d2log_mvn;
	dm->ddlogP = _DistributionModel_ddlog_mvn;
	dm->sample = _gaussian_copula_gamma_sample;
	dm->need_update = true;
	
	gsl_gaussian_copula_wrapper_t* wrapper = (gsl_gaussian_copula_wrapper_t*)malloc(sizeof(gsl_gaussian_copula_wrapper_t));
	size_t paramCount = Parameters_count(x);
	wrapper->mu = gsl_vector_calloc(paramCount);
	wrapper->L = gsl_matrix_calloc(paramCount, paramCount);
	wrapper->cor = gsl_matrix_calloc(paramCount, paramCount);
	wrapper->x = gsl_vector_calloc(paramCount);
	wrapper->work = gsl_vector_calloc(paramCount);
	wrapper->p = gsl_permutation_alloc(paramCount);
	size_t offset = paramCount*2;
	
	for (int i = 0; i < paramCount; i++) {
		gsl_vector_set(wrapper->mu, i, 0);
		for (int j = 0; j < paramCount; j++) {
			gsl_matrix_set(wrapper->cor, i, j, Parameters_value(dm->parameters, offset+i*paramCount+j));
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
	dm->ddlogP = _DistributionModel_ddlog_0;
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
	return 0;
}

static void _dist_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free distribution model %s\n", self->name);
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
	clone->sample = self->sample;
	clone->samplable = self->samplable;
	return clone;
}

static void _dist_model_get_free_parameters(Model* model, Parameters* parameters){
	DistributionModel* dm = (DistributionModel*)model->obj;
	
	Model* msimplex = (Model*)model->data;
}

static void _dist_model_sample(Model* model, double* samples, double* logP){
	DistributionModel* dm = (DistributionModel*)model->obj;
	dm->sample(dm, samples);
	if(logP != NULL){
		*logP = dm->logP_with_values(dm, samples);
	}
}

Model* new_DistributionModel2(const char* name, DistributionModel* dm){
	Model *model = new_Model("distribution",name, dm);
	model->logP = _dist_model_logP;
	model->dlogP = _dist_model_dlogP;
	model->d2logP = _dist_model_d2logP;
	model->ddlogP = _dist_model_ddlogP;
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
	model->ddlogP = _dist_model_ddlogP;
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
	model->ddlogP = _dist_model_ddlogP;
    model->free = _dist_model_free;
    model->clone = _dist_model_clone;
    model->get_free_parameters = _dist_model_get_free_parameters;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
    return model;
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
	Parameters* x = new_Parameters(1);
	DistributionModel* dm = NULL;
	Model* model = NULL;
	
	char* file = get_json_node_value_string(node, "file");
	Vector** samples = NULL;
	if (file != NULL) {
		get_parameters_references2(node, hash, x, "x");
		size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
		samples = read_log_for_parameters_t(file, burnin, x);
	}
	
	if (strcasecmp(d_string, "exponential") == 0) {
		parameters = new_Parameters(1);
		get_parameters_references(node, hash, parameters);
		
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
			get_parameters_references2(node, hash, x, "x");
		}
		dm = new_IndependantExpDistributionModel_with_parameters(parameters, x);
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->samplable = true;
	}
	else if (strcasecmp(d_string, "dirichlet") == 0) {
		parameters = new_Parameters(1);
		get_parameters_references(node, hash, parameters);
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
		model = new_DistributionModel3(id, dm, simplex);
		model->sample = _dist_model_sample;
		model->samplable = true;
    }
    else if (strcasecmp(d_string, "gamma") == 0) {
        parameters = new_Parameters(1);
		// empirical
		if (samples != NULL) {
			size_t paramCount = Parameters_count(x);
			
			for (int i = 0; i < paramCount; i++) {
				double* vec = Vector_data(samples[i]);
				double m = mean(vec, Vector_length(samples[i]));
				double v = variance(vec, Vector_length(samples[i]), m);
				Parameters_move(parameters, new_Parameter("alpha", m*m/v, NULL));
				Parameters_move(parameters, new_Parameter("beta", m/v, NULL));
			}
		}
		else if(get_json_node(node, "parameters") == NULL){
			get_parameters_references2(node, hash, x, "x");
			for (int i = 0; i < Parameters_count(x); i++) {
				Parameters_move(parameters, new_Parameter("alpha", 1, new_Constraint(0, INFINITY)));
				Parameters_move(parameters, new_Parameter("beta", 1, new_Constraint(0, INFINITY)));
			}
		}
		else{
			get_parameters_references2(node, hash, x, "x");
			get_parameters_references(node, hash, parameters);
		}
		
        dm = new_IndependantGammaDistributionModel_with_parameters(parameters, x);
        model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->samplable = true;
    }
	else if (strcasecmp(d_string, "lognormal") == 0) {
		parameters = new_Parameters(1);
		// empirical
		if (samples != NULL) {
			size_t paramCount = Parameters_count(x);
			
			for (int i = 0; i < paramCount; i++) {
				double* vec = Vector_data(samples[i]);
				double m = mean(vec, Vector_length(samples[i]));
				double v = variance(vec, Vector_length(samples[i]), m);
				Parameters_move(parameters, new_Parameter("mu", log(m*m/sqrt(v + m*m)), NULL));
				Parameters_move(parameters, new_Parameter("sigma", sqrt(log(1.0 + v/(m*m))), NULL));
			}
		}
		else if(get_json_node(node, "parameters") == NULL){
			get_parameters_references2(node, hash, x, "x");
			for (int i = 0; i < Parameters_count(x); i++) {
				Parameters_move(parameters, new_Parameter("mu", 1, new_Constraint(0, INFINITY)));
				Parameters_move(parameters, new_Parameter("sigma", 1, new_Constraint(0, INFINITY)));
			}
		}
		else{
			get_parameters_references2(node, hash, x, "x");
			get_parameters_references(node, hash, parameters);
		}
		
		dm = new_IndependantLognormalDistributionModel_with_parameters(parameters, x);
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->samplable = true;
	}
    else if (strcasecmp(d_string, "topology") == 0) {
        char* ref = get_json_node_value_string(node, "tree");
        Model* mtree = Hashtable_get(hash, ref+1);
        dm = new_UniformTreeDistribution(mtree);
        model = new_TreeDistributionModel(id, dm, mtree);
    }
	else if(strcasecmp(d_string, "normal") == 0){
		char* file = get_json_node_value_string(node, "file");
		parameters = new_Parameters(1);
		get_parameters_references2(node, hash, x, "x");
		
		if (file != NULL) {
			int n = Vector_length(samples[0]);
			size_t paramCount = Parameters_count(x);
			double* mu = dvector(paramCount);
			double* sigma = dvector(paramCount*paramCount);
			
			for (int i = 0; i < paramCount; i++) {
				double* vv = Vector_data(samples[i]);
				for (int j = 0; j < n; j++) {
					vv[j] = log(vv[j]);
				}
				mu[i] = mean(vv, n);
			}
			
			// Calculate covariance matrix
			for (int i = 0; i < paramCount; i++) {
				double* pp = Vector_data(samples[i]);
				sigma[i*paramCount+i] = variance(pp, n, mu[i]);
				for (int j = i+1; j < paramCount; j++) {
					double* pp2 = Vector_data(samples[j]);
					sigma[i*paramCount+j] = sigma[j*paramCount+i] = covariance(pp, pp2, mu[i], mu[j], n);
				}
			}
			
			dm = new_MultivariateNormalDistributionModel_with_parameters(mu, sigma, x);
			
			free(mu);
			free(sigma);
		}
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->samplable = true;
		
		gsl_multivariate_normal_wrapper_t* wrapper = dm->data;
		wrapper->transform = true;
	}
	else if(strcasecmp(d_string, "copula") == 0){
		char* margin_str = get_json_node_value_string(node, "margin");
		parameters = new_Parameters(1);
		
		if (samples != NULL) {
			int n = Vector_length(samples[0]);
			size_t paramCount = Parameters_count(x);
			double* means = dvector(paramCount);
			double* variances = dvector(paramCount);
			
			for (int i = 0; i < paramCount; i++) {
				double* vec = Vector_data(samples[i]);
				means[i] = mean(vec, n);
				variances[i] = variance(vec, n, means[i]);
			}
			
			// calculate margins parameters
			if(strcasecmp(margin_str, "gamma") == 0){
				for (int i = 0; i < paramCount; i++) {
					Parameters_move(parameters, new_Parameter("alpha", means[i]*means[i]/variances[i], NULL));
					Parameters_move(parameters, new_Parameter("beta", means[i]/variances[i], NULL));
				}
			}
			else if(strcasecmp(margin_str, "lognormal") == 0){
				for (int i = 0; i < paramCount; i++) {
					double mean2 = means[i]*means[i];
					Parameters_move(parameters, new_Parameter("mu", log(mean2/sqrt(variances[i] + mean2)), NULL));
					Parameters_move(parameters, new_Parameter("s", sqrt(log(1.0 + variances[i]/mean2)), NULL));
				}
			}
			
			// Calculate sample covariance matrix
			for (int i = 0; i < paramCount; i++) {
				double* pp = Vector_data(samples[i]);
				for (int j = 0; j < paramCount; j++) {
					double* pp2 = Vector_data(samples[j]);
					double covij = covariance(pp, pp2, means[i], means[j], n)/(sqrt(variances[i])*sqrt(variances[j]));
					Parameters_move(parameters, new_Parameter("ij", covij, NULL));
				}
			}
			
			// sample correlation matrix
			for (int i = 0; i < paramCount; i++) {
				double covii = sqrt(Parameters_value(parameters, 2*paramCount + i*paramCount + i));
				for (int j = 0; j < paramCount; j++) {
					double covjj = sqrt(Parameters_value(parameters, 2*paramCount + j*paramCount + j));
					double covij = Parameters_value(parameters, 2*paramCount + i*paramCount + j);
					Parameters_set_value(parameters, 2*paramCount + i*paramCount + j, covij/(covii*covjj));
				}
			}
			
			dm = new_CopulaDistributionModel_with_parameters(parameters, x);
			
			if(strcasecmp(margin_str, "lognormal") == 0){
				dm->logP = _gaussian_copula_lognormal_logP;
				dm->logP_with_values = _gaussian_copula_lognormal_logP_with_values;
				dm->sample = _gaussian_copula_lognormal_sample;
			}
			
			free(means);
			free(variances);
		}
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->samplable = true;
	}
	else{
		printf("%s\n", d_string);
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
	free_Parameters(parameters);
	free_Parameters(x);
	
	return model;
}
