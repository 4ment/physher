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
#include <tgmath.h>

#include "matrix.h"
#include "tree.h"

#include "exponential.h"
#include "gamma.h"
#include "dirichlet.h"

#include "filereader.h"
#include "statistics.h"
#include "parametersio.h"
#include "utilsio.h"

#include "distexp.h"
#include "distgamma.h"
#include "gmrf.h"

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
	if(dm->parameters != NULL) free_Parameters(dm->parameters);
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	// freeing data is left to the user
	free(dm);
}

static void _free_full_distribution(DistributionModel*dm){
	if(dm->x != NULL) free_Parameters(dm->x);
	if(dm->parameters != NULL) free_Parameters(dm->parameters);
	if(dm->tempx != NULL) free(dm->tempx);
	if(dm->tempp != NULL) free(dm->tempp);
	if(dm->simplex != NULL) free_Simplex(dm->simplex);
	if(dm->tree != NULL) free_Tree(dm->tree);
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
	clone->sample_evaluate = dm->sample_evaluate;
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
	clone->type = dm->type;
	return clone;
}

DistributionModel* clone_DistributionModel_with_parameters(DistributionModel* dm, const Parameters* params, const Parameters* x, Simplex* simplex){
	DistributionModel* clone = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(clone);
	clone->type = dm->type;
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
	clone->shift = dm->shift;
	clone->type = dm->type;
	clone->parameterization = dm->parameterization;
	clone->rng = dm->rng;
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
	dm->tree = NULL;
	dm->logP = NULL;
	dm->dlogP = NULL;
	dm->d2logP = NULL;
	dm->ddlogP = NULL;
	dm->sample = _DistributionModel_error_sample;
	dm->sample_evaluate = _DistributionModel_error_sample_evaluate;
	dm->free = _free_full_distribution;
	dm->clone = _clone_dist;
	dm->data = NULL;
	dm->tempx = NULL;
	dm->tempp = NULL;
	dm->need_update = true;
	return dm;
}

DistributionModel* new_DistributionModelSimplex(Parameters* p, Simplex* simplex){
	DistributionModel* dm = (DistributionModel*)malloc(sizeof(DistributionModel));
	assert(dm);
	dm->parameters = p;
	dm->x = NULL;
	dm->simplex = simplex;
	dm->tree = NULL;
	dm->logP = NULL;
	dm->dlogP = NULL;
	dm->d2logP = NULL;
	dm->ddlogP = NULL;
	dm->sample = _DistributionModel_error_sample;
	dm->sample_evaluate = _DistributionModel_error_sample_evaluate;
	dm->free = _free_full_distribution;
	dm->clone = _clone_dist;
	dm->data = NULL;
	dm->tempx = NULL;
	dm->tempp = NULL;
	dm->need_update = true;
	return dm;
}

double DistributionModel_log_one_on_x(DistributionModel* dm){
	return -log(Parameters_value(dm->x, 0));
}

double DistributionModel_log_one_on_x_with_values(DistributionModel* dm, const double* values){
	return -log(values[0]);
}

double DistributionModel_dlog_one_on_x(DistributionModel* dm, const Parameter* p){
	if (strcmp(Parameter_name(p), Parameters_name(dm->x, 0)) == 0) {
		return -1.0/Parameters_value(dm->x, 0);
	}
	return 0;
}

double DistributionModel_d2log_one_on_x(DistributionModel* dm, const Parameter* p){
	if (strcmp(Parameter_name(p), Parameters_name(dm->x, 0)) == 0) {
		return 1.0/Parameters_value(dm->x, 0)/Parameters_value(dm->x, 0);
	}
	return 0;
}

double DistributionModel_log_beta(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			double x = Parameters_value(dm->x, i);
			logP += log(gsl_ran_beta_pdf(x, alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i);
			logP += log(gsl_ran_beta_pdf(x, alpha, beta));
		}
	}
	return logP;
}

double DistributionModel_log_beta_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double alpha = Parameters_value(dm->parameters, i*2);
			double beta = Parameters_value(dm->parameters, i*2+1);
			logP += log(gsl_ran_beta_pdf(values[i], alpha, beta));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(gsl_ran_beta_pdf(values[i], alpha, beta));
		}
	}
	return logP;
}

double DistributionModel_dlog_beta(DistributionModel* dm, const Parameter* p){
	error("DistributionModel_dlog_beta not implemented\n");
	return 0;
}

double DistributionModel_d2log_beta(DistributionModel* dm, const Parameter* p){
	error("DistributionModel_d2log_beta not implemented\n");
	return 0;
}

static void DistributionModel_beta_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_beta(dm->rng, Parameters_value(dm->parameters, i*2), Parameters_value(dm->parameters, i*2+1));
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_beta(dm->rng, alpha, beta);
		}
	}
}

static double DistributionModel_beta_sample_evaluate(DistributionModel* dm){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_beta(dm->rng, Parameters_value(dm->parameters, i*2), Parameters_value(dm->parameters, i*2+1));
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		double alpha = Parameters_value(dm->parameters, 0);
		double beta = Parameters_value(dm->parameters, 1);
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_beta(dm->rng, alpha, beta);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_log_beta(dm);
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


static double DistributionModel_lognormal_sample_evaluate(DistributionModel* dm){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters, i*2), Parameters_value(dm->parameters, i*2+1));
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sample = gsl_ran_lognormal(dm->rng, Parameters_value(dm->parameters, 0), Parameters_value(dm->parameters, 1));
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_lognormal_logP(dm);
}

double DistributionModel_normal_logP(DistributionModel* dm){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double mu = Parameters_value(dm->parameters, i*2);
			double sigma = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
				sigma = sqrt(1.0/sigma);
			}
			double x = Parameters_value(dm->x, i);
			logP += log(gsl_ran_gaussian_pdf(x - mu, sigma));
		}
	}
	else{
		double mu = Parameters_value(dm->parameters, 0);
		double sigma = Parameters_value(dm->parameters, 1);
		if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
			sigma = sqrt(1.0/sigma);
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double x = Parameters_value(dm->x, i);
			logP += log(gsl_ran_gaussian_pdf(x - mu, sigma));
		}
	}
	return logP;
}

double DistributionModel_normal_logP_with_values(DistributionModel* dm, const double* values){
	double logP = 0;
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double mu = Parameters_value(dm->parameters, i*2);
			double sigma = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
				sigma = sqrt(1.0/sigma);
			}
			logP += log(gsl_ran_gaussian_pdf(values[i] - mu, sigma));
		}
	}
	else{
		double mu = Parameters_value(dm->parameters, 0);
		double sigma = Parameters_value(dm->parameters, 1);
		if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
			sigma = sqrt(1.0/sigma);
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			logP += log(gsl_ran_gaussian_pdf(values[i] - mu, sigma));
		}
	}
	return logP;
}

double DistributionModel_normal_dlogP(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
		if (strcmp(Parameter_name(p), Parameters_name(dm->x,i)) == 0) {
			double mu = Parameters_value(dm->parameters, 0);
			double sigma = Parameters_value(dm->parameters, 1);
			if(Parameters_count(dm->parameters) > 2){
				mu = Parameters_value(dm->parameters, i*2);
				sigma = Parameters_value(dm->parameters, i*2+1);
			}
			error("DistributionModel_normal_dlogP not yet implemented");
			return INFINITY;
			
		}
	}
	return 0;
}

double DistributionModel_normal_d2logP(DistributionModel* dm, const Parameter* p){
	error("DistributionModel_normal_d2logP not yet implemented");
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

static void DistributionModel_normal_sample(DistributionModel* dm, double* samples){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sigma = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
				sigma = sqrt(1.0/sigma);
			}
			samples[i] = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters, i*2);
		}
	}
	else{
		double mean = Parameters_value(dm->parameters, 0);
		double sigma = Parameters_value(dm->parameters, 1);
		if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
			sigma = sqrt(1.0/sigma);
		}
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			samples[i] = gsl_ran_gaussian(dm->rng, sigma) + mean;
		}
	}
}


static double DistributionModel_normal_sample_evaluate(DistributionModel* dm){
	if(Parameters_count(dm->parameters) > 2){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double sigma = Parameters_value(dm->parameters, i*2+1);
			if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
				sigma = sqrt(1.0/sigma);
			}
			double sample = gsl_ran_gaussian(dm->rng, sigma) + Parameters_value(dm->parameters, i*2);
			Parameters_set_value(dm->x, i, sample);
		}
	}
	else{
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			double mean = Parameters_value(dm->parameters, 0);
			double sigma = Parameters_value(dm->parameters, 1);
			if (dm->parameterization == DISTRIBUTION_NORMAL_MEAN_TAU) {
				sigma = sqrt(1.0/sigma);
			}
			double sample = gsl_ran_gaussian(dm->rng, sigma) + mean;
			Parameters_set_value(dm->x, i, sample);
		}
	}
	return DistributionModel_lognormal_logP(dm);
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
		return gsl_ran_dirichlet_lnpdf(dm->simplex->K, dm->tempp, dm->simplex->get_values(dm->simplex));
	}
	else{
		assert(0);
	}
	return log(ddirchlet(dm->tempx, Parameters_count(dm->x), dm->tempp));
}

double DistributionModel_log_dirichlet_with_values(DistributionModel* dm, const double* values){
	return gsl_ran_dirichlet_lnpdf(dm->simplex->K, dm->tempp, dm->simplex->get_values(dm->simplex));
}
static void DistributionModel_dirichlet_sample(DistributionModel* dm, double* samples){
	gsl_ran_dirichlet(dm->rng, dm->simplex->K, dm->tempp, samples);
//	rdirichlet(samples, dm->simplex->K, dm->tempp);
}


static double DistributionModel_dirichlet_sample_evaluate(DistributionModel* dm){
	double* samples = dvector(dm->simplex->K);
	gsl_ran_dirichlet(dm->rng, dm->simplex->K, dm->tempp, samples);
	double logP = DistributionModel_log_dirichlet_with_values(dm, samples);
	dm->simplex->set_values(dm->simplex, samples);
	free(samples);
	//	rdirichlet(samples, dm->simplex->K, dm->tempp);
	return logP;
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

DistributionModel* new_OneOnXDistributionModel(const Parameters* x){
	DistributionModel* dm = new_DistributionModel(NULL, x);
	dm->type = DISTRIBUTION_ONE_ON_X;
	dm->logP = DistributionModel_log_one_on_x;
	dm->logP_with_values = DistributionModel_log_one_on_x_with_values;
	dm->dlogP = DistributionModel_dlog_one_on_x;
	dm->d2logP = DistributionModel_d2log_one_on_x;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	return dm;
}

DistributionModel* new_IndependantBetaDistributionModel(const double alpha, const double beta, const Parameters* x){
	Parameters* ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter("beta.alpha", alpha, NULL));
	Parameters_move(ps, new_Parameter("beta.beta", beta, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->type = DISTRIBUTION_BETA;
	dm->logP = DistributionModel_log_beta;
	dm->logP_with_values = DistributionModel_log_beta_with_values;
	dm->dlogP = DistributionModel_dlog_beta;
	dm->d2logP = DistributionModel_d2log_beta;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_beta_sample;
	dm->sample_evaluate = DistributionModel_beta_sample_evaluate;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantBetaDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_BETA;
	dm->logP = DistributionModel_log_beta;
	dm->logP_with_values = DistributionModel_log_beta_with_values;
	dm->dlogP = DistributionModel_dlog_beta;
	dm->d2logP = DistributionModel_d2log_beta;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_beta_sample;
	dm->sample_evaluate = DistributionModel_beta_sample_evaluate;
	return dm;
}


DistributionModel* new_IndependantNormalDistributionModel(const double mean, const double sigma, const Parameters* x){
	Parameters* ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter("normal.mean", mean, NULL));
	Parameters_move(ps, new_Parameter("normal.sigma", sigma, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->type = DISTRIBUTION_NORMAL;
	dm->logP = DistributionModel_normal_logP;
	dm->logP_with_values = DistributionModel_normal_logP_with_values;
	dm->dlogP = DistributionModel_normal_dlogP;
	dm->d2logP = DistributionModel_normal_d2logP;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_normal_sample;
	dm->sample_evaluate = DistributionModel_normal_sample_evaluate;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantNormalDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_NORMAL;
	dm->logP = DistributionModel_normal_logP;
	dm->logP_with_values = DistributionModel_normal_logP_with_values;
	dm->dlogP = DistributionModel_normal_dlogP;
	dm->d2logP = DistributionModel_normal_d2logP;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_normal_sample;
	dm->sample_evaluate = DistributionModel_normal_sample_evaluate;
	return dm;
}

DistributionModel* new_IndependantLognormalDistributionModel(const double shape, const double rate, const Parameters* x){
	Parameters* ps = new_Parameters(2);
	Parameters_move(ps, new_Parameter("gamma.shape", shape, NULL));
	Parameters_move(ps, new_Parameter("gamma.rate", rate, NULL));
	DistributionModel* dm = new_DistributionModel(ps, x);
	dm->type = DISTRIBUTION_LOGNORMAL;
	dm->logP = DistributionModel_lognormal_logP;
	dm->logP_with_values = DistributionModel_lognormal_logP_with_values;
	dm->dlogP = DistributionModel_lognormal_dlogP;
	dm->d2logP = DistributionModel_lognormal_d2logP;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_lognormal_sample;
	dm->sample_evaluate = DistributionModel_lognormal_sample_evaluate;
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_IndependantLognormalDistributionModel_with_parameters(Parameters* parameters, const Parameters* x){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_LOGNORMAL;
	dm->logP = DistributionModel_lognormal_logP;
	dm->logP_with_values = DistributionModel_lognormal_logP_with_values;
	dm->dlogP = DistributionModel_lognormal_dlogP;
	dm->d2logP = DistributionModel_lognormal_d2logP;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_lognormal_sample;
	dm->sample_evaluate = DistributionModel_lognormal_sample_evaluate;
	return dm;
}

DistributionModel* new_FlatDirichletDistributionModel(Simplex* simplex){
	DistributionModel* dm = new_DistributionModelSimplex(NULL, simplex);
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_flat_dirichlet;
	dm->logP_with_values = DistributionModel_log_flat_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_flat_dirichlet;
	dm->d2logP = DistributionModel_d2log_flat_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
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
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_dirichlet;
	dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
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
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_dirichlet;
	dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->clone = _clone_dist;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
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
	dm->type = DISTRIBUTION_NORMAL_MULTIVARIATE;
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
	dm->tree = NULL;
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
    dm->free = _free_full_distribution;
    dm->clone = _clone_dist;
    dm->tree = tree;
    dm->tempx = NULL;
    dm->tempp = NULL;
	dm->logP = DistributionModel_log_uniform_tree;
	dm->dlogP = DistributionModel_dlog_0;
	dm->d2logP = DistributionModel_d2log_0;
	dm->ddlogP = DistributionModel_ddlog_0;
    dm->clone = _clone_dist;
	dm->need_update = true;
    return dm;
}

void _dist_model_handle_change( Model *self, Model *model, int index ){
	DistributionModel* dm = self->obj;
	dm->need_update = true;
	self->listeners->fire( self->listeners, self, index );
}

void _dist_model_handle_restore( Model *self, Model *model, int index ){
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _dist_model_store(Model* self){
	self->storedLogP = self->lp;
	DistributionModel* dm = self->obj;
	dm->stored_lp = dm->lp;
	Parameters_store(dm->parameters);
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
	for (int i = 0; i < Parameters_count(dm->parameters); i++) {
		p = Parameters_at(dm->parameters, i);
		if (Parameter_changed(p)) {
			changed = true;
			Parameter_restore_quietly(p);
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
	for (int i = 0; i < Parameters_count(cm->parameters); i++) {
		if(Parameters_at(cm->parameters, i) == p){
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
	for (int i = 0; i < Parameters_count(cm->parameters); i++) {
		if(Parameters_at(cm->parameters, i) == p){
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
	
	found1 = false;
	found2 = false;
	for (int i = 0; i < Parameters_count(cm->parameters); i++) {
		if(Parameters_at(cm->parameters, i) == p1){
			found1 = true;
		}
		else if(Parameters_at(cm->parameters, i) == p2){
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
		// tree or simplex
		if(self->data != NULL){
			Model* model = self->data;
			model->free(model);
		}
		_free_partial_distribution(cm);
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
	model->get_free_parameters = _dist_model_get_free_parameters;
	model->store = _dist_model_store;
	model->restore = _dist_model_restore;
	model->update = _dist_model_handle_change;
	model->handle_restore = _dist_model_handle_restore;
	model->sample = _dist_model_sample;
	model->sample_evaluate = _dist_model_sample_evaluate;
	model->samplable = false;
	
	for ( int i = 0; i < Parameters_count(dm->parameters); i++ ) {
		Parameters_at(dm->parameters, i)->listeners->add( Parameters_at(dm->parameters, i)->listeners, model );
	}
	for ( int i = 0; i < Parameters_count(dm->x); i++ ) {
		Parameters_at(dm->x, i)->listeners->add( Parameters_at(dm->x, i)->listeners, model );
	}
	return model;
}

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* simplex){
	Model *model = new_Model(MODEL_DISTRIBUTION,name, dm);
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
	model->update = _dist_model_handle_change;
	model->handle_restore = _dist_model_handle_restore;
	model->sample = _dist_model_sample;
	model->sample_evaluate = _dist_model_sample_evaluate;
	model->samplable = false;
	
	for ( int i = 0; i < Parameters_count(dm->parameters); i++ ) {
		Parameters_at(dm->parameters, i)->listeners->add( Parameters_at(dm->parameters, i)->listeners, model );
	}
	simplex->listeners->add(simplex->listeners, model);
	return model;
}

Model* get_simplex(json_node* node, Hashtable* hash){
	json_node* x_node = get_json_node(node, "x");
	char* ref = (char*)x_node->value;
	return Hashtable_get(hash, ref+1);
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
	
	if (strcasecmp(d_string, "exponential") == 0) {
		return new_ExponentialDistributionModel_from_json(node, hash);
	}
	else if (strcasecmp(d_string, "gamma") == 0) {
		return new_GammaDistributionModel_from_json(node, hash);
	}
	else if (strcasecmp(d_string, "gmrf") == 0) {
		return new_GMRFModel_from_json(node, hash);
	}
	
	char* id = get_json_node_value_string(node, "id");
	json_node* tree_node = get_json_node(node, "tree");
	
	Parameters* parameters = NULL;
	Parameters* x = new_Parameters(1);
	DistributionModel* dm = NULL;
	Model* model = NULL;
	
	char* file = get_json_node_value_string(node, "file");
	Vector** samples = NULL;
	if (file != NULL && strcasecmp(d_string, "dirichlet") != 0){
		get_parameters_references2(node, hash, x, "x");
		size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
		samples = read_log_for_parameters_t(file, burnin, x);
	}
	
	if (strcasecmp(d_string, "dirichlet") == 0) {
		parameters = new_Parameters(1);
		Model* msimplex = get_simplex(node, hash);
		Simplex* simplex = msimplex->obj;

		if (file != NULL) {
			char* name = msimplex->name;
			char** names = malloc(sizeof(char*)*simplex->K);
			StringBuffer* buffer = new_StringBuffer(10);
			for(int i = 0; i < simplex->K; i++){
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "%s%d", name, i+1);
				names[i] = StringBuffer_tochar(buffer);
			}
			size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
			samples = read_log_for_names_t(file, burnin, names, simplex->K);
			
			size_t paramCount = simplex->K;
			double* means = dvector(paramCount);
			double* variances = dvector(paramCount);
			for (int i = 0; i < paramCount; i++) {
				means[i] = mean(Vector_data(samples[i]), Vector_length(samples[i]));
				variances[i] = variance(Vector_data(samples[i]), Vector_length(samples[i]), means[i]);
			}
			double num = 0;
			double denom = 0;
			for (int i = 0; i < paramCount; i++) {
				num += means[i]*pow(1.0 - means[i], 2);
				denom += means[i]*variances[i]*(1.0 - means[i]);
			}
			double mhat = num/denom - 1.0;
			for (int i = 0; i < paramCount; i++) {
				Parameters_move(parameters, new_Parameter("alpha", mhat*means[i], NULL));
//				printf("%s %f %f %f\n", names[i], mhat*means[i], mhat, means[i]);
				free(names[i]);
			}
//			printf("\n");
			dm = new_DirichletDistributionModel_with_parameters(parameters, simplex);
			
			free(means);
			free(variances);
			free(names);
			free_StringBuffer(buffer);
		}
		else{
			get_parameters_references(node, hash, parameters);
			int i = 0;
			for ( ; i < Parameters_count(parameters); i++) {
				if(Parameters_value(parameters, i) != 1) break;
			}
			if(i == Parameters_count(parameters)){
				dm = new_FlatDirichletDistributionModel(simplex);
			}
			else{
				dm = new_DirichletDistributionModel_with_parameters(parameters, simplex);
			}
		}
		
		model = new_DistributionModel3(id, dm, msimplex);
		model->sample = _dist_model_sample;
		model->sample_evaluate = _dist_model_sample_evaluate;
		model->samplable = true;
    }
	else if (strcasecmp(d_string, "beta") == 0) {
		parameters = new_Parameters(1);
		// empirical
		if (samples != NULL) {
			size_t paramCount = Parameters_count(x);
			
			for (int i = 0; i < paramCount; i++) {
				double* vec = Vector_data(samples[i]);
				double m = mean(vec, Vector_length(samples[i]));
				double v = variance(vec, Vector_length(samples[i]), m);
				double alpha = m*(m*(1.0 - m)/v - 1.0);
				double beta = (1.0 - m)*(m*(1.0 - m)/v - 1.0);
				Parameters_move(parameters, new_Parameter("alpha", alpha, NULL));
				Parameters_move(parameters, new_Parameter("beta", beta, NULL));
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
			json_node* x_node = get_json_node(node, "parameters");
			get_parameters_references2(node, hash, x, "x");
			get_parameters_references(node, hash, parameters);
			if (strcasecmp(x_node->children[0]->key, "alpha") != 0) {
				Parameters_swap_index(parameters, 0, 1);
			}
			for (int i = 0; i < Parameters_count(parameters); i++) {
				Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
			}
		}
		
		dm = new_IndependantBetaDistributionModel_with_parameters(parameters, x);
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->sample_evaluate = _dist_model_sample_evaluate;
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
			json_node* x_node = get_json_node(node, "parameters");
			get_parameters_references2(node, hash, x, "x");
			get_parameters_references(node, hash, parameters);
			if (strcasecmp(x_node->children[0]->key, "mu") != 0) {
				Parameters_swap_index(parameters, 0, 1);
			}
			for (int i = 0; i < Parameters_count(parameters); i++) {
				Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
			}
		}
		
		dm = new_IndependantLognormalDistributionModel_with_parameters(parameters, x);
		dm->shift = get_json_node_value_double(node, "shift", 0);
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->sample_evaluate = _dist_model_sample_evaluate;
		model->samplable = true;
	}
	else if (strcasecmp(d_string, "normal") == 0) {
		parameters = new_Parameters(2);
		bool tau = false;
		// empirical
		if (samples != NULL) {
			size_t paramCount = Parameters_count(x);
			
			for (int i = 0; i < paramCount; i++) {
				double* vec = Vector_data(samples[i]);
				double m = mean(vec, Vector_length(samples[i]));
				double v = variance(vec, Vector_length(samples[i]), m);
				Parameters_move(parameters, new_Parameter("mu", m, NULL));
				Parameters_move(parameters, new_Parameter("sigma", sqrt(v), NULL));
			}
		}
		else if(get_json_node(node, "parameters") == NULL){
			get_parameters_references2(node, hash, x, "x");
			for (int i = 0; i < Parameters_count(x); i++) {
				Parameters_move(parameters, new_Parameter("mu", 1, new_Constraint(-INFINITY, INFINITY)));
				Parameters_move(parameters, new_Parameter("sigma", 1, new_Constraint(0, INFINITY)));
			}
		}
		else{
			get_parameters_references2(node, hash, x, "x");
			json_node* x_node = get_json_node(node, "parameters");
			for (int i = 0; i < x_node->child_count; i++) {
				if (strcasecmp(x_node->children[i]->key, "tau") == 0) {
					tau = true;
				}
				else if (strcasecmp(x_node->children[i]->key, "mean") != 0 && strcasecmp(x_node->children[i]->key, "sigma") != 0) {
					fprintf(stderr, "Normal distribution should be parametrized with mean and (sigma ot tau)\n");
					exit(13);
				}
			}
			get_parameters_references(node, hash, parameters);
			if (strcasecmp(x_node->children[0]->key, "mean") != 0) {
				Parameters_swap_index(parameters, 0, 1);
			}
			for (int i = 0; i < Parameters_count(parameters); i++) {
				Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
			}
		}
		dm = new_IndependantNormalDistributionModel_with_parameters(parameters, x);
		
		dm->parameterization = tau ? DISTRIBUTION_NORMAL_MEAN_TAU: DISTRIBUTION_NORMAL_MEAN_SIGMA;
		dm->shift = get_json_node_value_double(node, "shift", -INFINITY);
		
		model = new_DistributionModel2(id, dm);
		model->sample = _dist_model_sample;
		model->sample_evaluate = _dist_model_sample_evaluate;
		model->samplable = true;
	}
    else if (strcasecmp(d_string, "topology") == 0) {
        char* ref = get_json_node_value_string(node, "tree");
        Model* mtree = Hashtable_get(hash, ref+1);
        dm = new_UniformTreeDistribution(mtree->obj);
        model = new_DistributionModel3(id, dm, mtree);
    }
	else if(strcasecmp(d_string, "multivariatenormal") == 0){
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
		model->sample_evaluate = _dist_model_sample_evaluate;
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
			
			// calculate margins parameters
			if(from_node == NULL){
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
			}
			else{
				char* ref = (char*)from_node->value;
				Model* simpleModel = Hashtable_get(hash, ref+1);
				DistributionModel* simpleDM = simpleModel->obj;
				for (int i = 0; i < paramCount*2; i++) {
					Parameters_move(parameters, clone_Parameter(Parameters_at(simpleDM->parameters, i)));
				}
			}
		
			double* cov = dvector(paramCount*paramCount);
			// Calculate sample covariance matrix
			for (int i = 0; i < paramCount; i++) {
				double* pp = Vector_data(samples[i]);
				cov[i*paramCount+i] = covariance(pp, pp, means[i], means[i], n);
				for (int j = i+1; j < paramCount; j++) {
					double* pp2 = Vector_data(samples[j]);
					cov[i*paramCount+j] = cov[j*paramCount+i] = covariance(pp, pp2, means[i], means[j], n);
				}
			}
			
			// sample correlation matrix
			for (int i = 0; i < paramCount; i++) {
				double covii = sqrt(cov[i*paramCount + i]);
				for (int j = 0; j < paramCount; j++) {
					double covjj = sqrt(cov[j*paramCount + j]);
					double covij = cov[i*paramCount + j];
					Parameters_move(parameters, new_Parameter("ij", covij/(covii*covjj), NULL));
				}
			}
			free(cov);
			
			dm = new_CopulaDistributionModel_with_parameters(parameters, x);
			
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
			parameters = new_Parameters(1);
			
			size_t paramCount = Parameters_count(x);
			
			json_node* from_node = get_json_node(node, "from");
			char* ref = (char*)from_node->value;
			Model* simpleModel = Hashtable_get(hash, ref+1);
			DistributionModel* simpleDM = simpleModel->obj;
			for (int i = 0; i < paramCount*2; i++) {
				Parameters_move(parameters, clone_Parameter(Parameters_at(simpleDM->parameters, i)));
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
					Parameters_move(parameters, new_Parameter("ij", covij/(covii*covjj), NULL));
				}
			}
			gsl_permutation_free(perm);
			gsl_matrix_free(H);
			gsl_matrix_free(L);
			
			dm = new_CopulaDistributionModel_with_parameters(parameters, x);
			
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
	free_Parameters(parameters);
	free_Parameters(x);
	
	return model;
}
