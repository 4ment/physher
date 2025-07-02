//
//  distmodel.h
//  physher
//
//  Created by Mathieu Fourment on 3/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef distmodel_h
#define distmodel_h

#include <stdio.h>

#include "parameters.h"
#include "tree.h"
#include "mjson.h"

#ifndef GSL_DISABLED
#include <gsl/gsl_rng.h>
#endif

struct _DistributionModel;
typedef struct _DistributionModel DistributionModel;

typedef enum distribution_parameterization{
	DISTRIBUTION_EXPONENTIAL_MEAN,
	DISTRIBUTION_EXPONENTIAL_RATE,
	DISTRIBUTION_GAMMA_SHAPE_RATE,
	DISTRIBUTION_GAMMA_SHAPE_SCALE,
    DISTRIBUTION_HALFNORMAL_MEAN_SIGMA,
    DISTRIBUTION_HALFNORMAL_MEAN_TAU,
    DISTRIBUTION_NORMAL_MEAN_SIGMA,
    DISTRIBUTION_NORMAL_MEAN_TAU
}distribution_parameterization;

typedef enum distribution_t{
	DISTRIBUTION_BETA = 0,
	DISTRIBUTION_BETA_PRIME,
	DISTRIBUTION_CAUCHY,
	DISTRIBUTION_CTMC_SCALE,
	DISTRIBUTION_DIRICHLET,
	DISTRIBUTION_DISCRETE,
	DISTRIBUTION_EXPONENTIAL,
	DISTRIBUTION_GAMMA,
	DISTRIBUTION_GMRF,
    DISTRIBUTION_HALFNORMAL,
	DISTRIBUTION_KUMARASWAMY,
	DISTRIBUTION_LOGNORMAL,
	DISTRIBUTION_NORMAL,
	DISTRIBUTION_NORMAL_MULTIVARIATE,
	DISTRIBUTION_ONE_ON_X,
	DISTRIBUTION_UNIFORM,
	DISTRIBUTION_WEIBULL
}distribution_t;

char* DISTRIBUTION_NAME[] = {
	"beta",
	"beta prime",
	"cauchy",
	"ctmc scale",
	"dirichlet",
	"discrete",
	"exponential",
	"gamma",
	"gmrf",
    "half normal",
	"Kumaraswamy",
	"lognormal",
	"normal",
	"multivariate normal",
	"one on x",
	"uniform",
	"Weibull"
};

struct _DistributionModel{
	distribution_t type;
	Parameters* parameters;
	Parameters* x;
	Tree* tree;
	double* tempx; // array to pass to multivariate distributions and sampling in general
	double* tempp;
	double (*log_prob)(DistributionModel*);
	double (*log_prob_grad)(DistributionModel*, const Parameters*);
	void (*log_prob_hessian_diag)(DistributionModel*, const Parameters*);
	void (*reparam_backprop)(DistributionModel*);
	double (*dlogP)(DistributionModel*, const Parameter*);
	double (*d2logP)(DistributionModel*, const Parameter*);
	double (*ddlogP)(DistributionModel*, const Parameter*, const Parameter*);
	void (*sample)(DistributionModel*);
	void (*rsample)(DistributionModel*);
	double (*entropy)(DistributionModel*);
	void (*entropy_grad)(DistributionModel*, const Parameters*);
	void (*free)(DistributionModel*);
	DistributionModel* (*clone)(DistributionModel*);
	void* data;
	double lp;
	double stored_lp;
	bool need_update;
	distribution_parameterization parameterization;
#ifndef GSL_DISABLED
	gsl_rng* rng;
#endif
	double shift;
	
	int prepared_gradient;
	double* gradient;
	size_t gradient_length;
	bool need_update_gradient;
	double support[2];
};


DistributionModel* new_UniformTreeDistribution(Tree* tree);

DistributionModel* new_DistributionModel(Parameters* p, Parameters* x);

DistributionModel* clone_DistributionModel_with_parameters(DistributionModel* dm, Parameters* params, Parameters* x);

Model* new_DistributionModel2(const char* name, DistributionModel* dm);

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* amodel);

double DistributionModel_dlog_0(DistributionModel* dm, const Parameter* p);

double DistributionModel_d2log_0(DistributionModel* dm, const Parameter* p);

double DistributionModel_ddlog_0(DistributionModel* dm, const Parameter* p1, const Parameter* p2);

void distmodel_get_parameters(json_node* parameters_node, Hashtable* hash, Parameters* parameters);

Parameters* distmodel_get_x(const char* who, json_node* node, Hashtable* hash);

Parameter* distmodel_parse_parameter(json_node* parameter_node, Hashtable* hash, const char* id, double lower, double upper);

void distmodel_expand_2parameters(Parameters* x, Parameters* parameters, double** aValues, double** bValues);

void distmodel_expand_parameter(Parameters* x, Parameters* parameters, double** aValues);

#endif /* distmodel_h */
