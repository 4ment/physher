// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
	DISTRIBUTION_LOGNORMAL_MEAN_STDEV,
	DISTRIBUTION_LOGNORMAL_MU_SIGMA,
    DISTRIBUTION_NORMAL_MEAN_SIGMA,
    DISTRIBUTION_NORMAL_MEAN_TAU,
    DISTRIBUTION_WEIBULL_SCALE_SHAPE
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
	DISTRIBUTION_GAR,
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

// char* DISTRIBUTION_NAME[] = {
// 	"beta",
// 	"beta prime",
// 	"cauchy",
// 	"ctmc scale",
// 	"dirichlet",
// 	"discrete",
// 	"exponential",
// 	"gamma",
// 	"gar",
// 	"gmrf",
//     "half normal",
// 	"KUMARASWAMY",
// 	"lognormal",
// 	"normal",
// 	"multivariate normal",
// 	"one on x",
// 	"uniform",
// 	"weibull"
// };

struct _DistributionModel{
	distribution_t type;
	Parameters* parameters;
	Parameters* x;
	// Simplex* simplex;
	Tree* tree;
	double* tempx; // array to pass to multivariate distributions and sampling in general
	double* tempp;
	double (*logP)(DistributionModel*);
	void (*gradient)(DistributionModel*, Parameters*);
	void (*rgradient)(DistributionModel*);
	void (*sample)(DistributionModel*);
	void (*rsample)(DistributionModel*);
	double (*entropy)(DistributionModel*);
	void (*gradient_entropy)(DistributionModel*, const Parameters*);
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
	double support[2];
};


DistributionModel* new_UniformTreeDistribution(Tree* tree);

DistributionModel* new_DistributionModel(Parameters* p, Parameters* x);

DistributionModel* clone_DistributionModel_with_parameters(DistributionModel* dm, Parameters* params, Parameters* x);

Model* new_DistributionModel2(const char* name, DistributionModel* dm);

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* amodel);

void distmodel_get_parameters(json_node* parameters_node, Hashtable* hash, Parameters* parameters);

Parameters* distmodel_get_x(const char* who, json_node* node, Hashtable* hash);

Parameter* distmodel_parse_parameter(json_node* parameter_node, Hashtable* hash, const char* id, double lower, double upper);

#endif /* distmodel_h */
