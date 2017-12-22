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

#include "model.h"
#include "parameters.h"
#include "simplex.h"
#include "mjson.h"


struct _DistributionModel;
typedef struct _DistributionModel DistributionModel;

struct _DistributionModel{
	Parameters* parameters;
	Parameters* x;
	Simplex* simplex;
	double* tempx; // array to pass to multivariate distributions
	double* tempp;
	double (*logP)(DistributionModel*);
	double (*dlogP)(DistributionModel*, const Parameter*);
	void (*free)(DistributionModel*);
	DistributionModel* (*clone)(DistributionModel*);
	void* data;
	double lp;
	bool need_update;
};


DistributionModel* new_IndependantGammaDistributionModel(const double shape, const double rate, const Parameters* x);

DistributionModel* new_IndependantExpDistributionModel(const double lambda, const Parameters* x);

DistributionModel* new_FlatDirichletDistributionModel(Simplex* simplex);

DistributionModel* new_UniformTreeDistribution(Model* tree);


Model* new_DistributionModel2(const char* name, DistributionModel* dm);

Model* new_DistributionModel3(const char* name, DistributionModel* dm, Model* simplex);

Model* new_TreeDistributionModel(const char* name, DistributionModel* dm, Model* tree);

Model* new_DistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distmodel_h */
