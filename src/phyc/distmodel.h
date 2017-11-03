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


struct _DistributionModel;
typedef struct _DistributionModel DistributionModel;

struct _DistributionModel{
	Parameters* parameters;
	Parameters* x;
	double* tempx; // array to pass to multivariate distributions
	double* tempp;
	double (*logP)(DistributionModel*);
	double (*dlogP)(DistributionModel*, const Parameter*);
	void (*free)(DistributionModel*);
	DistributionModel* (*clone)(DistributionModel*);
	void* data;
};


DistributionModel* new_IndependantGammaDistributionModel(const double shape, const double rate, const Parameters* x);

DistributionModel* new_IndependantExpDistributionModel(const double lambda, const Parameters* x);

DistributionModel* new_FlatDirichletDistributionModel(const Parameters* x);

Model* new_DistributionModel2(const char* name, DistributionModel* cm);

#endif /* distmodel_h */
