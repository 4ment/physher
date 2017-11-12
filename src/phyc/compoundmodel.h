//
//  compoundmodel.h
//  physher
//
//  Created by Mathieu Fourment on 1/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef compoundmodel_h
#define compoundmodel_h

#include "model.h"
#include "parameters.h"
#include "mjson.h"

struct _CompoundModel;
typedef struct _CompoundModel CompoundModel;

struct _CompoundModel{
	Model** models;
	int count;
	double (*logP)(CompoundModel*);
	double (*dlogP)(CompoundModel*, const Parameter*);
	void (*free)(CompoundModel*);
	void(*add)(CompoundModel*, Model*);
	void(*move)(CompoundModel*, Model*);
	void(*remove)(CompoundModel*, Model*);
	void(*removeAll)(CompoundModel*);
};

CompoundModel* new_CompoundModel();

Model* new_CompoundModel2(const char* name, CompoundModel* cm);

#endif /* compoundmodel_h */
