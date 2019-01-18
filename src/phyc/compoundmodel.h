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
#include "simplex.h"
#include "mjson.h"

struct _CompoundModel;
typedef struct _CompoundModel CompoundModel;

struct _CompoundModel{
	Model** models;
	int count;
	Model* weights;
	double (*logP)(CompoundModel*);
	double (*dlogP)(CompoundModel*, const Parameter*);
	double (*d2logP)(CompoundModel*, const Parameter*);
	double (*ddlogP)(CompoundModel*, const Parameter*, const Parameter*);
	void (*free)(CompoundModel*);
	void(*add)(CompoundModel*, Model*);
	void(*move)(CompoundModel*, Model*);
	void(*remove)(CompoundModel*, Model*);
	void(*removeAll)(CompoundModel*);
};

CompoundModel* new_CompoundModel();

Model* new_CompoundModel2(const char* name, CompoundModel* cm);

Model* new_CompoundModel_from_json(json_node*node, Hashtable*hash);

#endif /* compoundmodel_h */
