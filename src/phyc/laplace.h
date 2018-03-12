//
//  laplace.h
//  physher
//
//  Created by Mathieu Fourment on 17/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef laplace_h
#define laplace_h

#include <stdio.h>

#include "parameters.h"

typedef struct Laplace{
	Model* model;
	Parameters* parameters;
	Model* refdist;
	double(*calculate)(struct Laplace*);
	void(*free)(struct Laplace*);
	// something print
}Laplace;

Laplace* new_Laplace_from_json(json_node* node, Hashtable* hash);

#endif /* laplace_h */
