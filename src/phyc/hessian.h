//
//  hessian.h
//  physher
//
//  Created by Mathieu Fourment on 5/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef hessian_h
#define hessian_h

#include <stdio.h>

#include "parameters.h"
#include "mjson.h"
#include "hashtable.h"

typedef struct Hessian{
	Model* likelihood;
	Parameters* parameters;
	double* hessian;
	void(*calculate)(struct Hessian*);
	// something print
}Hessian;

Hessian* new_Hessian_from_json(json_node* node, Hashtable* hash);

void print_hessian(Hessian* hessian);

#endif /* hessian_h */
