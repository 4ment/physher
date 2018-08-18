//
//  is.h
//  physher
//
//  Created by Mathieu Fourment on 19/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef is_h
#define is_h

#include <stdio.h>
#include "parameters.h"

typedef struct ImportanceSampler{
	Model* model; // distribution of interest
	Model** distribution; // importance distribution or mixture of distribution
	double* weights;
	size_t distribution_count;
	size_t samples;
	Parameters* parameters;
	bool normalize;
	double(*calculate)(struct ImportanceSampler*);
	void(*free)(struct ImportanceSampler*);
}ImportanceSampler;

ImportanceSampler* new_ImportanceSampler_from_json(json_node* node, Hashtable* hash);

#endif /* is_h */
