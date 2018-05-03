//
//  mc.h
//  physher
//
//  Created by Mathieu Fourment on 3/5/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef mc_h
#define mc_h

#include <stdio.h>

#include "mjson.h"
#include "parameters.h"

typedef struct MC{
	Model* likelihood;
	Model* prior;
	Parameters* parameters;
	size_t samples;
	double(*calculate)(struct MC*);
	void(*free)(struct MC*);
}MC;

MC* new_MonteCarlo_from_json(json_node* node, Hashtable* hash);

#endif /* mc_h */
