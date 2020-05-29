//
//  onlineopt.h
//  physher
//
//  Created by mathieu on 28/5/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef onlineopt_h
#define onlineopt_h

#include "parameters.h"

typedef struct OnlineOptimizer{
	Model* model;
	Model* treelikelihood;
	size_t* indexes; // index in sitepattern
	size_t indexes_count;
	double new_taxon_bl;
	
	double (*optimize)( struct OnlineOptimizer * );
	void (*free)( struct OnlineOptimizer * );
}OnlineOptimizer;

OnlineOptimizer* new_OnlineOptimizer_from_json(json_node* node, Hashtable* hash);

#endif /* onlineopt_h */
