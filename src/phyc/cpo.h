//
//  cpo.h
//  physher
//
//  Created by Mathieu Fourment on 22/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef cpo_h
#define cpo_h

#include "mjson.h"
#include "hashtable.h"

typedef struct CPO{
	char* filename;
	size_t burnin;
	void(*calculate)(struct CPO*);
	void(*free)(struct CPO*);
	// something print
}CPO;

CPO* new_CPO_from_json(json_node* node, Hashtable* hash);

//void print_hessian(Hessian* hessian);

#endif /* cpo_h */
