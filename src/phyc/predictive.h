//
//  predictive.h
//  physher
//
//  Created by Mathieu Fourment on 16/02/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef predictive_h
#define predictive_h

#include <stdio.h>

#include "mjson.h"
#include "hashtable.h"

#define JSON_PREDICTIVE "predictive"

struct _Predictive;
typedef struct _Predictive Predictive;

struct _Predictive{
	char* filename;
	size_t burnin;
	void(*calculate)(Predictive*);
	void(*free)(Predictive*);
};

Predictive* new_Predictive_from_json(json_node* node, Hashtable* hash);

#endif /* predictive_h */
