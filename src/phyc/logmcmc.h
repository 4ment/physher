//
//  logmcmc.h
//  physher
//
//  Created by Mathieu Fourment on 30/01/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef logmcmc_h
#define logmcmc_h

#include <stdio.h>

#include "parameters.h"

typedef struct Log{
	Parameters* x;
	Model** models;
	size_t model_count;
	Model** simplexes;
	size_t simplex_count;
	FILE* file;
	char* filename;
	size_t every;
	bool append;
	void(*initialize)(struct Log* logger);
	void(*finalize)(struct Log* logger);
	void(*write)(struct Log* logger, size_t);
	void(*write_header)(struct Log* logger, bool);
	void(*write_with)(struct Log* logger, size_t, const char*);
	void(*free)(struct Log*);
	bool cpo;
}Log;
	
 Log* new_Log_from_json(json_node* node, Hashtable* hash);

#endif /* logmcmc_h */
