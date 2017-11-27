//
//  logger.h
//  physher
//
//  Created by Mathieu Fourment on 26/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef logger_h
#define logger_h

#include <stdio.h>

#include "parameters.h"
#include "simplex.h"

struct Logger{
	Parameters* parameters;
	Model** simplexes;
	size_t simplexCount;
	void (*log)(struct Logger*);
	FILE* file;
	char* filename;
	char* format;
};


struct Logger* new_logger_from_json(json_node* node, Hashtable* hash);

void free_Logger(struct Logger* logger);

#endif /* logger_h */
