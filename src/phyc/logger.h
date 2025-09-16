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

typedef struct Logger{
	Parameters* parameters;
	Model** models;
	size_t model_count;
	void (*log)(struct Logger*);
	FILE* file;
	char* filename;
	char* format;
	bool tree;
	bool internal; // show internal node name
    char sep; //separator
}Logger;

struct Dumper{
    Parameters** parameters;
    size_t parameter_count;
    Model** models;
    size_t model_count;
    void (*dump)(struct Dumper*);
    void (*free)(struct Dumper*);
    FILE* file;
    char* filename;
};


struct Logger* new_logger_from_json(json_node* node, Hashtable* hash);

void free_Logger(struct Logger* logger);


struct Dumper* new_Dumper_from_json(json_node* node, Hashtable* hash);

void free_Dumper(struct Dumper* logger);

#endif /* logger_h */
