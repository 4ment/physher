/*
 *  args.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 6/7/11.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef _ARGS_H_
#define _ARGS_H_

#include "utils.h"
#include "hashtable.h"

struct argsparser_option;
struct argsparser2_option;

enum argsparser_option_type {
    ARGS_OPTION_FLAG,
    ARGS_OPTION_BOOLEAN,
    ARGS_OPTION_INTEGER,
    ARGS_OPTION_LONG,
    ARGS_OPTION_FLOAT,
    ARGS_OPTION_DOUBLE,
    ARGS_OPTION_STRING,
};

struct argsparser_option {
    enum argsparser_option_type type;
    const char short_name;
    const char* long_name;
    const char* config_name;
    void* value;
    const char *help;
};

struct argsparser2_option {
	enum argsparser_option_type type;
	const char short_name;
	const char* long_name;
	char* value;
	const char *help;
};

typedef struct args_parser {
    struct argsparser_option *options;
    int option_count;
}args_parser;

typedef struct args_parser2 {
	struct argsparser2_option *options;
	int option_count;
}args_parser2;

args_parser2* argsparser2_init(struct argsparser2_option options[], int count);

Hashtable * argsparser2_parse(args_parser2* args, const char* argv[], int argc);

void argsparser2_free(args_parser2* args);

args_parser* argsparser_init(struct argsparser_option options[], int count);

void argsparser_parse(args_parser* args, char* argv[], int argc);

void argsparser_parse_file(args_parser* args, const char* filename);

void argsparser_help(args_parser* args, int level);

void argsparser2_help(args_parser2* args);


int args_get_index( int argc, char* argv[], const char flag[] );

bool args_contains( int argc, char* argv[], const char *flag );


char * args_get_string( int argc, char* argv[], const char flag[] );

int * args_get_pint( int argc, char* argv[], const char flag[] );

int args_get_int( int argc, char* argv[], const char flag[], bool *success );

double * args_get_pdouble( int argc, char* argv[], const char flag[] );

double args_get_double( int argc, char* argv[], const char flag[], bool *success );

// should be between double quotes if the separator is a space
int * args_get_int_array( int argc, char* argv[], const char flag[], const char sep, int *count );

// should be between double quotes if the separator is a space
double * args_get_double_array( int argc, char* argv[], const char flag[], const char sep, int *count );

bool args_get_boolean( int argc, char* argv[], const char flag[] );

const char* get_program_name(const char* argv[]);


#endif
