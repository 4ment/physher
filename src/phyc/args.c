/*
 *  args.c
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

#include "args.h"

#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <assert.h>

#include "utils.h"
#include "matrix.h"
#include "mstring.h"
#include "filereader.h"
#include "hashtable.h"

args_parser* argsparser_init(struct argsparser_option options[], int count){
    struct args_parser* args = malloc(sizeof(args_parser));
    assert(args);
    args->option_count = count;
    args->options = options;
    assert(args->options);
    
    return args;
}

args_parser2* argsparser2_init(struct argsparser2_option options[], int count){
	struct args_parser2* args = malloc(sizeof(args_parser2));
	assert(args);
	args->option_count = count;
	args->options = options;
	assert(args->options);
	
	return args;
}

void argsparser_help(args_parser* args, int level){
    int max = 0;
    for(int j = 0; j < args->option_count; j++){
        if(args->options[j].long_name != NULL){
            max = dmax(max, strlen(args->options[j].long_name));
        }
    }
    
    int max2 = 0;
    if(level > 0){
        for(int j = 0; j < args->option_count; j++){
            if(args->options[j].config_name != NULL){
                max2 = dmax(max2, strlen(args->options[j].config_name));
            }
        }
    }
    
    for(int j = 0; j < args->option_count; j++){
        if(args->options[j].short_name != 0){
            fprintf(stdout, "-%c ", args->options[j].short_name);
        }
        else{
            fprintf(stdout, "   ");
        }
        
        if(args->options[j].long_name != NULL){
            fprintf(stdout, "--%-*s ", max+2, args->options[j].long_name);
        }
        else{
            fprintf(stdout, " ");
        }
        
        if(level > 0){
            fprintf(stdout, "%-*s ", max2+1, args->options[j].config_name);
        }
        
        fprintf(stdout, "%s ", args->options[j].help);
        
        if(args->options[j].type == ARGS_OPTION_STRING && *(char**)args->options[j].value != NULL){
            fprintf(stdout, " [default: %s] ", *(char**)args->options[j].value);
        }
        if(args->options[j].type == ARGS_OPTION_INTEGER){
            fprintf(stdout, " [default: %d] ", *(int*)args->options[j].value);
        }
        if(args->options[j].type == ARGS_OPTION_LONG){
            fprintf(stdout, " [default: %ld] ", *(long*)args->options[j].value);
        }
        if(args->options[j].type == ARGS_OPTION_BOOLEAN){
            fprintf(stdout, " [default: %s] ", (*(bool*)args->options[j].value ? "true": "false"));
        }
        if(args->options[j].type == ARGS_OPTION_DOUBLE){
            fprintf(stdout, " [default: %g] ", *(double*)args->options[j].value);
        }
        if(args->options[j].type == ARGS_OPTION_FLOAT){
            fprintf(stdout, " [default: %g] ", *(float*)args->options[j].value);
        }
        fprintf(stdout, "\n");
    }
}

void argsparser2_help(args_parser2* args){
	int max = 0;
	for(int j = 0; j < args->option_count; j++){
		if(args->options[j].long_name != NULL){
			max = dmax(max, strlen(args->options[j].long_name));
		}
	}
	
	for(int j = 0; j < args->option_count; j++){
		if(args->options[j].short_name != 0){
			fprintf(stdout, "-%c ", args->options[j].short_name);
		}
		else{
			fprintf(stdout, "   ");
		}
		
		if(args->options[j].long_name != NULL){
			fprintf(stdout, "--%-*s ", max+2, args->options[j].long_name);
		}
		else{
			fprintf(stdout, " ");
		}
		
		fprintf(stdout, "%s ", args->options[j].help);
		
		if(args->options[j].value != NULL && args->options[j].type != ARGS_OPTION_FLAG){
			fprintf(stdout, " [default: %s] ", args->options[j].value);
		}
		fprintf(stdout, "\n");
	}
}

void argsparser_check(args_parser* args, char* argv[], int argc){
    Hashtable *hash = new_Hashtable_string(10);
    hashtable_set_key_ownership( hash, false );
    hashtable_set_value_ownership( hash, false );
    char* a = "a";
    for(int i = 1; i < argc; i++){
        if (Hashtable_exists(hash, argv[i])){
            fprintf(stderr, "Duplicate option %s\n", argv[i]);
            exit(1);
        }
        Hashtable_add(hash, argv[i], a);
    }
    free_Hashtable(hash);
}

void argsparser2_check(args_parser2* args, char* argv[], int argc){
	Hashtable *hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	char* a = "a";
	for(int i = 1; i < argc; i++){
		if (Hashtable_exists(hash, argv[i])){
			fprintf(stderr, "Duplicate option %s\n", argv[i]);
			exit(1);
		}
		Hashtable_add(hash, argv[i], a);
	}
	free_Hashtable(hash);
}

void argsparser_parse(args_parser* args, char* argv[], int argc){
    
    argsparser_check(args, argv, argc);
    
    int i = 0;
    for ( ; i < argc; i++) {
        int j = 0;
        bool option_flag = false;
        for(; j < args->option_count; j++){
            if(argv[i][0] == '-') option_flag = true;
            if ( strlen(argv[i]) > 2 && argv[i][0] == '-' && argv[i][1] == '-' && strcmp(argv[i]+2, args->options[j].long_name) == 0 ) {
                break;
            }
            else if (strlen(argv[i]) == 2 && argv[i][0] == '-' && args->options[j].short_name == argv[i][1]){
                break;
            }
        }
        
        if(j != args->option_count && option_flag == true){
            switch(args->options[j].type){
                case ARGS_OPTION_FLAG:{
                    *((bool*)args->options[j].value) = true;
                    break;
                }
                case ARGS_OPTION_BOOLEAN:{
                    bool v;
                    if(strcasecmp(argv[i+1], "true") == 0 || strcmp(argv[i+1], "1") == 0){
                        v = true;
                    }
                    else if(strcasecmp(argv[i+1], "false") == 0 || strcmp(argv[i+1], "0") == 0){
                        v = false;
                    }
                    else{
                        fprintf(stderr, "Option %s %s is not valid. Expect a boolean argument", argv[i], argv[i+1]);
                        exit(1);
                    }
                    *((bool*)args->options[j].value) = v;
                    break;
                }
                case ARGS_OPTION_INTEGER:{
                    *((int*)args->options[j].value) = atoi(argv[i+1]);
                    break;
                }
                case ARGS_OPTION_LONG:{
                    *((long*)args->options[j].value) = atol(argv[i+1]);
                    break;
                }
                case ARGS_OPTION_STRING:{
                    *((char **)args->options[j].value) = String_clone(argv[i+1]);
                    break;
                }
                case ARGS_OPTION_DOUBLE:{
                    *((double*)args->options[j].value) = atof(argv[i+1]);
                    break;
                }
                case ARGS_OPTION_FLOAT:{
                    *((float*)args->options[j].value) = atof(argv[i+1]);
                    break;
                }
            }
            if (args->options[j].type != ARGS_OPTION_FLAG) {
                i++;
            }
        }
        else if(option_flag == true){
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            exit(1);
        }
    }
}

void argsparser_parse_file(args_parser* args, const char* filename){

    FileReader *reader = new_FileReader(filename, 1000);
    while ( reader->read_line(reader) ) {
        StringBuffer_trim(reader->buffer);
        if ( reader->buffer->length == 0 || String_start_with(reader->line, "#", true)) continue;
        char *p = reader->buffer->c;
        while ( *p != '\0') {
            if ( *p == '#') {
                *p = '\0';
                break;
            }
            p++;
        }
        reader->buffer->length = strlen(reader->buffer->c);
        StringBuffer_trim(reader->buffer);
        
        
        int l = 0;
        char **temp = String_split_char( reader->buffer->c, '=', &l );
        int index1 = 0;
        int index2 = 0;
        char* line = reader->buffer->c;
        while (*line == ' '){ line++;index1++;index2++;}
        while (*line != '='){ line++;index2++;}
        line++;index2++;
        while (*line == ' '){ line++;index2++;}

        if ( l != 2 ) {
            fprintf(stderr, "Could not parse: %s . There should be only 1 '='(%d)\n", reader->line, l);
            continue;
        }
        char *key   = String_rtrim( temp[0] );
        char *value = String_trim( temp[1] );
        free(temp);
        
        int j = 0;
        for(; j < args->option_count; j++){
            if ( strcasecmp(key, args->options[j].config_name) == 0 ) {
                break;
            }
        }
        
        if(j != args->option_count){
            switch(args->options[j].type){
                case ARGS_OPTION_FLAG:{
                    bool v;
                    if(strcasecmp(value, "true") == 0 || strcmp(value, "1") == 0){
                        v = true;
                    }
                    else if(strcasecmp(value, "false") == 0 || strcmp(value, "0") == 0){
                        v = false;
                    }
                    else{
                        fprintf(stderr, "Option %s %s is not valid. Expect a boolean argument", key, value);
                        exit(1);
                    }
                    *((bool*)args->options[j].value) = v;
                    break;
                }
                case ARGS_OPTION_BOOLEAN:{
                    bool v;
                    if(strcasecmp(value, "true") == 0 || strcmp(value, "1") == 0){
                        v = true;
                    }
                    else if(strcasecmp(value, "false") == 0 || strcmp(value, "0") == 0){
                        v = false;
                    }
                    else{
                        fprintf(stderr, "Option %s %s is not valid. Expect a boolean argument", key, value);
                        exit(1);
                    }
                    *((bool*)args->options[j].value) = v;
                    break;
                }
                case ARGS_OPTION_INTEGER:{
                    *((int*)args->options[j].value) = atoi(value);
                    break;
                }
                case ARGS_OPTION_LONG:{
                    *((long*)args->options[j].value) = atol(value);
                    break;
                }
                case ARGS_OPTION_STRING:{
                    char* option = *((char **)args->options[j].value);
                    if(option != NULL){
                        free(option);
                    }
                    *((char **)args->options[j].value) = String_clone(value);
                    break;
                }
                case ARGS_OPTION_DOUBLE:{
                    *((double*)args->options[j].value) = atof(value);
                    break;
                }
                case ARGS_OPTION_FLOAT:{
                    *((float*)args->options[j].value) = atof(value);
                    break;
                }
            }
        }
    }
}



int args_get_index( int argc, char* argv[], const char flag[] ){
    int i = 0;
    for (; i < argc; i++) {
        if ( strncmp(argv[i], flag, strlen(flag)) == 0 ) return i;
	}
	
	return -1;
}

bool args_contains( int argc, char* argv[], const char *flag ){
    int i = 0;
	for ( ; i < argc; i++) {
        if ( strncmp(argv[i], flag, strlen(flag)) == 0 ) return true;
	}
	
	return false;
}

char * args_get_string( int argc, char* argv[], const char flag[] ){
	char *option = NULL;
    char *str = NULL;
    int i = 0;
    
	for (; i < argc; i++) {
        if ( strncmp(argv[i], flag, strlen(flag)) == 0 ) {
            str = NULL;
            if( strlen(argv[i]) > strlen(flag) ){
                str = argv[i]+strlen(flag);
            }
            else if( i+1 < argc ){
                str = argv[i+1];
            }
            
            if( str != NULL ){
                option = String_clone(argv[i+1]);
                return option;
            }
            break;
        }
	}
	
	return option;
}


int * args_get_pint( int argc, char* argv[], const char flag[] ){
    int *option = NULL;
    char *str = NULL;
    int i = 0;
    
	for ( ; i < argc; i++) {
        if ( strncmp(argv[i], flag, strlen(flag)) == 0 ) {
            str = NULL;
            if( strlen(argv[i]) > strlen(flag) ){
                str = argv[i]+strlen(flag);
            }
            else if( i+1 < argc ){
                str = argv[i+1];
            }
            
            if( str != NULL && isInt(str) ){
                option = (int*)malloc( sizeof(int) );
                assert(option);
                *option = atoi( str );
                return option;
            }
            break;
        }
	}
	
	return option;
}

int args_get_int( int argc, char* argv[], const char flag[], bool *success ){
    int i = 0;
    char *str = NULL;
    
    for ( ; i < argc; i++) {
        if ( strncmp(argv[i], flag, strlen(flag)) == 0  ) {
            str = NULL;
            if( strlen(argv[i]) > strlen(flag) ){
                str = argv[i]+strlen(flag);
            }
            else if( i+1 < argc ){
                str = argv[i+1];
            }
            
            if( str != NULL && isInt(str) ){
                *success = true;
                return atoi( str );
            }
            break;
        }
        
	}
	*success = false;
	return 0;
}

double * args_get_pdouble( int argc, char* argv[], const char flag[] ){
	double *option = NULL;
    char *str = NULL;
    int i = 0;
    
    for ( ; i < argc; i++) {
        
        if ( strncmp(argv[i], flag, strlen(flag)) == 0 ) {
            str = NULL;
            if( strlen(argv[i]) > strlen(flag) ){
                str = argv[i]+strlen(flag);
            }
            else if( i+1 < argc ){
                str = argv[i+1];
            }
            
            if( str!= NULL && isFloat2(str) ){
                option = (double*)malloc( sizeof(double) );
                assert(option);
                *option = atof( str );
                return option;
            }
            break;
        }
		
	}
	
	return option;
}

double args_get_double( int argc, char* argv[], const char flag[], bool *success ){
    int i = 0;
    char *str = NULL;
    
    for ( ; i < argc; i++) {
			
        if ( strncmp(argv[i], flag, strlen(flag)) == 0 ) {
            str = NULL;
            if( strlen(argv[i]) > strlen(flag) ){
                str = argv[i]+strlen(flag);
            }
            else if( i+1 < argc ){
                str = argv[i+1];
            }
            
            if( str != NULL && isFloat2(str) ){
                *success = true;
                return atof( str );
            }
            break;
        }
		
	}
	*success = false;
	return 0;
}

// should be between double quotes if the separator is a space
int * args_get_double_int( int argc, char* argv[], const char flag[], const char sep, int *count ){
	int *option = NULL;
	char *pstr = NULL;
    char *pstr2 = NULL;
    int i;
	StringBuffer *buffer = new_StringBuffer(10);
	*count = 0;
	
	for ( i = 0; i < argc; i++) {
		pstr = argv[i];
		if ( *pstr == '-' ) {
			while ( *pstr == '-' ) {
				pstr++;
			}
			
			if ( strcmp(pstr, flag) == 0 ) {
				pstr2 = argv[i+1];
				while ( *pstr2 != '\0') {
					if ( *pstr2 == sep ){
                        (*count)++;
                    }
                    pstr2++;
				}
                
				pstr2 = argv[i+1];
				option = ivector( *count+1 );
				*count = 0;
				
				
				while ( *pstr2 != '\0' ) {
					if( *pstr2 == sep && buffer->length != 0 ){
						option[(*count)++] = atoi( buffer->c );
						StringBuffer_empty( buffer );
					}
					else {
						StringBuffer_append_char(buffer, *pstr2);
					}
					pstr2++;
				}
				if ( buffer->length != 0) {
					option[(*count)++] = atoi( buffer->c );
				}
				break;
			}
		}
	}
	
	free_StringBuffer(buffer);
	return option;
}

// should be between double quotes if the separator is a space
double * args_get_double_array( int argc, char* argv[], const char flag[], const char sep, int *count ){
	double *option = NULL;
    char *pstr = NULL;
    char *pstr2 = NULL;
    int i;
	StringBuffer *buffer = new_StringBuffer(10);
	*count = 0;
	
    for ( i = 0; i < argc; i++ ) {
			
			if ( strncmp(argv[i], flag, strlen(flag)) == 0 ) {
                pstr = NULL;
                pstr2 = NULL;
                
                if( strlen(argv[i]) > strlen(flag) ){
                    pstr = pstr2 = argv[i]+strlen(flag);
                }
                else if( i+1 < argc ) {
                    pstr = pstr2 = argv[i+1];
                }
                else {
                    break;
                }
				int size = 1;
				while ( *pstr2 != '\0') {
					if ( *pstr2 == sep ){
                        size++;
                    }
                    pstr2++;
				}
                
				option = dvector( size );
				*count = 0;
                
				while ( *pstr != '\0' ) {
					if( *pstr == sep && buffer->length != 0 ){
                        if ( !isFloat2(buffer->c) ) {
                            free(option);
                            option = NULL;
                            break;
                        }
						option[(*count)++] = atof( buffer->c );
						StringBuffer_empty( buffer );
					}
					else {
						StringBuffer_append_char(buffer, *pstr);
					}
					pstr++;
				}
                if( option != NULL ){
                    if ( buffer->length != 0) {
                        if ( !isFloat2(buffer->c) ) {
                            free(option);
                            option = NULL;
                        }
                        else {
                            option[(*count)++] = atof( buffer->c );
                        }
                    }
                }
				break;
			}
		
	}
	
	free_StringBuffer(buffer);
	return option;
}

bool args_get_boolean( int argc, char* argv[], const char flag[] ){
    int i;
	for ( i = 0; i < argc; i++) {
        if ( strcmp(argv[i], flag) == 0 ) return true;		
	}
	return false;
}

char* get_program_name(char* argv[]){
	char *name = argv[0]+strlen(argv[0]);
	while( name != argv[0] ){
		if ( *name == '/' || *name == '\\' ) {
			name++;
			break;
		}
		name--;
	}
	return name;
}

void argsparser2_free(args_parser2* args){
	for (int i = 0; i < args->option_count; i++) {
		if(args->options[i].value != NULL){
			free(args->options[i].value);
		}
	}
	free(args);
}

Hashtable * argsparser2_parse(args_parser2* args, char* argv[], int argc){
	
	argsparser2_check(args, argv, argc);
	
	int i = 0;
	for ( ; i < argc; i++) {
		int j = 0;
		bool option_flag = false;
		for(; j < args->option_count; j++){
			if(argv[i][0] == '-') option_flag = true;
			if ( strlen(argv[i]) > 2 && argv[i][0] == '-' && argv[i][1] == '-' && strcmp(argv[i]+2, args->options[j].long_name) == 0 ) {
				break;
			}
			else if (strlen(argv[i]) == 2 && argv[i][0] == '-' && args->options[j].short_name == argv[i][1]){
				break;
			}
		}
		
		if(j != args->option_count && option_flag == true){
			if(args->options[j].value != NULL) free(args->options[j].value);
			
			switch(args->options[j].type){
				case ARGS_OPTION_FLAG:{
					args->options[j].value = String_clone("1");
					break;
				}
				case ARGS_OPTION_BOOLEAN:{
					if(strcasecmp(argv[i+1], "true") == 0 || strcmp(argv[i+1], "1") == 0){
						args->options[j].value = String_clone("1");
					}
					else if(strcasecmp(argv[i+1], "false") == 0 || strcmp(argv[i+1], "0") == 0){
						args->options[j].value = String_clone("0");
					}
					else{
						fprintf(stderr, "Option %s %s is not valid. Expect a boolean argument", argv[i], argv[i+1]);
						exit(1);
					}
					break;
				}
				case ARGS_OPTION_INTEGER:{
					if (!isInt(argv[i+1])) {
						fprintf(stderr, "Option %s %s is not an integer", argv[i], argv[i+1]);
					}
					args->options[j].value = String_clone(argv[i+1]);
					break;
				}
				case ARGS_OPTION_LONG:{
					if (!isInt(argv[i+1])) {
						fprintf(stderr, "Option %s %s is not a long", argv[i], argv[i+1]);
					}
					args->options[j].value = String_clone(argv[i+1]);
					break;
				}
				case ARGS_OPTION_STRING:{
					args->options[j].value = String_clone(argv[i+1]);
					break;
				}
				case ARGS_OPTION_DOUBLE:{
					if (!isFloat(argv[i+1])) {
						fprintf(stderr, "Option %s %s is not a double", argv[i], argv[i+1]);
					}
					args->options[j].value = String_clone(argv[i+1]);
					break;
				}
				case ARGS_OPTION_FLOAT:{
					if (!isFloat(argv[i+1])) {
						fprintf(stderr, "Option %s %s is not a float", argv[i], argv[i+1]);
					}
					args->options[j].value = String_clone(argv[i+1]);
					break;
				}
			}
			if (args->options[j].type != ARGS_OPTION_FLAG) {
				i++;
			}
		}
		else if(option_flag == true){
			fprintf(stderr, "Unknown option: %s\n", argv[i]);
			exit(1);
		}
	}
	
	
	Hashtable* hash = new_Hashtable_string(10);
	for(int i = 0; i < args->option_count; i++){
		char* name = args->options[i].long_name;
		char* value = args->options[i].value;
		if (value == NULL) {
			value = "";
		}
		Hashtable_add(hash, String_clone(name), String_clone(value));
	}
	return hash;
}
