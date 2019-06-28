//
//  mjson.h
//  physher
//
//  Created by Mathieu Fourment on 11/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef mjson_h
#define mjson_h

#include <stdio.h>
#include <stdbool.h>

typedef enum {
	MJSON_ARRAY = 0,
	MJSON_STRING = 1,
	MJSON_PRIMITIVE = 2,
	MJSON_OBJECT = 3,
	MJSON_UNDEFINED = 4,
} json_node_t;

typedef struct json_node {
	struct json_node* parent;
	struct json_node** children;
	char* key;
	void*value;
	size_t child_count;
	json_node_t node_type;
}json_node;

json_node* create_json_node(json_node* parent);
json_node* create_json_node_object(json_node* parent, const char* name);
json_node* create_json_node_array(json_node* parent, char* name);
json_node* create_json_node_parameter(json_node* parent, const char* name, double value, double lower, double upper);
json_node* create_json_node_parameters(json_node* parent, const char* name, size_t dimension, double lower, double upper);
json_node* create_json_node_parameters2(json_node* parent, const char* name, size_t dimension, const double * values, double lower, double upper);
json_node* create_json_node_simplex(json_node* parent, const char* name, size_t dimension);
json_node* create_json_node_simplex2(json_node* parent, const char* name, size_t dimension, const double* values);
void add_json_node(json_node* parent, json_node* child);
json_node* add_json_node_bool(json_node* parent, const char* key, bool value);
json_node* add_json_node_string(json_node* parent, const char* key, const char* value);
json_node* add_json_node_size_t(json_node* parent, const char* key, size_t value);
json_node* add_json_node_double(json_node* parent, const char* key, double value);
json_node* add_json_node_array_double(json_node* parent, const char* key, double* values, size_t dim);
json_node* add_json_node_array_unsigned(json_node* parent, const char* key, unsigned* values, size_t dim);

json_node* create_json_tree(const char* json);
json_node* get_json_node(json_node* node, const char* key);
char* get_json_node_value_string(json_node* node, const char* key);
bool get_json_node_value_bool(json_node* node, const char* key, bool defaultv);
double get_json_node_value_double(json_node* node, const char* key, double defaultv);
int get_json_node_value_int(json_node* node, const char* key, int defaultv);
size_t get_json_node_value_size_t(json_node* node, const char* key, size_t defaultv);
void json_tree_to_string(json_node* node);
bool json_prune_ignored(json_node* node);
void json_tree_print(json_node* node);
void json_free_tree(json_node* node);

void json_check_allowed(json_node* node, char** allowed, int length);
void json_check_required(json_node* node, char** required, int length);
#endif /* mjson_h */
