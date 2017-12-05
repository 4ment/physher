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
	char* id;
	char* type;
	char* key;
	void*value;
	size_t child_count;
	json_node_t node_type;
}json_node;

json_node* create_json_tree(const char* json);
json_node* get_json_node(json_node* node, const char* key);
char* get_json_node_value_string(json_node* node, const char* key);
bool get_json_node_value_bool(json_node* node, const char* key, bool defaultv);
double get_json_node_value_double(json_node* node, const char* key, double defaultv);
int get_json_node_value_int(json_node* node, const char* key, int defaultv);
size_t get_json_node_value_size_t(json_node* node, const char* key, size_t defaultv);
void json_tree_to_string(json_node* node);
void json_free_tree(json_node* node);
#endif /* mjson_h */
