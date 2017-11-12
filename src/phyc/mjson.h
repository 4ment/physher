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


typedef enum {
	MJSON_ARRAY = 0,
	MJSON_STRING = 1,
	MJSON_PRIMITIVE = 2,
	MJSON_OBJECT = 3
} json_node_t;

typedef struct json_node {
	struct json_node* parent;
	struct json_node** children;
	char* id;
	char* type;
	char* key;
	void*value;
	size_t start;
	size_t end;
	size_t child_count;
	json_node_t node_type;
}json_node;

json_node* create_json_tree(const char* json);
json_node* get_json_node(json_node* node, const char* key);
void json_tree_to_string(json_node* node);
#endif /* mjson_h */
