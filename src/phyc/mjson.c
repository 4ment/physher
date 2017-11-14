//
//  mjson.c
//  physher
//
//  Created by Mathieu Fourment on 11/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "mjson.h"

#include <string.h>

#include "mstring.h"

json_node* create_json_node(json_node* parent){
	json_node* node = malloc(sizeof(json_node));
	node->parent = parent;
	node->children = NULL;
	node->child_count = 0;
	node->id = NULL;
	node->type = NULL;
	node->key = NULL;
	node->value = NULL;
	return node;
}

void add_json_node(json_node* parent, json_node* child){
	if(parent->child_count == 0){
		parent->children = calloc(1, sizeof(json_node));
		parent->children[0] = child;
	}
	parent->child_count++;
	parent->children = realloc(parent->children, sizeof(json_node)*parent->child_count);
	parent->children[parent->child_count-1] = child;
}

json_node* get_json_node(json_node* node, const char* key){
	for (int i = 0; i < node->child_count; i++) {
		if(strcmp(key, node->children[i]->key) == 0){
			return node->children[i];
		}
	}
	
	return NULL;
}

char* get_json_node_value_string(json_node* node, const char* key){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		return (char*)n->value;
	}
	return NULL;
}

json_node* create_json_tree(const char* json){
	json_node* current = create_json_node(NULL);
	json_node* root = current;
	current->node_type = MJSON_OBJECT;
	size_t len = strlen(json);
	StringBuffer* buffer = new_StringBuffer(100);
	size_t i = 1;
	while (json[i] != '{') i++;
	i++;
	while (json[len-1] != '}') {
		len--;
	}
	for(size_t i = 1; i < len-1; i++){
		if (json[i] == '{') {
			if(current->node_type == MJSON_ARRAY){
				json_node* n = create_json_node(current);
				add_json_node(current, n);
				printf("Add object to array %s\n", current->key);
				current = n;
				n->node_type = MJSON_OBJECT;
			}
			else{
				printf("set %s as object\n", current->key);
				current->node_type = MJSON_OBJECT;
			}
		}
		// key
		else if(json[i] == '"' && current->node_type == MJSON_OBJECT){
			json_node* n = create_json_node(current);
			add_json_node(current, n);
			i++;
			StringBuffer_empty(buffer);
			while (json[i] != '"') {
				StringBuffer_append_char(buffer, json[i]);
				i++;
			}
			
			n->key = StringBuffer_tochar(buffer);
			n->node_type = MJSON_UNDEFINED;
			printf("key: %s\n", buffer->c);
			current = n;
		}
		// string value
		else if(json[i] == '"'){
			i++;
			StringBuffer_empty(buffer);
			while (json[i] != '"') {
				StringBuffer_append_char(buffer, json[i]);
				i++;
			}

			// Value is a string
			if (current->node_type == MJSON_UNDEFINED) {
				
				printf("Add string %s to %s\n", buffer->c, current->key);
				current->node_type = MJSON_STRING;
				current->value = StringBuffer_tochar(buffer);
			}
			// Value is part of an array
			else if(current->node_type == MJSON_ARRAY){
				
				printf("Add string %s to array %s\n", buffer->c, current->key);
				json_node* n = create_json_node(current);
				add_json_node(current, n);
				n->value = StringBuffer_tochar(buffer);
				n->node_type = MJSON_STRING;
				current = n;
			}
			else exit(1);
		}
		// array
		else if(json[i] == '['){
			printf("%s is an array\n", current->key);
			current->node_type = MJSON_ARRAY;
		}
		else if(json[i] == ']'){
			printf("close array %s\n", current->key);
			// allow empty array
			if(current->child_count > 0) current = current->parent;
			current = current->parent;
		}
		else if(json[i] == ','){
			if(current->parent != NULL)current = current->parent;
		}
		else if(json[i] == '}'){
			if(current == NULL){
				free_StringBuffer(buffer);
				return root;
			}
			printf("close object %s  %d\n", current->key, (current->parent == NULL));
			if(current->parent != NULL)
				current = current->parent;
		}
		else if((json[i] >=48 && json[i] <= 57) || json[i] == '.'){
			StringBuffer_empty(buffer);
			while (json[i] != ',' && json[i] != ']' && json[i] != '}') {
				StringBuffer_append_char(buffer, json[i]);
				i++;
			}
			StringBuffer_trim(buffer);
			i--;
			
			if (current->node_type == MJSON_ARRAY) {
				printf("Add primitive %s to array %s\n", buffer->c,current->key);
				json_node* n = create_json_node(current);
				add_json_node(current, n);
				n->value = StringBuffer_tochar(buffer);
				n->node_type = MJSON_PRIMITIVE;
				current = n;
			}
			else{
				printf("Add primitive %s to %s\n", buffer->c,current->key);
				current->value = StringBuffer_tochar(buffer);
				current->node_type = MJSON_PRIMITIVE;
			}
		}
	}
	free_StringBuffer(buffer);
	return root;
}


void json_tree_to_string(json_node* node){
	//printf("key: %s %s %s %zu %zu %zu %d\n", node->key, node->id, node->type, node->start, node->end, node->child_count, node->node_type);
	if(node->node_type == MJSON_STRING){
		printf("key %s value %s\n", node->key, (char*)node->value);
	}
	else if(node->node_type == MJSON_PRIMITIVE){
		printf("key %s value %s*\n", node->key, (char*)node->value);
	}
	else{
		printf("key %s\n", node->key);
	}
	for (int i = 0; i < node->child_count; i++) {
		json_tree_to_string(node->children[i]);
	}
}

void json_free_tree(json_node* node){
	for (int i = 0; i < node->child_count; i++) {
		json_free_tree(node->children[i]);
	}
	free(node->key);
	if(node->child_count > 0)free(node->children);
	if(node->value!= NULL) free(node->value);
	free(node);
}
