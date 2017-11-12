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

void free_json_tree(json_node* node){
	for (int i = 0; i < node->child_count; i++) {
		free_json_tree(node->parent);
	}
	if(node->child_count > 0){
		free(node->children);
	}
	free(node);
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

json_node* create_json_tree2(const char* json){
	json_node* current = NULL;
	size_t len = strlen(json);
	StringBuffer* buffer = new_StringBuffer(100);
	bool flag = true;
	for(size_t i = 0; i < len; i++){
		if (json[i] == '{') {
			if(current == NULL){
				current = create_json_node(NULL);
			}
			else{
				json_node* n = create_json_node(current);
				add_json_node(current, n);
				current = n;
			}
			current->start = i;
		}
		// inside a key
		else if(json[i] == '"' && flag){
			StringBuffer_empty(buffer);
			i++;
			while (json[i] != '"') {
				StringBuffer_append_char(buffer, json[i]);
				i++;
			}
			
			//printf("=%s= ", buffer->c);
			if(strcasecmp("id", buffer->c) == 0){
				i++;
				while (json[i] != '"') i++;
				i++;
				StringBuffer_empty(buffer);
				while (json[i] != '"') {
					StringBuffer_append_char(buffer, json[i]);
					i++;
				}
				current->id = StringBuffer_tochar(buffer);
				//printf(" %s", buffer->c);
			}
			else if(strcasecmp("type", buffer->c) == 0){
				i++;
				while (json[i] != '"') i++;
				i++;
				StringBuffer_empty(buffer);
				while (json[i] != '"') {
					StringBuffer_append_char(buffer, json[i]);
					i++;
				}
				current->type = StringBuffer_tochar(buffer);
				//printf(" %s", buffer->c);
			}
			//printf("\n");
		}
		// string value
		else if(json[i] == '"'){
			i++;
			while (json[i] != '"') i++;
		}
		// array
		else if(json[i] == '['){
			flag = false;
		}
		else if(json[i] == ':'){
			flag = false;
		}
		if(json[i] == ']' || json[i] == '{'){
			flag = true;
		}
		else if(json[i] == '}'){
			current->end = i;
			if(current->id == NULL){
				json_node* parent = current->parent;
				
			}
			if(current->parent != NULL) current = current->parent;
		}
	}
	return current;
}

json_node* create_json_tree(const char* json){
	json_node* current = NULL;
	size_t len = strlen(json);
	StringBuffer* buffer = new_StringBuffer(100);
	bool flag = true;
	
	for(size_t i = 0; i < len; i++){
		if (json[i] == '{') {
			if(current != NULL && current->node_type == MJSON_ARRAY){
				json_node* n = create_json_node(current);
				add_json_node(current, n);
				current = n;
//				current->start = i;
				n->node_type = MJSON_OBJECT;
			}
			else if(current == NULL){
				current = create_json_node(NULL);
				current->node_type = MJSON_ARRAY;
			}
			flag = true;
		}
		// inside a key
		else if(json[i] == '"' && flag && (current == NULL || current->node_type != MJSON_ARRAY)){
			json_node* n = create_json_node(current);
			add_json_node(current, n);
			i++;
			StringBuffer_empty(buffer);
			while (json[i] != '"') {
				StringBuffer_append_char(buffer, json[i]);
				i++;
			}
			
			n->key = StringBuffer_tochar(buffer);
			n->node_type = MJSON_OBJECT;
			printf("key: %s\n", buffer->c);
			flag = false;
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
			printf("value: %s (%d)\n", buffer->c, current->node_type);
			if (current->node_type == MJSON_ARRAY) {
				json_node* n = create_json_node(current);
				add_json_node(current, n);
				n->value = StringBuffer_tochar(buffer);
				n->node_type = MJSON_STRING;
				current = n;
			}
			else{
				current->value = StringBuffer_tochar(buffer);
				current->node_type = MJSON_STRING;
			}
		}
		// array
		else if(json[i] == '['){
			current->node_type = MJSON_ARRAY;
		}
		else if(json[i] == ']'){
			flag = true;
			printf("%d\n", current->child_count);
			// allow empty array
			if(current->child_count > 0) current = current->parent;
			current = current->parent;
		}
		else if(json[i] == ':'){
			flag = false;
		}
		else if(json[i] == ','){
			flag = true;
			current = current->parent;
		}
		else if(json[i] == '}'){
			current->end = i;
			if(current->parent != NULL)
				current = current->parent;
		}
		else if((json[i] >=48 && json[i] <= 57) || json[i] == '.'){
			StringBuffer_empty(buffer);
			while (json[i] != ' ' && json[i] != ',' && json[i] != ']' && json[i] != '}') {
				StringBuffer_append_char(buffer, json[i]);
				i++;
			}
			json_node* n = create_json_node(current);
			add_json_node(current, n);
			n->value = StringBuffer_tochar(buffer);
			n->node_type = MJSON_PRIMITIVE;
		}
	}
	return current;
}


void json_tree_to_string(json_node* node){
	//printf("key: %s %s %s %zu %zu %zu %d\n", node->key, node->id, node->type, node->start, node->end, node->child_count, node->node_type);
	if(node->node_type == MJSON_STRING){
		printf("key %s value %s\n", node->key, (char*)node->value);
	}
	else{
		printf("key %s\n", node->key);
	}
	for (int i = 0; i < node->child_count; i++) {
		json_tree_to_string(node->children[i]);
	}
}