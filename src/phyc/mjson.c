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
#include "ctype.h"

json_node* create_json_node(json_node* parent){
	json_node* node = malloc(sizeof(json_node));
	node->parent = parent;
	node->children = NULL;
	node->child_count = 0;
	node->key = NULL;
	node->value = NULL;
	node->node_type = MJSON_UNDEFINED;
	return node;
}

json_node* create_json_node_object(json_node* parent, const char* name){
	json_node* node = malloc(sizeof(json_node));
	node->parent = parent;
	node->children = NULL;
	node->child_count = 0;
	node->key = (name == NULL ? NULL : String_clone(name));
	node->value = NULL;
	node->node_type = MJSON_OBJECT;
	return node;
}

json_node* create_json_node_array(json_node* parent, char* name){
	json_node* node = malloc(sizeof(json_node));
	node->parent = parent;
	node->children = NULL;
	node->child_count = 0;
	node->key = String_clone(name);
	node->value = NULL;
	node->node_type = MJSON_ARRAY;
	return node;
}

json_node* create_json_node_parameter(json_node* parent, const char* name, double value, double lower, double upper){
	json_node* jnode = create_json_node_object(parent, name);
	add_json_node(parent, jnode);
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "parameter");
	add_json_node_double(jnode, "value", value);
	add_json_node_double(jnode, "lower", lower);
	add_json_node_double(jnode, "upper", upper);
	return jnode;
}

json_node* add_json_node_aux(json_node* parent, char* key, char* value, json_node_t type){
//	if(parent->node_type != MJSON_UNDEFINED && parent->node_type != MJSON_OBJECT){
//		error("Can only add key and value to an object (add_json_node_aux)\n");
//	}
	json_node* new = create_json_node(parent);
	//new->id = id;
	new->node_type = type;
	new->key = key;
	new->value = value;
	add_json_node(parent, new);
	return new;
}

json_node* add_json_node_bool(json_node* parent, const char* key, bool value){
	char* nkey = String_clone(key);
	char* nvalue = NULL;
	if(value) nvalue = String_clone("1");
	else nvalue = String_clone("0");
	return add_json_node_aux(parent, nkey, nvalue, MJSON_STRING);
}

json_node* add_json_node_string(json_node* parent, const char* key, const char* value){
	StringBuffer* buffer = new_StringBuffer(10);
	StringBuffer_append_format(buffer, "%s", key);
	char* nkey = StringBuffer_tochar(buffer);
	StringBuffer_empty(buffer);
	StringBuffer_append_format(buffer, "%s", value);
	char* nvalue = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return add_json_node_aux(parent, nkey, nvalue, MJSON_STRING);
}

json_node* add_json_node_size_t(json_node* parent, const char* key, size_t value){
	StringBuffer* buffer = new_StringBuffer(10);
	StringBuffer_append_format(buffer, "%s", key);
	char* nkey = StringBuffer_tochar(buffer);
	StringBuffer_empty(buffer);
	StringBuffer_append_format(buffer, "%zu", value);
	char* nvalue = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return add_json_node_aux(parent, nkey, nvalue, MJSON_PRIMITIVE);
}


json_node* add_json_node_unsigned(json_node* parent, const char* key, unsigned value){
	StringBuffer* buffer = new_StringBuffer(10);
	StringBuffer_append_format(buffer, "%s", key);
	char* nkey = StringBuffer_tochar(buffer);
	StringBuffer_empty(buffer);
	StringBuffer_append_format(buffer, "%u", value);
	char* nvalue = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return add_json_node_aux(parent, nkey, nvalue, MJSON_PRIMITIVE);
}

json_node* add_json_node_double(json_node* parent, const char* key, double value){
	StringBuffer* buffer = new_StringBuffer(10);
	char* nkey = NULL;
	if(key != NULL){
		StringBuffer_append_format(buffer, "%s", key);
		nkey = StringBuffer_tochar(buffer);
		StringBuffer_empty(buffer);
	}
	if(isinf(value)){
		if (value < 0) {
			StringBuffer_append_format(buffer, "-infinity", value);
		}
		else{
			StringBuffer_append_format(buffer, "infinity", value);
		}
	}
	else if(isnan(value)){
		StringBuffer_append_format(buffer, "nan", value);
	}
	else{
		StringBuffer_append_format(buffer, "%f", value);
	}
	char* nvalue = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return add_json_node_aux(parent, nkey, nvalue, MJSON_PRIMITIVE);
}

json_node* add_json_node_array_double(json_node* parent, const char* key, double* values, size_t dim){
	json_node* new = create_json_node(parent);
	add_json_node(parent, new);
	new->node_type = MJSON_ARRAY;
	
	StringBuffer* buffer = new_StringBuffer(10);
	StringBuffer_append_format(buffer, "%s", key);
	new->key = StringBuffer_tochar(buffer);
	
	for(size_t i = 0; i < dim; i++){
		add_json_node_double(new, NULL, values[i]);
	}
	
	free_StringBuffer(buffer);
	return new;
}

json_node* add_json_node_array_unsigned(json_node* parent, const char* key, unsigned* values, size_t dim){
	json_node* new = create_json_node(parent);
	add_json_node(parent, new);
	new->node_type = MJSON_ARRAY;
	
	StringBuffer* buffer = new_StringBuffer(10);
	StringBuffer_append_format(buffer, "%s", key);
	new->key = StringBuffer_tochar(buffer);
	
	for(size_t i = 0; i < dim; i++){
		add_json_node_unsigned(new, NULL, values[i]);
	}
	
	free_StringBuffer(buffer);
	return new;
}

void add_json_node(json_node* parent, json_node* child){
	parent->child_count++;
	if(parent->child_count == 0){
		parent->children = calloc(1, sizeof(json_node*));
	}
	else{
		parent->children = realloc(parent->children, sizeof(json_node*)*parent->child_count);
	}
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

bool get_json_node_value_bool(json_node* node, const char* key, bool defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		return atoi((char*)n->value);
	}
	return defaultv;
}

double get_json_node_value_double(json_node* node, const char* key, double defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		return atof((char*)n->value);
	}
	return defaultv;
}

size_t get_json_node_value_size_t(json_node* node, const char* key, size_t defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		size_t v = 0;
		int result = sscanf((char*)n->value, "%zu", &v);
		return v;
	}
	return defaultv;
}

int get_json_node_value_int(json_node* node, const char* key, int defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		return atoi((char*)n->value);
	}
	return defaultv;
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
				//printf("Add object to array %s\n", current->key);
				current = n;
				n->node_type = MJSON_OBJECT;
			}
			else{
				//printf("set %s as object\n", current->key);
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
			//printf("key: %s\n", buffer->c);
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
				
				//printf("Add string %s to %s\n", buffer->c, current->key);
				current->node_type = MJSON_STRING;
				current->value = StringBuffer_tochar(buffer);
			}
			// Value is part of an array
			else if(current->node_type == MJSON_ARRAY){
				
				//printf("Add string %s to array %s\n", buffer->c, current->key);
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
			//printf("%s is an array\n", current->key);
			current->node_type = MJSON_ARRAY;
		}
		else if(json[i] == ']'){
			//printf("close array %s\n", current->key);
			// allow empty array
//			if(current->child_count > 0) current = current->parent;
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
			//printf("close object %s  %d\n", current->key, (current->parent == NULL));
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
				//printf("Add primitive %s to array %s\n", buffer->c,current->key);
				json_node* n = create_json_node(current);
				add_json_node(current, n);
				n->value = StringBuffer_tochar(buffer);
				n->node_type = MJSON_PRIMITIVE;
				current = n;
			}
			else{
				//printf("Add primitive %s to %s\n", buffer->c,current->key);
				current->value = StringBuffer_tochar(buffer);
				current->node_type = MJSON_PRIMITIVE;
			}
		}
		else if( tolower(json[i]) == 't' ){
			current->value = String_clone("1");
			current->node_type = MJSON_PRIMITIVE;
			i += 3;
		}
		else if( tolower(json[i]) == 'f' ){
			current->value = String_clone("0");
			current->node_type = MJSON_PRIMITIVE;
			i+= 4;
		}
	}
	free_StringBuffer(buffer);
	return root;
}


void json_tree_to_string(json_node* node){
	//printf("key: %s %s %s %zu %zu %zu %d\n", node->key, node->id, node->type, node->start, node->end, node->child_count, node->node_type);
	if(node->node_type == MJSON_STRING){
		//printf("key %s value %s\n", node->key, (char*)node->value);
	}
	else if(node->node_type == MJSON_PRIMITIVE){
		//printf("key %s value %s*\n", node->key, (char*)node->value);
	}
	else{
		//printf("key %s\n", node->key);
	}
	for (int i = 0; i < node->child_count; i++) {
		json_tree_to_string(node->children[i]);
	}
}

void json_tree_print_aux(json_node* node, size_t level){
	for(size_t i = 0; i < level; i++) printf("  ");
	if(node->node_type == MJSON_STRING && node->key != NULL){
		printf("\"%s\":\"%s\"", (char*)node->key, (char*)node->value);
	}
	else if(node->node_type == MJSON_PRIMITIVE && node->key != NULL){
		printf("\"%s\":%s", (char*)node->key, (char*)node->value);
	}
	else if(node->node_type == MJSON_STRING){
		printf("\"%s\"", (char*)node->value);
	}
	else if(node->node_type == MJSON_PRIMITIVE){
		printf("%s", (char*)node->value);
	}
	else if(node->node_type == MJSON_ARRAY){
		printf("\"%s\": [\n", (char*)node->key);
	}
	else if(node->node_type == MJSON_OBJECT){
		// root node and anonymous object in arrays
		if(node->parent != NULL && node->key != NULL) printf("\"%s\": ", (char*)node->key);
		printf("{\n");
	}
	else{
		error("error json_tree_print_aux");
	}
	level++;
	for (int i = 0; i < node->child_count; i++) {
		json_tree_print_aux(node->children[i], level);
		if(i != node->child_count-1) printf(",\n");
		else printf("\n");
	}
	if(node->node_type == MJSON_ARRAY){
		for(size_t i = 0; i < level-1; i++) printf("  ");
		printf("]");
	}
	else if(node->node_type == MJSON_OBJECT){
		for(size_t i = 0; i < level-1; i++) printf("  ");
		printf("}");
	}
}

void json_tree_print(json_node* node){
	json_tree_print_aux(node, 0);
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
