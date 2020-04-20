//
//  mjson.c
//  physher
//
//  Created by Mathieu Fourment on 11/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "mjson.h"

#include <string.h>
#include <strings.h>

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

json_node* create_json_node_parameters(json_node* parent, const char* name, size_t dimension, double lower, double upper){
	json_node* jnode = create_json_node_object(parent, name);
	add_json_node(parent, jnode);
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "parameter");
	add_json_node_size_t(jnode, "dimension", dimension);
	add_json_node_double(jnode, "lower", lower);
	add_json_node_double(jnode, "upper", upper);
	return jnode;
}

json_node* create_json_node_parameters2(json_node* parent, const char* name, size_t dimension, const double* values, double lower, double upper){
	json_node* jnode = create_json_node_object(parent, NULL);
	add_json_node(parent, jnode);
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "parameter");
	add_json_node_size_t(jnode, "dimension", dimension);
	add_json_node_array_double(jnode, "values", values, dimension);
	add_json_node_double(jnode, "lower", lower);
	add_json_node_double(jnode, "upper", upper);
	return jnode;
}

json_node* create_json_node_simplex(json_node* parent, const char* name, size_t dimension){
	json_node* jnode = create_json_node_object(parent, name);
	add_json_node(parent, jnode);
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "simplex");
	add_json_node_size_t(jnode, "dimension", dimension);
	return jnode;
}

json_node* create_json_node_simplex2(json_node* parent, const char* name, size_t dimension, const double* values){
	json_node* jnode = create_json_node_object(parent, name);
	add_json_node(parent, jnode);
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "simplex");
	add_json_node_size_t(jnode, "dimension", dimension);
	add_json_node_array_double(jnode, "values", values, dimension);
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
	return add_json_node_aux(parent, nkey, nvalue, MJSON_PRIMITIVE);
}

json_node* add_json_node_string(json_node* parent, const char* key, const char* value){
	StringBuffer* buffer = new_StringBuffer(10);
	char* nkey = NULL;
	if (key != NULL) {
		StringBuffer_append_format(buffer, "%s", key);
		nkey = StringBuffer_tochar(buffer);
		StringBuffer_empty(buffer);
	}

	StringBuffer_append_format(buffer, "%s", value);
	char* nvalue = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return add_json_node_aux(parent, nkey, nvalue, MJSON_STRING);
}

json_node* add_json_node_size_t(json_node* parent, const char* key, size_t value){
	StringBuffer* buffer = new_StringBuffer(10);
	char* nkey = NULL;
	if(key != NULL){
		StringBuffer_append_format(buffer, "%s", key);
		nkey = StringBuffer_tochar(buffer);
		StringBuffer_empty(buffer);
	}
	StringBuffer_append_format(buffer, "%zu", value);
	char* nvalue = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	return add_json_node_aux(parent, nkey, nvalue, MJSON_PRIMITIVE);
}


json_node* add_json_node_unsigned(json_node* parent, const char* key, unsigned value){
	StringBuffer* buffer = new_StringBuffer(10);
	char* nkey = NULL;
	if(key != NULL){
		StringBuffer_append_format(buffer, "%s", key);
		nkey = StringBuffer_tochar(buffer);
		StringBuffer_empty(buffer);
	}
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
			StringBuffer_append_format(buffer, "\"-infinity\"", value);
		}
		else{
			StringBuffer_append_format(buffer, "\"infinity\"", value);
		}
	}
	else if(isnan(value)){
		StringBuffer_append_format(buffer, "\"nan\"", value);
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
	new->key = String_clone(key);
	
	for(size_t i = 0; i < dim; i++){
		add_json_node_unsigned(new, NULL, values[i]);
	}
	return new;
}

json_node* add_json_node_array_string(json_node* parent, const char* key, char** values, size_t dim){
	json_node* new = create_json_node(parent);
	add_json_node(parent, new);
	new->node_type = MJSON_ARRAY;
	new->key = String_clone(key);
	
	for(size_t i = 0; i < dim; i++){
		add_json_node_string(new, NULL, values[i]);
	}
	return new;
}

void add_json_node(json_node* parent, json_node* child){
	child->parent = parent;
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
//		if(n->node_type != MJSON_STRING){
//			fprintf(stderr, "value (%s) for key %s is not a string\n", (char*)n->value, key);
//		}
		return (char*)n->value;
	}
	return NULL;
}

bool get_json_node_value_bool(json_node* node, const char* key, bool defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		if(n->node_type != MJSON_PRIMITIVE){
			fprintf(stderr, "value (%s) for key %s is not a number\n", (char*)n->value, key);
		}
		return atoi((char*)n->value);
	}
	return defaultv;
}

double get_json_node_value_double(json_node* node, const char* key, double defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		if(n->node_type != MJSON_PRIMITIVE){
			fprintf(stderr, "value (%s) for key %s is not a number\n", (char*)n->value, key);
		}
		return atof((char*)n->value);
	}
	return defaultv;
}

size_t get_json_node_value_size_t(json_node* node, const char* key, size_t defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		if(n->node_type != MJSON_PRIMITIVE){
			fprintf(stderr, "value (%s) for key %s is not a number\n", (char*)n->value, key);
		}
		size_t v = 0;
		int result = sscanf((char*)n->value, "%zu", &v);
		return v;
	}
	return defaultv;
}

int get_json_node_value_int(json_node* node, const char* key, int defaultv){
	json_node* n = get_json_node(node, key);
	if(n != NULL){
		if(n->node_type != MJSON_PRIMITIVE){
			fprintf(stderr, "value (%s) for key %s is not a number\n", (char*)n->value, key);
		}
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
		else if((json[i] == '"' || json[i] == '\'')  && current->node_type == MJSON_OBJECT){
			json_node* n = create_json_node(current);
			add_json_node(current, n);
			i++;
			StringBuffer_empty(buffer);
			while (json[i] != '"' && json[i] != '\'') {
				StringBuffer_append_char(buffer, json[i]);
				i++;
			}
			
			n->key = StringBuffer_tochar(buffer);
			n->node_type = MJSON_UNDEFINED;
			//printf("key: %s\n", buffer->c);
			current = n;
		}
		// string value
		else if(json[i] == '"' || json[i] == '\''){
			i++;
			StringBuffer_empty(buffer);
			while (json[i] != '"' && json[i] != '\'') {
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
		else if((json[i] >=48 && json[i] <= 57) || json[i] == '.' || json[i] == '+' || json[i] == '-'){
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

bool json_prune_ignored(json_node* node){
	for (int i = 0; i < node->child_count; i++) {
		if (node->children[i]->node_type == MJSON_PRIMITIVE && node->children[i]->key != NULL &&
			strcasecmp(node->children[i]->key, "ignore") == 0 && strcasecmp((char*)node->children[i]->value, "1") == 0) {
			return true;
		}
	}
	
	for (int i = 0; i < node->child_count; i++) {
		bool remove = json_prune_ignored(node->children[i]);
		if (remove) {
			json_free_tree(node->children[i]);
			for (int j = i; j < node->child_count-1; j++) {
				node->children[j] = node->children[j+1];
			}
			node->child_count--;
			i--;
		}
	}
	return false;
}

void json_tree_print_aux(json_node* node, size_t level, FILE* file){
	for(size_t i = 0; i < level; i++) fprintf(file, "  ");
	if(node->node_type == MJSON_STRING && node->key != NULL){
		fprintf(file, "\"%s\":\"%s\"", (char*)node->key, (char*)node->value);
	}
	else if(node->node_type == MJSON_PRIMITIVE && node->key != NULL){
		fprintf(file, "\"%s\":%s", (char*)node->key, (char*)node->value);
	}
	else if(node->node_type == MJSON_STRING){
		fprintf(file, "\"%s\"", (char*)node->value);
	}
	else if(node->node_type == MJSON_PRIMITIVE){
		fprintf(file, "%s", (char*)node->value);
	}
	else if(node->node_type == MJSON_ARRAY){
		fprintf(file, "\"%s\": [\n", (char*)node->key);
	}
	else if(node->node_type == MJSON_OBJECT){
		// root node and anonymous object in arrays
		if(node->parent != NULL && node->key != NULL) fprintf(file, "\"%s\": ", (char*)node->key);
		fprintf(file, "{\n");
	}
	else{
		error("error json_tree_print_aux");
	}
	level++;
	for (int i = 0; i < node->child_count; i++) {
		json_tree_print_aux(node->children[i], level, file);
		if(i != node->child_count-1) fprintf(file, ",\n");
		else fprintf(file, "\n");
	}
	if(node->node_type == MJSON_ARRAY){
		for(size_t i = 0; i < level-1; i++) fprintf(file, "  ");
		fprintf(file, "]");
	}
	else if(node->node_type == MJSON_OBJECT){
		for(size_t i = 0; i < level-1; i++) fprintf(file, "  ");
		fprintf(file, "}");
	}
}

void json_tree_print(json_node* node){
	json_tree_print_aux(node, 0, stdout);
}

void json_tree_fprint(json_node* node, FILE* file){
    json_tree_print_aux(node, 0, file);
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

void json_check_allowed(json_node* node, char** allowed, int length){
	int id = -1;
	int type = -1;
	for (int i = 0; i < node->child_count; i++) {
		if (strcasecmp(node->children[i]->key, "id") == 0) {
			id = i;
		}
		else if (strcasecmp(node->children[i]->key, "type") == 0) {
			type = i;
		}
	}
	if (node->node_type == MJSON_OBJECT && (id == -1 || type == -1)) {
		fprintf(stderr, "Missing id or type in:\n");
		json_tree_print(node);
		exit(12);
	}
	
	for (int i = 0; i < node->child_count; i++) {
		if (node->children[i]->key[0] == '_' || i == id || i == type) continue; // keys starting with _ are comments
		
		int j = 0;
		for ( ; j < length; j++) {
			if (strcasecmp(node->children[i]->key, allowed[j]) == 0) {
				break;
			}
		}
		if (j == length) {
			fprintf(stderr, "Key not recognised: %s in %s of type %s\n", node->children[i]->key, node->children[id]->value, node->children[type]->value);
			fprintf(stderr, "Possible keys:\n");
			for (int j = 0; j < length; j++) {
				fprintf(stderr, " %s\n", allowed[j]);
			}
			fprintf(stderr, "\n");
			exit(12);
		}
	}
}

void json_check_required(json_node* node, char** required, int length){
	int id = -1;
	int type = -1;
	for (int i = 0; i < node->child_count; i++) {
		if (strcasecmp(node->children[i]->key, "id") == 0) {
			id = i;
		}
		else if (strcasecmp(node->children[i]->key, "type") == 0) {
			type = i;
		}
	}
	if (node->node_type == MJSON_OBJECT && (id == -1 || type == -1)) {
		fprintf(stderr, "Missing id or type in:\n");
		json_tree_print(node);
		exit(12);
	}
	bool* found = malloc(length*sizeof(bool));
	for (int i = 0; i < node->child_count; i++) {
		if (node->children[i]->key[0] == '_' || i == id || i == type) continue; // keys starting with _ are comments
		
		for (int j = 0 ; j < length; j++) {
			if (strcasecmp(node->children[i]->key, required[j]) == 0) {
				found[j] = true;
				break;
			}
		}
	}
	for (int j = 0 ; j < length; j++) {
		if(!found[j]){
			fprintf(stderr, "Key not found: %s in %s of type %s\n", required[j], node->children[id]->value, node->children[type]->value);
			exit(12);
		}
	}
	free(found);
}
