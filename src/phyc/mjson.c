// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "mjson.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "mstring.h"
#include "ctype.h"
#include "transforms.h"

json_node* clone_json_node(json_node* parent, json_node* node) {
    json_node* clone = NULL;
    if (node->node_type == MJSON_OBJECT) {
        clone = create_json_node_object(parent, node->key);
    } else if (node->node_type == MJSON_ARRAY) {
        clone = create_json_node_array(parent, node->key);
    } else if (node->node_type == MJSON_STRING || node->node_type == MJSON_PRIMITIVE) {
        clone = create_json_node(parent);
        clone->key = node->key == NULL ? NULL : String_clone(node->key);
        clone->value = String_clone(node->value);
        clone->node_type = node->node_type;
    }

    if (node->child_count > 0) {
        clone->children = calloc(node->child_count, sizeof(json_node*));
        for (size_t i = 0; i < node->child_count; i++) {
            clone->children[i] = clone_json_node(clone, node->children[i]);
            clone->child_count++;
        }
    }
    return clone;
}

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
	return create_json_node_parameter_n(parent, name, &value, 1, lower, upper);
}

json_node* create_json_node_parameter_n(json_node* jnode, const char* name, const double* value, size_t dimension, double lower, double upper){
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "parameter");
	if(dimension > 1){
		add_json_node_size_t(jnode, "dimension", dimension);
		add_json_node_array_double(jnode, "value", value, dimension);
	}
	else{
		add_json_node_double(jnode, "value", *value);
	}
	if(isinf(upper)){
		add_json_node_string(jnode, "upper", "infinity");
	}
	else{
		add_json_node_double(jnode, "upper", upper);
	}
	if(isinf(lower) && lower < 0){
		add_json_node_string(jnode, "lower", "-infinity");
	}
	else{
		add_json_node_double(jnode, "lower", lower);
	}
	return jnode;
}

json_node* create_json_node_parameters(json_node* jnode, const char* name, size_t dimension, double lower, double upper){
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "parameter");
	add_json_node_size_t(jnode, "dimension", dimension);

	if(isinf(upper)){
		add_json_node_string(jnode, "upper", "infinity");
	}
	else{
		add_json_node_double(jnode, "upper", upper);
	}
	if(isinf(lower) && lower < 0){
		add_json_node_string(jnode, "lower", "-infinity");
	}
	else{
		add_json_node_double(jnode, "lower", lower);
	}
	return jnode;
}

json_node* create_json_node_parameter_full(json_node* jnode, const char* name, double value, size_t dimension, double lower, double upper){
	add_json_node_string(jnode, "id", name);
	add_json_node_string(jnode, "type", "parameter");
	add_json_node_size_t(jnode, "dimension", dimension);
	add_json_node_double(jnode, "value", value);

	if(isinf(upper)){
		add_json_node_string(jnode, "upper", "infinity");
	}
	else{
		add_json_node_double(jnode, "upper", upper);
	}
	if(isinf(lower) && lower < 0){
		add_json_node_string(jnode, "lower", "-infinity");
	}
	else{
		add_json_node_double(jnode, "lower", lower);
	}
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

// Fill an existing (already parent-attached) node as a simplex parameter in the
// transformed-parameter form: a "simplex" parameter wrapping an unconstrained
// "<id>.unconstrained" parameter of dimension K-1, seeded from the constrained
// `values` via the stick-breaking transform.
json_node* fill_json_node_simplex(json_node* jnode, const char* id, size_t dimension, const double* values){
	add_json_node_string(jnode, "id", id);
	add_json_node_string(jnode, "type", "simplex");

	double* unconstrained = malloc(sizeof(double) * (dimension - 1));
	_transform_simplex_stan(values, unconstrained, dimension - 1, 0, 1);

	json_node* junc = create_json_node_object(jnode, "x");
	add_json_node(jnode, junc);
	StringBuffer* buffer = new_StringBuffer(10);
	StringBuffer_append_format(buffer, "%s.unconstrained", id);
	add_json_node_string(junc, "id", buffer->c);
	free_StringBuffer(buffer);
	add_json_node_string(junc, "type", "parameter");
	add_json_node_array_double(junc, "x", unconstrained, dimension - 1);
	free(unconstrained);
	return jnode;
}

json_node* create_json_node_simplex2(json_node* parent, const char* name, size_t dimension, const double* values){
	json_node* jnode = create_json_node_object(parent, name);
	add_json_node(parent, jnode);
	return fill_json_node_simplex(jnode, name, dimension, values);
}

json_node* create_json_node_simplex(json_node* parent, const char* name, size_t dimension){
	double* values = malloc(sizeof(double) * dimension);
	for (size_t i = 0; i < dimension; i++) values[i] = 1.0 / dimension;
	json_node* jnode = create_json_node_simplex2(parent, name, dimension, values);
	free(values);
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

json_node* add_json_node_array_double(json_node* parent, const char* key, const double* values, size_t dim){
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

json_node* add_json_node_array_unsigned(json_node* parent, const char* key, const unsigned* values, size_t dim){
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

void get_json_node_value_array_double(json_node* node, const char* key,
                                      double* values) {
    json_node* n = get_json_node(node, key);
    for (size_t i = 0; i < n->child_count; i++) {
        values[i] = atof((char*)n->children[i]->value);
    }
}

void get_json_node_value_array_int(json_node* node, const char* key, int* values) {
    json_node* n = get_json_node(node, key);
    for (size_t i = 0; i < n->child_count; i++) {
        values[i] = atoi((char*)n->children[i]->value);
    }
}

json_node* create_json_tree(const char* json){
	json_node* current = create_json_node(NULL);
	json_node* root = current;
	current->node_type = MJSON_OBJECT;
	size_t len = strlen(json);
	StringBuffer* buffer = new_StringBuffer(100);
	size_t i = 0;
	size_t lineNbr = 0;
	while (i < len && json[i] != '{') {
		if(json[i] == '\n')
			lineNbr++;
		i++;
	}
	i++;
	while (json[len-1] != '}') {
		len--;
	}
	for(size_t i = 1; i < len-1; i++){
		if(json[i] == '\n'){
			lineNbr++;
		}
		// support comment but breaks JSON standard
		// only use for debugging purpose
		else if (json[i] == '/') {
			i++;
			if (json[i] == '/') {
				i++;
				// printf("%zu = ", lineNbr);
				while (json[i] != '\n') {
					// printf("%c", json[i]);
					i++;
				}
				i--;
				// printf("\n");
			}
			else if(json[i] == '*'){
				size_t startComment = lineNbr;
				i++;
				while(json[i] != '*' && json[i+1] != '/' && i < len-2) {
					i++;
					if(json[i] == '\n'){
						lineNbr++;
					}
				}
				if(i == len-2){
					fprintf(stderr, "Runnaway comment starting at %zu\n", startComment);
					exit(1);
				}
				i++;
			}
			else{
				fprintf(stderr, "Error parsing JSON: unexpected character after /\n");
				exit(1);
			}

		}
		else if (json[i] == '{') {
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
			while (json[i] != ',' && json[i] != ']' && json[i] != '}' && json[i] != '\n') {
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


bool json_prune_underscored(json_node* node){
    for (int i = 0; i < node->child_count; i++) {
        if(node->children[i]->key != NULL && node->children[i]->key[0] == '_'){
            json_free_tree(node->children[i]);
            for (int j = i; j < node->child_count-1; j++) {
                node->children[j] = node->children[j+1];
            }
            node->child_count--;
            i--;
        }
        else{
            json_prune_underscored(node->children[i]);
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

void json_free_node(json_node* node) {
    free(node->key);
    if (node->child_count > 0) free(node->children);
    if (node->value != NULL) free(node->value);
    free(node);
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

// --- required getters: die instead of silently returning a default --------

char* get_json_node_value_string_required(json_node* node, const char* key) {
	json_node* n = get_json_node(node, key);
	if (n == NULL) json_die(node, "required key \"%s\" is missing", key);
	if (n->node_type != MJSON_STRING) {
		json_die(node, "value for key \"%s\" must be a string", key);
	}
	return (char*)n->value;
}

double get_json_node_value_double_required(json_node* node, const char* key) {
	json_node* n = get_json_node(node, key);
	if (n == NULL) json_die(node, "required key \"%s\" is missing", key);
	if (n->node_type != MJSON_PRIMITIVE) {
		json_die(node, "value for key \"%s\" must be a number", key);
	}
	return atof((char*)n->value);
}

int get_json_node_value_int_required(json_node* node, const char* key) {
	json_node* n = get_json_node(node, key);
	if (n == NULL) json_die(node, "required key \"%s\" is missing", key);
	if (n->node_type != MJSON_PRIMITIVE) {
		json_die(node, "value for key \"%s\" must be a number", key);
	}
	return atoi((char*)n->value);
}

size_t get_json_node_value_size_t_required(json_node* node, const char* key) {
	json_node* n = get_json_node(node, key);
	if (n == NULL) json_die(node, "required key \"%s\" is missing", key);
	if (n->node_type != MJSON_PRIMITIVE) {
		json_die(node, "value for key \"%s\" must be a number", key);
	}
	size_t v = 0;
	sscanf((char*)n->value, "%zu", &v);
	return v;
}

bool get_json_node_value_bool_required(json_node* node, const char* key) {
	json_node* n = get_json_node(node, key);
	if (n == NULL) json_die(node, "required key \"%s\" is missing", key);
	if (n->node_type != MJSON_PRIMITIVE) {
		json_die(node, "value for key \"%s\" must be a boolean", key);
	}
	return atoi((char*)n->value);
}

// --- schema validation -----------------------------------------------------

// Locate the "id"/"type" children (case-insensitive); -1 if absent.
static void json_find_id_type(json_node* node, int* id, int* type) {
	*id = -1;
	*type = -1;
	for (int i = 0; i < node->child_count; i++) {
		if (strcasecmp(node->children[i]->key, "id") == 0) *id = i;
		else if (strcasecmp(node->children[i]->key, "type") == 0) *type = i;
	}
}

void json_die(json_node* node, const char* fmt, ...) {
	int id = -1, type = -1;
	json_find_id_type(node, &id, &type);
	fprintf(stderr, "physher: error in node");
	if (id != -1) fprintf(stderr, " \"%s\"", (char*)node->children[id]->value);
	if (type != -1) {
		fprintf(stderr, " (type \"%s\")", (char*)node->children[type]->value);
	}
	fprintf(stderr, ":\n  ");
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, "\n");
	exit(12);
}

// Levenshtein distance, used to suggest a likely-intended key on typos.
static int json_edit_distance(const char* a, const char* b) {
	size_t la = strlen(a), lb = strlen(b);
	int* prev = malloc((lb + 1) * sizeof(int));
	int* cur = malloc((lb + 1) * sizeof(int));
	for (size_t j = 0; j <= lb; j++) prev[j] = (int)j;
	for (size_t i = 1; i <= la; i++) {
		cur[0] = (int)i;
		for (size_t j = 1; j <= lb; j++) {
			int cost = (tolower(a[i - 1]) == tolower(b[j - 1])) ? 0 : 1;
			int del = prev[j] + 1;
			int ins = cur[j - 1] + 1;
			int sub = prev[j - 1] + cost;
			int m = del < ins ? del : ins;
			cur[j] = m < sub ? m : sub;
		}
		int* tmp = prev; prev = cur; cur = tmp;
	}
	int d = prev[lb];
	free(prev);
	free(cur);
	return d;
}

// True if the node type satisfies the declared field type.
static bool json_type_matches(json_node_t got, json_field_type want) {
	if (want == JSON_ANY) return true;
	json_field_type got_flag;
	switch (got) {
		case MJSON_STRING: got_flag = JSON_STRING; break;
		// The parser does not distinguish numbers, booleans and null: they are
		// all primitives, so a primitive satisfies either JSON_NUMBER or JSON_BOOL.
		case MJSON_PRIMITIVE: got_flag = JSON_NUMBER | JSON_BOOL; break;
		case MJSON_OBJECT: got_flag = JSON_OBJECT; break;
		case MJSON_ARRAY: got_flag = JSON_ARRAY; break;
		default: return false;
	}
	return (got_flag & want) != 0;
}

static const char* json_type_name(json_field_type want) {
	static const struct {
		json_field_type bit;
		const char* name;
	} names[] = {
		{JSON_STRING, "a string"}, {JSON_NUMBER, "a number"},
		{JSON_BOOL, "a boolean"},  {JSON_OBJECT, "an object"},
		{JSON_ARRAY, "an array"},
	};
	static char buf[128];
	buf[0] = '\0';
	int count = 0;
	for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); i++) {
		if ((want & names[i].bit) == 0) continue;
		if (count++ > 0) strcat(buf, " or ");
		strcat(buf, names[i].name);
	}
	return count > 0 ? buf : "valid";
}

void json_validate(json_node* node, const json_field* schema, size_t n) {
	int id = -1, type = -1;
	json_find_id_type(node, &id, &type);
	if (node->node_type == MJSON_OBJECT && (id == -1 || type == -1)) {
		json_die(node, "object is missing an \"id\" or \"type\"");
	}

	// Check every present key against the schema.
	for (int i = 0; i < node->child_count; i++) {
		const char* key = node->children[i]->key;
		if (key[0] == '_' || i == id || i == type) continue;

		const json_field* match = NULL;
		for (size_t j = 0; j < n; j++) {
			if (strcasecmp(key, schema[j].key) == 0) {
				match = &schema[j];
				break;
			}
		}

		if (match == NULL) {
			// Suggest the closest schema key for typos.
			const char* best = NULL;
			int best_d = 0;
			for (size_t j = 0; j < n; j++) {
				int d = json_edit_distance(key, schema[j].key);
				if (best == NULL || d < best_d) { best = schema[j].key; best_d = d; }
			}
			if (best != NULL && best_d <= 2) {
				json_die(node, "unknown key \"%s\" — did you mean \"%s\"?", key, best);
			}
			json_die(node, "unknown key \"%s\"", key);
		}

		if (match->req == JSON_FORBIDDEN) {
			if (match->hint != NULL) {
				json_die(node, "key \"%s\" is no longer allowed — %s", key,
				         match->hint);
			}
			json_die(node, "key \"%s\" is no longer allowed", key);
		}

		if (!json_type_matches(node->children[i]->node_type, match->type)) {
			json_die(node, "value for key \"%s\" must be %s", key,
			         json_type_name(match->type));
		}
	}

	// Ensure every required key is present.
	for (size_t j = 0; j < n; j++) {
		if (schema[j].req != JSON_REQUIRED) continue;
		if (get_json_node(node, schema[j].key) == NULL) {
			json_die(node, "required key \"%s\" is missing", schema[j].key);
		}
	}
}

void json_validate_xor(json_node* node, ...) {
	va_list ap;
	const char* key;

	// First pass: count how many keys are present (and how many in total).
	size_t found = 0, total = 0;
	va_start(ap, node);
	while ((key = va_arg(ap, const char*)) != NULL) {
		total++;
		if (get_json_node(node, key) != NULL) found++;
	}
	va_end(ap);
	if (found == 1) return;

	// Second pass: build a "a", "b" or "c" list for the diagnostic.
	StringBuffer* buf = new_StringBuffer(64);
	size_t j = 0;
	va_start(ap, node);
	while ((key = va_arg(ap, const char*)) != NULL) {
		if (j > 0) StringBuffer_append_string(buf, j + 1 == total ? " or " : ", ");
		StringBuffer_append_format(buf, "\"%s\"", key);
		j++;
	}
	va_end(ap);

	if (found == 0) {
		json_die(node, "exactly one of %s must be defined", buf->c);
	} else {
		json_die(node, "%s are mutually exclusive; define exactly one", buf->c);
	}
	free_StringBuffer(buf);  // unreachable: json_die exits
}

void json_validate_co_required(json_node* node, ...) {
	va_list ap;
	const char* key;

	// First pass: count how many keys are present (and how many in total).
	size_t found = 0, total = 0;
	va_start(ap, node);
	while ((key = va_arg(ap, const char*)) != NULL) {
		total++;
		if (get_json_node(node, key) != NULL) found++;
	}
	va_end(ap);
	// Co-required group: valid when all are present, or none are.
	if (found == 0 || found == total) return;

	// Second pass: build a "a", "b" or "c" list for the diagnostic.
	StringBuffer* buf = new_StringBuffer(64);
	size_t j = 0;
	va_start(ap, node);
	while ((key = va_arg(ap, const char*)) != NULL) {
		if (j > 0) StringBuffer_append_string(buf, j + 1 == total ? " or " : ", ");
		StringBuffer_append_format(buf, "\"%s\"", key);
		j++;
	}
	va_end(ap);

	json_die(node, "%s must be defined together; define all or none", buf->c);
	free_StringBuffer(buf);  // unreachable: json_die exits
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
			fprintf(stderr, "Key not recognised: %s in %s of type %s\n", node->children[i]->key, (char*)node->children[id]->value, (char*)node->children[type]->value);
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
			fprintf(stderr, "Key not found: %s in %s of type %s\n", required[j], (char*)node->children[id]->value, (char*)node->children[type]->value);
			exit(12);
		}
	}
	free(found);
}
