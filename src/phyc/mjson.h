// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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

json_node* clone_json_node(json_node* parent, json_node* node);
json_node* create_json_node(json_node* parent);
json_node* create_json_node_object(json_node* parent, const char* name);
json_node* create_json_node_array(json_node* parent, char* name);
json_node* create_json_node_parameter(json_node* parent, const char* name, double value, double lower, double upper);
json_node* create_json_node_parameter_n(json_node* parent, const char* name, const double* value, size_t dimension, double lower, double upper);
json_node* create_json_node_parameters(json_node* parent, const char* name, size_t dimension, double lower, double upper);
json_node* create_json_node_parameters2(json_node* parent, const char* name, size_t dimension, const double * values, double lower, double upper);
json_node* create_json_node_parameter_full(json_node* parent, const char* name, double value, size_t dimension, double lower, double upper);
json_node* create_json_node_simplex(json_node* parent, const char* name, size_t dimension);
json_node* create_json_node_simplex2(json_node* parent, const char* name, size_t dimension, const double* values);
json_node* fill_json_node_simplex(json_node* jnode, const char* id, size_t dimension, const double* values);
void add_json_node(json_node* parent, json_node* child);
json_node* add_json_node_bool(json_node* parent, const char* key, bool value);
json_node* add_json_node_string(json_node* parent, const char* key, const char* value);
json_node* add_json_node_size_t(json_node* parent, const char* key, size_t value);
json_node* add_json_node_double(json_node* parent, const char* key, double value);
json_node* add_json_node_array_double(json_node* parent, const char* key, const double* values, size_t dim);
json_node* add_json_node_array_unsigned(json_node* parent, const char* key, const unsigned* values, size_t dim);

json_node* create_json_tree(const char* json);
json_node* get_json_node(json_node* node, const char* key);
char* get_json_node_value_string(json_node* node, const char* key);
bool get_json_node_value_bool(json_node* node, const char* key, bool defaultv);
double get_json_node_value_double(json_node* node, const char* key, double defaultv);
int get_json_node_value_int(json_node* node, const char* key, int defaultv);
size_t get_json_node_value_size_t(json_node* node, const char* key, size_t defaultv);
void get_json_node_value_array_double(json_node* node, const char* key, double* values);
void get_json_node_value_array_int(json_node* node, const char* key, int* values);

void json_tree_to_string(json_node* node);
bool json_prune_ignored(json_node* node);
bool json_prune_underscored(json_node* node);
void json_tree_print(json_node* node);
void json_tree_fprint(json_node* node, FILE* file);
void json_free_tree(json_node* node);

char* get_json_node_value_string_required(json_node* node, const char* key);
double get_json_node_value_double_required(json_node* node, const char* key);
int get_json_node_value_int_required(json_node* node, const char* key);
size_t get_json_node_value_size_t_required(json_node* node, const char* key);
bool get_json_node_value_bool_required(json_node* node, const char* key);

void json_check_allowed(json_node* node, char** allowed, int length);
void json_check_required(json_node* node, char** required, int length);

// Whether a key may, must, or must not appear in a node.
typedef enum {
    JSON_OPTIONAL = 0,  // may be present (default)
    JSON_REQUIRED,      // must be present; missing => peaceful exit
    JSON_FORBIDDEN,     // must not be present (e.g. removed/renamed keys)
} json_field_req;

// Expected JSON type of a field's value, as a bitmask: OR the base types
// together to accept more than one shape, e.g. JSON_OBJECT | JSON_STRING for
// the common "reference or definition" pattern (an inline object or a string
// id-reference). JSON_ANY (0) skips the type check entirely.
typedef enum {
    JSON_ANY = 0,
    JSON_STRING = 1 << 0,
    JSON_NUMBER = 1 << 1,
    JSON_BOOL = 1 << 2,
    JSON_OBJECT = 1 << 3,
    JSON_ARRAY = 1 << 4,
} json_field_type;

// One row of a node's schema. Designated/partial init is intended, e.g.
//   {"epsilon", JSON_OPTIONAL, JSON_NUMBER}
//   {"model", JSON_REQUIRED, JSON_OBJECT | JSON_STRING}
//   {"substitutionmodel", JSON_FORBIDDEN, JSON_ANY, "renamed to 'model'"}
// "id" and "type" are implicitly allowed/required on object nodes and must not
// be listed. Keys starting with '_' are treated as comments and ignored.
typedef struct {
    const char* key;
    json_field_req req;
    json_field_type type;
    const char* hint;  // optional extra context shown on error; may be NULL
} json_field;

// Report a configuration error against `node` (printing its id/type when
// available) and exit cleanly. Single chokepoint for all validation failures.
void json_die(json_node* node, const char* fmt, ...);

// Validate `node` against `schema` (n rows): enforces id/type presence on
// object nodes, rejects unknown keys (with a "did you mean" suggestion),
// rejects JSON_FORBIDDEN keys, requires JSON_REQUIRED keys, and checks value
// types where a concrete json_field_type is given. Replaces json_check_allowed
// and json_check_required. Dies via json_die on the first violation.
void json_validate(json_node* node, const json_field* schema, size_t n);

// Enforce that exactly one of the given keys is defined on `node` — the
// "either/or, but not both, and not neither" pattern. The key list is variadic
// and NULL-terminated, e.g.
//   json_validate_xor(node, "heights", "branch_lengths", NULL);
// Dies via json_die if none or more than one is present. The keys should still
// be declared (as JSON_OPTIONAL) in the node's schema so json_validate accepts
// them; this adds the exclusivity rule.
void json_validate_xor(json_node* node, ...);

// Enforce that the given keys are co-required — the "all together, or none at
// all" pattern. The key list is variadic and NULL-terminated, e.g.
//   json_validate_co_required(node, "mu", "sigma", NULL);
// Dies via json_die if some (but not all) of the keys are present. The keys
// should still be declared (as JSON_OPTIONAL) in the node's schema so
// json_validate accepts them; this adds the grouping rule.
void json_validate_co_required(json_node* node, ...);
#endif /* mjson_h */
