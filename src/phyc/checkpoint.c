// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "checkpoint.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "filereader.h"  // load_file
#include "mjson.h"

// On-disk checkpoint format. The file is a JSON object:
//
//   {
//     "format": "physher-checkpoint",
//     "version": 1,
//     "parameters": {
//       "<parameter name>": [v0, v1, ...],
//       ...
//     }
//   }
//
// Only the "parameters" section is written/read for now, but the envelope is
// deliberately a top-level object so future state can be added as sibling
// sections without breaking older files, e.g. an "optimizer" section (learning
// rate, iteration count, moment estimates) or an "operators" section (MCMC
// operator tuning parameters and accept/reject counts). A reader ignores
// sections it does not understand, so adding one is backward compatible; bump
// "version" only for incompatible changes to an existing section.
#define CHECKPOINT_FORMAT "physher-checkpoint"
#define CHECKPOINT_VERSION 1

// Look up a direct child by key, tolerating children with a NULL key (which a
// malformed/foreign file may contain) rather than dereferencing them like
// get_json_node would.
static json_node* checkpoint_child(json_node* node, const char* key) {
    for (size_t i = 0; i < node->child_count; i++) {
        const char* child_key = node->children[i]->key;
        if (child_key != NULL && strcmp(child_key, key) == 0) {
            return node->children[i];
        }
    }
    return NULL;
}

// Write one parameter value as a JSON number at full double precision so a
// restart round-trips exactly. Non-finite values are emitted as JSON strings so
// the file stays valid JSON; atof() parses them back on read.
static void checkpoint_write_value(FILE* file, double value) {
    if (isfinite(value)) {
        fprintf(file, "%.17g", value);
    } else if (isnan(value)) {
        fprintf(file, "\"nan\"");
    } else {
        fprintf(file, value < 0 ? "\"-inf\"" : "\"inf\"");
    }
}

// Write the "parameters" section (indented under the root object, no trailing
// newline) mapping each parameter name to the JSON array of its values.
static void checkpoint_write_parameters(FILE* file, Parameters* parameters) {
    fprintf(file, "  \"parameters\": {\n");
    size_t count = Parameters_count(parameters);
    for (size_t i = 0; i < count; i++) {
        Parameter* p = Parameters_at(parameters, i);
        fprintf(file, "    \"%s\": [", Parameter_name(p));
        size_t dim = Parameter_size(p);
        for (size_t j = 0; j < dim; j++) {
            if (j != 0) fprintf(file, ", ");
            checkpoint_write_value(file, Parameter_value_at(p, j));
        }
        fprintf(file, "]%s\n", (i + 1 < count) ? "," : "");
    }
    fprintf(file, "  }");
}

void checkpoint_save(const char* file_path, Parameters* parameters) {
    // Write to a temporary file first, then atomically rename it over the
    // destination. A crash mid-write therefore never corrupts an existing
    // checkpoint: on disk we always have either the intact previous file or the
    // complete new one.
    StringBuffer* tmp = new_StringBuffer(10);
    StringBuffer_set_string(tmp, file_path);
    StringBuffer_append_string(tmp, ".tmp");

    FILE* file = fopen(tmp->c, "w");
    if (file == NULL) {
        fprintf(stderr, "checkpoint_save: cannot open '%s' for writing\n", tmp->c);
        free_StringBuffer(tmp);
        return;
    }

    fprintf(file, "{\n");
    fprintf(file, "  \"format\": \"%s\",\n", CHECKPOINT_FORMAT);
    fprintf(file, "  \"version\": %d,\n", CHECKPOINT_VERSION);
    // Future sections (optimizer state, MCMC operator tuning, ...) go here as
    // sibling keys, each followed by a comma, before "parameters".
    checkpoint_write_parameters(file, parameters);
    fprintf(file, "\n}\n");

    if (fclose(file) != 0) {
        fprintf(stderr, "checkpoint_save: error writing '%s'\n", tmp->c);
        remove(tmp->c);
        free_StringBuffer(tmp);
        return;
    }

    if (rename(tmp->c, file_path) != 0) {
        fprintf(stderr, "checkpoint_save: cannot move '%s' to '%s'\n", tmp->c,
                file_path);
        remove(tmp->c);
    }
    free_StringBuffer(tmp);
}

void checkpoint_apply(const char* file_path, Parameters* parameters) {
    char* content = load_file(file_path);  // exits if the file cannot be read
    json_node* root = create_json_tree(content);
    free(content);
    if (root == NULL) {
        fprintf(stderr, "checkpoint_apply: '%s' is not valid JSON\n", file_path);
        return;
    }

    int version = get_json_node_value_int(root, "version", 0);
    if (version != 0 && version != CHECKPOINT_VERSION) {
        fprintf(stderr,
                "checkpoint_apply: '%s' has version %d (expected %d); "
                "attempting to read anyway\n",
                file_path, version, CHECKPOINT_VERSION);
    }

    json_node* params = checkpoint_child(root, "parameters");
    if (params == NULL) {
        fprintf(stderr, "checkpoint_apply: '%s' has no \"parameters\" section\n",
                file_path);
        json_free_tree(root);
        return;
    }

    size_t restored = 0;
    size_t count = Parameters_count(parameters);
    for (size_t i = 0; i < count; i++) {
        Parameter* p = Parameters_at(parameters, i);
        const char* name = Parameter_name(p);
        json_node* entry = checkpoint_child(params, name);
        if (entry == NULL) {
            continue;  // not in the checkpoint; leave the parameter untouched
        }
        if (entry->node_type != MJSON_ARRAY) {
            fprintf(stderr, "checkpoint_apply: \"%s\" is not an array; skipping\n",
                    name);
            continue;
        }
        size_t dim = Parameter_size(p);
        if (entry->child_count != dim) {
            fprintf(stderr,
                    "checkpoint_apply: \"%s\" has %zu value(s) but parameter has "
                    "dimension %zu; skipping\n",
                    name, entry->child_count, dim);
            continue;
        }
        double* values = malloc(sizeof(double) * dim);
        for (size_t j = 0; j < dim; j++) {
            values[j] = atof((char*)entry->children[j]->value);
        }
        Parameter_set_values(p, values);
        free(values);
        restored++;
    }

    printf("checkpoint: restored %zu of %zu parameters from %s\n", restored, count,
           file_path);
    json_free_tree(root);
}
