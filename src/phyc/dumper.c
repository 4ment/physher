// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "dumper.h"

#include <string.h>

#include "mjson.h"
#include "mstring.h"

void dumper_get_references(json_node* node, Hashtable* hash, struct Dumper* dumper){
    json_node* models_node = get_json_node(node, "models");
    dumper->model_count = 0;
    if (models_node != NULL) {
        if(models_node->node_type == MJSON_ARRAY){
            for (int i = 0; i < models_node->child_count; i++) {
                json_node* child = models_node->children[i];
                char* ref = (char*)child->value;
                // it's a ref
                if (child->node_type == MJSON_STRING && ref[0] == '&') {
                    if (dumper->model_count == 0) {
                        dumper->models = malloc(sizeof(Model*));
                    }
                    else{
                        dumper->models = realloc(dumper->models, sizeof(Model*)*(dumper->model_count+1));
                    }
                    dumper->models[dumper->model_count] = Hashtable_get(hash, ref+1);
                    dumper->model_count++;
                }
                else{
                    fprintf(stderr, "Dumper only accepts models: `%s'\n", ref);
                    exit(12);
                }
            }
        }
        // it's a ref
        else if(models_node->node_type == MJSON_STRING){
            char* ref = (char*)models_node->value;
            if (ref[0] == '&') {
                dumper->models = malloc(sizeof(Model*));
                dumper->models[dumper->model_count] = Hashtable_get(hash, ref+1);
                dumper->model_count++;
            }
        }
        else{
            fprintf(stderr, "Dmuper %s cannot read value of key %s\n", get_json_node_value_string(node, "id"), models_node->key);
            exit(2);
        }
    }
}

static void _dump(struct Dumper* dumper){
    json_node* jroot = create_json_node(NULL);
    jroot->node_type = MJSON_OBJECT;
    for (int i = 0; i < dumper->model_count; i++) {
        Model* model = dumper->models[i];
        model->jsonize(model, jroot);
    }
    json_tree_fprint(jroot, dumper->file);
    json_free_tree(jroot);
}

static void _free_Dumper(struct Dumper* dumper){
    if(dumper->filename != NULL){
        free(dumper->filename);
        fclose(dumper->file);
    }
    if(dumper->model_count > 0) free(dumper->models);
    free(dumper);
}

struct Dumper* new_Dumper_from_json(json_node* node, Hashtable* hash){
    static const json_field schema[] = {
        {"file", JSON_OPTIONAL, JSON_ANY},
        {"models", JSON_OPTIONAL, JSON_ANY},
        {"parameters", JSON_OPTIONAL, JSON_ANY},
    };
    json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

    struct Dumper* dumper = malloc(sizeof(struct Dumper));
    dumper->parameters = NULL;
    dumper->parameter_count = 0;
    dumper->model_count = 0;
    dumper->models = NULL;

    dumper_get_references(node, hash, dumper);

    json_node* file_node = get_json_node(node, "file");
    dumper->file = stdout;
    dumper->filename = NULL;
    if(file_node != NULL){
        char* filename = file_node->value;
        if (strcmp(filename, "stderr") != 0 && strcmp(filename, "stdout") != 0) {
            dumper->filename = String_clone(filename);
            dumper->file = fopen(dumper->filename, "w");
        }
        else if(strcmp(filename, "stderr") == 0){
            dumper->file = stderr;
        }
        else{
            dumper->file = stdout;
        }
    }
    dumper->dump = _dump;
    dumper->free = _free_Dumper;

    return dumper;
}
