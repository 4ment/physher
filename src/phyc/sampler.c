//
//  sampler.c
//  physher
//
//  Created by Mathieu Fourment.
//  Copyright © 2022 Mathieu Fourment. All rights reserved.
//

#include "sampler.h"

static void _sampler_sample(Sampler* sampler) {
    Model* model = sampler->model;
    size_t samples = sampler->samples;
    for (size_t i = 0; i < samples; i++) {
        model->sample(model);
        double logQ = model->logP(model);
        for (size_t j = 0; j < sampler->logger_count; j++) {
            sampler->loggers[j]->write(sampler->loggers[j], i);
        }
    }
}

void _sampler_initialize(Sampler* sampler) {
    for (size_t i = 0; i < sampler->logger_count; i++) {
        sampler->loggers[i]->initialize(sampler->loggers[i]);
    }
}

void _sampler_finalize(Sampler* sampler) {
    for (size_t i = 0; i < sampler->logger_count; i++) {
        sampler->loggers[i]->finalize(sampler->loggers[i]);
    }
}

void _free_Sampler(Sampler* sampler) {
    sampler->model->free(sampler->model);
    for (int i = 0; i < sampler->logger_count; i++) {
        sampler->loggers[i]->free(sampler->loggers[i]);
    }
    free(sampler->loggers);
    free(sampler);
}

Sampler* new_Sampler_from_json(json_node* node, Hashtable* hash) {
    char* allowed[] = {"model", "loggers", "samples"};
    json_check_allowed(node, allowed, sizeof(allowed) / sizeof(allowed[0]));

    char* ref = get_json_node_value_string(node, "model");
    json_node* loggers_node = get_json_node(node, "loggers");

    Sampler* sampler = malloc(sizeof(Sampler));
    sampler->model = Hashtable_get(hash, ref + 1);
    sampler->model->ref_count++;
    sampler->loggers = malloc(loggers_node->child_count * sizeof(Log*));
    sampler->logger_count = loggers_node->child_count;
    for (size_t i = 0; i < loggers_node->child_count; i++) {
        json_node* child = loggers_node->children[i];
        sampler->loggers[i] = new_Log_from_json(child, hash);
    }
    sampler->samples = get_json_node_value_size_t(node, "samples", 1000);
    sampler->sample = _sampler_sample;
    sampler->initialize = _sampler_initialize;
    sampler->finalize = _sampler_finalize;
    sampler->free = _free_Sampler;
    return sampler;
}
