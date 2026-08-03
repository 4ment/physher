// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "bootstrap.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#ifdef GSL_DISABLED
#include "random.h"
#else
#include <gsl/gsl_randist.h>
#endif

#include "matrix.h"
#include "treelikelihood.h"

#ifdef GSL_DISABLED
// Multinomial(n, w/sum(w)) on the Mersenne Twister of random.c, the stand-in for
// gsl_ran_multinomial in a GSL-free build. Drawing each of the n sites
// independently and adding it to the pattern whose cumulative-weight interval
// contains it is the definition of the resampling, at O(n log K); a pattern of
// weight 0 spans an empty interval and can never be drawn.
// `cumulative` is scratch of length count.
static void _multinomial(size_t n, size_t count, const double* w, double* cumulative,
                         unsigned int* draws) {
    double total = 0;
    for (size_t j = 0; j < count; j++) {
        total += w[j];
        cumulative[j] = total;
        draws[j] = 0;
    }
    for (size_t i = 0; i < n; i++) {
        // genrand_real2 is [0,1)
        double u = total * genrand_real2();
        // First interval whose upper bound is above u.
        size_t low = 0;
        size_t high = count - 1;
        while (low < high) {
            size_t mid = low + (high - low) / 2;
            if (u < cumulative[mid]) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        draws[low]++;
    }
}
#endif

// Draw the multiplicities of one replicate. Resampling nsites sites with
// replacement, each site landing on the pattern it belongs to, makes the vector
// of pattern counts Multinomial(nsites, weights/nsites) -- exact, not an
// approximation. Neither implementation needs the weights normalised.
static void _bootstrap_draw(Bootstrap* bootstrap, size_t p) {
    const SitePattern* sp = bootstrap->patterns[p];
    unsigned int* draws = bootstrap->draws[p];
    double* counts = bootstrap->counts[p];

#ifdef GSL_DISABLED
    // counts is the cumulative-weight scratch until it receives the draw below.
    _multinomial(bootstrap->nsites[p], sp->count, bootstrap->weights[p], counts, draws);
#else
    gsl_ran_multinomial(bootstrap->rng, sp->count, bootstrap->nsites[p],
                        bootstrap->weights[p], draws);
#endif
    for (int j = 0; j < sp->count; j++) {
        counts[j] = draws[j];
    }
}

// Install `counts` (or the observed weights when counts is NULL) as the pattern
// multiplicities of pattern set p, and invalidate every likelihood reading it.
static void _bootstrap_apply(Bootstrap* bootstrap, size_t p, const double* counts) {
    SitePattern* sp = bootstrap->patterns[p];
    const double* multiplicities = counts != NULL ? counts : bootstrap->weights[p];

    if (bootstrap->compact) {
        SitePattern_compact_into(sp, multiplicities, bootstrap->views[p]);
    } else {
        memcpy(sp->weights, multiplicities, sizeof(double) * sp->count);
    }

    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        if (bootstrap->pattern_index[i] != p) continue;
        SingleTreeLikelihood* tlk = bootstrap->treelikelihoods[i]->obj;
        if (bootstrap->compact) {
            SingleTreeLikelihood_set_sitepattern(tlk, bootstrap->views[p]);
        } else {
            SingleTreeLikelihood_update_weights(tlk);
        }
    }
}

// Put the observed alignment back, so anything running after the bootstrap
// (asr, ppsites, a final logger) sees the real data.
static void _bootstrap_restore_observed(Bootstrap* bootstrap) {
    for (size_t p = 0; p < bootstrap->pattern_count; p++) {
        SitePattern* sp = bootstrap->patterns[p];
        memcpy(sp->weights, bootstrap->weights[p], sizeof(double) * sp->count);
        for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
            if (bootstrap->pattern_index[i] != p) continue;
            SingleTreeLikelihood* tlk = bootstrap->treelikelihoods[i]->obj;
            if (bootstrap->compact) {
                SingleTreeLikelihood_set_sitepattern(tlk, sp);
            } else {
                SingleTreeLikelihood_update_weights(tlk);
            }
        }
    }
}

static void _bootstrap_run(Bootstrap* bootstrap) {
    Parameters_store_value(bootstrap->parameters, bootstrap->estimate);

    for (size_t i = 0; i < bootstrap->logger_count; i++) {
        bootstrap->loggers[i]->initialize(bootstrap->loggers[i]);
    }

    for (size_t r = 0; r < bootstrap->replicates; r++) {
        if (bootstrap->reset) {
            Parameters_restore_value(bootstrap->parameters, bootstrap->estimate);
        }

        bool observed = bootstrap->include_observed && r == 0;
        for (size_t p = 0; p < bootstrap->pattern_count; p++) {
            if (observed) {
                _bootstrap_apply(bootstrap, p, NULL);
            } else {
                // Partitions are resampled independently of each other.
                _bootstrap_draw(bootstrap, p);
                _bootstrap_apply(bootstrap, p, bootstrap->counts[p]);
            }
        }

        for (size_t i = 0; i < bootstrap->optimizer_count; i++) {
            double logP;
            opt_optimize(bootstrap->optimizers[i], &logP);
        }

        for (size_t i = 0; i < bootstrap->logger_count; i++) {
            bootstrap->loggers[i]->write(bootstrap->loggers[i], r);
        }
    }

    for (size_t i = 0; i < bootstrap->logger_count; i++) {
        bootstrap->loggers[i]->finalize(bootstrap->loggers[i]);
    }

    _bootstrap_restore_observed(bootstrap);
    Parameters_restore_value(bootstrap->parameters, bootstrap->estimate);
}

static void _free_Bootstrap(Bootstrap* bootstrap) {
    for (size_t p = 0; p < bootstrap->pattern_count; p++) {
        if (bootstrap->views != NULL) free_SitePattern_view(bootstrap->views[p]);
        free(bootstrap->weights[p]);
        free(bootstrap->counts[p]);
        free(bootstrap->draws[p]);
    }
    free(bootstrap->views);
    free(bootstrap->weights);
    free(bootstrap->counts);
    free(bootstrap->draws);
    free(bootstrap->nsites);
    free(bootstrap->patterns);
    free(bootstrap->pattern_index);

    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        bootstrap->treelikelihoods[i]->free(bootstrap->treelikelihoods[i]);
    }
    free(bootstrap->treelikelihoods);

    for (size_t i = 0; i < bootstrap->optimizer_count; i++) {
        free_Optimizer(bootstrap->optimizers[i]);
    }
    free(bootstrap->optimizers);

    for (size_t i = 0; i < bootstrap->logger_count; i++) {
        bootstrap->loggers[i]->free(bootstrap->loggers[i]);
    }
    free(bootstrap->loggers);

    free_Parameters(bootstrap->parameters);
    free(bootstrap->estimate);
    free(bootstrap);
}

// Collect the tree likelihoods named by "treelikelihood" (a ref or an array of
// refs) and the distinct SitePatterns they read.
static void _bootstrap_collect_models(Bootstrap* bootstrap, json_node* node,
                                      Hashtable* hash, const char* id) {
    json_node* tlk_node = get_json_node(node, "treelikelihood");
    size_t count = tlk_node->node_type == MJSON_ARRAY ? tlk_node->child_count : 1;

    bootstrap->treelikelihoods = malloc(sizeof(Model*) * count);
    bootstrap->pattern_index = malloc(sizeof(size_t) * count);
    bootstrap->treelikelihood_count = count;
    bootstrap->patterns = malloc(sizeof(SitePattern*) * count);
    bootstrap->pattern_count = 0;

    for (size_t i = 0; i < count; i++) {
        json_node* child =
            tlk_node->node_type == MJSON_ARRAY ? tlk_node->children[i] : tlk_node;
        if (child->node_type != MJSON_STRING) {
            json_die(node,
                     "%s - \"treelikelihood\" must be a reference such as "
                     "\"@treelikelihood\"",
                     id);
        }
        Model* model = safe_get_reference_model((char*)child->value, hash, id);
        if (model->type != MODEL_TREELIKELIHOOD) {
            json_die(node,
                     "%s - \"treelikelihood\" must reference a tree "
                     "likelihood, got a %s ('%s')",
                     id, model_type_strings[model->type], (char*)child->value);
        }
        model->ref_count++;
        bootstrap->treelikelihoods[i] = model;

        // Two partitions can share one alignment: resample it once and swap it
        // into both, or the two would disagree on the replicate.
        SitePattern* sp = ((SingleTreeLikelihood*)model->obj)->sp;
        size_t p = 0;
        while (p < bootstrap->pattern_count && bootstrap->patterns[p] != sp) p++;
        if (p == bootstrap->pattern_count) {
            bootstrap->patterns[p] = sp;
            bootstrap->pattern_count++;
        }
        bootstrap->pattern_index[i] = p;
    }
}

static void _bootstrap_allocate_patterns(Bootstrap* bootstrap) {
    size_t n = bootstrap->pattern_count;
    bootstrap->weights = malloc(sizeof(double*) * n);
    bootstrap->counts = malloc(sizeof(double*) * n);
    bootstrap->draws = malloc(sizeof(unsigned int*) * n);
    bootstrap->nsites = malloc(sizeof(size_t) * n);
    bootstrap->views = bootstrap->compact ? malloc(sizeof(SitePattern*) * n) : NULL;

    for (size_t p = 0; p < n; p++) {
        const SitePattern* sp = bootstrap->patterns[p];
        bootstrap->weights[p] = clone_dvector(sp->weights, sp->count);
        bootstrap->counts[p] = dvector(sp->count);
        bootstrap->draws[p] = uivector(sp->count);
        // The weights are integer site counts stored as doubles; sum them
        // rather than trusting sp->nsites, which some pattern constructors set
        // from the alignment length instead.
        double total = 0;
        for (int j = 0; j < sp->count; j++) total += sp->weights[j];
        bootstrap->nsites[p] = (size_t)(total + 0.5);
        if (bootstrap->compact) {
            bootstrap->views[p] = new_SitePattern_view(sp);
        }
    }
}

Bootstrap* new_Bootstrap_from_json(json_node* node, Hashtable* hash) {
    static const json_field schema[] = {
        {"compact", JSON_OPTIONAL, JSON_BOOL},
        {"include_observed", JSON_OPTIONAL, JSON_BOOL},
        {"loggers", JSON_REQUIRED, JSON_ARRAY},
        {"optimizer", JSON_REQUIRED, JSON_OBJECT | JSON_ARRAY},
        {"parameters", JSON_OPTIONAL, JSON_ANY},
        {"replicates", JSON_OPTIONAL, JSON_NUMBER},
        {"reset", JSON_OPTIONAL, JSON_BOOL},
        {"treelikelihood", JSON_REQUIRED, JSON_STRING | JSON_ARRAY},
    };
    json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

    const char* id = get_json_node_value_string(node, "id");
    Bootstrap* bootstrap = malloc(sizeof(Bootstrap));

    bootstrap->replicates = get_json_node_value_size_t(node, "replicates", 100);
    bootstrap->compact = get_json_node_value_bool(node, "compact", true);
    bootstrap->reset = get_json_node_value_bool(node, "reset", true);
    bootstrap->include_observed =
        get_json_node_value_bool(node, "include_observed", false);
#ifndef GSL_DISABLED
    bootstrap->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
#endif

    _bootstrap_collect_models(bootstrap, node, hash, id);
    _bootstrap_allocate_patterns(bootstrap);

    json_node* opt_node = get_json_node(node, "optimizer");
    bootstrap->optimizer_count =
        opt_node->node_type == MJSON_ARRAY ? opt_node->child_count : 1;
    bootstrap->optimizers = malloc(sizeof(Optimizer*) * bootstrap->optimizer_count);
    for (size_t i = 0; i < bootstrap->optimizer_count; i++) {
        json_node* child =
            opt_node->node_type == MJSON_ARRAY ? opt_node->children[i] : opt_node;
        bootstrap->optimizers[i] = new_Optimizer_from_json(child, hash);
    }

    json_node* loggers_node = get_json_node(node, "loggers");
    bootstrap->logger_count = loggers_node->child_count;
    bootstrap->loggers = malloc(sizeof(Trace*) * bootstrap->logger_count);
    for (size_t i = 0; i < bootstrap->logger_count; i++) {
        bootstrap->loggers[i] = new_Trace_from_json(loggers_node->children[i], hash);
    }

    // What gets snapshotted at the point estimate and restored before each
    // replicate. Defaults to everything the optimizers move.
    bootstrap->parameters = new_Parameters(1);
    if (get_json_node(node, "parameters") != NULL) {
        get_parameters_references(node, hash, bootstrap->parameters);
    } else {
        for (size_t i = 0; i < bootstrap->optimizer_count; i++) {
            Parameters* parameters = opt_parameters(bootstrap->optimizers[i]);
            for (size_t j = 0; parameters != NULL && j < Parameters_count(parameters);
                 j++) {
                Parameter* parameter = Parameters_at(parameters, j);
                bool seen = false;
                for (size_t k = 0; k < Parameters_count(bootstrap->parameters); k++) {
                    seen |= Parameters_at(bootstrap->parameters, k) == parameter;
                }
                if (!seen) Parameters_add(bootstrap->parameters, parameter);
            }
        }
    }
    if (Parameters_count(bootstrap->parameters) == 0) {
        json_die(node,
                 "%s - no parameters to restart replicates from: the "
                 "optimizers declare none, so \"parameters\" is required",
                 id);
    }

    bootstrap->estimate_size = 0;
    for (size_t i = 0; i < Parameters_count(bootstrap->parameters); i++) {
        bootstrap->estimate_size +=
            Parameter_size(Parameters_at(bootstrap->parameters, i));
    }
    bootstrap->estimate = dvector(bootstrap->estimate_size);

    bootstrap->run = _bootstrap_run;
    bootstrap->free = _free_Bootstrap;
    return bootstrap;
}
