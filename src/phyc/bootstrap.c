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
#include "mstring.h"
#include "treelikelihood.h"

#if defined(PTHREAD_ENABLED)
#include <pthread.h>
#endif

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

// Take the per-pattern site-model state every replicate has to be rebuilt from.
// Only the empirical CAT site model keeps any; see the fields in bootstrap.h.
static void _bootstrap_store_assignments(Bootstrap* bootstrap) {
    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        if (bootstrap->site_categories[i] != NULL) free(bootstrap->site_categories[i]);
        bootstrap->site_categories[i] = NULL;
        bootstrap->site_patterns[i] = NULL;

        SingleTreeLikelihood* tlk = bootstrap->treelikelihoods[i]->obj;
        if (tlk->sm == NULL || tlk->sm->site_category == NULL) continue;

        const SitePattern* sp = bootstrap->patterns[bootstrap->pattern_index[i]];
        bootstrap->site_categories[i] = clone_ivector(tlk->sm->site_category, sp->count);
        bootstrap->site_patterns[i] = tlk->sm->sp;
    }
}

// Reinstate the per-pattern site-model state of tree likelihood `i` for the
// pattern set that has just been installed. Under compaction the surviving
// patterns are packed densely, so position k of the replicate is the k-th
// original pattern with a nonzero multiplicity and the assignment is gathered in
// that order; carrying it over positionally would score every pattern with
// another pattern's rate category. The site model also holds its own reference
// to the pattern set, which is what _cat_update normalises the category rates
// against and which SingleTreeLikelihood_set_sitepattern does not move.
//
// Gathering from the stored assignment overrides any carry-over between
// replicates under "reset": false. There is nothing to carry over: consecutive
// replicates compact to different pattern sets, so the previous assignment is
// indexed in a space this one does not share.
static void _bootstrap_apply_sitemodel(Bootstrap* bootstrap, Model* mtlk, size_t i,
                                       size_t p, SitePattern** views,
                                       const double* multiplicities) {
    SingleTreeLikelihood* tlk = mtlk->obj;
    SiteModel* sm = tlk->sm;
    if (sm == NULL || sm->site_category == NULL) return;

    if (bootstrap->compact) {
        const SitePattern* sp = bootstrap->patterns[p];
        const int* stored = bootstrap->site_categories[i];
        size_t k = 0;
        for (int j = 0; j < sp->count; j++) {
            if (multiplicities[j] != 0.0) sm->site_category[k++] = stored[j];
        }
        sm->sp = views[p];
    }
    // The CAT normaliser is a function of the weights, and nothing else marks it
    // stale when they move underneath it, as _cv_install in crossvalidation.c
    // also has to do.
    sm->need_update = true;
}

// Install `counts` (or the observed weights when counts is NULL) as the pattern
// multiplicities of pattern set p, and invalidate every likelihood reading it.
// `treelikelihoods` and `views` are those of the model graph being set up: the
// bootstrap's own, or a worker's (see _BootstrapGraph).
static void _bootstrap_apply(Bootstrap* bootstrap, Model** treelikelihoods,
                             SitePattern** views, size_t p, const double* counts) {
    SitePattern* sp = bootstrap->patterns[p];
    const double* multiplicities = counts != NULL ? counts : bootstrap->weights[p];

    if (bootstrap->compact) {
        SitePattern_compact_into(sp, multiplicities, views[p]);
    } else {
        memcpy(sp->weights, multiplicities, sizeof(double) * sp->count);
    }

    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        if (bootstrap->pattern_index[i] != p) continue;
        SingleTreeLikelihood* tlk = treelikelihoods[i]->obj;
        if (bootstrap->compact) {
            SingleTreeLikelihood_set_sitepattern(tlk, views[p]);
        } else {
            SingleTreeLikelihood_update_weights(tlk);
        }
        _bootstrap_apply_sitemodel(bootstrap, treelikelihoods[i], i, p, views,
                                   multiplicities);
    }
}

// Put the observed alignment back into the tree likelihoods of a model graph,
// so anything running after the bootstrap (asr, ppsites, a final logger) sees
// the real data, and so nothing frees a replicate pattern set as its own.
static void _bootstrap_restore_observed(Bootstrap* bootstrap, Model** treelikelihoods) {
    for (size_t p = 0; p < bootstrap->pattern_count; p++) {
        SitePattern* sp = bootstrap->patterns[p];
        if (!bootstrap->compact) {
            memcpy(sp->weights, bootstrap->weights[p], sizeof(double) * sp->count);
        }
        for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
            if (bootstrap->pattern_index[i] != p) continue;
            SingleTreeLikelihood* tlk = treelikelihoods[i]->obj;
            if (bootstrap->compact) {
                SingleTreeLikelihood_set_sitepattern(tlk, sp);
            } else {
                SingleTreeLikelihood_update_weights(tlk);
            }
            // The observed weights are back in place, so the assignment and the
            // site model's own reference to the pattern set go back with them:
            // the site model outlives the bootstrap and frees that reference.
            SiteModel* sm = tlk->sm;
            if (sm != NULL && sm->site_category != NULL &&
                bootstrap->site_categories[i] != NULL) {
                memcpy(sm->site_category, bootstrap->site_categories[i],
                       sizeof(int) * sp->count);
                sm->sp = bootstrap->site_patterns[i];
                sm->need_update = true;
            }
        }
    }
}

static void _bootstrap_optimize(Optimizer** optimizers, size_t count) {
    for (size_t i = 0; i < count; i++) {
        double logP;
        opt_optimize(optimizers[i], &logP);
    }
}

static void _bootstrap_run_serial(Bootstrap* bootstrap, size_t rounds) {
    for (size_t r = 0; r < rounds; r++) {
        if (bootstrap->reset) {
            Parameters_restore_value(bootstrap->parameters, bootstrap->estimate);
        }

        bool observed = bootstrap->include_observed && r == 0;
        for (size_t p = 0; p < bootstrap->pattern_count; p++) {
            if (observed) {
                _bootstrap_apply(bootstrap, bootstrap->treelikelihoods,
                                 bootstrap->views, p, NULL);
            } else {
                // Partitions are resampled independently of each other.
                _bootstrap_draw(bootstrap, p);
                _bootstrap_apply(bootstrap, bootstrap->treelikelihoods,
                                 bootstrap->views, p, bootstrap->counts[p]);
            }
        }

        _bootstrap_optimize(bootstrap->optimizers, bootstrap->optimizer_count);

        for (size_t i = 0; i < bootstrap->logger_count; i++) {
            bootstrap->loggers[i]->write(bootstrap->loggers[i], r);
        }
    }
}

#if defined(PTHREAD_ENABLED)

#define BOOTSTRAP_RNG_KEY "RANDOM_GENERATOR!@"

// A worker's copy of the model graph: every model and parameter of the
// bootstrap's hashtable, cloned into a hashtable of its own, with optimizers
// built against it. Only the pattern sets themselves are shared, and only read.
typedef struct _BootstrapGraph {
    Hashtable* hash;
    Model** models;  // references taken by Model_clone_into
    size_t model_count;
    Parameters* orphans;  // parameters registered on their own, owned
    Parameters** lists;   // parameter lists registered on their own, owned
    size_t list_count;
    Model** treelikelihoods;  // borrowed from `models`, as bootstrap->treelikelihoods
    SitePattern** views;      // owned, one per pattern set
    Optimizer** optimizers;   // owned, as bootstrap->optimizers
    Parameters* state;        // as bootstrap->state
#ifndef GSL_DISABLED
    gsl_rng* rng;  // owned, so that nothing a worker samples races on the original
#endif
} _BootstrapGraph;

// Copy of an optimizer node with every optimizer it holds silenced: workers run
// concurrently and their progress would interleave on stdout.
static void _bootstrap_silence(json_node* node) {
    if (node->node_type == MJSON_OBJECT) {
        char* type = get_json_node_value_string(node, "type");
        if (type != NULL && strcasecmp(type, "optimizer") == 0) {
            json_node* verbosity = get_json_node(node, "verbosity");
            if (verbosity == NULL) {
                add_json_node_size_t(node, "verbosity", 0);
            } else {
                free(verbosity->value);
                verbosity->value = String_clone("0");
                verbosity->node_type = MJSON_PRIMITIVE;
            }
        }
    }
    for (size_t i = 0; i < node->child_count; i++) {
        _bootstrap_silence(node->children[i]);
    }
}

static void _bootstrap_graph_keep_model(_BootstrapGraph* graph, Model* model) {
    graph->models = realloc(graph->models, sizeof(Model*) * (graph->model_count + 1));
    graph->models[graph->model_count++] = model;
}

// Clone every entry of the bootstrap's hashtable of the given kind, and make it
// resolvable under the key it was registered with there.
static void _bootstrap_graph_clone_entries(Bootstrap* bootstrap, _BootstrapGraph* graph,
                                           phyc_tag_t tag, bool trees) {
    Hashtable_init_iterator(bootstrap->hash);
    HashEntry* entry;
    while ((entry = Hashtable_next(bootstrap->hash)) != NULL) {
        const char* key = HashEntry_key(entry);
        void* value = HashEntry_value(entry);
        if (phyc_tag_of(value) != tag || Hashtable_exists(graph->hash, key)) continue;
        void* clone = NULL;
        if (tag == PHYC_TAG_MODEL) {
            Model* model = value;
            // Trees first: they own their node heights, which other models (a
            // prior on the root height) may read.
            if ((model->type == MODEL_TREE) != trees) continue;
            clone = Model_clone_into(model, graph->hash);
            _bootstrap_graph_keep_model(graph, clone);
        } else if (tag == PHYC_TAG_PARAMETER) {
            clone = clone_Parameter_into(value, graph->hash);
            Parameters_move(graph->orphans, clone);
        } else {
            clone = clone_Parameters_into(value, graph->hash);
            graph->lists =
                realloc(graph->lists, sizeof(Parameters*) * (graph->list_count + 1));
            graph->lists[graph->list_count++] = clone;
        }
        if (!Hashtable_exists(graph->hash, key)) Hashtable_add(graph->hash, key, clone);
    }
}

static _BootstrapGraph* _new_BootstrapGraph(Bootstrap* bootstrap) {
    _BootstrapGraph* graph = calloc(1, sizeof(_BootstrapGraph));
    assert(graph);
    graph->hash = new_Hashtable_string(100);
    hashtable_set_key_ownership(graph->hash, false);
    hashtable_set_value_ownership(graph->hash, false);
#ifndef GSL_DISABLED
    // A copy rather than a draw: seeding from the bootstrap's generator would
    // shift every replicate drawn after it away from those of a serial run.
    graph->rng = gsl_rng_clone(bootstrap->rng);
    Hashtable_add(graph->hash, BOOTSTRAP_RNG_KEY, graph->rng);
#endif

    // Everything the optimizers might reference, not only the tree likelihoods.
    graph->orphans = new_Parameters(1);
    _bootstrap_graph_clone_entries(bootstrap, graph, PHYC_TAG_MODEL, true);
    _bootstrap_graph_clone_entries(bootstrap, graph, PHYC_TAG_MODEL, false);
    _bootstrap_graph_clone_entries(bootstrap, graph, PHYC_TAG_PARAMETER, false);
    _bootstrap_graph_clone_entries(bootstrap, graph, PHYC_TAG_PARAMETERS, false);

    graph->treelikelihoods = malloc(sizeof(Model*) * bootstrap->treelikelihood_count);
    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        graph->treelikelihoods[i] =
            Hashtable_get(graph->hash, bootstrap->treelikelihoods[i]->name);
        assert(graph->treelikelihoods[i] != NULL);
    }

    graph->views = malloc(sizeof(SitePattern*) * bootstrap->pattern_count);
    for (size_t p = 0; p < bootstrap->pattern_count; p++) {
        graph->views[p] = new_SitePattern_view(bootstrap->patterns[p]);
    }

    json_node* node = bootstrap->optimizer_node;
    graph->optimizers = malloc(sizeof(Optimizer*) * bootstrap->optimizer_count);
    for (size_t i = 0; i < bootstrap->optimizer_count; i++) {
        json_node* child = node->node_type == MJSON_ARRAY ? node->children[i] : node;
        json_node* silent = clone_json_node(NULL, child);
        _bootstrap_silence(silent);
        graph->optimizers[i] = new_Optimizer_from_json(silent, graph->hash);
        json_free_tree(silent);
    }

    graph->state = new_Parameters(Parameters_count(bootstrap->state));
    for (size_t i = 0; i < Parameters_count(bootstrap->state); i++) {
        const char* name = Parameters_name(bootstrap->state, i);
        Parameter* p = Hashtable_get(graph->hash, name);
        if (p == NULL || phyc_tag_of(p) != PHYC_TAG_PARAMETER) {
            fprintf(stderr,
                    "bootstrap - parameter %s has no counterpart in the clone of "
                    "the model\n",
                    name);
            exit(1);
        }
        Parameters_add(graph->state, p);
    }
    return graph;
}

static void _free_BootstrapGraph(Bootstrap* bootstrap, _BootstrapGraph* graph) {
    // The tree likelihoods free the pattern set they read: it must be the
    // original again, not one of the views about to be freed.
    _bootstrap_restore_observed(bootstrap, graph->treelikelihoods);
    for (size_t i = 0; i < bootstrap->optimizer_count; i++) {
        free_Optimizer(graph->optimizers[i]);
    }
    free(graph->optimizers);
    for (size_t p = 0; p < bootstrap->pattern_count; p++) {
        free_SitePattern_view(graph->views[p]);
    }
    free(graph->views);
    free(graph->treelikelihoods);
    free_Parameters(graph->state);
    for (size_t i = 0; i < graph->list_count; i++) free_Parameters(graph->lists[i]);
    free(graph->lists);
    free_Parameters(graph->orphans);
    for (size_t i = 0; i < graph->model_count; i++) {
        graph->models[i]->free(graph->models[i]);
    }
    free(graph->models);
    free_Hashtable(graph->hash);
#ifndef GSL_DISABLED
    gsl_rng_free(graph->rng);
#endif
    free(graph);
}

// The category assignment of every empirical CAT site model, which a CAT
// optimizer rewrites without moving any parameter. Under compaction it is
// indexed by the replicate pattern set, which is the same in every graph.
static int** _bootstrap_copy_categories(Bootstrap* bootstrap, Model** treelikelihoods) {
    int** categories = NULL;
    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        SingleTreeLikelihood* tlk = treelikelihoods[i]->obj;
        if (tlk->sm == NULL || tlk->sm->site_category == NULL) continue;
        if (categories == NULL) {
            categories = calloc(bootstrap->treelikelihood_count, sizeof(int*));
        }
        categories[i] = clone_ivector(tlk->sm->site_category, tlk->sp->count);
    }
    return categories;
}

// Install, and release, assignments taken by _bootstrap_copy_categories.
static void _bootstrap_install_categories(Bootstrap* bootstrap, int** categories) {
    if (categories == NULL) return;
    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        if (categories[i] == NULL) continue;
        SingleTreeLikelihood* tlk = bootstrap->treelikelihoods[i]->obj;
        memcpy(tlk->sm->site_category, categories[i], sizeof(int) * tlk->sp->count);
        tlk->sm->need_update = true;
        SingleTreeLikelihood_update_all_nodes(tlk);
        free(categories[i]);
    }
    free(categories);
}

// State shared by the workers and the thread logging their results.
typedef struct _BootstrapQueue {
    Bootstrap* bootstrap;
    size_t rounds;
    size_t next;        // next round to hand out
    double*** counts;   // [round][pattern set], NULL for the observed round
    double* estimates;  // [round * state_size], the state at the optimum
    int*** categories;  // [round][tree likelihood], CAT assignments, NULL otherwise
    bool* done;         // [round]
    pthread_mutex_t lock;
    pthread_cond_t finished;
} _BootstrapQueue;

typedef struct _BootstrapWorker {
    _BootstrapQueue* queue;
    _BootstrapGraph* graph;
    pthread_t thread;
} _BootstrapWorker;

// Rounds are handed out and drawn under one lock, so replicate r is the r-th
// draw from the bootstrap's generator whichever worker takes it.
static void* _bootstrap_worker(void* data) {
    _BootstrapWorker* worker = data;
    _BootstrapQueue* queue = worker->queue;
    _BootstrapGraph* graph = worker->graph;
    Bootstrap* bootstrap = queue->bootstrap;

    while (true) {
        pthread_mutex_lock(&queue->lock);
        if (queue->next == queue->rounds) {
            pthread_mutex_unlock(&queue->lock);
            break;
        }
        size_t r = queue->next++;
        bool observed = bootstrap->include_observed && r == 0;
        if (!observed) {
            queue->counts[r] = malloc(sizeof(double*) * bootstrap->pattern_count);
            for (size_t p = 0; p < bootstrap->pattern_count; p++) {
                _bootstrap_draw(bootstrap, p);
                queue->counts[r][p] =
                    clone_dvector(bootstrap->counts[p], bootstrap->patterns[p]->count);
            }
        }
        pthread_mutex_unlock(&queue->lock);

        // The whole state, not only `parameters`: whatever a worker's previous
        // replicate left behind would otherwise depend on the scheduling.
        Parameters_restore_value(graph->state, bootstrap->state_estimate);
        for (size_t p = 0; p < bootstrap->pattern_count; p++) {
            _bootstrap_apply(bootstrap, graph->treelikelihoods, graph->views, p,
                             observed ? NULL : queue->counts[r][p]);
        }
        _bootstrap_optimize(graph->optimizers, bootstrap->optimizer_count);

        int** categories = _bootstrap_copy_categories(bootstrap, graph->treelikelihoods);
        pthread_mutex_lock(&queue->lock);
        Parameters_store_value(graph->state, queue->estimates + r * bootstrap->state_size);
        queue->categories[r] = categories;
        queue->done[r] = true;
        pthread_cond_broadcast(&queue->finished);
        pthread_mutex_unlock(&queue->lock);
    }
    return NULL;
}

// Fit the replicates on clones of the model graph and log each one, in order,
// from the original graph: its pattern weights and the fitted parameters are
// installed there before the loggers read it.
static void _bootstrap_run_threaded(Bootstrap* bootstrap, size_t rounds) {
    size_t thread_count = bootstrap->threads < rounds ? bootstrap->threads : rounds;

    _BootstrapQueue queue;
    queue.bootstrap = bootstrap;
    queue.rounds = rounds;
    queue.next = 0;
    queue.counts = calloc(rounds, sizeof(double**));
    queue.estimates = dvector(rounds * bootstrap->state_size);
    queue.categories = calloc(rounds, sizeof(int**));
    queue.done = calloc(rounds, sizeof(bool));
    pthread_mutex_init(&queue.lock, NULL);
    pthread_cond_init(&queue.finished, NULL);

    // Cloned here, on one thread, from the point estimate in place.
    _BootstrapWorker* workers = malloc(sizeof(_BootstrapWorker) * thread_count);
    for (size_t w = 0; w < thread_count; w++) {
        workers[w].queue = &queue;
        workers[w].graph = _new_BootstrapGraph(bootstrap);
    }
    for (size_t w = 0; w < thread_count; w++) {
        pthread_create(&workers[w].thread, NULL, _bootstrap_worker, &workers[w]);
    }

    for (size_t r = 0; r < rounds; r++) {
        pthread_mutex_lock(&queue.lock);
        while (!queue.done[r]) pthread_cond_wait(&queue.finished, &queue.lock);
        double** counts = queue.counts[r];
        int** categories = queue.categories[r];
        pthread_mutex_unlock(&queue.lock);

        for (size_t p = 0; p < bootstrap->pattern_count; p++) {
            _bootstrap_apply(bootstrap, bootstrap->treelikelihoods, bootstrap->views, p,
                             counts == NULL ? NULL : counts[p]);
        }
        Parameters_restore_value(bootstrap->state,
                                 queue.estimates + r * bootstrap->state_size);
        _bootstrap_install_categories(bootstrap, categories);
        for (size_t i = 0; i < bootstrap->logger_count; i++) {
            bootstrap->loggers[i]->write(bootstrap->loggers[i], r);
        }
        if (counts != NULL) {
            for (size_t p = 0; p < bootstrap->pattern_count; p++) free(counts[p]);
            free(counts);
        }
    }

    for (size_t w = 0; w < thread_count; w++) {
        pthread_join(workers[w].thread, NULL);
        _free_BootstrapGraph(bootstrap, workers[w].graph);
    }
    free(workers);
    pthread_mutex_destroy(&queue.lock);
    pthread_cond_destroy(&queue.finished);
    free(queue.counts);
    free(queue.estimates);
    free(queue.categories);
    free(queue.done);
}

#endif  // PTHREAD_ENABLED

static void _bootstrap_run(Bootstrap* bootstrap) {
    Parameters_store_value(bootstrap->parameters, bootstrap->estimate);
    Parameters_store_value(bootstrap->state, bootstrap->state_estimate);
    _bootstrap_store_assignments(bootstrap);

    for (size_t i = 0; i < bootstrap->logger_count; i++) {
        bootstrap->loggers[i]->initialize(bootstrap->loggers[i]);
    }

    // include_observed prepends a fit on the observed data; it is extra, so the
    // number of resampled replicates is `replicates` either way.
    size_t rounds = bootstrap->replicates + (bootstrap->include_observed ? 1 : 0);

#if defined(PTHREAD_ENABLED)
    if (bootstrap->threads > 1 && rounds > 1) {
        _bootstrap_run_threaded(bootstrap, rounds);
    } else {
        _bootstrap_run_serial(bootstrap, rounds);
    }
#else
    _bootstrap_run_serial(bootstrap, rounds);
#endif

    for (size_t i = 0; i < bootstrap->logger_count; i++) {
        bootstrap->loggers[i]->finalize(bootstrap->loggers[i]);
    }

    _bootstrap_restore_observed(bootstrap, bootstrap->treelikelihoods);
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
        free(bootstrap->site_categories[i]);
        bootstrap->treelikelihoods[i]->free(bootstrap->treelikelihoods[i]);
    }
    free(bootstrap->site_categories);
    free(bootstrap->site_patterns);
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
    free_Parameters(bootstrap->state);
    free(bootstrap->state_estimate);
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
    // Filled by _bootstrap_store_assignments when the run starts, which is the
    // first moment the point estimate, and so the assignment, is in place.
    bootstrap->site_categories = calloc(count, sizeof(int*));
    bootstrap->site_patterns = calloc(count, sizeof(SitePattern*));
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
        {"threads", JSON_OPTIONAL, JSON_NUMBER},
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
    bootstrap->hash = hash;
    bootstrap->optimizer_node = get_json_node(node, "optimizer");
    bootstrap->threads = get_json_node_value_size_t(node, "threads", 1);
    if (bootstrap->threads == 0) {
        json_die(node, "%s - \"threads\" must be at least 1", id);
    }
    if (bootstrap->threads > 1) {
        // Each worker starts its replicates from the point estimate: with
        // "reset": false a replicate starts where the previous one ended, which
        // is a sequence, not something to share out.
        if (!bootstrap->reset) {
            json_die(node, "%s - \"threads\" > 1 requires \"reset\": true", id);
        }
        // Without compaction a replicate rewrites the weights of the pattern
        // set itself, which every worker reads.
        if (!bootstrap->compact) {
            json_die(node, "%s - \"threads\" > 1 requires \"compact\": true", id);
        }
#if !defined(PTHREAD_ENABLED)
        fprintf(stderr, "%s - built without pthreads, \"threads\" is ignored\n", id);
        bootstrap->threads = 1;
#endif
    }

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

    bootstrap->state = new_Parameters(1);
    for (size_t i = 0; i < bootstrap->optimizer_count; i++) {
        opt_collect_parameters(bootstrap->optimizers[i], bootstrap->state);
    }
    for (size_t i = 0; i < bootstrap->treelikelihood_count; i++) {
        // tree, substitution, site and branch models (see new_TreeLikelihoodModel)
        Model** submodels = bootstrap->treelikelihoods[i]->data;
        // The node heights of a reparameterized time tree are computed from the
        // parameters of its transform, which are what is fitted.
        bool reparameterized = submodels[0]->data != NULL;  // its transform model
        for (size_t k = 0; k < 4; k++) {
            if (submodels[k] == NULL) continue;
            Parameters* ps = submodels[k]->parameters;
            for (size_t j = 0; j < Parameters_count(ps); j++) {
                Parameter* p = Parameters_at(ps, j);
                if (k == 0 && reparameterized && p->model == MODEL_TREE) continue;
                // A transformed parameter follows the one it is computed from. One
                // without a name is internal (the unused branch lengths of a time
                // tree) and cannot be matched to its clone.
                if (p->transform != NULL || !p->estimate) continue;
                if (Parameter_name(p)[0] == '\0') continue;
                if (!Parameters_contains(bootstrap->state, p)) {
                    Parameters_add(bootstrap->state, p);
                }
            }
        }
    }

    // What gets snapshotted at the point estimate and restored before each
    // replicate. Defaults to everything a replicate fit can move.
    bootstrap->parameters = new_Parameters(1);
    if (get_json_node(node, "parameters") != NULL) {
        get_parameters_references(node, hash, bootstrap->parameters);
        for (size_t i = 0; i < Parameters_count(bootstrap->parameters); i++) {
            Parameter* p = Parameters_at(bootstrap->parameters, i);
            if (!Parameters_contains(bootstrap->state, p)) Parameters_add(bootstrap->state, p);
        }
    } else {
        Parameters_add_parameters(bootstrap->parameters, bootstrap->state);
    }
    if (Parameters_count(bootstrap->parameters) == 0) {
        json_die(node,
                 "%s - no parameters to restart replicates from: neither the "
                 "optimizers nor the tree likelihoods have any, so \"parameters\" "
                 "is required",
                 id);
    }

    bootstrap->estimate_size = 0;
    for (size_t i = 0; i < Parameters_count(bootstrap->parameters); i++) {
        bootstrap->estimate_size +=
            Parameter_size(Parameters_at(bootstrap->parameters, i));
    }
    bootstrap->estimate = dvector(bootstrap->estimate_size);
    bootstrap->state_size = Parameters_size(bootstrap->state);
    bootstrap->state_estimate = dvector(bootstrap->state_size);

    bootstrap->run = _bootstrap_run;
    bootstrap->free = _free_Bootstrap;
    return bootstrap;
}
