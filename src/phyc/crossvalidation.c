// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "crossvalidation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "cat.h"
#include "matrix.h"
#include "treelikelihood.h"

// splitmix64: a self-contained generator, so the split does not depend on how
// much of the global stream the actions before it happened to consume.
static uint64_t _cv_next(CrossValidation* cv) {
    uint64_t z = (cv->rng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// A uniform integer in [0, n), by rejection so the range is exact.
static size_t _cv_uniform(CrossValidation* cv, size_t n) {
    uint64_t limit = UINT64_MAX - (UINT64_MAX % n) - 1;
    uint64_t draw;
    do {
        draw = _cv_next(cv);
    } while (draw > limit);
    return (size_t)(draw % n);
}

// Deal every column of pattern set p across the folds, writing the held-out
// multiplicities into cv->test[p] as a folds x count matrix.
//
// The labels are built balanced (fold j gets every N/folds-th column) and then
// shuffled, so the folds are equal in *sites* to within one rather than equal in
// expectation. A pattern of multiplicity 300 therefore contributes about 300/F
// columns to each fold, and a singleton contributes its one column to exactly one
// fold and none to the rest -- which is the point: the singletons are the
// informative tail, and they have to be held out somewhere.
static void _cv_split(CrossValidation* cv, size_t p) {
    const SitePattern* sp = cv->patterns[p];
    const double* weights = cv->weights[p];
    double* test = cv->test[p];
    memset(test, 0, sizeof(double) * cv->folds * sp->count);

    size_t total = 0;
    for (int i = 0; i < sp->count; i++) total += (size_t)(weights[i] + 0.5);

    size_t* labels = malloc(sizeof(size_t) * total);
    for (size_t j = 0; j < total; j++) labels[j] = j % cv->folds;
    // Fisher-Yates, downward, so every permutation is equally likely.
    for (size_t j = total; j > 1; j--) {
        size_t k = _cv_uniform(cv, j);
        size_t swap = labels[j - 1];
        labels[j - 1] = labels[k];
        labels[k] = swap;
    }

    size_t pos = 0;
    for (int i = 0; i < sp->count; i++) {
        size_t n = (size_t)(weights[i] + 0.5);
        for (size_t j = 0; j < n; j++) {
            test[labels[pos] * sp->count + i] += 1.0;
            pos++;
        }
    }
    free(labels);

    // The whole criterion rests on every site being held out exactly once: if a
    // site were held out twice it would be scored twice, and if it were held out
    // never the elpd would silently be over a shorter alignment than the one
    // reported. Both are cheap to rule out and neither is visible in the result.
    for (int i = 0; i < sp->count; i++) {
        double dealt = 0;
        for (size_t f = 0; f < cv->folds; f++) dealt += test[f * sp->count + i];
        if (fabs(dealt - weights[i]) > 1.e-9) {
            fprintf(stderr,
                    "crossvalidation: pattern %d has weight %g but %g sites were "
                    "dealt to folds\n",
                    i, weights[i], dealt);
            exit(1);
        }
    }
}

// Install the multiplicities of pattern set p and invalidate everything reading
// it. sp->nsites is part of the installation, not bookkeeping: the empirical CAT
// site model divides its category rates by their mean over sp->weights *and*
// sp->nsites, so leaving nsites at the observed total while the weights sum to
// less would put the mean rate somewhere other than one and rescale the tree.
static void _cv_install(CrossValidation* cv, size_t p, const double* multiplicities) {
    SitePattern* sp = cv->patterns[p];
    double total = 0;
    for (int i = 0; i < sp->count; i++) total += multiplicities[i];
    memcpy(sp->weights, multiplicities, sizeof(double) * sp->count);
    sp->nsites = (int)(total + 0.5);

    for (size_t i = 0; i < cv->treelikelihood_count; i++) {
        if (cv->pattern_index[i] != p) continue;
        SingleTreeLikelihood* tlk = cv->treelikelihoods[i]->obj;
        SingleTreeLikelihood_update_weights(tlk);
        // The CAT normaliser is a function of the weights, and nothing else marks
        // it stale when they move underneath it.
        if (tlk->sm != NULL && tlk->sm->site_category != NULL) {
            tlk->sm->need_update = true;
        }
    }
}

// Per-pattern log predictive density of tree likelihood `t`, taken at the fit the
// training fold left behind and with the training weights still installed.
static void _cv_predictive(CrossValidation* cv, size_t t) {
    SingleTreeLikelihood* tlk = cv->treelikelihoods[t]->obj;
    double* out = cv->logpred[t];

    // The optimizer leaves the parameters at the optimum but the partials --
    // and so pattern_lk -- at whichever point its last line-search evaluated,
    // and a clean model returns its cached lk without refilling the array. The
    // per-pattern densities are the whole output here, so force the traversal.
    SingleTreeLikelihood_update_all_nodes(tlk);

    if (tlk->sm != NULL && tlk->sm->site_category != NULL) {
        // K + 2 traversals, and it puts the assignment, the rates and the
        // partials back as it found them.
        cat_mixture_patterns(tlk, 0, out);
        return;
    }
    // Everything else already reports a per-pattern density marginal over its
    // rate categories; tlk->lk is this array dotted with the weights.
    tlk->calculate(tlk);
    memcpy(out, tlk->pattern_lk, sizeof(double) * tlk->sp->count);
}

static void _cv_run(CrossValidation* cv) {
    Parameters_store_value(cv->parameters, cv->estimate);

    for (size_t p = 0; p < cv->pattern_count; p++) _cv_split(cv, p);
    for (size_t t = 0; t < cv->treelikelihood_count; t++) cv->elpd[t] = 0;

    FILE* file = NULL;
    if (cv->filename != NULL) {
        file = fopen(cv->filename, "w");
        if (file == NULL) {
            fprintf(stderr, "crossvalidation: cannot write '%s'\n", cv->filename);
            exit(1);
        }
        fprintf(file, "fold\tmodel\tpattern\ttest_weight\tlog_predictive\n");
    }

    for (size_t f = 0; f < cv->folds; f++) {
        Parameters_restore_value(cv->parameters, cv->estimate);

        for (size_t p = 0; p < cv->pattern_count; p++) {
            const SitePattern* sp = cv->patterns[p];
            const double* test = cv->test[p] + f * sp->count;
            for (int i = 0; i < sp->count; i++) {
                cv->train[p][i] = cv->weights[p][i] - test[i];
            }
            _cv_install(cv, p, cv->train[p]);
        }

        for (size_t i = 0; i < cv->optimizer_count; i++) {
            double logP;
            opt_optimize(cv->optimizers[i], &logP);
        }

        // Scored before anything touches the weights again: the held-out
        // multiplicities are a reduction vector here, never an installed state.
        for (size_t t = 0; t < cv->treelikelihood_count; t++) {
            _cv_predictive(cv, t);
            const SitePattern* sp = cv->patterns[cv->pattern_index[t]];
            const double* test = cv->test[cv->pattern_index[t]] + f * sp->count;
            double fold_elpd = 0;
            for (int i = 0; i < sp->count; i++) {
                fold_elpd += test[i] * cv->logpred[t][i];
                if (file != NULL && test[i] > 0) {
                    fprintf(file, "%zu\t%s\t%d\t%g\t%.10e\n", f,
                            cv->treelikelihoods[t]->name, i, test[i],
                            cv->logpred[t][i]);
                }
            }
            cv->elpd[t] += fold_elpd;
            if (cv->verbosity > 0) {
                printf(
                    "Cross-validation fold %zu/%zu: %s held-out log-likelihood "
                    "%f\n",
                    f + 1, cv->folds, cv->treelikelihoods[t]->name, fold_elpd);
            }
        }
    }

    if (file != NULL) fclose(file);

    // Put the observed alignment back, so anything running after this sees the
    // real data.
    for (size_t p = 0; p < cv->pattern_count; p++) {
        _cv_install(cv, p, cv->weights[p]);
        cv->patterns[p]->nsites = cv->observed_nsites[p];
    }
    Parameters_restore_value(cv->parameters, cv->estimate);

    double total = 0;
    double sites = 0;
    for (size_t p = 0; p < cv->pattern_count; p++) {
        for (int i = 0; i < cv->patterns[p]->count; i++) sites += cv->weights[p][i];
    }
    // Reported as "key: value" so the same log scrapers that read every other
    // scalar physher prints pick these up too.
    printf("Cross-validation over %zu folds, %.0f sites\n", cv->folds, sites);
    for (size_t t = 0; t < cv->treelikelihood_count; t++) {
        total += cv->elpd[t];
        printf("%s.elpd: %f\n", cv->treelikelihoods[t]->name, cv->elpd[t]);
    }
    if (cv->treelikelihood_count > 1) {
        printf("crossvalidation.elpd: %f\n", total);
    }
    // Every site is held out exactly once, so the denominator is the alignment
    // length and this number is comparable across models and across alignments.
    printf("crossvalidation.elpd_per_site: %f\n", total / sites);
}

static void _free_CrossValidation(CrossValidation* cv) {
    for (size_t p = 0; p < cv->pattern_count; p++) {
        free(cv->weights[p]);
        free(cv->train[p]);
        free(cv->test[p]);
    }
    free(cv->weights);
    free(cv->train);
    free(cv->test);
    free(cv->observed_nsites);
    free(cv->patterns);

    for (size_t t = 0; t < cv->treelikelihood_count; t++) {
        cv->treelikelihoods[t]->free(cv->treelikelihoods[t]);
        free(cv->logpred[t]);
    }
    free(cv->treelikelihoods);
    free(cv->logpred);
    free(cv->elpd);
    free(cv->pattern_index);

    for (size_t i = 0; i < cv->optimizer_count; i++) {
        free_Optimizer(cv->optimizers[i]);
    }
    free(cv->optimizers);

    free_Parameters(cv->parameters);
    free(cv->estimate);
    free(cv->filename);
    free(cv);
}

// Collect the tree likelihoods named by "treelikelihood" (a ref or an array of
// refs) and the distinct SitePatterns they read.
static void _cv_collect_models(CrossValidation* cv, json_node* node, Hashtable* hash,
                               const char* id) {
    json_node* tlk_node = get_json_node(node, "treelikelihood");
    size_t count = tlk_node->node_type == MJSON_ARRAY ? tlk_node->child_count : 1;

    cv->treelikelihoods = malloc(sizeof(Model*) * count);
    cv->pattern_index = malloc(sizeof(size_t) * count);
    cv->logpred = malloc(sizeof(double*) * count);
    cv->elpd = dvector(count);
    cv->treelikelihood_count = count;
    cv->patterns = malloc(sizeof(SitePattern*) * count);
    cv->pattern_count = 0;

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
                     "%s - \"treelikelihood\" must reference a tree likelihood, "
                     "got a %s ('%s')",
                     id, model_type_strings[model->type], (char*)child->value);
        }
        model->ref_count++;
        cv->treelikelihoods[i] = model;

        SingleTreeLikelihood* tlk = model->obj;
        cv->logpred[i] = dvector(tlk->sp->count);

        // Two partitions can share one alignment: split it once and install it
        // into both, or the two would disagree on the fold.
        SitePattern* sp = tlk->sp;
        size_t p = 0;
        while (p < cv->pattern_count && cv->patterns[p] != sp) p++;
        if (p == cv->pattern_count) {
            cv->patterns[p] = sp;
            cv->pattern_count++;
        }
        cv->pattern_index[i] = p;
    }
}

static void _cv_allocate_patterns(CrossValidation* cv, json_node* node,
                                  const char* id) {
    size_t n = cv->pattern_count;
    cv->weights = malloc(sizeof(double*) * n);
    cv->train = malloc(sizeof(double*) * n);
    cv->test = malloc(sizeof(double*) * n);
    cv->observed_nsites = malloc(sizeof(int) * n);

    for (size_t p = 0; p < n; p++) {
        const SitePattern* sp = cv->patterns[p];
        cv->weights[p] = clone_dvector(sp->weights, sp->count);
        cv->train[p] = dvector(sp->count);
        cv->test[p] = dvector((size_t)sp->count * cv->folds);
        cv->observed_nsites[p] = sp->nsites;

        // The weights are integer site counts stored as doubles; sum them rather
        // than trusting sp->nsites, which some pattern constructors set from the
        // alignment length instead.
        double total = 0;
        for (int i = 0; i < sp->count; i++) total += sp->weights[i];
        if (total < (double)cv->folds) {
            json_die(node,
                     "%s - \"folds\" is %zu but the alignment has only %.0f "
                     "sites to deal out",
                     id, cv->folds, total);
        }
    }
}

CrossValidation* new_CrossValidation_from_json(json_node* node, Hashtable* hash) {
    static const json_field schema[] = {
        {"file", JSON_OPTIONAL, JSON_STRING},
        {"folds", JSON_OPTIONAL, JSON_NUMBER},
        {"optimizer", JSON_REQUIRED, JSON_OBJECT | JSON_ARRAY},
        {"parameters", JSON_OPTIONAL, JSON_ANY},
        {"seed", JSON_OPTIONAL, JSON_NUMBER},
        {"treelikelihood", JSON_REQUIRED, JSON_STRING | JSON_ARRAY},
        {"verbosity", JSON_OPTIONAL, JSON_NUMBER},
    };
    json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));

    const char* id = get_json_node_value_string(node, "id");
    CrossValidation* cv = malloc(sizeof(CrossValidation));

    cv->folds = get_json_node_value_size_t(node, "folds", 10);
    if (cv->folds < 2) {
        json_die(node, "%s - \"folds\" must be at least 2", id);
    }
    cv->verbosity = get_json_node_value_int(node, "verbosity", 0);
    const char* filename = get_json_node_value_string(node, "file");
    cv->filename = filename == NULL ? NULL : String_clone(filename);
    cv->rng_state = (uint64_t)get_json_node_value_size_t(node, "seed", 1);

    _cv_collect_models(cv, node, hash, id);
    _cv_allocate_patterns(cv, node, id);

    json_node* opt_node = get_json_node(node, "optimizer");
    cv->optimizer_count =
        opt_node->node_type == MJSON_ARRAY ? opt_node->child_count : 1;
    cv->optimizers = malloc(sizeof(Optimizer*) * cv->optimizer_count);
    for (size_t i = 0; i < cv->optimizer_count; i++) {
        json_node* child =
            opt_node->node_type == MJSON_ARRAY ? opt_node->children[i] : opt_node;
        cv->optimizers[i] = new_Optimizer_from_json(child, hash);
    }

    // What gets snapshotted and restored before each fold, so no fold starts from
    // the fit of the one before it. Defaults to everything the optimizers move.
    cv->parameters = new_Parameters(1);
    if (get_json_node(node, "parameters") != NULL) {
        get_parameters_references(node, hash, cv->parameters);
    } else {
        for (size_t i = 0; i < cv->optimizer_count; i++) {
            Parameters* parameters = opt_parameters(cv->optimizers[i]);
            for (size_t j = 0; parameters != NULL && j < Parameters_count(parameters);
                 j++) {
                Parameter* parameter = Parameters_at(parameters, j);
                bool seen = false;
                for (size_t k = 0; k < Parameters_count(cv->parameters); k++) {
                    seen |= Parameters_at(cv->parameters, k) == parameter;
                }
                if (!seen) Parameters_add(cv->parameters, parameter);
            }
        }
    }
    if (Parameters_count(cv->parameters) == 0) {
        json_die(node,
                 "%s - no parameters to restart folds from: the optimizers "
                 "declare none, so \"parameters\" is required",
                 id);
    }

    cv->estimate_size = 0;
    for (size_t i = 0; i < Parameters_count(cv->parameters); i++) {
        cv->estimate_size += Parameter_size(Parameters_at(cv->parameters, i));
    }
    cv->estimate = dvector(cv->estimate_size);

    cv->run = _cv_run;
    cv->free = _free_CrossValidation;
    return cv;
}
