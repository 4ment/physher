// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "minunit.h"

#include "phyc/filereader.h"
#include "phyc/hashtable.h"
#include "phyc/mjson.h"
#include "phyc/modelfactory.h"
#include "phyc/optimizer.h"
#include "phyc/parameters.h"

#define MAX_MODELS 8

static Hashtable* _new_hash(void) {
    Hashtable* hash = new_Hashtable_string(10);
    hashtable_set_key_ownership(hash, false);
    hashtable_set_value_ownership(hash, false);
    return hash;
}

static Parameter* _find(Parameters* ps, const char* name) {
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        if (strcmp(Parameters_name(ps, i), name) == 0) return Parameters_at(ps, i);
    }
    return NULL;
}

// Scale every untransformed parameter of `ps` by `factor`, within its bounds.
static void _scale(Parameters* ps, double factor) {
    for (size_t i = 0; i < Parameters_count(ps); i++) {
        Parameter* p = Parameters_at(ps, i);
        if (p->transform != NULL || p->simplex) continue;
        for (size_t j = 0; j < Parameter_size(p); j++) {
            double v = Parameter_value_at(p, j) * factor;
            if (v > 0 && v >= Parameter_lower(p) && v <= Parameter_upper(p)) {
                Parameter_set_value_at(p, v, j);
            }
        }
    }
}

static Optimizer* _new_optimizer(Model* model, Parameters* parameters) {
    Optimizer* opt = new_Optimizer(OPT_SERIAL_BRENT);
    opt_set_data(opt, model);
    opt_set_objective_function(opt, model_negative_logP);
    opt_set_parameters(opt, parameters);
    opt_set_max_iteration(opt, 2);
    opt_set_verbosity(opt, 0);
    return opt;
}

// Clone every model of `file` into a second hashtable (Model_clone_into) and
// check the clone is an equivalent and independent model graph: same logP, its
// own parameters registered under the same ids, unaffected by changes to the
// original and vice versa, and optimized to the same point by an optimizer built
// against its hashtable.
static char* _check_clone_json(const char* content) {
    Hashtable* hash = _new_hash();
    Hashtable* hash2 = _new_hash();

    json_node* json = create_json_tree((char*)content);

    Model* models[MAX_MODELS];
    Model* clones[MAX_MODELS];
    double logP[MAX_MODELS];
    size_t count = 0;
    for (size_t i = 0; i < json->child_count; i++) {
        json_node* child = json->children[i];
        if (strcasecmp(child->key, "physher") == 0) continue;
        mu_assert(count < MAX_MODELS, "too many models");
        models[count] = model_factory_from_json(child, hash);
        Hashtable_add(hash, models[count]->name, models[count]);
        count++;
    }

    // Last first: a model is then cloned before the models it reads from, whose
    // parameters it must end up sharing.
    for (size_t m = count; m > 0; m--) {
        clones[m - 1] = Model_clone_into(models[m - 1], hash2);
    }

    for (size_t m = 0; m < count; m++) {
        logP[m] = models[m]->full_logP(models[m]);
        mu_assert(clones[m] != models[m], "clone is the original");
        mu_assert(clones[m]->full_logP(clones[m]) == logP[m], "logP not matching");
        mu_assert(Hashtable_get(hash2, models[m]->name) == clones[m],
                  "clone not registered");

        Parameters* ps = models[m]->parameters;
        Parameters* qs = clones[m]->parameters;
        mu_assert(Parameters_count(ps) == Parameters_count(qs), "parameter count");
        for (size_t i = 0; i < Parameters_count(ps); i++) {
            Parameter* p = Parameters_at(ps, i);
            Parameter* q = _find(qs, Parameter_name(p));
            mu_assert(q != NULL && q != p, "parameter not cloned");
            mu_assert(memcmp(Parameter_values(p), Parameter_values(q),
                             sizeof(double) * Parameter_size(p)) == 0,
                      "parameter values not matching");
            if (Hashtable_exists(hash, Parameter_name(p))) {
                mu_assert(Hashtable_get(hash2, Parameter_name(p)) == q,
                          "parameter not registered");
            }
        }
    }

    // independence, both ways
    for (size_t m = 0; m < count; m++) {
        _scale(clones[m]->parameters, 1.1);
        double perturbed = clones[m]->logP(clones[m]);
        for (size_t k = 0; k < count; k++) {
            mu_assert(models[k]->logP(models[k]) == logP[k], "original moved");
        }
        _scale(models[m]->parameters, 1.1);
        mu_assert(fabs(models[m]->logP(models[m]) - perturbed) < 1e-9,
                  "perturbed logP not matching");
        _scale(clones[m]->parameters, 1.0 / 1.1);
        _scale(models[m]->parameters, 1.0 / 1.1);
        mu_assert(fabs(clones[m]->logP(clones[m]) - logP[m]) < 1e-9,
                  "clone logP not restored");
        mu_assert(fabs(models[m]->logP(models[m]) - logP[m]) < 1e-9,
                  "original logP not restored");
    }

    // The same optimization on both graphs, over every free scalar parameter of
    // the first model, must reach the same point.
    Parameters* free_parameters = new_Parameters(1);
    Parameters* free_parameters2 = new_Parameters(1);
    for (size_t i = 0; i < Parameters_count(models[0]->parameters); i++) {
        Parameter* p = Parameters_at(models[0]->parameters, i);
        if (p->transform != NULL || p->simplex || !p->estimate) continue;
        Parameters_add(free_parameters, p);
        Parameters_add(free_parameters2, _find(clones[0]->parameters, Parameter_name(p)));
    }
    Optimizer* opt = _new_optimizer(models[0], free_parameters);
    Optimizer* opt2 = _new_optimizer(clones[0], free_parameters2);
    double fmin;
    double fmin2;
    opt_optimize(opt, &fmin);
    for (size_t m = 0; m < count; m++) {
        mu_assert(fabs(clones[m]->logP(clones[m]) - logP[m]) < 1e-9,
                  "clone moved by the optimization of the original");
    }
    opt_optimize(opt2, &fmin2);
    for (size_t m = 0; m < count; m++) {
        mu_assert(fabs(models[m]->logP(models[m]) - clones[m]->logP(clones[m])) < 1e-8,
                  "optima not matching");
    }
    free_Optimizer(opt);
    free_Optimizer(opt2);
    free_Parameters(free_parameters);
    free_Parameters(free_parameters2);

    for (size_t m = 0; m < count; m++) clones[m]->free(clones[m]);
    for (size_t m = 0; m < count; m++) {
        mu_assert(isfinite(models[m]->logP(models[m])), "original freed with clone");
    }
    for (size_t m = count; m > 0; m--) models[m - 1]->free(models[m - 1]);
    free_Hashtable(hash);
    free_Hashtable(hash2);
    json_free_tree(json);
    return NULL;
}

static char* _check_clone(const char* file) {
    char* content = load_file(file);
    char* result = _check_clone_json(content);
    free(content);
    return result;
}

char* test_jc69_clone() { return _check_clone("jc69.json"); }

char* test_distance_clone() { return _check_clone("jc69-distance.json"); }

// two likelihoods of reparameterized time trees
char* test_time_clone() { return _check_clone("jc69-time.json"); }

// GTR+G with simplex rates and frequencies, and priors reading the parameters of
// the likelihood
char* test_gtr_priors_clone() {
    return _check_clone_json(
        "{"
        "\"model\": {"
        "  \"id\": \"treelikelihood\", \"type\": \"treelikelihood\","
        "  \"sitepattern\": {\"id\": \"patterns\", \"type\": \"sitepattern\","
        "    \"datatype\": \"nucleotide\","
        "    \"alignment\": {\"id\": \"seqs\", \"type\": \"alignment\","
        "      \"file\": \"tiny.fa\"}},"
        "  \"sitemodel\": {\"id\": \"sitemodel\", \"type\": \"sitemodel\","
        "    \"distribution\": \"gamma\", \"categories\": 4,"
        "    \"shape\": {\"id\": \"alpha\", \"type\": \"parameter\","
        "      \"x\": 0.5, \"lower\": 0}},"
        "  \"substitutionmodel\": {\"id\": \"sm\", \"type\": \"substitutionmodel\","
        "    \"kind\": \"gtr\", \"datatype\": \"nucleotide\","
        "    \"frequencies\": {\"id\": \"freqs\", \"type\": \"Simplex\","
        "      \"x\": {\"id\": \"freqs.unconstrained\", \"type\": \"parameter\","
        "        \"x\": [0.1, -0.2, 0.3]}},"
        "    \"rates\": {\"id\": \"rates\", \"type\": \"Simplex\","
        "      \"x\": {\"id\": \"rates.unconstrained\", \"type\": \"parameter\","
        "        \"x\": [0.2, 0.1, -0.1, 0.3, -0.3]}}},"
        "  \"tree\": {\"id\": \"tree\", \"type\": \"tree\","
        "    \"branch_lengths\": {\"id\": \"tree.distances\", \"type\": \"parameter\","
        "      \"dimension\": 17, \"x\": 0.1},"
        "    \"initializer\": {\"id\": \"nj\", \"algorithm\": \"nj\","
        "      \"model\": \"jc69\", \"sitepattern\": \"&patterns\"}}"
        "},"
        "\"prior\": {\"id\": \"prior\", \"type\": \"compound\", \"distributions\": ["
        "  {\"id\": \"priorrates\", \"type\": \"distribution\","
        "    \"distribution\": \"dirichlet\", \"x\": \"&rates\"},"
        "  {\"id\": \"prioralpha\", \"type\": \"distribution\","
        "    \"distribution\": \"exponential\", \"x\": \"&alpha\","
        "    \"parameters\": {\"lambda\": 1.0}},"
        "  {\"id\": \"priorblens\", \"type\": \"distribution\","
        "    \"distribution\": \"exponential\", \"x\": \"&tree.distances\","
        "    \"parameters\": {\"lambda\": 10.0}}]},"
        "\"joint\": {\"id\": \"joint\", \"type\": \"compound\","
        "  \"distributions\": [\"@treelikelihood\", \"@prior\"]}"
        "}");
}

char* all_tests() {
    mu_suite_start();
    mu_run_test(test_jc69_clone);
    mu_run_test(test_distance_clone);
    mu_run_test(test_time_clone);
    mu_run_test(test_gtr_priors_clone);
    return NULL;
}

RUN_TESTS(all_tests);
