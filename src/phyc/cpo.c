// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "cpo.h"

#include <float.h>
#include <math.h>

#include "matrix.h"
#include "mstring.h"
#include "patterntrace.h"
#include "treelikelihood.h"

void _cpo_calculate(struct CPO* cpo){
	SingleTreeLikelihood* tlk = cpo->model->obj;
	const int pattern_count = tlk->sp->count;
	const double* weights = tlk->sp->weights;

	size_t count = 0;
	double** trace = read_pattern_log_likelihoods(cpo->model, cpo->filename,
	                                              cpo->burnin, "cpo", &count);

	double logCPO = 0;
	for (int i = 0; i < pattern_count; i++) {
		const double* lnl = trace[i];
		// The CPO of a pattern is the harmonic mean of its likelihood over the
		// sample; shifting by the smallest log-likelihood keeps exp() in range.
		double min = lnl[0];
		for (size_t j = 1; j < count; j++) {
			min = fmin(lnl[j], min);
		}
		double sum = -DBL_MAX;
		for (size_t j = 0; j < count; j++) {
			sum = logaddexp(sum, min-lnl[j]);
		}
		logCPO += (log(count) + min - sum)*weights[i];
	}
	printf("logCPO: %f\n", logCPO);

	free_dmatrix(trace, pattern_count);
}

void _free_cpo(CPO* cpo){
	free(cpo->filename);
	free(cpo);
}


CPO* new_CPO_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"burnin", JSON_OPTIONAL, JSON_NUMBER},
	    {"filename", JSON_REQUIRED, JSON_STRING},
	    {"model", JSON_REQUIRED, JSON_STRING},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	
	char* filename = get_json_node_value_string(node, "filename");
	char* ref = get_json_node_value_string(node, "model");
	// The likelihood the trace was logged from: it supplies the pattern weights
	// and the column names to look for, so the file needs to carry neither.
	Model* model = Hashtable_get(hash, ref + 1);
	if (model == NULL || model->type != MODEL_TREELIKELIHOOD) {
		json_die(node, "\"model\" must reference a tree likelihood (e.g. \"@treelikelihood\"): '%s'", ref);
	}

	CPO* cpo = malloc(sizeof(CPO));
	cpo->filename = String_clone(filename);
	cpo->model = model;
	cpo->burnin = get_json_node_value_size_t(node, "burnin", 0);
	cpo->calculate = _cpo_calculate;
	cpo->free = _free_cpo;
	return cpo;
}
