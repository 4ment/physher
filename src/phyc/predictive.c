// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "predictive.h"

#include <math.h>

#include "matrix.h"
#include "mstring.h"
#include "patterntrace.h"
#include "statistics.h"
#include "treelikelihood.h"

// Log pointwise predictive density and the effective number of parameters, the
// two halves of the WAIC, from a trace of per-pattern log-likelihoods.
void _predictive_calculate(Predictive* predictive){
	SingleTreeLikelihood* tlk = predictive->model->obj;
	const int pattern_count = tlk->sp->count;
	const double* weights = tlk->sp->weights;

	size_t count = 0;
	double** trace = read_pattern_log_likelihoods(predictive->model,
	                                              predictive->filename,
	                                              predictive->burnin,
	                                              "predictive", &count);

	double lppd = 0;
	double pwaic = 0;
	for (int i = 0; i < pattern_count; i++) {
		const double* lnl = trace[i];
		double sum = lnl[0];
		for (size_t j = 1; j < count; j++) {
			sum = logaddexp(sum, lnl[j]);
		}
		lppd += (sum - log(count))*weights[i];
		pwaic += weights[i]*variance(lnl, count, mean(lnl, count));
	}
	printf("lppd: %f\n", lppd);
	printf("pwaic: %f\n", pwaic);

	free_dmatrix(trace, pattern_count);
}

void _free_predictive(Predictive* predictive){
	free(predictive->filename);
	free(predictive);
}


Predictive* new_Predictive_from_json(json_node* node, Hashtable* hash){
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

	Predictive* predictive = malloc(sizeof(Predictive));
	predictive->filename = String_clone(filename);
	predictive->model = model;
	predictive->burnin = get_json_node_value_size_t(node, "burnin", 0);
	predictive->calculate = _predictive_calculate;
	predictive->free = _free_predictive;
	return predictive;
}
