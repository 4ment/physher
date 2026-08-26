// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "marginal.h"

#include <float.h>
#include <string.h>
#include <strings.h>

#include "utils.h"
#include "utilsio.h"
#include "matrix.h"
#include "descriptivestats.h"
#include "statistics.h"
#include "mstring.h"
#include "beta.h"

// A trace with no usable samples silently turns every estimator into nan, and
// the usual cause is a burnin larger than the number of rows in the file --
// burnin counts logged samples, not MCMC iterations.
static void check_samples(const char* filename, size_t size, size_t burnin){
	if(size == 0){
		fprintf(stderr, "No sample left in %s after discarding %zu of them as burnin.\n", filename, burnin);
		fprintf(stderr, "`burnin` counts logged samples (rows), not MCMC iterations.\n\n");
		exit(1);
	}
}

double log_arithmetic_mean(const Vector* vecvalues){
	double sum = -DBL_MAX;
	const double* values = Vector_data(vecvalues);
	size_t size = Vector_length(vecvalues);
	for(size_t i = 0; i < size; i++){
		sum = logaddexp(sum, values[i]);
	}
	return sum - log(size);
}

double log_harmonic_mean(const Vector* vecvalues){
	double sum = 0;
	const double* values = Vector_data(vecvalues);
	size_t size = Vector_length(vecvalues);
	for(size_t i = 0; i < size; i++){
		sum += values[i];
	}
	
	double denominator = -DBL_MAX;
	
	for(size_t i = 0; i < size; i++){
		denominator = logaddexp(denominator, sum - values[i]);
	}
	return sum - denominator + log(size);
}

double log_smoothed_harmonic_mean(double logP, const Vector* values, double delta){
	double ldelta = log(delta);
	double l1_delta = log(1.0-delta);
	
	size_t n = Vector_length(values);
	double num = log(n) + ldelta - l1_delta + logP;
	double denom = log(n) + ldelta - l1_delta;
	
	for(size_t i = 0; i < n; i++){
		double value = Vector_at(values, i);
		double norm = -logaddexp(ldelta, l1_delta + value - logP);
		num = logaddexp(num, norm + value);
		denom = logaddexp(denom, norm);
	}
	return num - denom;
}

double log_stablilized_harmonic_mean(double guess, const Vector* values, double delta){
	double logP = guess;
	double logPp = 10;
	for(size_t i = 0; i < 10000; i++){
		logP = log_smoothed_harmonic_mean(logP, values, delta);
		if(fabs(logP - logPp) < 1.e-7) break;
		logPp = logP;
	}
	return logP;
}

// temperatures should be sorted in increasing order
double log_marginal_stepping_stone(const Vector** values, size_t temp_count, const double* temperatures, double* lrssk){
	double lrss = 0;
	for(size_t i = 1; i < temp_count; i++){
		double tempdiff = temperatures[i]-temperatures[i-1];
		const double* previousll = Vector_data(values[i-1]);
		size_t size = Vector_length(values[i-1]);
		double logmaxll = dmax_vector(previousll, size);
		double temp = 0;
		for (size_t j = 0; j < size; j++) {
			temp += exp(tempdiff * (previousll[j] - logmaxll));
		}
		lrssk[i-1] = tempdiff*logmaxll + log(temp/size);
		lrss += lrssk[i-1];
	}
	return lrss;
}

// temperatures should be sorted in increasing order
double log_marginal_path_sampling(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk){
	double lrps = 0;
	double* means = dvector(temp_count);
	for(size_t i = 0; i < temp_count; i++){
		means[i] = dmean(Vector_data(values[i]), Vector_length(values[i]));
	}
	
	for(size_t i = 1; i < temp_count; i++){
		double weight = temperatures[i]-temperatures[i-1];
		lrpsk[i-1] = weight * (means[i]+means[i-1])/2.0;
		lrps += lrpsk[i-1];
	}
	free(means);
	return lrps;
}

double log_marginal_path_sampling_modified(const Vector** values, size_t temp_count, const double* temperatures, double* lrpsk){
	double lrps = 0;
	double* means = dvector(temp_count);
	double* variances = dvector(temp_count);
	for(size_t i = 0; i < temp_count; i++){
		means[i] = dmean(Vector_data(values[i]), Vector_length(values[i]));
		variances[i] = variance(Vector_data(values[i]), Vector_length(values[i]), means[i]);
	}
	
	for(size_t i = 1; i < temp_count; i++){
		double weight = temperatures[i]-temperatures[i-1];
		lrpsk[i-1] = weight * (means[i]+means[i-1])/2.0 - ((weight*weight)/12)*(variances[i]-variances[i-1]);
		lrps += lrpsk[i-1];
	}
	free(means);
	free(variances);
	return lrps;
}

// temperatures should be sorted in increasing order
double log_bf_stepping_stone(const Vector** values, const Vector** values1, size_t temp_count, const double* temperatures, double* lrssk){
	size_t maxsize = 0;
	for(size_t i = 0; i < temp_count; i++){
		maxsize = fmax(maxsize, Vector_length(values[i]));
	}
	double* diffll = dvector(maxsize);
	double lrss = 0;
	for(size_t i = 1; i < temp_count; i++){
		size_t size = Vector_length(values[i-1]);
		double tempdiff = temperatures[i]-temperatures[i-1];
		const double* v = Vector_data(values[i-1]);
		const double* v1 = Vector_data(values1[i-1]);
		double logmaxll = -INFINITY;
		for (size_t j = 0; j < size; j++) {
			diffll[j] = v1[j] - v[j];
			logmaxll = fmax(logmaxll, diffll[j]);
		}
		double temp = 0;
		for (size_t j = 0; j < size; j++) {
			temp += exp(tempdiff * (diffll[j] - logmaxll));
		}
		lrssk[i-1] = tempdiff*logmaxll + log(temp/size);
		lrss += lrssk[i-1];
	}
	free(diffll);
	return lrss;
}

// temperatures should be sorted in increasing order
double log_bf_path_sampling(const Vector** values, const Vector** values1, size_t temp_count, const double* temperatures, double* lrpsk){
	double lrps = 0;
	double* means = dvector(temp_count);
	for(size_t i = 0; i < temp_count; i++){
		size_t size = Vector_length(values[i]);
		const double* v = Vector_data(values[i]);
		const double* v1 = Vector_data(values1[i]);
		double mean = 0;
		for (size_t j = 0; j < size; j++) {
			mean += v1[j] - v[j];
		}
		means[i] = mean/size;
	}
	
	for(size_t i = 1; i < temp_count; i++){
		double weight = temperatures[i]-temperatures[i-1];
		lrpsk[i-1] = weight * (means[i]+means[i-1])/2.0;
		lrps += lrpsk[i-1];
	}
	free(means);
	return lrps;
}

// code from BEAST
double ess(double* values, size_t samples, size_t stepSize) {
	
	size_t maxLag = samples - 1;//dmin(samples - 1, MAX_LAG);
	
	double* gammaStat = dvector(maxLag);
	double varStat = 0.0;
	double mean = dmean(values, samples);
	
	for (size_t lag = 0; lag < maxLag; lag++) {
		for (size_t j = 0; j < samples - lag; j++) {
			double del1 = values[j] - mean;
			double del2 = values[j + lag] - mean;
			gammaStat[lag] += (del1 * del2);
		}
		
		gammaStat[lag] /= ((double) (samples - lag));
		
		if (lag == 0) {
			varStat = gammaStat[0];
		} else if (lag % 2 == 0) {
			// fancy stopping criterion :)
			if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
				varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
			}
			// stop
			else
				maxLag = lag;
		}
	}
	
	double ACT;
	double ESS;
	
	// auto correlation time
	if (gammaStat[0]==0)
		ACT = 0;
	else
		ACT = stepSize * varStat / gammaStat[0];
	
	// effective sample size
	if (ACT==0)
		ESS=1;
	else
		ESS = (stepSize * samples) / ACT;
	
	free(gammaStat);
	return ESS;
}

void test_marg(const Vector** values, size_t temp_count, const double* temperatures){
	const Vector** lls = malloc(temp_count*sizeof(Vector*));
	double* temps = dvector(temp_count);
	double* means = dvector(temp_count);
	for (size_t i = 3; i < temp_count; i*=2) {
		double incr = ((double)temp_count-0.5)/(i-1);
		double dindex = 0;
		for (size_t j = 0; j < i; j++) {
			int index = dindex;
			lls[j] = values[index];
			temps[j] = temperatures[index];
			dindex += incr;
		}
		double* lrssk = malloc(temp_count*sizeof(double));
		
		double lrss = log_marginal_stepping_stone(lls, i, temps, lrssk);
		printf("%zu Stepping stone marginal likelihood: %f\n", i, lrss);
		//		for (int i = 0; i < temp_count-1; i++) {
		//			printf("%f: %f\n", temperatures[i], lrssk[i]);
		//		}
		
		double lrps = log_marginal_path_sampling(lls, i, temps, lrssk);
		printf("%zu Path sampling marginal likelihood: %f\n", i, lrps);
		//		for (int i = 0; i < temp_count-1; i++) {
		//			printf("%f: %f\n", temperatures[i], lrssk[i]);
		//		}
		lrps = log_marginal_path_sampling_modified(lls, i, temps, lrssk);
		printf("%zu Modified Path sampling marginal likelihood: %f\n", i, lrps);
		//		for (int i = 0; i < temp_count-1; i++) {
		//			printf("%f: %f\n", temperatures[i], lrssk[i]);
		//		}
		free(lrssk);
		
		for (size_t j = 0; j < i; j++) {
			means[j] = dmean(Vector_data(lls[j]), Vector_length(lls[j]));
		}
		double S = 0;
		for (size_t j = 0; j < i-1; j++) {
			S += (temps[j+1] - temps[j])*(means[j+1] - means[j]);
			printf("  %f\n", (temps[j+1] - temps[j])*(means[j+1] - means[j]));
		}
		printf("S: %f\n", S);
	}
	free(means);
	free(temps);
	free(lls);
}

void test_marg_gss(const Vector** values, size_t temp_count, const double* temperatures){
	const Vector** lls = malloc(temp_count*sizeof(Vector*));
	double* temps = dvector(temp_count);
	for (size_t i = 3; i < temp_count; i*=2) {
		double incr = ((double)temp_count-0.5)/(i-1);
		double dindex = 0;
		for (size_t j = 0; j < i; j++) {
			int index = dindex;
			lls[j] = values[index];
			temps[j] = temperatures[index];
			dindex += incr;
		}
		double* lrssk = malloc(temp_count*sizeof(double));
		
		double lrss = log_marginal_stepping_stone(lls, i, temps, lrssk);
		printf("%zu Generalized stepping stone marginal likelihood: %f\n", i, lrss);
		free(lrssk);
		
	}
	free(temps);
	free(lls);
}

static void _bayes_factor_run(MarginaLikelihood* margl){
	StringBuffer* buffer = new_StringBuffer(10);
	const Vector** lls = malloc(margl->temperature_count*sizeof(Vector*));
	const Vector** lls1 = malloc(margl->temperature_count*sizeof(Vector*));
	double* temperatures = malloc(margl->temperature_count*sizeof(double));
	
	int start_file = strlen(margl->file)-1;
	for (; start_file >= 0; start_file--) {
		if (margl->file[start_file] == '/') break;
	}
	
	// temperatures should be in decreasing order
	for (size_t i = 0; i < margl->temperature_count; i++) {
		// saved in increasing order
		StringBuffer_empty(buffer);
		if (start_file >= 0) {
			StringBuffer_append_substring(buffer, margl->file, start_file+1);
			StringBuffer_append_format(buffer, "%zu%s",i, margl->file+start_file+1);

		}
		else{
			StringBuffer_append_format(buffer, "%zu%s",i, margl->file);
		}
		printf("Temperature: %f - %s\n", margl->temperatures[i],buffer->c);
		lls[margl->temperature_count-i-1] = read_log_column_with_id(buffer->c, margl->burnin, margl->likelihood_tag[0]);
		lls1[margl->temperature_count-i-1] = read_log_column_with_id(buffer->c, margl->burnin, margl->likelihood_tag[1]);
		check_samples(buffer->c, Vector_length(lls[margl->temperature_count-i-1]), margl->burnin);
		temperatures[margl->temperature_count-i-1] = margl->temperatures[i];
	}

	double* lrssk = malloc(margl->temperature_count*sizeof(double));
	
	if(margl->ss){
		double lrss = log_bf_stepping_stone(lls, lls1, margl->temperature_count, temperatures, lrssk);
		printf("Stepping stone Bayes factor: %f\n", lrss);
		for (size_t i = 0; i < margl->temperature_count-1; i++) {
			printf("%f: %f\n", temperatures[i], lrssk[i]);
		}
	}
	if (margl->ps) {
		double lrps = log_bf_path_sampling(lls, lls1, margl->temperature_count, temperatures, lrssk);
		printf("Path sampling Bayes factor: %f\n", lrps);
		for (size_t i = 0; i < margl->temperature_count-1; i++) {
			printf("%f: %f\n", temperatures[i], lrssk[i]);
		}
	}
	
	free(lrssk);
	
	free_StringBuffer(buffer);
	free(temperatures);
	for (size_t i = 0; i < margl->temperature_count; i++) {
		free_Vector((Vector*)lls[i]);
		free_Vector((Vector*)lls1[i]);
	}
	free(lls);
	free(lls1);
}

static void _marginal_likelihood_run(MarginaLikelihood* margl){
	StringBuffer* buffer = new_StringBuffer(10);
	Vector** lls = malloc(margl->temperature_count*sizeof(Vector*));
	double* temperatures = malloc(margl->temperature_count*sizeof(double));
	
	if(margl->temperature_count == 1){
		lls[0] = read_log_column_with_id(margl->file, margl->burnin, margl->likelihood_tag[0]);
		check_samples(margl->file, Vector_length(lls[0]), margl->burnin);
	}
	// Generalized stepping stone
	else if(margl->refdist_tag != NULL){
		double* lrssk = malloc(margl->temperature_count*sizeof(double));
		const char* ids[2];
		ids[0] = margl->likelihood_tag[0];
		ids[1] = margl->refdist_tag;
		
		int start_file = strlen(margl->file)-1;
		for (; start_file >= 0; start_file--) {
			if (margl->file[start_file] == '/') break;
		}
		
		for (size_t i = 0; i < margl->temperature_count; i++) {
			// saved in increasing order
			StringBuffer_empty(buffer);
			if (start_file >= 0) {
				StringBuffer_append_substring(buffer, margl->file, start_file+1);
				StringBuffer_append_format(buffer, "%zu%s",i, margl->file+start_file+1);
    
			}
			else{
				StringBuffer_append_format(buffer, "%zu%s",i, margl->file);
			}
			temperatures[margl->temperature_count-i-1] = margl->temperatures[i];
			
			// The posterior (beta==1) chain is not run for GSS, only for GPS, so
			// its trace may be missing. Stepping stone never looks at that entry
			// (it only reads values[0..count-2]) but path sampling does, so leave
			// it NULL and let the GPS block below skip itself.
			if (!file_exists(buffer->c)) {
				printf("Temperature: %f - %s (not run)\n", margl->temperatures[i], buffer->c);
				lls[margl->temperature_count-i-1] = NULL;
				continue;
			}
			
			printf("Temperature: %f - %s\n", margl->temperatures[i], buffer->c);
			Vector** all = read_log_column_with_ids(buffer->c, margl->burnin, ids, 2);
			size_t size = Vector_length(all[0]);
			check_samples(buffer->c, size, margl->burnin);
			lls[margl->temperature_count-i-1] = new_Vector(size);
			for (size_t j = 0; j < size; j++) {
				Vector_push(lls[margl->temperature_count-i-1], Vector_at(all[0], j) - Vector_at(all[1], j));
			}
			free_Vector(all[0]);
			free_Vector(all[1]);
			free(all);
		}
		
		double lrss = log_marginal_stepping_stone((const Vector**)lls, margl->temperature_count, temperatures, lrssk);
		printf("Generalized stepping stone marginal likelihood: %f\n", lrss);
		for (size_t i = 0; i < margl->temperature_count-1; i++) {
			printf("%f: %f\n", temperatures[i], lrssk[i]);
		}
		
		if (lls[margl->temperature_count-1] != NULL) {
			double lrps = log_marginal_path_sampling((const Vector**)lls, margl->temperature_count, temperatures, lrssk);
			printf("Generalized path sampling marginal likelihood: %f\n", lrps);
			for (size_t i = 0; i < margl->temperature_count-1; i++) {
				printf("%f: %f\n", temperatures[i], lrssk[i]);
			}
			
			lrps = log_marginal_path_sampling_modified((const Vector**)lls, margl->temperature_count, temperatures, lrssk);
			printf("Modified Path sampling marginal likelihood: %f\n", lrps);
			for (size_t i = 0; i < margl->temperature_count-1; i++) {
				printf("%f: %f\n", temperatures[i], lrssk[i]);
			}
		}
		free(lrssk);
	}
	else{
		int start_file = strlen(margl->file)-1;
		for (; start_file >= 0; start_file--) {
			if (margl->file[start_file] == '/') break;
		}
		
		// temperatures should be in decreasing order
		for (size_t i = 0; i < margl->temperature_count; i++) {
			// saved in increasing order
			StringBuffer_empty(buffer);
			if (start_file >= 0) {
				StringBuffer_append_substring(buffer, margl->file, start_file+1);
				StringBuffer_append_format(buffer, "%zu%s",i, margl->file+start_file+1);
    
			}
			else{
				StringBuffer_append_format(buffer, "%zu%s",i, margl->file);
			}
			printf("Temperature: %f - %s\n", margl->temperatures[i],buffer->c);
			lls[margl->temperature_count-i-1] = read_log_column_with_id(buffer->c, margl->burnin, margl->likelihood_tag[0]);
			check_samples(buffer->c, Vector_length(lls[margl->temperature_count-i-1]), margl->burnin);
//			double* data = Vector_data(lls[margl->temperature_count-i-1]);
//			double n = Vector_length(lls[margl->temperature_count-i-1]);
//			double esss = ess(data, n, n);
//			printf("ESS: %f %f\n", esss, n);
			temperatures[margl->temperature_count-i-1] = margl->temperatures[i];
		}
		
		// sample from prior (temperature == 0)
//		double logNaive = log_arithmetic_mean(lls[0]);
//		printf("Naive arithmetic mean: %f\n", logNaive);
		
	}
	if(margl->temperature_count == 1 || margl->refdist_tag == NULL){
//		double logAM = log_arithmetic_mean(lls[margl->temperature_count-1]);
//		printf("Arithmetic mean: %f\n", logAM);
		
		if (margl->hm) {
			double logHM = log_harmonic_mean(lls[margl->temperature_count-1]);
			printf("Harmonic mean: %f\n", logHM);
		}
		
		if (margl->shm) {
			double logHM = log_harmonic_mean(lls[margl->temperature_count-1]);
			printf("Harmonic mean: %f\n", logHM);
			double logSHM = log_stablilized_harmonic_mean(logHM, lls[margl->temperature_count-1], 0.01);
			printf("Smoothed harmonic mean: %f\n", logSHM);
		}
		
		if(margl->temperature_count > 1){
			double* lrssk = malloc(margl->temperature_count*sizeof(double));
			
			if(margl->ss){
				double lrss = log_marginal_stepping_stone((const Vector**)lls, margl->temperature_count, temperatures, lrssk);
				printf("Stepping stone marginal likelihood: %f\n", lrss);
				for (size_t i = 0; i < margl->temperature_count-1; i++) {
					printf("%f: %f\n", temperatures[i], lrssk[i]);
				}
			}
			
			if(margl->ps){
				double lrps = log_marginal_path_sampling((const Vector**)lls, margl->temperature_count, temperatures, lrssk);
				printf("Path sampling marginal likelihood: %f\n", lrps);
				for (size_t i = 0; i < margl->temperature_count-1; i++) {
					printf("%f: %f\n", temperatures[i], lrssk[i]);
				}
			}
			
			if(margl->ps2){
				double lrps = log_marginal_path_sampling_modified((const Vector**)lls, margl->temperature_count, temperatures, lrssk);
				printf("Modified Path sampling marginal likelihood: %f\n", lrps);
				for (size_t i = 0; i < margl->temperature_count-1; i++) {
					printf("%f: %f\n", temperatures[i], lrssk[i]);
				}
			}
			free(lrssk);
			//test_marg(lls, margl->temperature_count, temperatures);
		}
	}
	else if(margl->refdist_tag != NULL){
		//test_marg_gss(lls, margl->temperature_count, temperatures);
	}
	
	free_StringBuffer(buffer);
	free(temperatures);
	for (size_t i = 0; i < margl->temperature_count; i++) {
		if(lls[i] != NULL)free_Vector(lls[i]);
	}
	free(lls);
}

static void _free_MarginaLikelihood(MarginaLikelihood* margl){
	free(margl->file);
	for (size_t i = 0; i < margl->likelihood_tag_count; i++) {
		free(margl->likelihood_tag[i]);
	}
	free(margl->likelihood_tag);
	if(margl->refdist_tag != NULL) free(margl->refdist_tag);
	if(margl->temperatures != NULL)free(margl->temperatures);
	free(margl);
}

MarginaLikelihood* new_MarginaLikelihood_from_json(json_node* node, Hashtable* hash){
	static const json_field schema[] = {
	    {"algorithm", JSON_OPTIONAL, JSON_ANY},
	    {"alpha", JSON_OPTIONAL, JSON_NUMBER},
	    {"beta", JSON_OPTIONAL, JSON_NUMBER},
	    {"burnin", JSON_OPTIONAL, JSON_ANY},
	    {"distribution", JSON_OPTIONAL, JSON_ANY},
	    {"file", JSON_OPTIONAL, JSON_ANY},
	    {"reference", JSON_OPTIONAL, JSON_ANY},
	    {"steps", JSON_OPTIONAL, JSON_ANY},
	    {"temperatures", JSON_OPTIONAL, JSON_ANY},
	    {"treelikelihood", JSON_OPTIONAL, JSON_ANY},
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	
	MarginaLikelihood* margl = malloc(sizeof(MarginaLikelihood));
	json_node* temp_node = get_json_node(node, "temperatures");
	json_node* steps_node = get_json_node(node, "steps");
	margl->burnin = get_json_node_value_double(node, "burnin", 0);
	json_node* likelihood_node = get_json_node(node, "treelikelihood");
	char* ref_dist_tag = get_json_node_value_string(node, "reference"); // reference distribution for GSS

	margl->refdist_tag = NULL;
	margl->bf = false;
	
	if(ref_dist_tag != NULL){
		margl->refdist_tag = String_clone(ref_dist_tag);
	}
	
	// BF
	if (likelihood_node != NULL && likelihood_node->node_type == MJSON_ARRAY) {
		if (likelihood_node->child_count != 2) {
			fprintf(stderr, "Exactly 2 IDs should be given in `treelikelihood` to calculate a Bayes factor (got %zu)\n\n", (size_t)likelihood_node->child_count);
			exit(1);
		}
		margl->bf = true;
		margl->likelihood_tag_count = likelihood_node->child_count;
		margl->likelihood_tag = malloc(margl->likelihood_tag_count*sizeof(char*));
		for (size_t i = 0; i < margl->likelihood_tag_count; i++) {
			margl->likelihood_tag[i] = String_clone((char*)likelihood_node->children[i]->value);
		}
	}
	else{
		char* likelihood_tag = get_json_node_value_string(node, "treelikelihood");
		margl->likelihood_tag_count = 1;
		margl->likelihood_tag = malloc(sizeof(char*));
		if(likelihood_tag == NULL){
			margl->likelihood_tag[0] = String_clone("treelikelihood");
		}
		else{
			margl->likelihood_tag[0] = String_clone(likelihood_tag);
		}
	}
	
	if (temp_node != NULL) {
		if (temp_node->node_type != MJSON_ARRAY) {
			fprintf(stderr, "attribute `temperatures` should be an array\n\n");
			exit(1);
		}
		margl->temperature_count = temp_node->child_count;
		if (margl->temperature_count == 0) {
			fprintf(stderr, "attribute `temperatures` should not be empty\n\n");
			exit(1);
		}
		margl->temperatures = dvector(margl->temperature_count);
		for (size_t i = 0; i < (size_t)temp_node->child_count; i++) {
			margl->temperatures[i] = atof((char*)temp_node->children[i]->value);
		}
	}
	else if (steps_node != NULL){
		char* dist_string = get_json_node_value_string(node, "distribution");
		margl->temperature_count = get_json_node_value_size_t(node, "steps", 100);
		if (margl->temperature_count < 2) {
			fprintf(stderr, "Attribute `steps` should be at least 2\n\n");
			exit(1);
		}
		margl->temperatures = dvector(margl->temperature_count);
		margl->temperatures[0] = 1;
		margl->temperatures[margl->temperature_count-1] = 0;
		
		// temperatures are in descreasing order
		if(dist_string == NULL){
			fprintf(stderr, "Attribute `distribution` should be specified (`uniform` or `beta`)\n\n");
			exit(1);
		}
		if(strcasecmp(dist_string, "beta") == 0){
			double alpha = get_json_node_value_double(node, "alpha", 0.3);
			double beta = get_json_node_value_double(node, "beta", 1.0);
			double value = 1;
			double incr = 1.0/(margl->temperature_count-1);
			for (size_t i = 1; i < margl->temperature_count-1; i++) {
				value -= incr;
				margl->temperatures[i] = invbetai(value, alpha, beta);
			}
		}
		else if(strcasecmp(dist_string, "uniform") == 0){
			double incr = 1.0/(margl->temperature_count-1);
			for (size_t i = 1; i < margl->temperature_count-1; i++) {
				margl->temperatures[i] = margl->temperatures[i-1] - incr;
			}
		}
		else{
			fprintf(stderr, "Attribute `distribution` should be specified (`uniform` or `beta`)\n\n");
			exit(1);
		}
	}
	else{
		margl->temperature_count = 1;
		margl->temperatures = NULL;
	}
	json_node* algorithm = get_json_node(node, "algorithm");
	margl->shm = margl->hm = margl->ss = margl->ps = margl->ps2 = true; // default is calculate everything
	if(algorithm != NULL){
		if (algorithm->node_type != MJSON_STRING) {
			fprintf(stderr, "Attribute `algorithm` should be a string (`ss`, `ps`, `ps2`, `hm` or `shm`)\n\n");
			exit(1);
		}
		char* meth_string = get_json_node_value_string(node, "algorithm");
		margl->ss = strcasecmp(meth_string, "ss") == 0;
		margl->ps = strcasecmp(meth_string, "ps") == 0;
		margl->ps2 = strcasecmp(meth_string, "ps2") == 0;
		margl->hm = strcasecmp(meth_string, "hm") == 0;
		margl->shm = strcasecmp(meth_string, "shm") == 0;
		if(!(margl->ss || margl->ps || margl->ps2 || margl->hm || margl->shm)){
			fprintf(stderr, "Unknown `algorithm` `%s` (should be `ss`, `ps`, `ps2`, `hm` or `shm`)\n\n", meth_string);
			exit(1);
		}
	}
	
	char* file = get_json_node_value_string(node, "file");
	if(file == NULL){
		fprintf(stderr, "Attribute `file` should be specified (the log file the chains were written to)\n\n");
		exit(1);
	}
	margl->file = String_clone(file);
	margl->run = _marginal_likelihood_run;
	margl->free = _free_MarginaLikelihood;
	if (margl->bf) {
		if (margl->temperature_count < 2) {
			fprintf(stderr, "Attribute `temperatures` or `steps` should be specified to calculate a Bayes factor\n\n");
			exit(1);
		}
		margl->run = _bayes_factor_run;
	}
	return margl;
}
