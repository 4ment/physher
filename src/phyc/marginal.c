//
//  marginal.c
//  physher
//
//  Created by Mathieu Fourment on 18/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "marginal.h"

#include <float.h>
#include <string.h>

#include "utils.h"
#include "matrix.h"
#include "descriptivestats.h"
#include "statistics.h"
#include "mstring.h"
#include "beta.h"

double log_arithmetic_mean(const Vector* vecvalues){
	double sum = -DBL_MAX;
	double* values = Vector_data(vecvalues);
	size_t size = Vector_length(vecvalues);
	for(int i = 0; i < size; i++){
		sum = logaddexp(sum, values[i]);
	}
	return sum - log(size);
}

double log_harmonic_mean(const Vector* vecvalues){
	double sum = 0;
	double* values = Vector_data(vecvalues);
	size_t size = Vector_length(vecvalues);
	for(int i = 0; i < size; i++){
		sum += values[i];
	}
	
	double denominator = -DBL_MAX;
	
	for(int i = 0; i < size; i++){
		denominator = logaddexp(denominator, sum - values[i]);
	}
	return sum - denominator + log(size);
}

double log_smoothed_harmonic_mean(double logP, const Vector* values, double delta){
	double ldelta = log(delta);
	double l1_delta = log(1.0-delta);
	
	int n = Vector_length(values);
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
	for(int i = 1; i < temp_count; i++){
		double tempdiff = temperatures[i]-temperatures[i-1];
		double* previousll = Vector_data(values[i-1]);
		int size = Vector_length(values[i-1]);
		double logmaxll = dmax_vector(previousll, size);
		double temp = 0;
		for (int j = 0; j < size; j++) {
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
	for(int i = 0; i < temp_count; i++){
		means[i] = dmean(Vector_data(values[i]), Vector_length(values[i]));
	}
	
	for(int i = 1; i < temp_count; i++){
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
	for(int i = 0; i < temp_count; i++){
		means[i] = dmean(Vector_data(values[i]), Vector_length(values[i]));
		variances[i] = variance(Vector_data(values[i]), Vector_length(values[i]), means[i]);
	}
	
	for(int i = 1; i < temp_count; i++){
		double weight = temperatures[i]-temperatures[i-1];
		lrpsk[i-1] = weight * (means[i]+means[i-1])/2.0 - ((weight*weight)/12)*(variances[i]-variances[i-1]);
		lrps += lrpsk[i-1];
	}
	free(means);
	free(variances);
	return lrps;
}


static void _marginal_likelihood_run(MarginaLikelihood* margl){
	StringBuffer* buffer = new_StringBuffer(10);
	Vector** lls = malloc(margl->temperature_count*sizeof(Vector*));
	double* temperatures = malloc(margl->temperature_count*sizeof(double));
	
	if(margl->temperature_count == 1){
		lls[0] = read_log_column_with_id(margl->file, margl->burnin, margl->likelihood_tag);
	}
	// Generalized stepping stone
	else if(margl->refdist_tag != NULL){
		double* lrssk = malloc(margl->temperature_count*sizeof(double));
		char* ids[2];
		ids[0] = margl->likelihood_tag;
		ids[1] = margl->refdist_tag;
		// we ignore temperature 1
		for (int i = 1; i < margl->temperature_count; i++) {
			// saved in increasing order
			StringBuffer_empty(buffer);
			StringBuffer_append_format(buffer, "%d%s",i, margl->file);
			printf("Temperature: %f - %s\n", margl->temperatures[i], buffer->c);
			Vector** all = read_log_column_with_ids(buffer->c, margl->burnin, ids, 2);
			size_t size = Vector_length(all[0]);
			lls[margl->temperature_count-i-1] = new_Vector(size);
			for (int j = 0; j < size; j++) {
				Vector_push(lls[margl->temperature_count-i-1], Vector_at(all[0], j) - Vector_at(all[1], j));
			}
			free_Vector(all[0]);
			free_Vector(all[1]);
			free(all);
			temperatures[margl->temperature_count-i-1] = margl->temperatures[i];
		}
		temperatures[margl->temperature_count-1] = margl->temperatures[0];
		lls[margl->temperature_count-1] = NULL;
		
		double lrss = log_marginal_stepping_stone(lls, margl->temperature_count, temperatures, lrssk);
		printf("Generalized stepping stone marginal likelihood: %f\n", lrss);
		for (int i = 0; i < margl->temperature_count-1; i++) {
			printf("%f: %f\n", temperatures[i], lrssk[i]);
		}
		free(lrssk);
	}
	else{
		// temperatures should be in decreasing order
		for (int i = 0; i < margl->temperature_count; i++) {
			printf("Temperature: %f - %s\n", margl->temperatures[i],buffer->c);
			// saved in increasing order
			StringBuffer_empty(buffer);
			StringBuffer_append_format(buffer, "%d%s",i, margl->file);
			lls[margl->temperature_count-i-1] = read_log_column_with_id(buffer->c, margl->burnin, margl->likelihood_tag);
			temperatures[margl->temperature_count-i-1] = margl->temperatures[i];
		}
		
		// sample from prior (temperature == 0)
		double logNaive = log_arithmetic_mean(lls[0]);
		printf("Naive arithmetic mean: %f\n", logNaive);
		
	}
	if(margl->temperature_count == 1 || margl->refdist_tag == NULL){
		double logAM = log_arithmetic_mean(lls[margl->temperature_count-1]);
		printf("Arithmetic mean: %f\n", logAM);
		
		double logHM = log_harmonic_mean(lls[margl->temperature_count-1]);
		printf("Harmonic mean: %f\n", logHM);
		
		if(margl->temperature_count > 1){
			double* lrssk = malloc(margl->temperature_count*sizeof(double));
			
			double lrss = log_marginal_stepping_stone(lls, margl->temperature_count, temperatures, lrssk);
			printf("Stepping stone marginal likelihood: %f\n", lrss);
			for (int i = 0; i < margl->temperature_count-1; i++) {
				printf("%f: %f\n", temperatures[i], lrssk[i]);
			}
			
			double lrps = log_marginal_path_sampling(lls, margl->temperature_count, temperatures, lrssk);
			printf("Path sampling marginal likelihood: %f\n", lrps);
			for (int i = 0; i < margl->temperature_count-1; i++) {
				printf("%f: %f\n", temperatures[i], lrssk[i]);
			}
			lrps = log_marginal_path_sampling_modified(lls, margl->temperature_count, temperatures, lrssk);
			printf("Modified Path sampling marginal likelihood: %f\n", lrps);
			for (int i = 0; i < margl->temperature_count-1; i++) {
				printf("%f: %f\n", temperatures[i], lrssk[i]);
			}
			free(lrssk);
		}
	}
	
	free_StringBuffer(buffer);
	free(temperatures);
	for (int i = 0; i < margl->temperature_count; i++) {
		if(lls[i] != NULL)free_Vector(lls[i]);
	}
	free(lls);
}

static void _free_MarginaLikelihood(MarginaLikelihood* margl){
	free(margl->file);
	free(margl->likelihood_tag);
	if(margl->refdist_tag != NULL) free(margl->refdist_tag);
	if(margl->temperatures != NULL)free(margl->temperatures);
	free(margl);
}

MarginaLikelihood* new_MarginaLikelihood_from_json(json_node* node, Hashtable* hash){
	MarginaLikelihood* margl = malloc(sizeof(MarginaLikelihood));
	json_node* temp_node = get_json_node(node, "temperatures");
	json_node* steps_node = get_json_node(node, "steps");
	margl->burnin = get_json_node_value_double(node, "burnin", 0);
	char* likelihood_tag = get_json_node_value_string(node, "treelikelihood");
	char* ref_dist_tag = get_json_node_value_string(node, "reference"); // reference distribution for GSS

	margl->refdist_tag = NULL;
	
	if(ref_dist_tag != NULL){
		margl->refdist_tag = String_clone(ref_dist_tag);
	}
	
	if(likelihood_tag == NULL){
		margl->likelihood_tag = String_clone("treelikelihood");
	}
	else{
		margl->likelihood_tag = String_clone(likelihood_tag);
	}
	
	if (temp_node != NULL) {
		if (temp_node->node_type != MJSON_ARRAY) {
			fprintf(stderr, "attribute `temperatures` should be an array\n\n");
			exit(1);
		}
		margl->temperature_count = temp_node->child_count;
		margl->temperatures = dvector(margl->temperature_count);
		for (int i = 0; i < temp_node->child_count; i++) {
			margl->temperatures[i] = atof((char*)temp_node->children[i]->value);
		}
	}
	else if (steps_node != NULL){
		char* dist_string = get_json_node_value_string(node, "distribution");
		margl->temperature_count = get_json_node_value_size_t(node, "steps", 100);
		margl->temperatures = dvector(margl->temperature_count);
		margl->temperatures[0] = 1;
		margl->temperatures[margl->temperature_count-1] = 0;
		
		// temperatures are in descreasing order
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
				margl->temperatures[i] = margl->temperatures[i+1] - incr;
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
	char* file = get_json_node_value_string(node, "file");
	margl->file = String_clone(file);
	margl->run = _marginal_likelihood_run;
	margl->free = _free_MarginaLikelihood;
	return margl;
}
