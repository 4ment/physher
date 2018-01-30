//
//  mmcmc.c
//  physher
//
//  Created by Mathieu Fourment on 17/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "mmcmc.h"

#include <string.h>

#include "matrix.h"
#include "filereader.h"
#include "marginal.h"
#include "beta.h"
#include "statistics.h"
#include "descriptivestats.h"
#include "random.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

Vector* read_log_last_column( const char *filename, size_t burnin ){
	int count = 0;
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	Vector* vec = new_Vector(1000);
	
	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);// discard header
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		
		if ( reader->buffer->length == 0){
			continue;
		}
		if ( count >= burnin){
			ptr = reader->line;
			l = 0;
			temp = String_split_char_double( ptr, '\t', &l );
			Vector_push(vec, temp[l-1]);
			free(temp);
		}
		count++;
	}
	free_FileReader(reader);
	Vector_pack(vec);
	return  vec;
}

Vector** read_log_except_last_column( const char *filename, size_t burnin, size_t* count ){
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	size_t capacity = 1000;
	*count = 0;
	Vector** vecs = malloc(capacity*sizeof(Vector*));
	size_t sample = 0;

	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);// discard header
	
	while ( reader->read_line(reader) ) {
		StringBuffer_trim(reader->buffer);
		
		if ( reader->buffer->length == 0){
			continue;
		}
		if ( sample >= burnin){
			if(*count == capacity){
				capacity *= 2;
				vecs = realloc(vecs, capacity*sizeof(Vector*));
			}
			ptr = reader->line;
			l = 0;
			temp = String_split_char_double( ptr, '\t', &l );
			vecs[*count] = new_Vector(l-4);
			for (int i = 3; i < l-1; i++) {
				Vector_push(vecs[*count], temp[i]);
			}
			free(temp);
			(*count)++;
		}
		sample++;
	}
	free_FileReader(reader);
	vecs = realloc(vecs, *count*sizeof(Vector*));
	return  vecs;
}

Vector** sample_multivariate_normal(gsl_rng *rng, size_t n, const double* means, const double* covar, size_t dim){
	Vector** values = malloc(n*sizeof(Vector*));
	gsl_vector * mu = gsl_vector_calloc(dim);
	gsl_matrix * Sigma = gsl_matrix_calloc(dim, dim);
	gsl_matrix * L = gsl_matrix_calloc(dim, dim);
	gsl_vector * samples = gsl_vector_calloc(dim);
	
	for (int i = 0; i < dim; i++) {
		gsl_vector_set(mu, i, means[i]);
		for (int j = 0; j < dim; j++) {
			gsl_matrix_set(Sigma, i, j, covar[i*dim+j]);
		}
	}
	
	gsl_matrix_memcpy(L, Sigma);
	gsl_linalg_cholesky_decomp1(L);
	
	for (int i = 0; i < n; i++) {
		gsl_ran_multivariate_gaussian(rng, mu, L, samples);
		values[i] = new_Vector(dim);
		for (int j = 0; j < dim; j++) {
			Vector_push(values[i], gsl_vector_get(samples, j));
		}
	}
	
	gsl_vector_free(mu);
	gsl_vector_free(samples);
	gsl_matrix_free(L);
	gsl_matrix_free(Sigma);
	return values;
}

double* multivariate_normal_logP(Vector** values, size_t n, const double* means, const double* covar, size_t dim){
	double* logPs = dvector(n);
	gsl_vector * mu = gsl_vector_calloc(dim);
	gsl_matrix * Sigma = gsl_matrix_calloc(dim, dim);
	gsl_matrix * L = gsl_matrix_calloc(dim, dim);
	gsl_vector * x = gsl_vector_calloc(dim);
	gsl_vector * work = gsl_vector_calloc(dim);
	
	for (int i = 0; i < dim; i++) {
		gsl_vector_set(mu, i, means[i]);
		for (int j = 0; j < dim; j++) {
			gsl_matrix_set(Sigma, i, j, covar[i*dim+j]);
		}
	}
	
	gsl_matrix_memcpy(L, Sigma);
	gsl_linalg_cholesky_decomp1(L);
	
	for (int i = 0; i < n; i++) {
		double* vv = Vector_data(values[i]);
		for (int j = 0; j < dim; j++) {
			gsl_vector_set(x, j,  vv[j]);
		}
		gsl_ran_multivariate_gaussian_log_pdf (x,  mu, L, &logPs[i], work);
	}
	
	gsl_vector_free(mu);
	gsl_matrix_free(L);
	gsl_matrix_free(Sigma);
	gsl_vector_free(x);
	gsl_vector_free(work);
	return logPs;
}

//double nested_sampling(Vector** lls, size_t count){
//
//	double logZ = -DBL_MAX;
//	size_t nLogL = Vector_length(lls[0]);
//	double* logL = clone_dvector(Vector_data(lls[0]), nLogL);
//	double N = nLogL;
//	int j = count - 1;
//	double logw = log(1.0 - exp(-1.0/N));
//	
//	for (int i = 1; i <= j; i++) {
//		double ii = i;
//		size_t worst = which_dmin(logL, nLogL);
//		double lw = logL[worst] + logw - (ii-1.0)/N;
//		logZ = logaddexp(logZ, lw);
//		printf("%f %lu %f\n", logZ, worst, logL[worst]);
//
//		while (1) {
//			size_t index = random_int(nLogL);
//			if (Vector_at(lls[i], index) > logL[worst]) {
//				logL[worst] = Vector_at(lls[i], index);
//				break;
//			}
//		}
//	}
//	double sum = -DBL_MAX;
//	for (int i = 0; i < nLogL; i++) {
//		sum = logaddexp(sum, logL[i]);
//	}
////	print_dvector(logL, nLogL); 
//	printf("%e\n", sum);
//	sum += -(double)j/N - log(N);
//	printf("%f %f\n", sum, logZ);
//	return logaddexp(sum, logZ);
//}

double bridge_sampling(MMCMC*mmcmc, const char* filename, double guess){
	size_t count = 0;
	Vector** params = read_log_except_last_column(filename, mmcmc->burnin, &count);

	size_t aN1 = count/2;
	size_t bN1 = count - aN1;
	size_t N2 = aN1;
	double* l1s = dvector(bN1);
	double* l2s = dvector(N2);
	size_t paramCount = Vector_length(params[0]);
	
	// First batch
	Vector** params_transformed = malloc(paramCount*sizeof(Vector*));
	double* means = dvector(paramCount);
	double* covariances = dvector(paramCount*paramCount);
	// Transform parameters and calculate mean for each vector of parameter
	// transposed log(params)
	for (int i = 0; i < paramCount; i++) {
		params_transformed[i] = new_Vector(aN1);
		for (int j = 0; j < aN1; j++) {
			Vector_push(params_transformed[i], log(Vector_at(params[j], i)));
		}
		means[i] = mean(Vector_data(params_transformed[i]), aN1);
	}
	
	
	// Calculate covariance matrix
	for (int i = 0; i < paramCount; i++) {
		double* pp = Vector_data(params_transformed[i]);
		covariances[i*paramCount+i] = variance(pp, aN1, means[i]);
		for (int j = i+1; j < paramCount; j++) {
			double* pp2 = Vector_data(params_transformed[j]);
			covariances[i*paramCount+j] = covariances[j*paramCount+i] = covariance(pp, pp2, means[i], means[j], aN1);
		}
	}

	// sample N2 samples from multivariate normal
	Vector** samples = sample_multivariate_normal(mmcmc->rng, N2, means, covariances, paramCount);
	
	// Calculate probability of samples
	double* proposal_logPs = multivariate_normal_logP(samples, N2, means, covariances, paramCount);
	Parameters* parameters = mmcmc->x;
	Model* posterior = mmcmc->mcmc->model;

	for (int i = 0; i < N2; i++) {
		l2s[i] = 0;
		for (int j = 0; j < paramCount; j++) {
			l2s[i] += Vector_at(samples[i], j); // change of variable
			Parameters_set_value(parameters, j, exp(Vector_at(samples[i], j)));
		}
		l2s[i] += posterior->logP(posterior) - proposal_logPs[i];
	}
	
	// Second batch
	
	Vector** params_transformed2 = malloc(bN1*sizeof(Vector*));
	// Transform parameters for each vector of parameters
	for (int i = 0; i < bN1; i++) {
		params_transformed2[i] = new_Vector(paramCount);
		for (int j = 0; j < paramCount; j++) {
			Vector_push(params_transformed2[i], log(Vector_at(params[i+aN1], j)));
		}
	}
	
	// Calculate probability of samples
	double* proposal_logPs2 = multivariate_normal_logP(params_transformed2, bN1, means, covariances, paramCount);

	for (int i = 0; i < bN1; i++) {
		l1s[i] = 0;
		for (int j = 0; j < Parameters_count(parameters); j++) {
			l1s[i] += Vector_at(params_transformed2[i], j);
//			Parameters_set_value(parameters, j, exp(Vector_at(params_transformed2[i], j)));
			Parameters_set_value(parameters, j, Vector_at(params[i+aN1], j));
		}
		l1s[i] += posterior->logP(posterior) - proposal_logPs2[i];
	}
	
	double s1 = (double)bN1/(bN1+N2);
	double s2 = (double)N2/(bN1+N2);

	double lstar = dmedian(l1s, aN1);
	double r = exp(guess-lstar);
	for(size_t i = 0; i < 10000; i++) {
		double rr = r;
		double num = 0;
		for (int i = 0; i < N2; i++) {
			num += exp(l2s[i]-lstar)/(s1*exp(l2s[i]-lstar) + s2*r);
		}
		num /= N2;
//		printf("num %e\n", num);
		double denom = 0;
		for (int i = 0; i < bN1; i++) {
			denom += 1.0/(s1*exp(l1s[i]-lstar) + s2*r);
		}
		denom /= bN1;
//		printf("denom %e\n", denom);
		r = num/denom;
//		printf("%f %f\n", marg, marg2);
		if(fabs(log(r)-log(rr)) < 1.0e-4){
			break;
		}
		rr = r;
		printf("marg %f %f %f\n", log(r)+lstar, lstar, r);
	}
	
	for (int i = 0; i < bN1; i++) {
		free_Vector(params_transformed2[i]);
	}
	for (int i = 0; i < paramCount; i++) {
		free_Vector(params_transformed[i]);
	}
	for (int i = 0; i < N2; i++) {
		free_Vector(samples[i]);
	}
	for (int i = 0; i < count; i++) {
		free_Vector(params[i]);
	}
	free(params_transformed);
	free(params_transformed2);
	free(proposal_logPs);
	free(proposal_logPs2);
	free(means);
	free(covariances);
	free(l1s);
	free(l2s);
	free(params);
	free(samples);
	
	return log(r)+lstar;
}

void mmcmc_run(MMCMC* mmcmc){
	StringBuffer* buffer = new_StringBuffer(10);
	MCMC* mcmc = mmcmc->mcmc;
	Vector** lls = malloc(mmcmc->temperature_count*sizeof(Vector*));
	double* lrssk = malloc(mmcmc->temperature_count*sizeof(double));
	double* temperatures = malloc(mmcmc->temperature_count*sizeof(double));
	char** filenames = malloc(sizeof(char*)*mcmc->log_count);
	
	for (int j = 0; j < mcmc->log_count; j++) {
		filenames[j] = NULL;
		if(mcmc->logs[j]->filename != NULL){
			filenames[j] = String_clone(mcmc->logs[j]->filename);
		}
	}
	
	// temperatures should be in decreasing order
	for (int i = 0; i < mmcmc->temperature_count; i++) {
		for (int j = 0; j < mcmc->log_count; j++) {
			if(filenames[j] != NULL){
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "%d%s",i, filenames[j]);
				free(mcmc->logs[j]->filename);
				mcmc->logs[j]->filename = StringBuffer_tochar(buffer);
			}
		}
		printf("Temperature: %f - %s\n", mmcmc->temperatures[i],buffer->c);
		mcmc->chain_temperature = mmcmc->temperatures[i];
		mcmc->run(mcmc);
		
		// saved in increasing order
		for (int j = 0; j < mcmc->log_count; j++) {
			if(filenames[j] != NULL && strcmp(mmcmc->log_file, filenames[j]) == 0){
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "%d%s",i, filenames[j]);
				lls[mmcmc->temperature_count-i-1] = read_log_last_column(buffer->c, mmcmc->burnin);
				temperatures[mmcmc->temperature_count-i-1] = mmcmc->temperatures[i];
			}
		}
	}

	// sample from prior (temperature == 0)
	double logNaive = log_arithmetic_mean(lls[0]);
	printf("Naive arithmetic mean: %f\n", logNaive);
	
	double logAM = log_arithmetic_mean(lls[mmcmc->temperature_count-1]);
	printf("Arithmetic mean: %f\n", logAM);
	
	double logHM = log_harmonic_mean(lls[mmcmc->temperature_count-1]);
	printf("Harmonic mean: %f\n", logHM);
	
	StringBuffer_empty(buffer);
	StringBuffer_append_format(buffer, "%d%s", 0, mmcmc->log_file);
	double logBS = bridge_sampling(mmcmc, buffer->c, logHM);
	printf("Bridge sampling: %f\n", logBS);
	
	double lrss = log_marginal_stepping_stone(lls, mmcmc->temperature_count, temperatures, lrssk);
	printf("Stepping stone marginal likelihood: %f\n", lrss);
	for (int i = 1; i < mmcmc->temperature_count; i++) {
		printf("%f: %f\n", temperatures[i-1], lrssk[i]);
	}
	
	double lrps = log_marginal_path_sampling(lls, mmcmc->temperature_count, temperatures, lrssk);
	printf("Path sampling marginal likelihood: %f\n", lrps);
	for (int i = 0; i < mmcmc->temperature_count; i++) {
		printf("%f: %f\n", temperatures[i], lrssk[i]);
	}
	lrps = log_marginal_path_sampling_modified(lls, mmcmc->temperature_count, temperatures, lrssk);
	printf("Modified Path sampling marginal likelihood: %f\n", lrps);
	for (int i = 0; i < mmcmc->temperature_count; i++) {
		printf("%f: %f\n", temperatures[i], lrssk[i]);
	}
	
	// Leave it as a standard mcmc with original loggers
	mcmc->chain_temperature = -1;
	for (int j = 0; j < mcmc->log_count; j++) {
		if(filenames != NULL){
			mcmc->logs[j]->filename = String_clone(filenames[j]);
			free(filenames[j]);
		}
	}
	
	free(filenames);
	free_StringBuffer(buffer);
	free(lrssk);
	free(temperatures);
	for (int i = 0; i < mmcmc->temperature_count; i++) {
		free_Vector(lls[i]);
	}
	free(lls);
}

MMCMC* new_MMCMC_from_json(json_node* node, Hashtable* hash){
	MMCMC* mmcmc = malloc(sizeof(MMCMC));
	json_node* mcmc_node = get_json_node(node, "mcmc");
	json_node* temp_node = get_json_node(node, "temperatures");
	json_node* steps_node = get_json_node(node, "steps");
	
	mmcmc->burnin = get_json_node_value_double(node, "burnin", 0);
	
	if (temp_node != NULL) {
		if (temp_node->node_type != MJSON_ARRAY) {
			fprintf(stderr, "attribute `temperatures` should be an array\n\n");
			exit(1);
		}
		mmcmc->temperature_count = temp_node->child_count;
		mmcmc->temperatures = dvector(mmcmc->temperature_count);
		for (int i = 0; i < temp_node->child_count; i++) {
			mmcmc->temperatures[i] = atof((char*)temp_node->children[i]->value);
		}
	}
	else if (steps_node != NULL){
		char* dist_string = get_json_node_value_string(node, "distribution");
		mmcmc->temperature_count = get_json_node_value_size_t(node, "steps", 100);
		mmcmc->temperatures = dvector(mmcmc->temperature_count);
		mmcmc->temperatures[0] = 1;
		mmcmc->temperatures[mmcmc->temperature_count-1] = 0;
		
		// temperatures are in descreasing order
		if(strcasecmp(dist_string, "beta") == 0){
			double alpha = get_json_node_value_double(node, "alpha", 0.3);
			double beta = get_json_node_value_double(node, "beta", 1.0);
			double value = 1;
			double incr = 1.0/(mmcmc->temperature_count-1);
			for (size_t i = 1; i < mmcmc->temperature_count-1; i++) {
				value -= incr;
				mmcmc->temperatures[i] = invbetai(value, alpha, beta);
			}
		}
		else if(strcasecmp(dist_string, "uniform") == 0){
			double incr = 1.0/(mmcmc->temperature_count-1);
			for (size_t i = 1; i < mmcmc->temperature_count-1; i++) {
				mmcmc->temperatures[i] = mmcmc->temperatures[i+1] - incr;
			}
		}
		else{
			fprintf(stderr, "Attribute `distribution` should be specified (`uniform` or `beta`)\n\n");
			exit(1);
		}
	}
	else{
		fprintf(stderr, "Attribute `temperatures` or `steps` should be specified\n\n");
		exit(1);
	}
	mmcmc->mcmc = new_MCMC_from_json(mcmc_node, hash);
	mmcmc->run = mmcmc_run;
	

	mmcmc->x = NULL;
	json_node* x_node = get_json_node(node, "x");
	char* lfile = get_json_node_value_string(node, "log_file");
	if (x_node != NULL) {
		 mmcmc->x = new_Parameters(1);
		get_parameters_references2(node, hash, mmcmc->x, "x");
	}
	if(lfile != NULL){
		mmcmc->log_file = String_clone(lfile);
	}
	mmcmc->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return mmcmc;
}

void free_MMCMC(MMCMC* mmcmc){
	mmcmc->mcmc->free(mmcmc->mcmc);
	free(mmcmc->temperatures);
	free_Parameters(mmcmc->x);
	if(mmcmc->log_file != NULL) free(mmcmc->log_file);
	free(mmcmc);
}
