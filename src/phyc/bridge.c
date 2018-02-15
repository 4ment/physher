//
//  bridge.c
//  physher
//
//  Created by Mathieu Fourment on 9/02/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "bridge.h"

#include <math.h>
#include <string.h>
#include <strings.h>

#include "matrix.h"
#include "statistics.h"
#include "descriptivestats.h"
#include "marginal.h"
#include "filereader.h"
#include "compoundmodel.h"
#include "treelikelihood.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

Vector** read_log_for_parameters( const char *filename, size_t burnin, size_t* count, Parameters* params ){
	char *ptr = NULL;
	double *temp = NULL;
	int l;
	size_t capacity = 1000;
	*count = 0;
	size_t paramCount = Parameters_count(params);
	Vector** vecs = malloc(capacity*sizeof(Vector*));
	size_t sample = 0;
	
	FileReader *reader = new_FileReader(filename, 1000);
	reader->read_line(reader);// discard header
	ptr = reader->line;
	l = 0;
	char** header = String_split_char(ptr, '\t', &l);
	bool* in = bvector(l);
	
	for (int j = 0; j < paramCount; j++) {
		for (int i = 0; i < l; i++) {
			if(strcmp(Parameters_name(params, j), header[i]) == 0){
				in[i] = true;
				break;
			}
		}
	}
	
	for (int i = 0; i < l; i++) {
		free(header[i]);
	}
	free(header);
	
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
			vecs[*count] = new_Vector(paramCount);
			for (int i = 0; i < l; i++) {
				if(in[i]){
					Vector_push(vecs[*count], temp[i]);
				}
			}
			free(temp);
			(*count)++;
		}
		sample++;
	}
	free_FileReader(reader);
	vecs = realloc(vecs, *count*sizeof(Vector*));
	free(in);
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

void bridge_sampling(BridgeSampling* bs){
	size_t count = 0;
	Vector** params = read_log_for_parameters(bs->file, bs->burnin, &count, bs->x);
	
	Vector* lls = read_log_column_with_id(bs->file, bs->burnin, bs->likelihood_tag);
	
	double guess = log_harmonic_mean(lls);
	printf("%f\n", guess);
	
	free_Vector(lls);
	
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
	Vector** samples = sample_multivariate_normal(bs->rng, N2, means, covariances, paramCount);
	
	// Calculate probability of samples
	double* proposal_logPs = multivariate_normal_logP(samples, N2, means, covariances, paramCount);
	Parameters* parameters = bs->x;
	Model* posterior = bs->model;
	
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
	
	fprintf(stdout, "Bridge sampling: %f\n", log(r)+lstar);
}

static void _free_BridgeSampling(BridgeSampling* bs){
	free_Parameters(bs->x);
	free_Model(bs->model);
	free(bs->likelihood_tag);
	free(bs->file);
	free(bs);
}

BridgeSampling* new_BridgeSampling_from_json(json_node* node, Hashtable* hash){
	BridgeSampling* bs = malloc(sizeof(BridgeSampling));
	bs->burnin = get_json_node_value_double(node, "burnin", 0);
	char* likelihood_tag = get_json_node_value_string(node, "treelikelihood");
	char* lfile = get_json_node_value_string(node, "file");
	json_node* model_node = get_json_node(node, "model");
	bs->x = new_Parameters(1);
	get_parameters_references2(node, hash, bs->x, "x");
	
	bs->file = String_clone(lfile);
	
	if(likelihood_tag == NULL){
		bs->likelihood_tag = String_clone("treelikelihood");
	}
	else{
		bs->likelihood_tag = String_clone(likelihood_tag);
	}
	
	if(model_node->node_type == MJSON_OBJECT){
		char* type = get_json_node_value_string(model_node, "type");
		char* id = get_json_node_value_string(model_node, "id");
		
		if (strcasecmp(type, "compound") == 0) {
			bs->model = new_CompoundModel_from_json(model_node, hash);
		}
		else if (strcasecmp(type, "treelikelihood") == 0) {
			bs->model = new_TreeLikelihoodModel_from_json(model_node, hash);
		}
		Hashtable_add(hash, id, bs->model);
		// could be a parametric distribution
	}
	// ref
	else if(model_node->node_type == MJSON_STRING){
		char* model_string = model_node->value;
		bs->model = Hashtable_get(hash, model_string+1);
		bs->model->ref_count++;
	}
	bs->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	bs->free = _free_BridgeSampling;
	bs->run = bridge_sampling;
	return bs;
}
