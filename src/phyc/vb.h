//
//  vb.h
//  physher
//
//  Created by Mathieu Fourment on 24/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef vb_h
#define vb_h

#include <stdio.h>

#include "mjson.h"
#include "parameters.h"
#include "simplex.h"
#include "matrix.h"

#include <gsl/gsl_rng.h>

typedef struct variational_block_t{
    Model* posterior;
    Model** simplices;
    size_t simplex_count;
    size_t simplex_parameter_count;
    Parameters* parameters; // parameters of the posterior
    Parameters** var_parameters; // parameters of variational distribution
    size_t var_parameters_count; // length of array var_parameters
    void (*sample1)(struct variational_block_t*, double* jacobian); // sample from q and set values in p
    void (*sample2)(struct variational_block_t*, const Parameters*); // sample from q and set values in p
    void (*sample)(struct variational_block_t*, double* values); // sample from q and set values in vector values
    double (*entropy)(struct variational_block_t*);
    void (*grad_elbo)(struct variational_block_t*, const Parameters*, double* grads);
    void (*grad_entropy)(struct variational_block_t*, const Parameters*, double*);
    double (*logP)(struct variational_block_t*, double* values);
    double (*logQ)(struct variational_block_t*, double* values);
    bool use_entropy;
    gsl_rng* rng;
    bool initialized;
    Vector* etas;
}variational_block_t;

typedef struct variational_t{
    Model* posterior;
    variational_block_t** blocks;
    size_t block_count;
	Model** simplices;
	size_t simplex_count;
    Parameters* parameters; // parameters of the posterior
    Parameters* var_parameters; // parameters of variational distribution
    double(*elbofn)(struct variational_t*);
    void(*grad_elbofn)(struct variational_t*, const Parameters*, double*);
    size_t elbo_samples;
    size_t grad_samples;
	bool (*sample)( struct variational_t*, double* values );
	bool (*sample_some)( struct variational_t*, const Parameters* parameters, double* values );
	double (*logP)( struct variational_t*, double* x );
	void (*finalize)( struct variational_t* );
    void (*print)( struct variational_t*, FILE* file );
	bool initialized;
	bool ready_to_sample;
	FILE* file;
	size_t log_samples;
	size_t iter;
	size_t elbo_multi;
	//char* filename;
    gsl_rng* rng;
    // transforms
}variational_t;

Model* new_Variational_from_json(json_node* node, Hashtable* hash);

Model* new_VariationalModel(const char* name, variational_t* var);

double elbo( Parameters *params, double *grad, void *data );

void grad_elbo( Parameters *params, double *grad, void *data );

void init_fullrank_normal(variational_t* var);

void init_meanfield_normal(variational_t* var);

#endif /* vb_h */
