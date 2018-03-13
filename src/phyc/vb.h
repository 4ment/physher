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

#include <gsl/gsl_rng.h>


typedef struct variational_t{
    Model* posterior;
    Parameters* parameters; // parameters of the posterior
    Parameters* var_parameters; // parameters of variational distribution
    double(*elbofn)(struct variational_t*);
    void(*grad_elbofn)(struct variational_t*, double*);
    size_t elbo_samples;
    size_t grad_samples;
	double (*f)( Parameters *x, double *gradient, void *data );
	void (*grad_f)( Parameters *x, double *gradient, void *data );
	bool (*sample)( struct variational_t*, double* values );
	bool (*sample_some)( struct variational_t*, const Parameters* parameters, double* values );
	double (*logP)( struct variational_t*, double* x );
    double (*parameters_logP)( struct variational_t*, const Parameters* parameters );
	void (*finalize)( struct variational_t* );
	bool initialized;
	bool ready_to_sample;
	FILE* file;
	size_t iter;
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
