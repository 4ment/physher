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

//#include <gsl/gsl_rng.h>

struct variational_t{
    Model* posterior;
    Parameters* parameters; // parameters of the posterior
    Parameters* var_parameters; // parameters of variational distribution
    int var_N; // number of variational parameters
    //size_t N; // number of parameters
    double(*elbofn)(struct variational_t*);
    void(*grad_elbofn)(struct variational_t*, double*);
    size_t elbo_samples;
    size_t grad_samples;
	double (*f)( Parameters *x, double *gradient, void *data );
	void (*grad_f)( Parameters *x, double *gradient, void *data );
	bool initialized;
	FILE* file;
	size_t iter;
	//char* filename;
    //gsl_rng* rng;
    // transforms
};

Model* new_Variational_from_json(json_node* node, Hashtable* hash);

Model* new_VariationalModel(const char* name, struct variational_t* var);

double elbo( Parameters *params, double *grad, void *data );

void grad_elbo( Parameters *params, double *grad, void *data );

#endif /* vb_h */
