//
//  gamvi.h
//  viphy
//
//  Created by Mathieu Fourment on 12/04/2017.
//  Copyright © 2017 University of Technology Sydney. All rights reserved.
//

#ifndef gamvi_h
#define gamvi_h

#include "vb.h"

void klqp_block_meanfield_gamma_sample1(variational_block_t* var, double* jacobian);
void klqp_block_meanfield_gamma_sample2(variational_block_t* var, const Parameters* parameters);

double klqp_block_meanfield_gamma_entropy(variational_block_t* var);

void klqp_block_meanfield_gamma_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads);
void klqp_block_meanfield_gamma_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads);

double klqp_block_meanfield_gamma_logP(variational_block_t* var, const double* values);
double klqp_block_meanfield_gamma_logQ(variational_block_t* var, const double* values);

void klqp_block_meanfield_gamma_sample(variational_block_t* var, double* values);
bool klqp_block_meanfield_gamma_sample_some(variational_block_t* var, const Parameters* parameters, double* values);


void grad_elbo_gamma_meanfield(variational_t* var, const Parameters* parameters, double* grads);

double elbo_gamma_meanfield(variational_t* var);

bool variational_sample_gamma_meanfield(variational_t* var, double* values);

double variational_gamma_meanfield_parameters_logP(variational_t* var, const Parameters* parameters);

double variational_gamma_meanfield_logP(variational_t* var, const double* values);

bool variational_sample_some_gamma_meanfield(variational_t* var, const Parameters* parameters, double* values);

void meanfield_gamma_log_samples(variational_t* var, FILE* file);

#endif /* gamvi_h */
