//
//  klqp.h
//  physher
//
//  Created by Mathieu Fourment on 28/3/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef klqp_h
#define klqp_h

#include <stdio.h>

#include "vb.h"


double variational_klqp_elbo(variational_t* var);

double variational_klqp_elbo_multi(variational_t* var);

void variational_klqp_grad_elbo(variational_t* var, const Parameters* parameters, double* grads);


//MARK: Meanfield block

void klqp_block_meanfield_normal_sample1(variational_block_t* var, double* jacobian);
void klqp_block_meanfield_normal_sample2(variational_block_t* var, const Parameters* parameters);

double klqp_block_meanfield_normal_entropy(variational_block_t* var);

void klqp_block_meanfield_normal_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads);
void klqp_block_meanfield_normal_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads);

double klqp_block_meanfield_normal_logP(variational_block_t* var, const double* values);
double klqp_block_meanfield_normal_logQ(variational_block_t* var, const double* values);

void klqp_block_meanfield_normal_sample(variational_block_t* var, double* values);

bool klqp_block_meanfield_normal_sample_some(variational_block_t* var, const Parameters* parameters, double* values);

void klqp_block_meanfield_normal_initialize(variational_block_t* var);

//MARK: Fullrank block

void klqp_block_fullrank_normal_sample1(variational_block_t* var, double* jacobian);
void klqp_block_fullrank_normal_sample2(variational_block_t* var, const Parameters* parameters);

double klqp_block_fullrank_normal_entropy(variational_block_t* var);

void klqp_block_fullrank_normal_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads);
void klqp_block_fullrank_normal_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads);

double klqp_block_fullrank_normal_logP(variational_block_t* var, const double* values);
double klqp_block_fullrank_normal_logQ(variational_block_t* var, const double* values);

void klqp_block_fullrank_normal_sample(variational_block_t* var, double* values);

bool klqp_block_fullrank_normal_sample_some(variational_block_t* var, const Parameters* parameters, double* values);

//MARK: Meanfield

void klqp_meanfield_normal_init(variational_t* var);

void klqp_meanfield_normal_finalize(variational_t* var);

double klqp_meanfield_normal_elbo(variational_t* var);

double klqp_meanfield_normal_elbo_multi(variational_t* var);

void klqp_meanfield_normal_grad_elbo(variational_t* var, const Parameters* parameters, double* grads);

double klqp_meanfield_normal_logP(variational_t* var, const double* values);

bool klqp_meanfield_normal_sample(variational_t* var, double* values);

bool klqp_meanfield_normal_sample_some(variational_t* var, const Parameters* parameters, double* values);

void klqp_meanfield_normal_log_samples(variational_t* var, FILE* file);

//MARK: lognormal meanfield

double klqp_meanfield_lognormal_elbo(variational_t* var);

void klqp_meanfield_lognormal_grad_elbo(variational_t* var, const Parameters* parameters, double* grads);

double klqp_meanfield_lognormal_logP(variational_t* var, const double* values);

bool klqp_meanfield_lognormal_sample(variational_t* var, double* values);

bool klqp_meanfield_lognormal_sample_some(variational_t* var, const Parameters* parameters, double* values);

void klqp_meanfield_lognormal_log_samples(variational_t* var, FILE* file);


//MARK: Fullrank

void klqp_fullrank_normal_init(variational_t* var);

double klqp_fullrank_normal_elbo(variational_t* var);

void klqp_fullrank_normal_grad_elbo(variational_t* var, const Parameters* parameters, double* grads);

double klqp_fullrank_normal_logP(variational_t* var, const double* values);

bool klqp_fullrank_normal_sample(variational_t* var, double* values);

void klqp_fullrank_log_samples(variational_t* var, FILE* file);

#endif /* klqp_h */
