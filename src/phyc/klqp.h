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

//MARK: Meanfield

void klqp_meanfield_normal_init(variational_t* var);

void klqp_meanfield_normal_finalize(variational_t* var);

double klqp_meanfield_normal_elbo(variational_t* var);

double klqp_meanfield_normal_elbo_multi(variational_t* var);

void klqp_meanfield_normal_grad_elbo(variational_t* var, double* grads);

double klqp_meanfield_normal_logP(variational_t* var, double* values);

double klqp_meanfield_normal_logP_parameters(variational_t* var, const Parameters* parameters);

bool klqp_meanfield_normal_sample(variational_t* var, double* values);

bool klqp_meanfield_normal_sample_some(variational_t* var, const Parameters* parameters, double* values);

//MARK: Fullrank

void klqp_fullrank_normal_init(variational_t* var);

double klqp_fullrank_normal_elbo(variational_t* var);

void klqp_fullrank_normal_grad_elbo(variational_t* var, double* grads);

double klqp_fullrank_normal_logP(variational_t* var, double* values);

double klqp_fullrank_normal_logP_parameters(variational_t* var, const Parameters* parameters);

bool klqp_fullrank_normal_sample(variational_t* var, double* values);

#endif /* klqp_h */
