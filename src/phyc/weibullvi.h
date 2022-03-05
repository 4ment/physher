//
//  weibullvi.h
//  physher
//
//  Created by mathieu on 21/5/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef weibullvi_h
#define weibullvi_h

#include "vb.h"

void klqp_block_meanfield_weibull_sample1(variational_block_t* var, double* jacobian);


double klqp_block_meanfield_weibull_entropy(variational_block_t* var);

void klqp_block_meanfield_weibull_sample2(variational_block_t* var, const Parameters* parameters);
void klqp_block_meanfield_weibull_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads);

void klqp_block_meanfield_weibull_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads);

double klqp_block_meanfield_weibull_logP(variational_block_t* var, const double* values);
double klqp_block_meanfield_weibull_logQ(variational_block_t* var, const double* values);
void klqp_block_meanfield_weibull_sample(variational_block_t* var, double* values);

bool klqp_block_meanfield_weibull_sample_some(variational_block_t* var, const Parameters* parameters, double* values);


#endif /* weibullvi_h */
