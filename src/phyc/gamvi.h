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

void grad_elbo_gamma_meanfield(variational_t* var, double* grads);

double elbo_gamma_meanfield(variational_t* var);

bool variational_sample_gamma_meanfield(variational_t* var, double* values);

double variational_gamma_meanfield_parameters_logP(variational_t* var, const Parameters* parameters);

double variational_gamma_meanfield_logP(variational_t* var, double* values);

bool variational_sample_some_gamma_meanfield(variational_t* var, const Parameters* parameters, double* values);

#endif /* gamvi_h */
