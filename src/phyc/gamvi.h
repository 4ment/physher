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

#endif /* gamvi_h */
