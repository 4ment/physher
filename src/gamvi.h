//
//  gamvi.h
//  viphy
//
//  Created by Mathieu Fourment on 12/04/2017.
//  Copyright © 2017 University of Technology Sydney. All rights reserved.
//

#ifndef gamvi_h
#define gamvi_h

//#include <stdio.h>

#include <gsl/gsl_rng.h>

#include "phyc/treelikelihood.h"

void grad_elbo_gamma_meanfield(SingleTreeLikelihood* tlk, Node **nodes, int nodeCount, int sampleCount, const double* params, double* grads,
                         const gsl_rng* r);

double elbo_gamma_meanfield(SingleTreeLikelihood* tlk, Node **nodes, int nodeCount, int sampleCount, const double* params, const gsl_rng* r);

#endif /* gamvi_h */
