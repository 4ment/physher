//
//  gradascent.h
//  physher
//
//  Created by Mathieu Fourment on 30/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef gradascent_h
#define gradascent_h

#include <stdio.h>
#include "optimizer.h"

opt_result optimize_stochastic_gradient(Parameters* parameters, opt_func f, opt_grad_func grad_f, double eta, void *data, OptStopCriterion *stop, int verbose, double *fmin);


opt_result optimize_stochastic_gradient_adapt(Parameters* parameters, opt_func f, opt_grad_func grad_f, void(*reset)(void*),
											  double* etas, size_t eta_count, void *data,
											  OptStopCriterion *stop, int verbose, double *best_eta, size_t nthreads);

#endif /* gradascent_h */
