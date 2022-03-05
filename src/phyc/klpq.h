//
//  klpq.h
//  physher
//
//  Created by Mathieu Fourment on 8/03/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef klpq_h
#define klpq_h

#include <stdio.h>

#include "vb.h"

void grad_klpq_normal_meanfield(variational_t* var, const Parameters* parameters, double* grads);

double klpq_normal_meanfield(variational_t* var);

#endif /* klpq_h */
