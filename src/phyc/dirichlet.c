// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "dirichlet.h"

#include <math.h>

double ddirchletln( const double *x, const size_t dim, const double *alphas ){
    double logp = 0;
    double sum = 0;
    for (size_t i = 0; i < dim; i++) {
        logp += (alphas[i]-1.0) * log(x[i]) - gammln(alphas[i]);;
		sum += alphas[i];
    }
    logp += gammln(sum);

    return logp;
}

double ddirchlet( const double *x, const size_t dim, const double *alphas ){
    return exp(ddirchletln(x, dim, alphas));
}

double ddirchlet_flat( const size_t dim ){
    return gamm(dim);
}

void rdirichlet(double*x, const size_t dim, const double* alphas){
	double sum = 0;
	for (int i = 0; i < dim; i++) {
		x[i] = rgamma(alphas[i]);
		sum += x[i];
	}
	for (int i = 0; i < dim; i++) {
		x[i] /= sum;
	}
}
