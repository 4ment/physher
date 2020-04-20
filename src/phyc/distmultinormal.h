//
//  distmultinormal.h
//  physher
//
//  Created by mathieu on 17/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef distmultinormal_h
#define distmultinormal_h

#include <stdlib.h>

#include <gsl/gsl_randist.h>

#include "parameters.h"
#include "mjson.h"
#include "hashtable.h"

typedef struct gsl_multivariate_normal_wrapper_t{
    gsl_vector* mu;
    gsl_matrix* L;
    gsl_vector* x; // used for sampling or pdf
    gsl_vector * work;
    gsl_rng* rng;
    bool transform;
}gsl_multivariate_normal_wrapper_t;

Model* new_MultivariateNormalDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distmultinormal_h */
