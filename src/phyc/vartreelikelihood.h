//
//  vartreelikelihood.h
//  physher
//
//  Created by mathieu on 16/4/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#ifndef vartreelikelihood_h
#define vartreelikelihood_h

#include "parameters.h"
#include "tree.h"
#include "branchmodel.h"
#include "distmodel.h"

typedef struct VariationalTreeLikelihood{
    Tree* tree;
    BranchModel* bm;
    DistributionModel* distribution; // most likely lognormal
    bool update;
    double* gradient;
    double* gradient_distance;
}VariationalTreeLikelihood;

Model * new_VariationalTreeLikelihoodModel_from_json(json_node*node, Hashtable*hash);

#endif /* vartreelikelihood_h */
