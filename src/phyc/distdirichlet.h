//
//  distdirichlet.h
//  physher
//
//  Created by Mathieu Fourment on 2/04/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef distdirichlet_h
#define distdirichlet_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_DirichletDistributionModel_with_parameters(Parameters** parameters, Simplex* simplex);

Model* new_DirichletDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distdirichlet_h */
