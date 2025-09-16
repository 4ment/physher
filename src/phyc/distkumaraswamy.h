//
//  distkumaraswamy.h
//  physher
//
//  Created by Mathieu Fourment on 19/8/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef distkumaraswamy_h
#define distkumaraswamy_h

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_KumaraswamyDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_KumaraswamyDistributionModel_from_json(json_node* node, Hashtable* hash);

double DistributionModel_kumaraswamy_inverse_CDF(double p, double a, double b);

#endif /* distkumaraswamy_h */
