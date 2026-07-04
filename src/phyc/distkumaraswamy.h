// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distkumaraswamy_h
#define distkumaraswamy_h

#include "hashtable.h"
#include "mjson.h"
#include "distmodel.h"

DistributionModel* new_KumaraswamyDistributionModel_with_parameters(Parameters* parameters, Parameters* x);

Model* new_KumaraswamyDistributionModel_from_json(json_node* node, Hashtable* hash);

double DistributionModel_kumaraswamy_inverse_CDF(double p, double a, double b);

#endif /* distkumaraswamy_h */
