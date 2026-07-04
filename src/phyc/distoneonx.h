// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef distoneonx_h
#define distoneonx_h

#include "distmodel.h"
#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"

DistributionModel* new_OneOnXDistributionModel(Parameters* x);

Model* new_OneOnXDistributionModel_from_json(json_node* node, Hashtable* hash);

#endif /* distoneonx_h */
