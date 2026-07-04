// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "mjson.h"
#include "parameters.h"

Model* new_JacobianTransformModel_from_json(json_node* node, Hashtable* hash);

#endif  // JACOBIAN_H