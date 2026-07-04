// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef gmrf_h
#define gmrf_h

#include <stdio.h>

#include "hashtable.h"
#include "mjson.h"
#include "parameters.h"

Model* new_GMRFModel_from_json(json_node* node, Hashtable* hash);

#endif /* gmrf_h */
