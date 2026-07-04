// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef klpq_h
#define klpq_h

#include <stdio.h>

#include "mjson.h"
#include "hashtable.h"
#include "parameters.h"

Model* new_KLpqBound_from_json(json_node* node, Hashtable* hash);

#endif /* klpq_h */
