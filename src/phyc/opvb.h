// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef opvb_h
#define opvb_h

#include <stdio.h>

#include "operator.h"

Operator* new_VariationalOperator_from_json(json_node* node, Hashtable* hash);

#endif /* opvb_h */
