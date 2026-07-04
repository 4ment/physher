// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef ophmc_h
#define ophmc_h

#include <stdio.h>

#include "operator.h"

Operator* new_HMCOperator_from_json(json_node* node, Hashtable* hash);

#endif /* ophmc_h */
