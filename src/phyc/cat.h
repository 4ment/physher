// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef cat_h
#define cat_h

#include <stdio.h>

#include "mjson.h"
#include "hashtable.h"

void cat_estimator_from_json(json_node* node, Hashtable* hash);

#endif /* cat_h */
