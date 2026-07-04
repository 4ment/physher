// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "parameters.h"
#include "mjson.h"
#include "hashtable.h"

Model* new_BoundModel_from_json(json_node* node, Hashtable* hash);