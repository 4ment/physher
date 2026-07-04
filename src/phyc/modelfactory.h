// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef MODELFACTORY_H
#define MODELFACTORY_H

#include "model.h"
#include "hashtable.h"

Model* model_factory_from_json(json_node* node, Hashtable* hash);

#endif // MODELFACTORY_H