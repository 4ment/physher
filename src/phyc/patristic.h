// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PHYC_PATRISTIC_H_
#define PHYC_PATRISTIC_H_

#include "tree.h"
#include "hashtable.h"

double ** calculate_patristic( Tree *tree );

double ** calculate_patristic2( Tree **trees, int nTree, Hashtable *hash );

#endif
