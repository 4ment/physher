// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_asr_h
#define PhyC_asr_h

#include "treelikelihood.h"

void asr_marginal( SingleTreeLikelihood *tlk );

void asr_calculator_from_json(json_node* node, Hashtable* hash);

#endif
