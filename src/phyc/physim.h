// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_physim_h
#define PhyC_physim_h

#include "sequence.h"
#include "branchmodel.h"
#include "sitemodel.h"
#include "substmodel.h"
#include "tree.h"

#define JSON_SIMULTRON "simultron"

Sequences * Sequence_simulate( Tree *tree, SubstitutionModel *m, SiteModel *sm, BranchModel *bm, DataType *datatype, unsigned len, bool keep_internal );

void SimulateSequences_from_json(json_node* node, Hashtable* hash);

#endif
