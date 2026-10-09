// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_physim_h
#define PhyC_physim_h

#include "sequence.h"
#include "branchmodel.h"
#include "sitemodel.h"
#include "substmodel.h"
#include "tree.h"

#define JSON_SIMULTRON "simulator"

// How sites are given a rate category when the site model is CAT.
typedef enum {
	// Every site of the alignment the assignment was made on is simulated in the
	// category its pattern was assigned to: the model CAT fitted, as a parametric
	// bootstrap. The length is the number of sites in that alignment.
	CAT_SIMULATION_CONDITIONAL,
	// Each site draws its category with the site-weighted assignment frequencies
	// as weights: the mixture cat_mixture scores. Any length.
	CAT_SIMULATION_MIXTURE
} cat_simulation_t;

unsigned Sequence_simulate_cat_length( const SiteModel *sm );

// cat_mode is ignored unless sm is a CAT site model.
Sequences * Sequence_simulate( Tree *tree, SubstitutionModel *m, SiteModel *sm, BranchModel *bm, DataType *datatype, unsigned len, bool keep_internal, cat_simulation_t cat_mode );

void SimulateSequences_from_json(json_node* node, Hashtable* hash);

#endif
