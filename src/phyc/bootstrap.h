// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef bootstrap_h
#define bootstrap_h

#ifndef GSL_DISABLED
#include <gsl/gsl_rng.h>
#endif

#include "optimizer.h"
#include "parameters.h"
#include "sitepattern.h"
#include "tracelogger.h"

// Nonparametric bootstrap over alignment sites, run from the "physher" section
// after the point estimate has been obtained: the parameter values in place
// when initialize() runs are the estimate every replicate is restarted from.
//
// A replicate resamples the sites of the alignment with replacement. Every
// drawn site is a copy of an existing site, so a replicate can never contain a
// pattern absent from the original alignment: its compressed form is the
// original pattern set with new multiplicities, drawn from
// Multinomial(nsites, weights/nsites). Nothing is recompressed and no alignment
// is materialised.
typedef struct Bootstrap {
    // Distinct SitePatterns to resample, and the replicate pattern sets built
    // from them (compact mode only). Both are indexed by pattern set.
    SitePattern** patterns;  // borrowed from the tree likelihoods
    SitePattern** views;     // owned, NULL unless compact
    double** weights;        // owned, the observed multiplicities
    double** counts;         // owned, the replicate multiplicities
    unsigned int** draws;    // owned, the multinomial draw
    size_t* nsites;          // owned, sites per pattern set
    size_t pattern_count;

    // Tree likelihoods reading those pattern sets, and the pattern set each one
    // reads (two partitions may share one alignment).
    Model** treelikelihoods;
    size_t* pattern_index;
    size_t treelikelihood_count;

    Optimizer** optimizers;
    size_t optimizer_count;
    Trace** loggers;
    size_t logger_count;

    Parameters* parameters;  // restored to `estimate` before each replicate
    double* estimate;        // owned
    size_t estimate_size;

    size_t replicates;      // resampled replicates, not counting the observed one
    bool compact;           // swap a compacted pattern set vs. overwrite weights
    bool reset;             // restart each replicate from the point estimate
    bool include_observed;  // prepend a round 0 fitted on the observed weights
#ifndef GSL_DISABLED
    gsl_rng* rng;
#endif

    void (*run)(struct Bootstrap*);
    void (*free)(struct Bootstrap*);
} Bootstrap;

Bootstrap* new_Bootstrap_from_json(json_node* node, Hashtable* hash);

#endif /* bootstrap_h */
