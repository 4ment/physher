// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_spropt_h
#define PhyC_spropt_h

#include "topologyopt.h"

double spr_optimize_bl_parsimony( struct TopologyOptimizer * opt );

double spr_optimize_bl_parsimony_only( struct TopologyOptimizer * opt );

double spr_optimize_bl_parsimony_only_openmp( struct TopologyOptimizer * opt );

double spr_optimize_bl_openmp( struct TopologyOptimizer * opt );

#endif
