// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef PhyC_nj_h
#define PhyC_nj_h

#include "sequence.h"
#include "tree.h"

struct _Tree * new_NJ( const char **taxa, size_t dim, double **matrix, Parameter* branchLengths );

struct _Tree* create_NJ_from_json( json_node* node, Hashtable* hash, Parameter* branchLengths );

#endif
