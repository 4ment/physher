// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef __PhyC__upgma__
#define __PhyC__upgma__

#include <stdio.h>

#include "tree.h"
#include "sequence.h"

struct _Tree * new_UPGMA( const char **taxa, size_t dim, double **matrix );

struct _Tree * new_UPGMA_float( const Sequences *sequences, float **_matrix );

struct _Tree* create_UPGMA_from_json( json_node* node, Hashtable* hash );

#endif /* defined(__PhyC__upgma__) */
