/*
 *  upgma.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 20/01/2015.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef __PhyC__upgma__
#define __PhyC__upgma__

#include <stdio.h>

#include "tree.h"
#include "sequence.h"

struct _Tree * new_UPGMA( const char **taxa, size_t dim, double **matrix );

struct _Tree * new_UPGMA_float( const Sequences *sequences, float **_matrix );

struct _Tree* create_UPGMA_from_json( json_node* node, Hashtable* hash );

#endif /* defined(__PhyC__upgma__) */
