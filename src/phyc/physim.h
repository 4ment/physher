/*
 *  physim.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 9/09/13.
 *  Copyright (C) 2010 Mathieu Fourment. All rights reserved.
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


#ifndef PhyC_physim_h
#define PhyC_physim_h

#include "sequence.h"
#include "branchmodel.h"
#include "sitemodel.h"
#include "tree.h"

#define JSON_SIMULTRON "simultron"

Sequences * Sequence_simulate( Tree *tree, SiteModel *sm, BranchModel *bm, DataType *datatype, unsigned len, bool keep_internal );

void SimulateSequences_from_json(json_node* node, Hashtable* hash);

#endif
