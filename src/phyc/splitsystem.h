/*
 *  splitsystem.h
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

#ifndef __PhyC__splitsystem__
#define __PhyC__splitsystem__

#include <stdio.h>

#include "hashtable.h"
#include "tree.h"

bool ** getSplitSystem( Hashtable *hash, Tree *tree );

bool ** getSplitSystemUnrooted( Hashtable *hash, Tree *tree );

bool ** getSplitSystemAll( Hashtable *hash, Tree *tree );

bool hasSplit( bool **splits, bool *split, int splitCount, int length );

int hasSplit2( bool **splits, bool *split, int splitCount, int length );

#endif /* defined(__PhyC__splitsystem__) */
