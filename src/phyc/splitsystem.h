// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

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
