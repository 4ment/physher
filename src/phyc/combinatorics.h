/*
 *  combinatorics.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 23/8/12.
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

#ifndef PhyC_combinatorics_h
#define PhyC_combinatorics_h

double bell_number( unsigned int n );

int choose ( unsigned int n, unsigned int k );

void combination_at_index( unsigned int n,  unsigned int k, long m, unsigned *ans );

#endif
