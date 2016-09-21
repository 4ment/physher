/*
 *  nj.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/08/13.
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

#ifndef PhyC_nj_h
#define PhyC_nj_h

#include "sequence.h"
#include "tree.h"


Tree * new_NJ_float( const Sequences *sequences, float **matrix );


Tree * new_NJ( const Sequences *sequences, double **matrix );

Tree * new_NJ2( const Sequences *sequences, double **_matrix );

#endif
