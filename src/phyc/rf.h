/*
 *  rf.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 14/6/12.
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


#ifndef Math_rf_h
#define Math_rf_h


double rf_distance( bool **s1, bool **s2, int splitCount, int length );

double rf_norm_distance( bool **s1, bool **s2, int splitCount, int length );

double Branch_score( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length );

double K_tree_score( bool **s1, bool **s2, Tree *t1, Tree *t2, int splitCount, int length );

#endif
