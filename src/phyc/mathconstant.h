/*
 *  mathconstant.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/23/1.
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

#include <float.h>


#ifndef _MATH_CONSTANT_H_
#define _MATH_CONSTANT_H_


#define GOLDEN_RATIO ((sqrt(5)+1)*0.5)

#define CGOLDEN_RATIO (2-GOLDEN_RATIO)

#define EPS DBL_EPSILON

#define SQRT_EPS 1.4901161193847656E-8

#define PI 3.141592653589793238462643383279502884197169399375

#define LOG_2PI (log(2*PI))

#define SQRT_2PI (sqrt(2*PI))


#define TINY 1.0e-25	//A small number.


#endif
