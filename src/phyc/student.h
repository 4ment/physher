/*
 *  student.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 10/20/11.
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

#ifndef _STUDENT_DISTRIBUTION_H_
#define _STUDENT_DISTRIBUTION_H_

/**
 * Returns the quantile for the two-tailed student's t distribution.
 * This is an implementation of the algorithm in
 * G. W. Hill. "Algorithm 396: Student's t-Quantiles." Communications
 * of the ACM 13(10):619--620.  ACM Press, October, 1970.
 */

double qt(double p, long n);

#endif
