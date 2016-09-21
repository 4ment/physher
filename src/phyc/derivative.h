/*
 *  derivative.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 6/20/11.
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



#ifndef _DERIVATIVE_H_
#define _DERIVATIVE_H_

#include "optimizer.h"
#include "parameters.h"

double dfridr( Parameters *x, opt_func func, void *data, double h, int index, double *err, bool quick );

double first_derivative( Parameters *x, opt_func func, void *data );

#endif
