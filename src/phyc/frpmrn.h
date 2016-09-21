/*
 *  frpmrn.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/15/10.
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

#ifndef _FRPMRN_H_
#define _FRPMRN_H_

#include "optimizer.h"
#include "parameters.h"


//typedef enum cg_algorithm {FLETCHER_REEVES, POLAK_RIBIERE, BEALE_SORENSON_HESTENES_STIEFEL} cg_algorithm;

opt_result frprmn_optimize( Parameters *x, opt_func f, void *data, OptStopCriterion stop, double *fmin, opt_algorithm algorithm );

#endif
