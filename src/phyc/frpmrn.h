// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _FRPMRN_H_
#define _FRPMRN_H_

#include "optimizer.h"
#include "parameters.h"


//typedef enum cg_algorithm {FLETCHER_REEVES, POLAK_RIBIERE, BEALE_SORENSON_HESTENES_STIEFEL} cg_algorithm;

opt_result frprmn_optimize( Parameters *x, opt_func f, opt_grad_func grad_f, void *data, OptStopCriterion stop, double *fmin, opt_algorithm algorithm );

#endif
