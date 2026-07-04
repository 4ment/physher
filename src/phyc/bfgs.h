// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _BFGS_H_
#define _BFGS_H_

#include "parameters.h"
#include "optimizer.h"

opt_result dfpmin_optimize( Parameters *p, opt_func f, opt_grad_func grad_f, void *data, OptStopCriterion stop, double *fmin, double alpha);

#endif
