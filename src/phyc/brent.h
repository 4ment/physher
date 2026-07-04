// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _BRENT_H_
#define _BRENT_H_

#include "optimizer.h"

opt_result serial_brent_optimize( Parameters *ps, opt_func f, void *data, OptStopCriterion *stop, double *fmin );

opt_result brent_optimize( Parameters *ps, opt_func f, void *data, OptStopCriterion *stop, double *fmin );

opt_result brent_optimize2( Parameter *parameter, size_t index, opt_func f, void *data, OptStopCriterion *stop, double *fminp );

#endif
