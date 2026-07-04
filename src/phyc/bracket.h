// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _BRACKET_H_
#define _BRACKET_H_

#include "parameters.h"
#include "optimizer.h"

void bracket2( Parameters *ps, double *ax, double *bx, double *cx, opt_func f, void *data );

#endif
