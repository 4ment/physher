// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _POWELL_H_
#define _POWELL_H_

#include "parameters.h"
#include "optimizer.h"

typedef struct Powell{
	double **xi;
	int n;
}Powell;

Powell * new_Powell( const int n );

void free_Powell( Powell *powell);



opt_result powell_optimize( Parameters *p, opt_func f, void *data, OptStopCriterion *stop, double *fmin, opt_update_data uf );



void powell_reset_xi(Powell *powell);


#endif
