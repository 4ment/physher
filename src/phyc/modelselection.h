// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _MODEL_SELECTION_H_
#define _MODEL_SELECTION_H_

typedef enum information_criterion{
	INFORMATION_CRITERION_AIC,
	INFORMATION_CRITERION_AICc,
	INFORMATION_CRITERION_BIC,
	INFORMATION_CRITERION_HQ
}information_criterion;

static const char * const INFORMATION_CRITERION[4] = {
	"AIC",
	"AICc",
	"BIC",
	"HQ"
};


double AIC( const double lk, const int k );

double AICc( const double lk, const int k, const int n );

double BIC( const double lk, const int k, const int n );

double HQ( const double lk, const int k, const int n );

double LRT( const double lk0, const double lk1, const unsigned n0, const unsigned n1  );

#endif
