// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SITE_MODEL_H_
#define _SITE_MODEL_H_

#include "parameters.h"
#include "mstring.h"
#include "distmodel.h"
#include "sitepattern.h"

// Lower bound on the Gamma/Weibull shape. The inverse CDFs are evaluated at
// quantiles clamped to [quantile_eps, 1-quantile_eps] (see _gamma_approx_quantile),
// and the unit-scale Weibull quantile exp(ln(-ln(1-p))/shape) overflows a double
// once shape < ln(-ln(quantile_eps))/ln(DBL_MAX) = 4.7e-3. gsl_cdf_gamma_Qinv is
// no better: below shape ~ 3.2e-3 it returns +Inf for the deepest quantile when
// there are >= 16 categories, and it is already quantitatively wrong (deep-tail
// quantiles floored around 5e-9) below ~1e-2. 0.01 clears both with margin and is
// far below any shape a real alignment supports.
#define SITEMODEL_ALPHA_MIN 0.01
#define SITEMODEL_ALPHA_MAX 100

typedef enum quadrature_t {
	QUADRATURE_BETA,
	QUADRATURE_DISCRETE,
	QUADRATURE_QUANTILE_MEDIAN,
	QUADRATURE_QUANTILE_MEAN,
	QUADRATURE_GAUSS_LAGUERRE,
	QUADRATURE_KUMARASWAMY
}quadrature_t;

// How the free rate parameters of a "discrete" (free-rates) site model map onto
// the category rates. Every parameterization lands on the same constraint surface
// sum_k p_k r_k = 1; they differ in conditioning and in the induced prior.
// See docs/models/sitemodel.md.
typedef enum rate_parameterization_t {
	// Inferred from the shape of the "rates" entry: a simplex selects the rate
	// shape, a plain vector the increments, an array of two parameters the ratios.
	// The mean-contribution simplex has the same shape as the rate shape and so
	// can only be reached by naming it.
	RATE_PARAMETERIZATION_AUTO,
	// Rate shape x is a simplex, normalised by the weighted mean:
	// r_k = x_k / sum_j p_j x_j.
	RATE_PARAMETERIZATION_RATE_SHAPE,
	// Mean-contribution simplex: s_k = p_k r_k is the free simplex, so the unit
	// mean holds by construction and r_k = s_k / p_k.
	RATE_PARAMETERIZATION_MEAN_CONTRIBUTION,
	// Ordered increments: the free vector holds the gaps between consecutive
	// rates, so the raw rates are its running sum and r is increasing.
	RATE_PARAMETERIZATION_RATE_INCREMENTS,
	// Ordered ratios: the free vector holds the ratio of each rate to the next
	// one up, applied to a free top rate, so the raw rates are running products.
	RATE_PARAMETERIZATION_RATE_RATIOS
}rate_parameterization_t;

typedef struct SiteModel{
	SitePattern* sp;

	distribution_t distribution; // parametric distribution
	bool invariant;
	quadrature_t quadrature;
	rate_parameterization_t rate_parameterization;

	bool need_update;
    
    void     (*set_rate)( struct SiteModel *, const int, const double );
    
	bool     (*update)( struct SiteModel * );    
    double   (*get_rate)( struct SiteModel *, const int );
	double   (*get_proportion)( struct SiteModel *, const int );
	double * (*get_proportions)( struct SiteModel * );
	
	int (*get_site_category)( struct SiteModel *, const int );
	void (*gradient)(struct SiteModel *, const double* ingrad, double* grad);
	double (*derivative)(struct SiteModel *, const double* ingrad, Parameter* p);
	
	unsigned cat_count;
	double *cat_rates;
	double *cat_proportions;
	
	int* site_category;
	
	bool integrate;
	
	// categories
	Parameters *rates;
	Parameter* proportions;
    Parameter *mu;
	
	// for finite difference approx of gamma site model gradient
	// double epsilon;
} SiteModel;

Model * new_SiteModel2( const char* name, SiteModel *sm );

Model* new_SiteModel_from_json(json_node*node, Hashtable*hash);

#pragma mark -

SiteModel * new_SiteModel_with_parameters( const Parameters *params, Parameter* proportions, const size_t cat_count, distribution_t distribution, bool invariant, quadrature_t quad, rate_parameterization_t rate_parameterization);

void free_SiteModel( SiteModel *sm );

SiteModel * clone_SiteModel( const SiteModel *sm );

SiteModel * clone_SiteModel_with( const SiteModel *sm );

SiteModel * clone_SiteModel_with_parameters( const SiteModel *sm, Parameter* props, const Parameters* params, Parameter* mu );

void SiteModel_set_mu(SiteModel *sm, Parameter* mu);

#endif
