// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SITE_MODEL_H_
#define _SITE_MODEL_H_

#include "parameters.h"
#include "mstring.h"
#include "distmodel.h"
#include "sitepattern.h"

#define SITEMODEL_ALPHA_MIN 0.001
#define SITEMODEL_ALPHA_MAX 100

typedef enum quadrature_t {
	QUADRATURE_BETA,
	QUADRATURE_DISCRETE,
	QUADRATURE_QUANTILE_MEDIAN,
	QUADRATURE_QUANTILE_MEAN,
	QUADRATURE_GAUSS_LAGUERRE,
	QUADRATURE_KUMARASWAMY
}quadrature_t;

typedef struct SiteModel{
	SitePattern* sp;
	
	distribution_t distribution; // parametric distribution
	bool invariant;
	quadrature_t quadrature;
	
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

SiteModel * new_SiteModel_with_parameters( const Parameters *params, Parameter* proportions, const size_t cat_count, distribution_t distribution, bool invariant, quadrature_t quad);

void free_SiteModel( SiteModel *sm );

SiteModel * clone_SiteModel( const SiteModel *sm );

SiteModel * clone_SiteModel_with( const SiteModel *sm );

SiteModel * clone_SiteModel_with_parameters( const SiteModel *sm, Parameter* props, const Parameters* params, Parameter* mu );

void SiteModel_set_mu(SiteModel *sm, Parameter* mu);

#endif
