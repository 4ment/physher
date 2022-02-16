/*
 *  sitemodel.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 12/13/10.
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
	Simplex* proportions;
    Parameter *mu;
} SiteModel;

Model * new_SiteModel2( const char* name, SiteModel *sm, Model* proportions );

Model* new_SiteModel_from_json(json_node*node, Hashtable*hash);

#pragma mark -

void free_SiteModel( SiteModel *sm );

SiteModel * clone_SiteModel( const SiteModel *sm );

SiteModel * clone_SiteModel_with( const SiteModel *sm );

SiteModel * clone_SiteModel_with_parameters( const SiteModel *sm, Simplex* props, const Parameters* params, Parameter* mu );

#endif
