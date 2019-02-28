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

#include "substmodel.h"
#include "parameters.h"
#include "mstring.h"

#define SITEMODEL_ALPHA_MIN 0.001
#define SITEMODEL_ALPHA_MAX 100

typedef enum sitemodel_heterogeneity {
	SITEMODEL_NONE,
	SITEMODEL_DISCRETE,
	SITEMODEL_GAMMA,
	SITEMODEL_GAMMA_MEAN,
	SITEMODEL_INV,
	SITEMODEL_GAMMAINV,
	SITEMODEL_GAMMAINV_MEAN,
	SITEMODEL_GAMMA_LAGUERRE,
	SITEMODEL_CAT,
	SITEMODEL_UNRESTRICTED
}sitemodel_heterogeneity;

typedef struct SiteModel{
	SubstitutionModel *m;
	int nstate;
	
	sitemodel_heterogeneity type;
	
	bool need_update;
    
    void     (*set_rate)( struct SiteModel *, const int, const double );
    
    double   (*get_rate)( struct SiteModel *, const int );
	double   (*get_proportion)( struct SiteModel *, const int );
	double * (*get_proportions)( struct SiteModel * );
	
	unsigned cat_count;
	double *cat_rates;
	double *cat_proportions;
	
	bool integrate;
	
	// categories
	Parameters *rates;
    
    // CAT
	unsigned *cat_assignment; //map pattern index to rate index
    size_t cat_assignment_length; //map pattern index to rate index
    
    Parameter *mu;
	
	// HMM?
} SiteModel;

Model * new_SiteModel2( const char* name, SiteModel *sm, Model *substmodel );

Model* new_SiteModel_from_json(json_node*node, Hashtable*hash);

#pragma mark -
// MARK: Gamma SiteModel


SiteModel * new_SiteModel( SubstitutionModel *m );

SiteModel * new_PinvSiteModel( SubstitutionModel *m, const double pinv );

SiteModel * new_PinvSiteModel_with_parameters( SubstitutionModel *m, const Parameters* pinv );

SiteModel * new_GammaSiteModel( SubstitutionModel *m, const double shape, const unsigned int cat_count );

SiteModel * new_GammaSiteModel_with_parameters( SubstitutionModel *m, const Parameters* shape, const size_t cat_count );

SiteModel * new_GammaPinvSiteModel( SubstitutionModel *m, const double pinv, const double shape, const unsigned int cat_count );

SiteModel * new_GammaPinvSiteModel_with_parameters( SubstitutionModel *m, const Parameters* params, const size_t cat_count );

SiteModel * new_GammaLaguerreSiteModel( SubstitutionModel *m, const double shape, const unsigned int cat_count );

SiteModel * new_GammaLaguerreSiteModel_with_parameters( SubstitutionModel *m, const Parameters* shape, const size_t cat_count );


void SiteModel_set_mu( SiteModel *sm, double mu );
double SiteModel_mu( const SiteModel *sm );

#pragma mark -
// MARK: Discrete SiteModel

SiteModel * new_DiscreteSiteModel( SubstitutionModel *m, int cat_count );

void sitemodel_set_discrete( SiteModel *sm, const int index, const double value );
double sitemodel_get_discrete( SiteModel *sm, const int index );

#pragma mark -

void free_SiteModel( SiteModel *sm );

SiteModel * clone_SiteModel( const SiteModel *sm );

SiteModel * clone_SiteModel_with( const SiteModel *sm, SubstitutionModel* m );

SiteModel * clone_SiteModel_with_parameters( const SiteModel *sm, SubstitutionModel* m, const Parameters* params, Parameter* mu );

#endif
