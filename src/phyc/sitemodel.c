/*
 *  sitemodel.c
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

#include "sitemodel.h"

#include <assert.h>

#include "substmodel.h"
#include "parameters.h"
#include "gamma.h"
#include "gausslaguerre.h"
#include "mstring.h"
#include "matrix.h"
#include "mathconstant.h"
#include "gaussian.h"

#include "model.h"

static void _gamma_approx_quantile( SiteModel *sm );
static void _invariant( SiteModel *sm );
static void _calculate_rates_discrete( SiteModel *sm );

static double _get_rate( SiteModel *sm, const int index );
static double _get_proportion( SiteModel *sm, const int index );
static double *_get_proportions( SiteModel *sm );

static double _get_rate_gamma( SiteModel *sm, const int index );
static double _get_proportion_gamma( SiteModel *sm, const int index );
static double *_get_proportions_gamma( SiteModel *sm );

static double _get_rate_invariant( SiteModel *sm, const int index );
static double _get_proportion_invariant( SiteModel *sm, const int index );
static double *_get_proportions_invariant( SiteModel *sm );

static double _get_rate_laguerre( SiteModel *sm, const int index );
static double _get_proportion_laguerre( SiteModel *sm, const int index );
static double * _get_proportions_laguerre( SiteModel *sm );


static double _get_rate_discrete( SiteModel *sm, const int index );
static double _get_proportion_discrete( SiteModel *sm, const int index );
static double *_get_proportions_discrete( SiteModel *sm );

static double _get_rate_cat( SiteModel *sm, const int index );

static SiteModel * _new_SiteModel_with_parameters( SubstitutionModel *m, Parameters *params, const unsigned int cat_count, sitemodel_heterogeneity type );

#ifdef LISTENERS
static void _site_model_handle_change( Model *self, Model *model, int index ){
	SiteModel *sm = (SiteModel*)self->obj;
	if( model == NULL )sm->need_update = true;
	ListenerList_fire( model->listeners, model, index );
}

// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SiteModel2( SiteModel *sm, Model *substmodel ){
	Model *model = new_Model("substmodel", sm, 1); // at least the sitemodel will listen to it

	if ( sm->rates != NULL ) {
		for ( int i = 0; i < Parameters_count(sm->rates); i++ ) {
			ListenerList_add( Parameters_at(sm->rates, i)->listeners, model );
		}
	}
	ListenerList_add( substmodel->listeners, model );
	
	model->update = _site_model_handle_change;
	return model;
}
#endif

#pragma mark -
// MARK: Gamma SiteModel

SiteModel * new_SiteModel( SubstitutionModel *m ){
	return _new_SiteModel_with_parameters( m, NULL, 0, SITEMODEL_NONE );
}

SiteModel * new_PinvSiteModel( SubstitutionModel *m, const double pinv ){
	assert( pinv>0 && pinv <1 );
	Parameters *params = new_Parameters(1);
	Parameter *inv = new_Parameter("sitemodel.pinv", pinv, new_Constraint(0+TINY, 1-TINY) );
	Parameters_add(params, inv);
	return _new_SiteModel_with_parameters( m, params, 0, SITEMODEL_INV );
}

SiteModel * new_GammaSiteModel( SubstitutionModel *m, const double shape, const unsigned int cat_count ){
	assert( shape >= SITEMODEL_ALPHA_MIN && shape <= SITEMODEL_ALPHA_MAX );
	Parameters *params = new_Parameters(1);
	Parameter *alpha = new_Parameter("sitemodel.alpha", shape, new_Constraint(SITEMODEL_ALPHA_MIN, SITEMODEL_ALPHA_MAX) );
	Parameters_add(params, alpha);
	return _new_SiteModel_with_parameters( m, params, cat_count, SITEMODEL_GAMMA );
}

SiteModel * new_GammaPinvSiteModel( SubstitutionModel *m, const double pinv, const double shape, const unsigned int cat_count ){

    if( pinv == 0 && cat_count == 1 ){
        return _new_SiteModel_with_parameters( m, NULL, 0, SITEMODEL_NONE );
    }
    
    Parameters *params = new_Parameters(2);
    
    if( cat_count > 1 ){
        Parameter *alpha = new_Parameter("sitemodel.alpha", shape, new_Constraint(SITEMODEL_ALPHA_MIN, SITEMODEL_ALPHA_MAX) );
        Parameters_add(params, alpha);
    }
    if( pinv > 0 ){
        Parameter *inv = new_Parameter("sitemodel.pinv", pinv, new_Constraint(0+TINY, 1-TINY) );
        Parameters_add(params, inv);
    }
    
    sitemodel_heterogeneity type = SITEMODEL_GAMMAINV;
    if( pinv == 0 && cat_count > 1 ){
        type = SITEMODEL_GAMMA;
    }
    else if ( pinv > 0.0 && cat_count == 1 ){
        type = SITEMODEL_INV;
    }
	return _new_SiteModel_with_parameters( m, params, cat_count, type );
}

SiteModel * new_GammaLaguerreSiteModel( SubstitutionModel *m, const double shape, const unsigned int cat_count ){
	assert( shape >= SITEMODEL_ALPHA_MIN && shape <= SITEMODEL_ALPHA_MAX );
	Parameters *params = new_Parameters(1);
	Parameter *alpha = new_Parameter("sitemodel.alpha", shape, new_Constraint(SITEMODEL_ALPHA_MIN, SITEMODEL_ALPHA_MAX) );
	Parameters_add(params, alpha);
	return _new_SiteModel_with_parameters( m, params, cat_count, SITEMODEL_GAMMA_LAGUERRE );
}

static void _set_rate(SiteModel* sm, const int index, const double value){
    Parameters_set_value(sm->rates, index, value);
    sm->need_update = true;
}

// (Gamma or/and Invariant) or one rate
// should not be used directly
SiteModel * _new_SiteModel_with_parameters( SubstitutionModel *m, Parameters *params, const unsigned int cat_count, sitemodel_heterogeneity type ){
	SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	sm->id = 0;
	sm->m = m;
	sm->nstate = m->nstate;
	
	sm->type = type;
	
	sm->cat_count = 1;
	
	sm->cat_rates       = NULL;
	sm->cat_proportions = NULL;
    
	sm->rates = NULL;
	sm->mu    = NULL;
	
    sm->cat_assignment = NULL;
    sm->cat_assignment_length = 0;
    
    sm->set_rate = _set_rate;
    
	switch ( type) {
		case SITEMODEL_NONE:			
			sm->get_rate        = _get_rate;
			sm->get_proportion  = _get_proportion;
			sm->get_proportions = _get_proportions;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
            sm->cat_rates[0] = sm->cat_proportions[0] = 1;
			break;
		case SITEMODEL_GAMMA:
            sm->rates = new_Parameters(1);
            Parameters_set(sm->rates, 0, Parameters_at(params, 0));
			sm->cat_count       = cat_count;
			sm->get_rate        = _get_rate_gamma;
			sm->get_proportion  = _get_proportion_gamma;
			sm->get_proportions = _get_proportions_gamma;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
			break;
		case SITEMODEL_GAMMAINV:
            sm->rates = new_Parameters(2);
            Parameters_set(sm->rates, 0, Parameters_at(params, 0));
            Parameters_set(sm->rates, 1, Parameters_at(params, 1));
            sm->cat_count       = cat_count+1;
			sm->get_rate        = _get_rate_gamma;
			sm->get_proportion  = _get_proportion_gamma;
			sm->get_proportions = _get_proportions_gamma;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
			break;
        case SITEMODEL_INV:
            sm->rates = new_Parameters(1);
            Parameters_set(sm->rates, 0, Parameters_at(params, 0));
			sm->cat_count       = 2;
			sm->get_rate        = _get_rate_invariant;
			sm->get_proportion  = _get_proportion_invariant;
			sm->get_proportions = _get_proportions_invariant;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
			break;
        case SITEMODEL_GAMMA_LAGUERRE:
            sm->rates = new_Parameters(1);
            Parameters_set(sm->rates, 0, Parameters_at(params, 0));
			sm->cat_count       = cat_count;
			sm->get_rate        = _get_rate_laguerre;
			sm->get_proportion  = _get_proportion_laguerre;
			sm->get_proportions = _get_proportions_laguerre;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
			break;
		default:
			error("SiteModel type unknown");
			break;
	}
	
	if ( params != NULL ) {
		free_Parameters_soft(params);
	}
    
    
	sm->integrate   = true;
	sm->need_update = true;
	
	return sm;
}



double _get_rate( SiteModel *sm, const int index ){
	return (sm->mu == NULL ? 1.0 : Parameter_value(sm->mu));
}

double _get_proportion( SiteModel *sm, const int index ){
	return 1;
}

double *_get_proportions( SiteModel *sm ){
	return sm->cat_proportions;
}

double _get_rate_invariant( SiteModel *sm, const int index ){
    if ( sm->need_update ) {
        _invariant(sm);
    }
    return sm->cat_rates[index] * (sm->mu == NULL ? 1.0 : Parameter_value(sm->mu));
}

double _get_rate_gamma( SiteModel *sm, const int index ){
    if ( sm->need_update ) {
        _gamma_approx_quantile(sm);
    }
    return sm->cat_rates[index] * (sm->mu == NULL ? 1.0 : Parameter_value(sm->mu));
}

double _get_proportion_gamma( SiteModel *sm, const int index ){
    if ( sm->need_update ) {
        _gamma_approx_quantile(sm);
    }
    return sm->cat_proportions[index];
}

double *_get_proportions_gamma( SiteModel *sm ){
    if ( sm->need_update ) {
        _gamma_approx_quantile(sm);
    }
    return sm->cat_proportions;
}

double _get_proportion_invariant( SiteModel *sm, const int index ){
    if ( sm->need_update ) {
        _invariant(sm);
    }
    return sm->cat_proportions[index];
}

double *_get_proportions_invariant( SiteModel *sm ){
    if ( sm->need_update ) {
        _invariant(sm);
    }
    return sm->cat_proportions;
}

void _gamma_approx_quantile( SiteModel *sm ) {
    
    double propVariable = 1.0;
    int cat = 0;
    
    // there is also invariant sites
    if ( Parameters_count(sm->rates) > 1 ) {
        sm->cat_rates[0] = 0.0;
        sm->cat_proportions[0] = Parameters_value(sm->rates, 1);
        
        propVariable = 1.0 - sm->cat_proportions[0];
        cat = 1;
    }
    
    double mean = 0.0;
    const int nCat = sm->cat_count - cat;
    
    const double alpha = Parameters_value(sm->rates, 0);
    for (int i = 0; i < nCat; i++) {
        sm->cat_rates[i + cat] = qgamma( (2.0 * i + 1.0) / (2.0 * nCat),alpha, 1.0 / alpha );
        //sm->cat_rates[i + cat] = qnorm( (2.0 * i + 1.0) / (2.0 * nCat), 1, alpha );
        mean += sm->cat_rates[i + cat];
        
        sm->cat_proportions[i + cat] = propVariable / nCat;
    }
    
    mean = (propVariable * mean) / nCat;
    
    for (int i = 0; i < nCat; i++) {
        sm->cat_rates[i + cat] /= mean;
    }
    
    sm->need_update = false;
}

void _invariant( SiteModel *sm ) {
    
    sm->cat_rates[0] = 0.0;
    sm->cat_proportions[0] = Parameters_value(sm->rates, 0);
    double propVariable = propVariable = 1.0 - sm->cat_proportions[0];
    
    sm->cat_rates[1] = 1.0 / propVariable;
    sm->cat_proportions[1] = propVariable;
    sm->need_update = false;
}

// Gamma distribution approximated using Laguerre quadrature
void _gamma_approx_laguerre( SiteModel *sm ){
    const double alpha = Parameters_value(sm->rates, 0);
    
	// calculate using alpha -1
	gaulag(sm->cat_rates, sm->cat_proportions, sm->cat_count, alpha-1);
	
	double sum = 0.0;
	double gamalpha = gamm(alpha+1);
	for ( int i = 0; i < sm->cat_count; i++ ) {
		sm->cat_rates[i] /= alpha+1;
		sm->cat_proportions[i] *= gamalpha;
		sum += sm->cat_proportions[i];
	}
	
	for ( int i = 0; i < sm->cat_count; i++ ) {
		sm->cat_proportions[i] /= sum;
	}
	
//	fprintf(stderr, "+++++++++++++++++++++++++++++\n");
//	print_dvector(sm->cat_proportions, sm->cat_count);
//	print_dvector(sm->cat_rates, sm->cat_count);
	sm->need_update = false;
}

double _get_rate_laguerre( SiteModel *sm, const int index ){
	if ( sm->need_update ) {
		_gamma_approx_laguerre(sm);		
	}
	return sm->cat_rates[index];
}

double _get_proportion_laguerre( SiteModel *sm, const int index ){
	if ( sm->need_update ) {
		_gamma_approx_laguerre(sm);
	}
	return sm->cat_proportions[index];
}

double *_get_proportions_laguerre( SiteModel *sm ){
	if ( sm->need_update ) {
		_gamma_approx_laguerre(sm);
	}
	return sm->cat_proportions;
}

void SiteModel_set_mu( SiteModel *sm, double mu ){
    sm->mu->value = mu;
	sm->need_update = true;
}

double SiteModel_mu( const SiteModel *sm ){
    return sm->mu->value;
}


#pragma mark -
// MARK: Discrete SiteModel

SiteModel * new_DiscreteSiteModel( SubstitutionModel *m, int cat_count ){
	SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	sm->id = 0;
	sm->m = m;
	sm->nstate = m->nstate;
	
    sm->set_rate        = _set_rate;
	sm->get_rate        = _get_rate_discrete;
	sm->get_proportion  = _get_proportion_discrete;
	sm->get_proportions = _get_proportions_discrete;
	
	sm->mu    = NULL;
	sm->cat_count = cat_count;
	
	sm->cat_rates       = dvector(sm->cat_count);
	sm->cat_proportions = dvector(sm->cat_count);
    sm->cat_assignment = NULL;
	
	sm->rates = new_Parameters( (cat_count * 2) -2 );
	
	StringBuffer *buffer = new_StringBuffer(20);
	double prop = 1.0/cat_count;
	for ( int i = 0; i < sm->cat_count; i++ ) {
		sm->cat_proportions[i] = prop;
	}
	
    StringBuffer_set_string(buffer, "sitemodel.p0");
	Parameters_set(sm->rates, 0, new_Parameter(buffer->c, prop, new_Constraint(0.00001, 0.999)) );
    int i = 1;
	for ( ; i < cat_count-1; i++ ) {
		StringBuffer_set_string(buffer, "sitemodel.p");
		StringBuffer_append_format(buffer, "%d", i);
		double denom = 1;
		for ( int j = 0; j < i; j++ ) {
			denom *= 1 - Parameters_value(sm->rates, j);
		}
		double r = sm->cat_proportions[i] / denom;
		Parameters_set(sm->rates, i, new_Parameter(buffer->c, r, new_Constraint(0.00001, 0.999)) );
	}
	
    for (int j = 0 ; j < cat_count-1; j++, i++ ) {
		StringBuffer_set_string(buffer, "sitemodel.r");
		StringBuffer_append_format(buffer, "%d", j);
		Parameters_set(sm->rates, i, new_Parameter(buffer->c, 1.0, new_Constraint(0.00001, 10)) );
	}

	free_StringBuffer(buffer);
	
	sm->integrate   = true;
	sm->need_update = true;
	
	sm->type = SITEMODEL_DISCRETE;
	return sm;
}


double _get_rate_discrete( SiteModel *sm, const int index ){
	if ( sm->need_update ) {
		_calculate_rates_discrete(sm);
	}
	return sm->cat_rates[index];
}

double _get_proportion_discrete( SiteModel *sm, const int index ){
	if ( sm->need_update ) {
		_calculate_rates_discrete(sm);
	}
	return sm->cat_proportions[index];
}

double *_get_proportions_discrete( SiteModel *sm ){
	if ( sm->need_update ) {
		_calculate_rates_discrete(sm);
	}
	return sm->cat_proportions;
}

void _calculate_rates_discrete( SiteModel *sm ) {
	// proportions
    sm->cat_proportions[0] = Parameters_value(sm->rates, 0);
    double p;
	for ( int i = 1; i < sm->cat_count; i++ ) {
		p = 1;
		for ( int j = 0; j < i; j++ ) {
			p *= 1 - Parameters_value(sm->rates, j);
		}
		if ( i != sm->cat_count-1 ) p *= Parameters_value(sm->rates, i);
		sm->cat_proportions[i] = p;
	}
    
    // rates
    int j = 0;
    int i = sm->cat_count-1;
    sm->cat_rates[0] = Parameters_value(sm->rates, i++);
    sm->cat_rates[1] = 1;
    
    double scaler = sm->cat_rates[0]*sm->cat_proportions[0] +  sm->cat_proportions[1];
    
    for ( j = 2; i < Parameters_count(sm->rates); i++,j++ ) {
        sm->cat_rates[j] = Parameters_value(sm->rates, i)*sm->cat_rates[j-1];
        scaler += sm->cat_rates[j]*sm->cat_proportions[j];
    }
    
    for ( i = 0; i < sm->cat_count; i++ ) {
        sm->cat_rates[i] /= scaler;
    }
//    int i = 0;
//    for ( int j = sm->cat_count-1; j < Parameters_count(sm->rates); i++,j++ ) {
//        sm->cat_rates[i] = Parameters_value(sm->rates, j);
//    }
	sm->need_update = false;
//    print_dvector(sm->cat_rates, sm->cat_count);
//    print_dvector(sm->cat_proportions, sm->cat_count);
//    Parameters_print(sm->rates);
//    printf("scaler %f\n", scaler);
}

void sitemodel_set_discrete( SiteModel *sm, const int index, const double value ){
	Parameters_set_value(sm->rates, index, value);
	sm->need_update = true;	
}

double sitemodel_get_discrete( SiteModel *sm, const int index ){
	return Parameters_value(sm->rates, index);
}

#pragma mark -
// MARK: CAT SiteModel

static double _get_rate_cat( SiteModel *sm, const int index );

SiteModel * new_CATSiteModel( SubstitutionModel *m, size_t cat_count, const unsigned npatterns ){
    
    if( cat_count == 1 ){
        return new_SiteModel(m);
    }
    SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	sm->id = 0;
	sm->m = m;
	sm->nstate = m->nstate;
	
	sm->type = SITEMODEL_CAT;
	
	sm->cat_count = cat_count;
	
    sm->cat_assignment = uivector(npatterns);
    sm->cat_assignment_length = npatterns;
    
	sm->cat_rates       = dvector(cat_count);
	sm->cat_proportions = NULL;
    
	sm->rates = NULL;
	sm->mu    = NULL;
    
    sm->get_rate = _get_rate_cat;
    
    sm->get_proportion  = NULL;
    sm->get_proportions = NULL;
    
    return sm;
}

// index is the index of a site/pattern
double _get_rate_cat( SiteModel *sm, const int index ){
	return sm->cat_rates[ sm->cat_assignment[index]];
}

#pragma mark -
// MARK: Unrestricted SiteModel

static double _get_rate_unres( SiteModel *sm, const int index );
static void _set_discrete_unres( SiteModel *sm, const int index, const double value );

SiteModel * new_UnrestrictedSiteModel( SubstitutionModel *m, size_t npatterns ){
    
    SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	sm->id = 0;
	sm->m = m;
	sm->nstate = m->nstate;
	
	sm->type = SITEMODEL_UNRESTRICTED;
	
	sm->cat_count = npatterns;
	
    sm->cat_assignment = NULL;
    sm->cat_assignment_length = 0;
    
	sm->cat_rates       = dvector(npatterns);
	sm->cat_proportions = NULL;
    
	sm->rates = new_Parameters(npatterns);
	sm->mu    = NULL;
    
    StringBuffer *buffer = new_StringBuffer(20);
    for ( int i = 0; i < npatterns; i++ ) {
		StringBuffer_set_string(buffer, "sitemodel.r");
		StringBuffer_append_format(buffer, "%d", i);
		Parameters_add(sm->rates, new_Parameter(buffer->c, 1., new_Constraint(0.00001, 100)) );
	}
	free_StringBuffer(buffer);
    
    sm->get_rate = _get_rate_unres;
    
    sm->get_proportion  = NULL;
    sm->get_proportions = NULL;
    
    return sm;
}

// index is the index of a site/pattern
double _get_rate_unres( SiteModel *sm, const int index ){
	return Parameters_value(sm->rates, index);
}

#pragma mark -

SiteModel * clone_SiteModel( const SiteModel *sm ){
	SiteModel *newsm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(newsm);
    newsm->m = sm->m->clone(sm->m);
	newsm->id = sm->id;
	newsm->nstate = sm->nstate;
	
	newsm->rates = NULL;
	
    newsm->cat_count = sm->cat_count;
    
	if ( sm->rates != NULL ){
		newsm->rates = clone_Parameters(sm->rates, true);
	}
    
    newsm->mu = NULL;
	if ( sm->mu != NULL ){
		newsm->mu = clone_Parameter(sm->mu, true);
	}
    
    newsm->cat_assignment_length = sm->cat_assignment_length;
    newsm->cat_assignment = NULL;
	if ( sm->cat_assignment != NULL ){
		newsm->cat_assignment = clone_uivector((unsigned *)sm->cat_assignment, sm->cat_assignment_length);
	}
	
	newsm->cat_rates = clone_dvector(sm->cat_rates, sm->cat_count);
	
    newsm->cat_proportions = NULL;
    if ( sm->cat_proportions != NULL ){
        newsm->cat_proportions = clone_dvector(sm->cat_proportions, sm->cat_count);
	}
	
	newsm->integrate = sm->integrate;
	newsm->need_update = false;
	
    newsm->set_rate = sm->set_rate;
    
	newsm->get_rate        = sm->get_rate;
	newsm->get_proportion  = sm->get_proportion;
	newsm->get_proportions = sm->get_proportions;
	
	newsm->type = sm->type;
	
	return newsm;
}

SiteModel * clone_SiteModel_share( const SiteModel *sm, bool share_gamma, bool share_pinv, bool share_freqs, bool share_rates ){
	SiteModel *newsm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(newsm);
	newsm->m = clone_substitution_model_share(sm->m, share_freqs, share_rates);
	newsm->id = sm->id;
	newsm->nstate = sm->nstate;
    
	if ( sm->rates != NULL ){
        newsm->rates = clone_Parameters(sm->rates, true);		
	}
    
    newsm->mu = NULL;
	if ( sm->mu != NULL ){
		newsm->mu = clone_Parameter(sm->mu, true);
	}
	
	newsm->cat_count = sm->cat_count;
	
	newsm->cat_rates = clone_dvector(sm->cat_rates, sm->cat_count);
	
	newsm->cat_proportions = clone_dvector(sm->cat_proportions, sm->cat_count);
	
	
	newsm->integrate = sm->integrate;
	newsm->need_update = false;
	
    newsm->set_rate = sm->set_rate;
    
	newsm->get_rate        = sm->get_rate;
	newsm->get_proportion  = sm->get_proportion;
	newsm->get_proportions = sm->get_proportions;
	
	newsm->type = sm->type;
	
	return newsm;
}


void free_SiteModel( SiteModel *sm ){
	if ( sm->rates != NULL ) free_Parameters(sm->rates);
	if ( sm->mu != NULL ) free_Parameter(sm->mu);
	if ( sm->cat_proportions != NULL ) free(sm->cat_proportions);
    if ( sm->cat_assignment != NULL ) free(sm->cat_assignment);
	free(sm->cat_rates);
	sm->m->free(sm->m);
	free(sm);
}

void free_SiteModel_share( SiteModel *sm, bool share_gamma, bool share_pinv, bool share_freqs, bool share_rates ){
	if ( sm->rates != NULL ) free_Parameters(sm->rates);
	if ( sm->mu != NULL ) free_Parameter(sm->mu);
	if ( sm->cat_proportions != NULL ) free(sm->cat_proportions);
    if ( sm->cat_assignment != NULL ) free(sm->cat_assignment);
	free(sm->cat_rates);
    free_SubstitutionModel_share(sm->m, share_freqs, share_rates);
	free(sm);
}

/*char * SiteModel_stringify( SiteModel *sm ){
	StringBuffer *buffer = new_StringBuffer(100);
	
	buffer = SiteModel_bufferize( buffer, sm );
	
	char *final = StringBuffer_tochar(buffer);
	free_StringBuffer(buffer);
	//fprintf(stderr, "%d\n", (int)strlen(final));
	return final;
}

StringBuffer * SiteModel_bufferize( StringBuffer *buffer, SiteModel *sm ){
	
	buffer = StringBuffer_append_format(buffer, "(SiteModel:\n (id:\"sm%d\")\n",sm->id);
	buffer = SubstitutionModel_bufferize( buffer, sm->m );
	buffer = StringBuffer_append_char(buffer,'\n');
	if( sm->pinv != NULL ){
		buffer = Parameter_bufferize(buffer, sm->pinv);
		buffer = StringBuffer_append_char(buffer,'\n');
	}
	
	if( sm->shape != NULL ){
		buffer = StringBuffer_append_format(buffer,"(categories: %d)\n", sm->cat_count);
		
		buffer = Parameter_bufferize(buffer, sm->shape);
		buffer = StringBuffer_append_char(buffer, '\n');
		
		buffer = StringBuffer_append_char(buffer, '\n');
	}
	buffer = StringBuffer_append_char(buffer, ')');
	
	return buffer;
}

void * SiteModel_SML_to_object( ObjectStore *store, SMLNode node ){
	fprintf(stderr, "SiteModel_SML_to_object\n");
	SiteModel *sm = NULL;
//	
//	SMLNode node_model = SML_get_element( node, "Model");
//	SubstitutionModel *m = SubstitutionModel_SML_to_object( store, node_model);
//	
//	SMLNode node_invariant = SML_get_element( node, "invariant");
//	Parameter *pinv = NULL;
//	if( node_invariant != NULL ){
//		
//	}
//	
//	
//	SMLNode node_gamma = SML_get_element( node, "gamma");
//	Parameter *shape = NULL;
//	unsigned int cat_count = 1;
//	if( node_gamma != NULL ){
//		
//	}
//	
//	sm = _new_SiteModel_with_parameters( m, pinv, shape, cat_count);
//	
//	char *id = SML_get_data_of_child( node, "id");
//	if ( id != NULL ) {
//		while( *id < '0' || *id > '9' ){
//			id++;
//		}
//		sm->id = atoi( id );
//	}
//	else sm->id = 0;
	return sm;
}*/

void compare_sitemodel( const SiteModel *sm1, const SiteModel *sm2){
	assert(sm1->id == sm2->id );
	assert( sm1->nstate == sm2->nstate );
	assert( sm1->integrate == sm2->integrate );
	
	//compare_model( sm1->m, sm2->m );
	
	assert( sm1->cat_count == sm2->cat_count );	
	assert( memcmp( sm1->cat_rates, sm2->cat_rates, sm1->cat_count *sizeof(double) ) == 0);
	assert( memcmp( sm1->cat_proportions, sm2->cat_proportions, sm1->cat_count *sizeof(double) ) == 0);
}
