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
#include <strings.h>

#include "substmodel.h"
#include "parameters.h"
#include "gamma.h"
#include "gausslaguerre.h"
#include "mstring.h"
#include "matrix.h"
#include "mathconstant.h"
#include "gaussian.h"

#include "model.h"

#include <gsl/gsl_cdf.h>

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

static SiteModel * _new_SiteModel_with_parameters( SubstitutionModel *m, const Parameters *params, const size_t cat_count, sitemodel_heterogeneity type );


static void _site_model_handle_change( Model *self, Model *model, int index ){
	SiteModel *sm = (SiteModel*)self->obj;
	sm->need_update = true;// one of the sitemodel parameters
	self->listeners->fire( self->listeners, self, index );
}

static void _site_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free site model %s\n", self->name);
		SiteModel *sm = (SiteModel*)self->obj;
		Model* mm = (Model*)self->data;
		mm->free(mm); // substitution model
	//	free_SiteModel(sm);
		
		if(sm->rates != NULL) free_Parameters(sm->rates);
		//TODO: deal with this
		if(sm->mu!=NULL)free_Parameter(sm->mu);
		if ( sm->cat_proportions != NULL ) free(sm->cat_proportions);
		if ( sm->cat_assignment != NULL ) free(sm->cat_assignment);
		free(sm->cat_rates);
		free(sm);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static Model* _site_model_clone( Model *self, Hashtable* hash ){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	Model* mm = (Model*)self->data;
	Model* mmclone = NULL;
	// Susbtitution model may have been parsed already
	if (Hashtable_exists(hash, mm->name)) {
		mmclone = Hashtable_get(hash, mm->name);
		mmclone->ref_count++; // it is decremented at the end using free
	}
	else{
		mmclone = mm->clone(mm, hash);
		Hashtable_add(hash, mmclone->name, mmclone);
	}
	
	SiteModel* sm = (SiteModel*)self->obj;
	Parameters* ps = new_Parameters(1);
	for (int i = 0; i < Parameters_count(sm->rates); i++) {
		char* name = Parameters_name(sm->rates, i);
		if (Hashtable_exists(hash, name)) {
			Parameters_add(ps, Hashtable_get(hash, name));
		}
		else{
			Parameter* p = clone_Parameter(Parameters_at(sm->rates, i));
			Parameters_move(ps, p);
			Hashtable_add(hash, name, p);
		}
	}
	SiteModel* smclone = clone_SiteModel_with_parameters(sm, (SubstitutionModel*)mmclone->obj, ps);
	free_Parameters(ps);
	Model* clone = new_SiteModel2(self->name, smclone, mmclone);
	Hashtable_add(hash, clone->name, clone);
	mmclone->free(mmclone);
	return clone;
}

static void _site_model_get_free_parameters(Model* model, Parameters* parameters){
	SiteModel* m = (SiteModel*)model->obj;
	if (m->rates != NULL) {
		Parameters_add_free_parameters(parameters, m->rates);
	}
	Model* msimplex = (Model*)model->data;
	msimplex->get_free_parameters(msimplex,parameters);
}

// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SiteModel2( const char* name, SiteModel *sm, Model *substmodel ){
	Model *model = new_Model(MODEL_SITEMODEL, name, sm);

	if ( sm->rates != NULL ) {
		for ( int i = 0; i < Parameters_count(sm->rates); i++ ) {
			Parameters_at(sm->rates, i)->listeners->add( Parameters_at(sm->rates, i)->listeners, model );
		}
	}
	// Listen to substitution model
	substmodel->listeners->add( substmodel->listeners, model );
	
	model->update = _site_model_handle_change;
	model->free = _site_model_free;
	model->clone = _site_model_clone;
	model->get_free_parameters = _site_model_get_free_parameters;
	model->data = substmodel;
	substmodel->ref_count++;
	return model;
}


Model* new_SiteModel_from_json(json_node*node, Hashtable*hash){
	char* allowed[] = {
		"categories",
		"discretization",
		"distribution",
		"model",
		"mu",
		"rates",
		"substitutionmodel"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* categories_node = get_json_node(node, "categories");
	json_node* model_node = get_json_node(node, "model");
	json_node* rates_node = get_json_node(node, "rates");
	json_node* distribution_node = get_json_node(node, "distribution");
	json_node* discretization_node = get_json_node(node, "discretization");
	json_node* m_node = get_json_node(node, "substitutionmodel");
	json_node* mu_node = get_json_node(node, "mu");

	Model* mm = NULL;
	if (m_node->node_type == MJSON_STRING && ((char*)m_node->value)[0] == '&') {
		char* ref = (char*)m_node->value;
		mm = Hashtable_get(hash, ref+1);
		mm->ref_count++;
	}
	else{
		mm = new_SubstitutionModel_from_json(m_node, hash);
	}
	SubstitutionModel* m = (SubstitutionModel*)mm->obj;
	
	size_t cat = 1;
	if (categories_node != NULL) {
		cat = atoi((char*)categories_node->value);
	}

	bool invariant = false;
	Parameters* rates = NULL;
	if (rates_node != NULL && rates_node->node_type == MJSON_ARRAY) {
		rates = new_Parameters_from_json(rates_node, hash);
		check_constraints(rates, 0, INFINITY, 0.001, 100);
	}
	else if (rates_node != NULL && rates_node->node_type == MJSON_OBJECT) {
		rates = new_Parameters(rates_node->child_count);
		for(int i = 0; i < rates_node->child_count; i++){
			json_node* p_node = rates_node->children[i];
			Parameter* p = new_Parameter_from_json(p_node, hash);
			Parameters_move(rates, p);
			if(strcasecmp(p_node->key, "invariant") == 0){
				invariant = true;
			}
		}
		check_constraints(rates, 0, INFINITY, 0.001, 100);
	}
	
	for (int i = 0; i < Parameters_count(rates); i++) {
		Hashtable_add(hash, Parameters_name(rates, i), Parameters_at(rates, i));
	}
	
	SiteModel* sm = NULL;
	// Distribution
	if ( cat > 1 ) {
		const char* average_str = get_json_node_value_string(node, "average");
		if (average_str != NULL) {
			char* average_allowed[] = {"mean", "median", "laguerre"};
			if(!array_of_string_contains(average_str, average_allowed, 3, false)){
				fprintf(stderr, "In %s of type %s\n", get_json_node_value_string(node, "id"), get_json_node_value_string(node, "type"));
				fprintf(stderr, "%s is not a valid value for key %s\n", average_str, "average");
				exit(13);
			}
			
		}
		// Just laguerre discretization
		if(discretization_node != NULL && strcasecmp("laguerre", (char*)discretization_node->value) == 0){
			sm = new_GammaLaguerreSiteModel_with_parameters(m, rates, cat);
		}
		// And invariant
		else if (invariant) {
			//swapsy
			if (strcasecmp(Parameters_name(rates, 0), "invariant") == 0) {
				Parameters_swap_index(rates, 0, 1);
			}
			sm = new_GammaPinvSiteModel_with_parameters(m, rates, cat);
			if (average_str != NULL && strcasecmp(average_str, "mean") == 0) {
				sm->type = SITEMODEL_GAMMAINV_MEAN;
			}
		}
		// Just distribution
		else{
			sm = new_GammaSiteModel_with_parameters(m, rates, cat);
			if (average_str != NULL && strcasecmp(average_str, "mean") == 0) {
				sm->type = SITEMODEL_GAMMA_MEAN;
			}
		}

	}
	else if(invariant){
		sm = new_PinvSiteModel_with_parameters(m, rates);
	}
	else{
		sm = new_SiteModel(m);
	}

	char* id_string = get_json_node_value_string(node, "id");

	Model* msm = new_SiteModel2(id_string, sm, mm);
	
	if (mu_node != NULL) {
		sm->mu = new_Parameter_from_json(mu_node, hash);
		check_constraint(sm->mu, 0, INFINITY, 0.001, 100);
		Hashtable_add(hash, Parameter_name(sm->mu), sm->mu);
		sm->mu->listeners->add( sm->mu->listeners, msm );
	}
	
	mm->free(mm);
	free_Parameters(rates);
	return msm;
}

#pragma mark -
// MARK: Gamma SiteModel

SiteModel * new_SiteModel( SubstitutionModel *m ){
	return _new_SiteModel_with_parameters( m, NULL, 0, SITEMODEL_NONE );
}

SiteModel * new_PinvSiteModel( SubstitutionModel *m, const double pinv ){
	assert( pinv>0 && pinv <1 );
	Parameters *params = new_Parameters(1);
	Parameters_move(params, new_Parameter("sitemodel.pinv", pinv, new_Constraint(0+TINY, 1-TINY) ));
	SiteModel* sm = _new_SiteModel_with_parameters( m, params, 0, SITEMODEL_INV );
	free_Parameters(params);
	return sm;
}

SiteModel * new_PinvSiteModel_with_parameters( SubstitutionModel *m, const Parameters* pinv ){
	assert( Parameters_value(pinv, 0)>0 && Parameters_value(pinv, 0) <1 );
	return _new_SiteModel_with_parameters( m, pinv, 0, SITEMODEL_INV );
}

SiteModel * new_GammaSiteModel( SubstitutionModel *m, const double shape, const unsigned int cat_count ){
	assert( shape >= SITEMODEL_ALPHA_MIN && shape <= SITEMODEL_ALPHA_MAX );
	Parameters *params = new_Parameters(1);
	Parameters_move(params, new_Parameter("sitemodel.alpha", shape, new_Constraint(SITEMODEL_ALPHA_MIN, SITEMODEL_ALPHA_MAX) ));
	SiteModel* sm = _new_SiteModel_with_parameters( m, params, cat_count, SITEMODEL_GAMMA );
	free_Parameters(params);
	return sm;
}

SiteModel * new_GammaSiteModel_with_parameters( SubstitutionModel *m, const Parameters* shape, const size_t cat_count ){
	assert( Parameters_value(shape, 0) >= SITEMODEL_ALPHA_MIN && Parameters_value(shape, 0) <= SITEMODEL_ALPHA_MAX );
	return _new_SiteModel_with_parameters( m, shape, cat_count, SITEMODEL_GAMMA );
}

SiteModel * new_GammaPinvSiteModel( SubstitutionModel *m, const double pinv, const double shape, const unsigned int cat_count ){

    if( pinv == 0 && cat_count == 1 ){
        return _new_SiteModel_with_parameters( m, NULL, 0, SITEMODEL_NONE );
    }
    
    Parameters *params = new_Parameters(2);
    
    if( cat_count > 1 ){
        Parameter *alpha = new_Parameter("sitemodel.alpha", shape, new_Constraint(SITEMODEL_ALPHA_MIN, SITEMODEL_ALPHA_MAX) );
        Parameters_move(params, alpha);
    }
    if( pinv > 0 ){
        Parameter *inv = new_Parameter("sitemodel.pinv", pinv, new_Constraint(0+TINY, 1-TINY) );
        Parameters_move(params, inv);
    }
    
    sitemodel_heterogeneity type = SITEMODEL_GAMMAINV;
    if( pinv == 0 && cat_count > 1 ){
        type = SITEMODEL_GAMMA;
    }
    else if ( pinv > 0.0 && cat_count == 1 ){
        type = SITEMODEL_INV;
    }
	SiteModel* sm = _new_SiteModel_with_parameters( m, params, cat_count, type );
	free_Parameters(params);
	return sm;
}

SiteModel * new_GammaPinvSiteModel_with_parameters( SubstitutionModel *m, const Parameters* params, const size_t cat_count ){
	sitemodel_heterogeneity type = SITEMODEL_GAMMAINV;
	if( Parameters_count(params) == 1 && cat_count > 1 ){
		type = SITEMODEL_GAMMA;
	}
	else if ( Parameters_count(params) == 1 && cat_count == 1 ){
		type = SITEMODEL_INV;
	}
	return _new_SiteModel_with_parameters( m, params, cat_count, type );
}

SiteModel * new_GammaLaguerreSiteModel( SubstitutionModel *m, const double shape, const unsigned int cat_count ){
	assert( shape >= SITEMODEL_ALPHA_MIN && shape <= SITEMODEL_ALPHA_MAX );
	Parameters *params = new_Parameters(1);
	Parameter *alpha = new_Parameter("sitemodel.alpha", shape, new_Constraint(SITEMODEL_ALPHA_MIN, SITEMODEL_ALPHA_MAX) );
	Parameters_move(params, alpha);
	SiteModel* sm = _new_SiteModel_with_parameters( m, params, cat_count, SITEMODEL_GAMMA_LAGUERRE );
	free_Parameters(params);
	return sm;
}

SiteModel * new_GammaLaguerreSiteModel_with_parameters( SubstitutionModel *m, const Parameters* shape, const size_t cat_count ){
	assert( Parameters_value(shape, 0) >= SITEMODEL_ALPHA_MIN && Parameters_value(shape, 0) <= SITEMODEL_ALPHA_MAX );
	return _new_SiteModel_with_parameters( m, shape, cat_count, SITEMODEL_GAMMA_LAGUERRE );
}

void set_rate(SiteModel* sm, const int index, const double value){
    Parameters_set_value(sm->rates, index, value);
    sm->need_update = true;
}

// (Gamma or/and Invariant) or one rate
// should not be used directly
SiteModel * _new_SiteModel_with_parameters( SubstitutionModel *m, const Parameters *params, const size_t cat_count, sitemodel_heterogeneity type ){
	SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
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
    
    sm->set_rate = set_rate;
    
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
            Parameters_add(sm->rates, Parameters_at(params, 0));
			sm->cat_count       = cat_count;
			sm->get_rate        = _get_rate_gamma;
			sm->get_proportion  = _get_proportion_gamma;
			sm->get_proportions = _get_proportions_gamma;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
			break;
		case SITEMODEL_GAMMAINV:
            sm->rates = new_Parameters(2);
            Parameters_add(sm->rates, Parameters_at(params, 0));
            Parameters_add(sm->rates, Parameters_at(params, 1));
            sm->cat_count       = cat_count+1;
			sm->get_rate        = _get_rate_gamma;
			sm->get_proportion  = _get_proportion_gamma;
			sm->get_proportions = _get_proportions_gamma;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
			break;
        case SITEMODEL_INV:
            sm->rates = new_Parameters(1);
            Parameters_add(sm->rates, Parameters_at(params, 0));
			sm->cat_count       = 2;
			sm->get_rate        = _get_rate_invariant;
			sm->get_proportion  = _get_proportion_invariant;
			sm->get_proportions = _get_proportions_invariant;
            sm->cat_rates       = dvector(sm->cat_count);
            sm->cat_proportions = dvector(sm->cat_count);
			break;
        case SITEMODEL_GAMMA_LAGUERRE:
            sm->rates = new_Parameters(1);
            Parameters_add(sm->rates, Parameters_at(params, 0));
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
    
	sm->integrate   = true;
	sm->need_update = true;
	
	return sm;
}



double _get_rate( SiteModel *sm, const int index ){
	return (sm->mu == NULL ? 1.0 : Parameter_value(sm->mu));
}

double _get_proportion( SiteModel *sm, const int index ){
	return 1.0;
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
	
	// median
	if(sm->type == SITEMODEL_GAMMA || sm->type == SITEMODEL_GAMMAINV){
		for (int i = 0; i < nCat; i++) {
			//sm->cat_rates[i + cat] = qgamma( (2.0 * i + 1.0) / (2.0 * nCat), alpha, alpha );
			sm->cat_rates[i + cat] = gsl_cdf_gamma_Qinv( (2.0 * i + 1.0) / (2.0 * nCat), alpha, 1.0/alpha );
			mean += sm->cat_rates[i + cat];
			
			sm->cat_proportions[i + cat] = propVariable / nCat;
		}
		
		mean = (propVariable * mean) / nCat;
		
		for (int i = 0; i < nCat; i++) {
			sm->cat_rates[i + cat] /= mean;
		}
	}
	// mean
	else{
		for (int i = 0; i < nCat - 1; i++){
			sm->cat_proportions[i+cat] = qgamma((i + 1.0) / nCat, alpha, alpha);
			sm->cat_proportions[i+cat] = gammp(alpha + 1, sm->cat_proportions[i+cat] * alpha);
		}
		
		sm->cat_rates[cat] = sm->cat_proportions[cat]*nCat;
		for (int i = cat+1; i < nCat - 1; i++){
			sm->cat_rates[i] = (sm->cat_proportions[i] - sm->cat_proportions[i - 1])*nCat;
		}
		sm->cat_rates[nCat - 1] = (1 - sm->cat_proportions[nCat - 2])*nCat;
		
		for (int i = 0; i < nCat; i++){
			sm->cat_proportions[i + cat] = propVariable/nCat;
		}
	}
    
    sm->need_update = false;
}

void _invariant( SiteModel *sm ) {
    
    double propVariable = 1.0;
    int cat = 0;
    
    // there is also invariant sites
    if ( Parameters_count(sm->rates) > 0 ){
        sm->cat_rates[0] = 0.0;
        sm->cat_proportions[0] = Parameters_value(sm->rates, 1);
        
        propVariable = 1.0 - sm->cat_proportions[0];
        cat = 1;
    }
    sm->cat_rates[cat] = 1.0 / propVariable;
    sm->cat_proportions[cat] = propVariable;
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
	sm->m = m;
	sm->nstate = m->nstate;
	
    sm->set_rate        = set_rate;
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
	Parameters_move(sm->rates, new_Parameter(buffer->c, prop, new_Constraint(0.00001, 0.999)) );
    int i = 1;
	for ( ; i < cat_count-1; i++ ) {
		StringBuffer_set_string(buffer, "sitemodel.p");
		StringBuffer_append_format(buffer, "%d", i);
		double denom = 1;
		for ( int j = 0; j < i; j++ ) {
			denom *= 1 - Parameters_value(sm->rates, j);
		}
		double r = sm->cat_proportions[i] / denom;
		Parameters_move(sm->rates, new_Parameter(buffer->c, r, new_Constraint(0.00001, 0.999)) );
	}
	
    for (int j = 0 ; j < cat_count-1; j++, i++ ) {
		StringBuffer_set_string(buffer, "sitemodel.r");
		StringBuffer_append_format(buffer, "%d", j);
		Parameters_move(sm->rates, new_Parameter(buffer->c, 1.0, new_Constraint(0.00001, 10)) );
	}
    
	//Parameters_print(sm->rates);
    
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
		Parameters_move(sm->rates, new_Parameter(buffer->c, 1., new_Constraint(0.00001, 100)) );
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
	return clone_SiteModel_with(sm, sm->m->clone(sm->m));
}

SiteModel * clone_SiteModel_with( const SiteModel *sm, SubstitutionModel* m ){
	SiteModel *newsm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(newsm);
	newsm->m = m;
	newsm->nstate = sm->nstate;
	
	newsm->rates = NULL;
	
	newsm->cat_count = sm->cat_count;
	
	if ( sm->rates != NULL ){
		newsm->rates = clone_Parameters(sm->rates);
	}
	
	newsm->mu = NULL;
	if ( sm->mu != NULL ){
		newsm->mu = clone_Parameter(sm->mu);
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

SiteModel * clone_SiteModel_with_parameters( const SiteModel *sm, SubstitutionModel* m, const Parameters* params ){
	SiteModel *newsm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(newsm);
	newsm->m = m;
	newsm->type = sm->type;
	newsm->nstate = sm->nstate;
	
	newsm->rates = NULL;
	
	newsm->cat_count = sm->cat_count;
	
	if ( sm->rates != NULL ){
		newsm->rates = new_Parameters(Parameters_count(sm->rates));
		for (int i = 0; i < Parameters_count(sm->rates); i++) {
			Parameters_add(newsm->rates, Parameters_at(params, i));
		}
	}
	
	newsm->mu = NULL;
	if ( sm->mu != NULL ){
		newsm->mu = clone_Parameter(sm->mu);
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


void free_SiteModel( SiteModel *sm ){
	if ( sm->rates != NULL ) free_Parameters(sm->rates);
	if ( sm->mu != NULL ) free_Parameter(sm->mu);
	if ( sm->cat_proportions != NULL ) free(sm->cat_proportions);
    if ( sm->cat_assignment != NULL ) free(sm->cat_assignment);
	free(sm->cat_rates);
	if(sm->m!=NULL)sm->m->free(sm->m);
	free(sm);
}
