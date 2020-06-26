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
#include <math.h>

#include "substmodel.h"
#include "parameters.h"
#include "gamma.h"
#include "gausslaguerre.h"
#include "mstring.h"
#include "matrix.h"
#include "mathconstant.h"
#include "gaussian.h"
#include "distkumaraswamy.h"

#include "model.h"

#include <gsl/gsl_cdf.h>

static void _gamma_approx_quantile( SiteModel *sm );
static void _calculate_rates_discrete( SiteModel *sm );

static double _get_rate( SiteModel *sm, const int index );
static double _get_proportion( SiteModel *sm, const int index );
static double *_get_proportions( SiteModel *sm );

static double _get_rate_gamma( SiteModel *sm, const int index );
static double _get_proportion_gamma( SiteModel *sm, const int index );
static double *_get_proportions_gamma( SiteModel *sm );

static double _get_rate_laguerre( SiteModel *sm, const int index );
static double _get_proportion_laguerre( SiteModel *sm, const int index );
static double * _get_proportions_laguerre( SiteModel *sm );


static double _get_rate_discrete( SiteModel *sm, const int index );
static double _get_proportion_discrete( SiteModel *sm, const int index );
static double *_get_proportions_discrete( SiteModel *sm );


static void _site_model_handle_change( Model *self, Model *model, int index ){
	SiteModel *sm = (SiteModel*)self->obj;
	sm->need_update = true;// one of the sitemodel parameters
	self->listeners->fire( self->listeners, self, index );
}

static void _site_model_store(Model* self){
	Model** models = (Model**)self->data;
	Model *msubst = models[0];
	Model *mprop = models[1];
	msubst->store(msubst); // substitution model
	if(mprop != NULL) mprop->store(mprop); // simplex proportion model
	SiteModel* sm = self->obj;
	if (Parameters_count(sm->rates) > 0) {
		Parameters_store(sm->rates);
	}
	if (sm->mu != NULL) {
		Parameter_store(sm->mu);
	}
}

static void _site_model_restore(Model* self){
	Model** models = (Model**)self->data;
	Model *msubst = models[0];
	Model *mprop = models[1];
	msubst->restore(msubst); // substitution model
	if(mprop != NULL) mprop->restore(mprop); // simplex proportion model
	SiteModel* sm = self->obj;
	if (Parameters_count(sm->rates) > 0) {
		bool changed = false;
		Parameter*p = NULL;
		for (int i = 0; i < Parameters_count(sm->rates); i++) {
			p = Parameters_at(sm->rates, i);
			if (Parameter_changed(p)) {
				changed = true;
				Parameter_restore_quietly(p);
			}
		}
		if (changed) {
			p->listeners->fire_restore(p->listeners, NULL, p->id);
		}
	}
	if (sm->mu != NULL) {
		Parameter_restore(sm->mu);
	}
}

static void _site_model_handle_restore( Model *self, Model *model, int index ){
	SiteModel *sm = (SiteModel*)self->obj;
	sm->need_update = true;// one of the sitemodel parameters
	self->listeners->fire_restore( self->listeners, self, index );
}

static void _site_model_free( Model *self ){
	if(self->ref_count == 1){
		//printf("Free site model %s\n", self->name);
		SiteModel *sm = (SiteModel*)self->obj;
		Model** models = (Model**)self->data;
		models[0]->free(models[0]); // substitution model
		if(models[1] != NULL) models[1]->free(models[1]);
		free(self->data);
		
		if(sm->rates != NULL) free_Parameters(sm->rates);
		//TODO: deal with this
		if(sm->mu!=NULL)free_Parameter(sm->mu);
		if ( sm->cat_proportions != NULL ) free(sm->cat_proportions);
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
	Model** models = (Model**)self->data;
	Model *mm = models[0];
	Model *mprop = models[1];
	Model* mmclone = NULL;
	Model* mpropclone = NULL;
    Simplex* propclone = NULL;
	// Susbtitution model may have been parsed already
	if (Hashtable_exists(hash, mm->name)) {
		mmclone = Hashtable_get(hash, mm->name);
		mmclone->ref_count++; // it is decremented at the end using free
	}
	else{
		mmclone = mm->clone(mm, hash);
		Hashtable_add(hash, mmclone->name, mmclone);
	}
	
    if(mprop != NULL){
        if (Hashtable_exists(hash, mprop->name)) {
            mpropclone = Hashtable_get(hash, mprop->name);
            mpropclone->ref_count++; // it is decremented at the end using free
        }
        else{
            mpropclone = mprop->clone(mprop, hash);
            Hashtable_add(hash, mpropclone->name, mpropclone);
        }
        propclone = mpropclone->obj;
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
	Parameter* mu = NULL;
	if (sm->mu != NULL) {
		char* name = Parameter_name(sm->mu);
		if (Hashtable_exists(hash, name)) {
			mu = Hashtable_get(hash, name);
			mu->refCount++;
		}
		else{
			mu = clone_Parameter(sm->mu);
			Hashtable_add(hash, name, mu);
		}
	}
	SiteModel* smclone = clone_SiteModel_with_parameters(sm, (SubstitutionModel*)mmclone->obj, propclone, ps, mu);
	free_Parameters(ps);
	free_Parameter(mu);
	Model* clone = new_SiteModel2(self->name, smclone, mmclone, mpropclone);
	Hashtable_add(hash, clone->name, clone);
	mmclone->free(mmclone);
	if(mpropclone != NULL)mpropclone->free(mpropclone);
	clone->print = self->print;
	return clone;
}

// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SiteModel2( const char* name, SiteModel *sm, Model *substmodel, Model* proportions ){
	Model *model = new_Model(MODEL_SITEMODEL, name, sm);
	Model** data = (Model**)calloc(2, sizeof(Model*));
	data[0] = substmodel;
	data[1] = NULL;
	if ( sm->rates != NULL ) {
		for ( int i = 0; i < Parameters_count(sm->rates); i++ ) {
			Parameters_at(sm->rates, i)->listeners->add( Parameters_at(sm->rates, i)->listeners, model );
		}
	}
	if ( sm->mu != NULL ) {
		sm->mu->listeners->add( sm->mu->listeners, model );
	}
	if (proportions != NULL) {
		data[1] = proportions;
		proportions->ref_count++;
		proportions->listeners->add( proportions->listeners, model );
	}
	
	// Listen to substitution model
	substmodel->listeners->add( substmodel->listeners, model );
	
	model->update = _site_model_handle_change;
	model->handle_restore = _site_model_handle_restore;
	model->store = _site_model_store;
	model->restore = _site_model_restore;
	model->free = _site_model_free;
	model->clone = _site_model_clone;
	model->data = data;
	substmodel->ref_count++;
	return model;
}

static void _SiteModel_print(Model* model, FILE* out){
	SiteModel* sm = model->obj;
	sm->get_proportions(sm); // make sure it is updated
	if (Parameters_count(sm->rates)) {
		for (int i = 0; i < Parameters_count(sm->rates); i++) {
			fprintf(out, "%s %f\n", Parameters_name(sm->rates, i),Parameters_value(sm->rates, i));
		}
		double sum = 0;
		fprintf(out, "proportion rate\n");
		for (int i = 0; i < sm->cat_count; i++){
			fprintf(out, "%f %f\n",sm->cat_proportions[i], sm->cat_rates[i]);
			sum += sm->cat_proportions[i]*sm->cat_rates[i];
		}
		fprintf(out, "sum %f\n", sum);
	}
	else if(sm->proportions != NULL){
		fprintf(out, "proportion rate %d\n", sm->proportions->K);
		for (int i = 0; i < sm->proportions->K; i++){
			fprintf(out, "%f %f\n",sm->cat_proportions[i], sm->cat_rates[i]);
		}
	}
	else if(sm->cat_count > 1){
		fprintf(out, "proportion rate %d\n", sm->cat_count);
		for (int i = 0; i < sm->cat_count; i++){
			fprintf(out, "%f %f\n",sm->cat_proportions[i], sm->cat_rates[i]);
		}
	}
}

void set_rate(SiteModel* sm, const int index, const double value){
    Parameters_set_value(sm->rates, index, value);
    sm->need_update = true;
}

int _get_site_category(SiteModel* sm, const int pattern){
	return 0;
}

void _no_gradient( SiteModel *sm, const double* ingrad, double* grad ){

}

void _weibull_gradient( SiteModel *sm, const double* ingrad, double* grad ){
	double derivsum = 0;
	double sum = 0;
	double cat = sm->cat_count;
	double* proportions = sm->get_proportions(sm);
	double shape = Parameters_value(sm->rates, 0);

	for (int i = 0; i < cat; i++) {
		double prob = (2.0*i + 1.0)/(2.0*cat);
		derivsum += -pow(-log(1.0 - prob), 1.0/shape) * log(-log(1.0 - prob))/shape/shape*proportions[i];
		sum += pow(-log(1.0 - prob), 1.0/shape)*proportions[i];
	}

	double shape_gradient = 0;
	for (int i = 0; i < cat; i++) {
		double prob = (2.0*i + 1.0)/(2.0*cat);
		double deriv = -pow(-log(1.0 - prob), 1.0/shape) * log(-log(1.0 - prob))/shape/shape;
		double deriv_rate = deriv/sum - (pow(-log(1.0 - prob), 1.0/shape)*derivsum)/sum/sum;
		shape_gradient += ingrad[i] * deriv_rate * proportions[i];
	}
	grad[0] = shape_gradient;
}

// (Gamma or/and Invariant) or one rate
// should not be used directly
SiteModel * new_SiteModel_with_parameters( SubstitutionModel *m, const Parameters *params, Simplex* proportions, const size_t cat_count, distribution_t distribution, bool invariant, quadrature_t quad){
	SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	sm->site_category = NULL;
	sm->sp = NULL;
	sm->get_site_category = _get_site_category;
	sm->m = m;
	sm->nstate = m->nstate;
	
	sm->distribution = distribution;
	sm->invariant = invariant;
	sm->quadrature = quad;
	
	sm->cat_count = cat_count;
	if (invariant) sm->cat_count++;
	
	sm->cat_rates       = NULL;
	sm->cat_proportions = NULL;
	sm->proportions = proportions;
    
	sm->rates = NULL;
	
	if(Parameters_count(params) > 0){
		sm->rates = new_Parameters(Parameters_count(params));
		Parameters_add_parameters(sm->rates, params);
		if (Parameters_name2(params) != NULL) {
			Parameters_set_name2(sm->rates, Parameters_name2(params));
		}
	}
	sm->mu    = NULL;
    
    sm->set_rate = set_rate;
	
	sm->cat_rates       = dvector(sm->cat_count);
	sm->cat_proportions = dvector(sm->cat_count);
	sm->gradient = _no_gradient;
	
	if (distribution == DISTRIBUTION_UNIFORM) {
		sm->get_rate        = _get_rate;
		sm->get_proportion  = _get_proportion;
		sm->get_proportions = _get_proportions;
		sm->cat_rates[0] = sm->cat_proportions[0] = 1;
	}
	else if (distribution == DISTRIBUTION_DISCRETE && Parameters_count(sm->rates) != 0) {
		sm->get_rate        = _get_rate_discrete;
		sm->get_proportion  = _get_proportion_discrete;
		sm->get_proportions = _get_proportions_discrete;
	}
	else if (quad == QUADRATURE_GAUSS_LAGUERRE){
		sm->get_rate        = _get_rate_laguerre;
		sm->get_proportion  = _get_proportion_laguerre;
		sm->get_proportions = _get_proportions_laguerre;
	}
	else{
		sm->get_rate        = _get_rate_gamma;
		sm->get_proportion  = _get_proportion_gamma;
		sm->get_proportions = _get_proportions_gamma;
	}
    
	if (distribution == DISTRIBUTION_WEIBULL) {
		sm->gradient = _weibull_gradient;
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

double icdf_weibull(double p, double lambda, double k){
	return lambda*pow(-log(1.0 - p), 1.0/k);
}

double icdf_weibull_1(double p, double k){
	return pow(-log(1.0 - p), 1.0/k);
}

void _gamma_approx_quantile( SiteModel *sm ) {
	
	double propVariable = 1.0;
	int cat = (sm->invariant ? 1 : 0);
	const int nCat = sm->cat_count - cat;
	double* quantiles = dvector(sm->cat_count); // cat_count includes invariant sites if specified

	// proportions are estimated from beta distribution
	if (sm->quadrature == QUADRATURE_BETA || sm->quadrature == QUADRATURE_KUMARASWAMY){
		double shape_alpha = Parameters_value(sm->rates, Parameters_count(sm->rates)-2);
		double shape_beta = Parameters_value(sm->rates, Parameters_count(sm->rates)-1);
		
		// proportion of invariant is estimated from beta
		// cat can be equal to 0 or 1
		if(sm->proportions == NULL){
			if(sm->quadrature == QUADRATURE_BETA){
				for (int i = 0; i < sm->cat_count; i++) {
					quantiles[i] = gsl_cdf_beta_Pinv((double)i/sm->cat_count, shape_alpha, shape_beta);
				}
			}
			else{
				for (int i = 0; i < sm->cat_count; i++) {
					quantiles[i] = DistributionModel_kumaraswamy_inverse_CDF((double)i/sm->cat_count, shape_alpha, shape_beta);
				}
			}
			
			for (int i = 0; i < sm->cat_count-1; i++) {
				// calculate cat proportions
				sm->cat_proportions[i] = quantiles[i + 1] - quantiles[i];
				// calculate quantiles for gamma distribution
				quantiles[i] = quantiles[i] + (quantiles[i+1] - quantiles[i])/2.0;
			}
			sm->cat_proportions[sm->cat_count-1] = 1.0 - quantiles[sm->cat_count-1];
			quantiles[sm->cat_count-1] = quantiles[sm->cat_count-1] + (1.0 - quantiles[sm->cat_count-1])/2.0;
		}
		// proportion of invariant come from a simplex
		// cat should be equal to 1
		else{
			for (int i = 0; i < sm->cat_count-cat; i++) {
				quantiles[i+cat] = gsl_cdf_beta_Pinv((double)i/(sm->cat_count-cat), shape_alpha, shape_beta);
			}
			
			propVariable = sm->proportions->get_value(sm->proportions, 1);
			sm->cat_proportions[0] = sm->proportions->get_value(sm->proportions, 0);
			for (int i = 0; i < sm->cat_count-1-cat; i++) {
				// calculate cat proportions
				sm->cat_proportions[i+cat] = (quantiles[i + 1 + cat] - quantiles[i+cat])*propVariable;
				// calculate quantiles for gamma distribution
				quantiles[i+cat] = quantiles[i+cat] + (quantiles[i+1+cat] - quantiles[i+cat])/2.0;
			}
			sm->cat_proportions[sm->cat_count-1] = (1.0 - quantiles[sm->cat_count-1])*propVariable;
			quantiles[sm->cat_count-1] = quantiles[sm->cat_count-1] + (1.0 - quantiles[sm->cat_count-1])/2.0;
		}
	}
	// That's +G+D or +G+I+D
	else if(sm->quadrature == QUADRATURE_DISCRETE){
		const double* values = sm->proportions->get_values(sm->proportions);
		double sum = 0;
		for (int i = 0; i < sm->cat_count-cat; i++) {
			quantiles[i+cat] = sum + values[i+cat]/2.0; // pick midpoint
			sum += values[i];
		}
		memcpy(sm->cat_proportions, values, sizeof(double)*sm->proportions->K);
	}
	// if the dimension of the simplex is then only pinv is estimated
	// and the remaining categories split 1.0-pinv equally
	// That's the traditional +G+I or +I
	else if (sm->proportions != NULL && sm->proportions->K == 2) {
		const double* values = sm->proportions->get_values(sm->proportions);
		// +I
		if(sm->distribution == DISTRIBUTION_DISCRETE){
			memcpy(sm->cat_proportions, values, sizeof(double)*2);
			sm->cat_rates[0] = 0;
			sm->cat_rates[1] = 1.0 / sm->cat_proportions[1];
			sm->need_update = false;
			return;
		}
		// +G+I
		else{
			cat = 1;
			sm->cat_proportions[0] = values[0];
//			printf("%f %f\n", values[0], values[1]);
			propVariable = values[1];
			int gammaCat = sm->cat_count - 1;
			for (int i = 0; i < gammaCat; i++) {
				quantiles[i+cat] = (2.0 * i + 1.0) / (2.0 * gammaCat);
				sm->cat_proportions[i+cat] = propVariable/gammaCat;
			}
//			exit(2);
		}
	}
	// distribution without invariant and proportions not estimated
	// cat_proportions corresponds to quantiles not proportions
	// That's the traditional +G
	else if(sm->proportions == NULL){
		for (int i = 0; i < nCat; i++) {
			quantiles[i] = (2.0 * i + 1.0) / (2.0 * nCat);
			sm->cat_proportions[i] = 1.0/nCat;
		}
	}
	
	double mean = 0.0;
	
	const double alpha = Parameters_value(sm->rates, 0);
	
	// median
	if(sm->quadrature == QUADRATURE_QUANTILE_MEDIAN || sm->quadrature == QUADRATURE_DISCRETE ||
	   sm->quadrature == QUADRATURE_BETA ||	sm->quadrature == QUADRATURE_KUMARASWAMY){
		sm->cat_rates[0] = 0;
		if(sm->distribution == DISTRIBUTION_GAMMA){
			for (int i = 0; i < sm->cat_count - cat; i++) {
				sm->cat_rates[i + cat] = gsl_cdf_gamma_Qinv( quantiles[i+cat], alpha, 1.0/alpha );
			}
		}
		else if(sm->distribution == DISTRIBUTION_WEIBULL){
			for (int i = 0; i < sm->cat_count - cat; i++) {
				// Unit mean Weibull
				// sm->cat_rates[i + cat] = gsl_cdf_weibull_Qinv( sm->cat_proportions[i + cat], alpha, 1.0/exp(gammln(1.0 + 1.0/alpha)) );
				// Fix lambda:=1
				sm->cat_rates[i + cat] = icdf_weibull_1( quantiles[i+cat], alpha);
			}
		}
		else if(sm->distribution == DISTRIBUTION_LOGNORMAL){
			for (int i = 0; i < sm->cat_count - cat; i++) {
				sm->cat_rates[i + cat] = gsl_cdf_lognormal_Qinv( quantiles[i+cat], -alpha*alpha/2, alpha );
			}
		}
		else if(sm->distribution == DISTRIBUTION_BETA){
			for (int i = 0; i < sm->cat_count - cat; i++) {
				sm->cat_rates[i + cat] = gsl_cdf_beta_Qinv( quantiles[i+cat], alpha, Parameters_value(sm->rates, 1) );
			}
		}
		
		if (sm->quadrature == QUADRATURE_BETA || sm->quadrature == QUADRATURE_KUMARASWAMY) {
			for (int i = 0; i < sm->cat_count - cat; i++) {
				mean += sm->cat_rates[i + cat]*sm->cat_proportions[i + cat];
			}
		}
		// QUADRATURE_DISCRETE
		else if ( (sm->proportions != NULL && sm->proportions->K == sm->cat_count)){
			for (int i = 0; i < nCat; i++) {
				mean += sm->cat_rates[i + cat]*sm->cat_proportions[i + cat];
			}
		}
		// +G or +G+I
		else if(sm->proportions == NULL || (sm->proportions != NULL && sm->proportions->K != sm->cat_count)){
			for (int i = 0; i < sm->cat_count - cat; i++) {
				mean += sm->cat_rates[i + cat];
			}
			mean = (propVariable * mean) / nCat;
		}
		
		for (int i = 0; i < sm->cat_count - cat; i++) {
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
	free(quantiles);
	sm->need_update = false;
}

// Gamma distribution approximated using Laguerre quadrature
void _gamma_approx_laguerre( SiteModel *sm ){
    const double alpha = Parameters_value(sm->rates, 0);
    
	// calculate using alpha -1
	gaulag(sm->cat_rates, sm->cat_proportions, sm->cat_count, alpha-1);
	
	double gamalpha = gamm(alpha);
	for ( int i = 0; i < sm->cat_count; i++ ) {
		sm->cat_rates[i] /= alpha;
		sm->cat_proportions[i] /= gamalpha;
	}
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

int _get_site_category_CAT(SiteModel* sm, const int pattern){
	return sm->site_category[pattern];
}

#pragma mark -
// MARK: Discrete SiteModel

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
	int j = 0;
	int cat = (sm->invariant ? 1 : 0);
	int cat_count = sm->proportions->K;
	memset(sm->cat_rates, 0, sizeof(double)*cat_count);
	double* cat_proportions = sm->proportions->get_values(sm->proportions);
	memcpy(sm->cat_proportions, cat_proportions, sizeof(double)*cat_count);
	
	if(Parameters_count(sm->rates) == 0){
		sm->cat_rates[1] = 1.0/sm->cat_proportions[1];
		sm->need_update = false;
		return;
	}
	
	sm->cat_rates[cat] = Parameters_value(sm->rates, j++);
	sm->cat_rates[cat+1] = 1;
	double sum = sm->cat_rates[cat]*sm->cat_proportions[cat] + sm->cat_proportions[cat+1];
	if(sm->invariant) sum += sm->cat_proportions[0];
	double prod = 1;
	for (int i = cat+2; i < cat_count; i++, j++ ) {
		prod *= Parameters_value(sm->rates, j);
		sm->cat_rates[i] = prod;
		sum += prod*sm->cat_proportions[i];
	}
	
	for (int i = 0; i < cat_count; i++ ) {
		sm->cat_rates[i] /= sum;
	}
	sm->need_update = false;
}

#pragma  mark -
#pragma mark CAT

void _cat_update(SiteModel* sm){
	double avg = 0;
	for (int i = 0; i < sm->sp->count; i++ ) {
		avg += Parameters_value(sm->rates, sm->get_site_category(sm, i)) * sm->sp->weights[i];
	}
	avg /= sm->sp->nsites;
	for (int i = 0; i < sm->cat_count; i++ ) {
		sm->cat_rates[i] = Parameters_value(sm->rates, i)/avg;
	}
	sm->need_update = false;
}

double _get_rate_cat( SiteModel *sm, const int index ){
	if ( sm->need_update ) {
		_cat_update(sm);
	}
	return sm->cat_rates[index] * (sm->mu == NULL ? 1.0 : Parameter_value(sm->mu));
}

SiteModel * new_CATSiteModel_with_parameters( SubstitutionModel *m, const Parameters *params,  const size_t cat_count, SitePattern* sp){
	SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	sm->sp = sp;
	sm->site_category = ivector(sp->count);
	sm->get_site_category = _get_site_category_CAT;
	sm->m = m;
	sm->nstate = m->nstate;
	
	sm->distribution = -1;
	sm->invariant = false;
	sm->quadrature = -1;
	
	sm->cat_count = cat_count;
	
	sm->cat_rates = dvector(sm->cat_count);
	sm->cat_proportions = NULL;
	sm->proportions = NULL;
	
	sm->rates = NULL;
	
	if(Parameters_count(params) > 0){
		sm->rates = new_Parameters(Parameters_count(params));
		Parameters_add_parameters(sm->rates, params);
		if (Parameters_name2(params) != NULL) {
			Parameters_set_name2(sm->rates, Parameters_name2(params));
		}
	}
	sm->mu    = NULL;
	
	sm->set_rate = set_rate;
	
	sm->get_rate        = _get_rate_cat;
	sm->get_proportion  = _get_proportion;
	sm->get_proportions = _get_proportions;
	
	sm->integrate   = false;
	sm->need_update = true;
	
	return sm;
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
    newsm->get_site_category = sm->get_site_category;
	
	newsm->distribution = sm->distribution;
	newsm->invariant = sm->invariant;
	newsm->quadrature = sm->quadrature;
	
	return newsm;
}

SiteModel * clone_SiteModel_with_parameters( const SiteModel *sm, SubstitutionModel* m, Simplex* props, const Parameters* params, Parameter* mu ){
	SiteModel *newsm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(newsm);
	newsm->m = m;
	newsm->proportions = props;
	newsm->distribution = sm->distribution;
	newsm->invariant = sm->invariant;
	newsm->quadrature = sm->quadrature;
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
	if ( mu != NULL ){
		newsm->mu = mu;
		mu->refCount++;
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
    newsm->get_site_category = sm->get_site_category;
	
	return newsm;
}


void free_SiteModel( SiteModel *sm ){
	if ( sm->rates != NULL ) free_Parameters(sm->rates);
	if ( sm->mu != NULL ) free_Parameter(sm->mu);
	if ( sm->cat_proportions != NULL ) free(sm->cat_proportions);
	free(sm->cat_rates);
	if(sm->m!=NULL)sm->m->free(sm->m);
	free(sm);
}

Model* new_SiteModel_from_json(json_node*node, Hashtable*hash){
	char* allowed[] = {
		"distribution",
		"invariant",
		"model",
		"mu",
		"rates",
		"sitepattern",
		"substitutionmodel"
	};
	//json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	json_node* distribution_node = get_json_node(node, "distribution");
	json_node* m_node = get_json_node(node, "substitutionmodel");
	json_node* mu_node = get_json_node(node, "mu");
	char* sp_ref = get_json_node_value_string(node, "sitepattern");
	json_node* proportions_node = NULL;
	
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
	
	Parameters* rates = new_Parameters(1);
	Simplex* props_simplex = NULL;
	Model* mprops_simplex = NULL;
	
	distribution_t distribution = DISTRIBUTION_UNIFORM;
	quadrature_t quad = QUADRATURE_QUANTILE_MEDIAN;
	bool invariant = false;
	
	if (distribution_node != NULL) {
		invariant = get_json_node_value_bool(distribution_node, "invariant", false);

		proportions_node = get_json_node(distribution_node, "proportions");
		cat = get_json_node_value_int(distribution_node, "categories", 4);

		char* distribution_name = get_json_node_value_string(distribution_node, "distribution");
		if (strcasecmp(distribution_name, "gamma") == 0){
			distribution = DISTRIBUTION_GAMMA;
		}
		else if (strcasecmp(distribution_name, "weibull") == 0){
			distribution = DISTRIBUTION_WEIBULL;
		}
		else if (strcasecmp(distribution_name, "lognormal") == 0){
			distribution = DISTRIBUTION_LOGNORMAL;
		}
		else if (strcasecmp(distribution_name, "discrete") == 0){
			distribution = DISTRIBUTION_DISCRETE;
		}
		else if (strcasecmp(distribution_name, "beta") == 0){
			distribution = DISTRIBUTION_BETA;
		}
		else{
			fprintf(stderr, "Cannot not recognize distribution %s\n", distribution_name);
			exit(13);
		}
		
		json_node* discretization_node = get_json_node(distribution_node, "quadrature");
		
		if (proportions_node != NULL) {
			mprops_simplex = new_SimplexModel_from_json(proportions_node, hash);
			Hashtable_add(hash, mprops_simplex->name, mprops_simplex);
			props_simplex = (Simplex*)mprops_simplex->obj;
			for (int i = 0; i < Parameters_count(props_simplex->parameters); i++) {
				Parameters_at(props_simplex->parameters, i)->model = MODEL_SITEMODEL;
			}
		}
		
		if(discretization_node != NULL && distribution != DISTRIBUTION_DISCRETE){
			char* method = discretization_node->value;
			if(strcasecmp("laguerre", method) == 0){
				quad = QUADRATURE_GAUSS_LAGUERRE;
				if (proportions_node != NULL) {
					fprintf(stderr, "Gauss-Laguerre quadrature does not need proportions to be specified (%s)\n", proportions_node->key);
					exit(13);
				}
			}
			else if(strcasecmp("median", method) == 0){
				quad = QUADRATURE_QUANTILE_MEDIAN;
				if (props_simplex != NULL && cat > props_simplex->K) {
					invariant = true;
				}
			}
			else if(strcasecmp("mean", method) == 0){
				quad = QUADRATURE_QUANTILE_MEAN;
				if (props_simplex != NULL && cat > props_simplex->K) {
					invariant = true;
				}
			}
			else if(strcasecmp("discrete", method) == 0){
				quad = QUADRATURE_DISCRETE;
				if (cat < props_simplex->K) {
					invariant = true;
				}
			}
			else if(strcasecmp("beta", method) == 0 || strcasecmp("kumaraswamy", method) == 0){
				if(strcasecmp("beta", method) == 0){
					quad = QUADRATURE_BETA;
				}
				else if(strcasecmp("kumaraswamy", method) == 0){
					quad = QUADRATURE_KUMARASWAMY;
				}
				else{
					fprintf(stderr, "Could not recognize quadrature type. It should be beta or kumaraswamy");
				}
				if (props_simplex != NULL) {
					if(props_simplex->K != 2){
						fprintf(stderr, "QUADRATURE_BETA: simplex should be of dimension 2 or no simplex at all\n");
						exit(2);
					}
					invariant = true;
				}
			}
			else{
				fprintf(stderr, "Cannot not recognize quadrature method %s\n", method);
				exit(13);
			}
		}
		
		if (distribution == DISTRIBUTION_DISCRETE) {
			json_node* parameters_node = get_json_node(distribution_node, "parameters");
			if(parameters_node!= NULL){
				get_parameters_references(parameters_node, hash, rates);
				
				int i = 0;
				Parameters_set_bounds(rates, i++, 1.e-8, 0.99);
				for ( ; i < Parameters_count(rates); i++) {
					Parameters_set_bounds(rates, i, 1, 100);
				}
				
				char* id_ps = get_json_node_value_string(parameters_node, "id");
				if(id_ps != NULL){
					Parameters_set_name2(rates, id_ps);
				}
			}
			if (parameters_node == NULL || Parameters_count(rates) == props_simplex->K) {
				invariant = true;
			}
		}
		else{
			get_parameters_references(distribution_node, hash, rates);
		}
		
		for (int i = 0; i < Parameters_count(rates); i++) {
			Parameter* p = Parameters_at(rates, i);
			p->model = MODEL_SITEMODEL;
			Hashtable_add(hash, Parameters_name(rates, i), p);
		}
	}
	
	SiteModel* sm = NULL;
	
	if (sp_ref != NULL) {
		int cat = get_json_node_value_int(node, "categories", 4);
		SitePattern* sp = Hashtable_get(hash, sp_ref+1);
		json_node* params_node = get_json_node(node, "parameters");
		get_parameter_list_from_node(params_node, rates);
//		get_parameters_references(node, hash, rates);
		for (int i = 0; i < Parameters_count(rates); i++) {
			Hashtable_add(hash, Parameters_name(rates, i), Parameters_at(rates, i));
		}
		
		sm = new_CATSiteModel_with_parameters(m, rates, cat, sp);
	}
	else {
		sm = new_SiteModel_with_parameters(m, rates, props_simplex, cat, distribution, invariant, quad);
	}
	
	
	if (sm->rates != NULL && Parameters_name2(sm->rates) != NULL) {
		Hashtable_add(hash, Parameters_name2(sm->rates), sm->rates);
	}
	char* id_string = get_json_node_value_string(node, "id");
	
	Model* msm = new_SiteModel2(id_string, sm, mm, mprops_simplex);
	
	if (mu_node != NULL) {
		sm->mu = new_Parameter_from_json(mu_node, hash);
		sm->mu->model = MODEL_SITEMODEL;
		check_constraint(sm->mu, 0, INFINITY, 0.001, 100);
		Hashtable_add(hash, Parameter_name(sm->mu), sm->mu);
	}
	
	mm->free(mm);
	if(mprops_simplex != NULL) mprops_simplex->free(mprops_simplex);
	msm->print = _SiteModel_print;
	free_Parameters(rates);
	return msm;
}
