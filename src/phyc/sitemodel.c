// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "sitemodel.h"

#include <assert.h>
#include <strings.h>
#include <math.h>

#include "parameters.h"
#include "gamma.h"
#include "gausslaguerre.h"
#include "mstring.h"
#include "matrix.h"
#include "mathconstant.h"
#include "gaussian.h"
#include "distkumaraswamy.h"
#include "transforms.h"

#ifndef GSL_DISABLED
#include <gsl/gsl_cdf.h>
#endif

// Default lower bound for the alpha/beta shapes of the Beta quadrature. Chosen just
// above where GSL's Beta inverse-CDF stops converging (it already fails at 0.02, and
// is reliable from ~0.05 upwards over the useful range of the other shape).
#define BETA_QUADRATURE_SHAPE_LOWER 0.1

static bool _gamma_approx_quantile( SiteModel *sm );
static void _calculate_rates_discrete( SiteModel *sm );
static void _calculate_rates_discrete_rate_shape( SiteModel *sm );
static void _calculate_rates_discrete_mean_contribution( SiteModel *sm );
static void _calculate_rates_discrete_increments( SiteModel *sm );
static void _calculate_rates_discrete_ratios( SiteModel *sm );

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

static double icdf_weibull_1(double p, double k);

static bool _cumulative_partial_means(distribution_t distribution, double alpha,
                                      size_t cat_count, double* masses);

// Indexed by rate_parameterization_t. Doubles as the set of JSON keys naming a
// parameterization (rate_parameterization_keys) and as what error messages quote
// back, so a new parameterization must be added to the enum and here together.
static const char* rate_parameterization_strings[] = {
	"auto", "rate_shape", "mean_contribution", "rate_increments", "rate_ratios"
};

const char* SiteModel_rate_parameterization_name(
	rate_parameterization_t parameterization){
	return rate_parameterization_strings[parameterization];
}

static void _site_model_handle_change( Model *self, Model *model, Parameter* parameter, int index ){
	SiteModel *sm = (SiteModel*)self->obj;
	// printf("_site_model_handle_change\n");
	sm->need_update = true;// one of the sitemodel parameters
	self->listeners->fire( self->listeners, self, parameter, index );
}

static void _site_model_store(Model* self){
	if(!self->stored){
		Model *mprop = (Model*)self->data;
		if(mprop != NULL) mprop->store(mprop); // simplex proportion model
		SiteModel* sm = self->obj;
		if (Parameters_count(sm->rates) > 0) {
			Parameters_store(sm->rates);
		}
		if (sm->mu != NULL) {
			Parameter_store(sm->mu);
		}
		if (sm->class_rate != NULL) {
			Parameter_store(sm->class_rate);
		}
		self->stored = true;
	}
}

static void _site_model_restore(Model* self){
	if(self->stored){
		Model *mprop = (Model*)self->data;
		if(mprop != NULL) mprop->restore(mprop); // simplex proportion model
		SiteModel* sm = self->obj;
		sm->need_update = true; // cat_rates and cat_proportions are not stored/restored
		if(Parameters_count(sm->rates) > 0) Parameters_restore(sm->rates);

		if (sm->mu != NULL) {
			Parameter_restore(sm->mu);
		}
		if (sm->class_rate != NULL) {
			Parameter_restore(sm->class_rate);
		}
		self->stored = false;
	}
}

static void _site_model_accept(Model* self){
	if(self->stored){
		Model *mprop = (Model*)self->data;
		if(mprop != NULL) mprop->accept(mprop); // simplex proportion model
		SiteModel* sm = self->obj;
		if(Parameters_count(sm->rates) > 0) Parameters_accept(sm->rates);
		if (sm->mu != NULL) {
			Parameter_accept(sm->mu);
		}
		if (sm->class_rate != NULL) {
			Parameter_accept(sm->class_rate);
		}
		self->stored = false;
	}
}

static void _site_model_free( Model *self ){
#ifdef DEBUG_REF_COUNTING
	printf("Free site model: %d\n", self->ref_count);
#endif
	if(self->ref_count == 1){
		//printf("Free site model %s\n", self->name);
		SiteModel *sm = (SiteModel*)self->obj;
		if(self->data != NULL){
			Model* model = (Model*)self->data;
			model->free(model);
		}
		
		// Delegated rather than open-coded: the two had already drifted apart, and
		// whatever the site model owns is the same either way.
		free_SiteModel(sm);
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
	Model* mprop = (Model*)self->data;
	Model* mpropclone = NULL;
    Parameter* propclone = NULL;

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
	Parameter* class_rate = NULL;
	if (sm->class_rate != NULL) {
		char* name = Parameter_name(sm->class_rate);
		if (Hashtable_exists(hash, name)) {
			class_rate = Hashtable_get(hash, name);
			class_rate->refCount++;
		}
		else{
			class_rate = clone_Parameter(sm->class_rate);
			Hashtable_add(hash, name, class_rate);
		}
	}
	SiteModel* smclone = clone_SiteModel_with_parameters(sm, propclone, ps, mu);
	// Before new_SiteModel2, which registers the model as a listener of it.
	if (class_rate != NULL) {
		SiteModel_set_class_rate(smclone, class_rate,
		                         sm->class_rate_parameterization);
	}
	free_Parameters(ps);
	free_Parameter(mu);
	free_Parameter(class_rate);
	Model* clone = new_SiteModel2(self->name, smclone);
	Hashtable_add(hash, clone->name, clone);
	if(mpropclone != NULL)mpropclone->free(mpropclone);
	clone->print = self->print;
	return clone;
}

static size_t _SiteModel_log_count(Model* model, const char* quantity);
static void _SiteModel_log_name(Model* model, const char* quantity, size_t i,
                                struct StringBuffer* out);
static void _SiteModel_log_value(Model* model, const char* quantity, size_t i,
                                 const char* format, struct StringBuffer* out);

// SubstitutionModel2 listen to the rate and freq parameters
Model * new_SiteModel2( const char* name, SiteModel *sm ){
	Model *model = new_Model(MODEL_SITEMODEL, name, sm);
	if ( sm->rates != NULL ) {
		Parameters_add_listener(sm->rates, model);
		Parameters_add_parameters_recursively(model->parameters, sm->rates);
	}
	if ( sm->mu != NULL ) {
		sm->mu->listeners->add( sm->mu->listeners, model );
		Parameters_add_recursively(model->parameters, sm->mu);
	}
	if ( sm->class_rate != NULL ) {
		sm->class_rate->listeners->add( sm->class_rate->listeners, model );
		Parameters_add_recursively(model->parameters, sm->class_rate);
	}
	if (sm->proportions != NULL) {
		sm->proportions->listeners->add(sm->proportions->listeners, model);
		Parameters_add_recursively(model->parameters, sm->proportions);
	}
	
	model->update = _site_model_handle_change;
	model->store = _site_model_store;
	model->restore = _site_model_restore;
	model->accept = _site_model_accept;
	model->free = _site_model_free;
	model->clone = _site_model_clone;
	model->log_count = _SiteModel_log_count;
	model->log_name = _SiteModel_log_name;
	model->log_value = _SiteModel_log_value;
	return model;
}

static void _SiteModel_print(Model* model, FILE* out){
	SiteModel* sm = model->obj;
	sm->get_proportions(sm); // make sure it is updated
	if (Parameters_count(sm->rates)) {
		for (int i = 0; i < Parameters_count(sm->rates); i++) {
			fprintf(out, "%s %f\n", Parameters_name(sm->rates, i),Parameters_value(sm->rates, i));
		}
		// CAT assigns each site to one category instead of mixing them, so there are
		// no weights to report and nothing that sums to one; what it does have is
		// the number of patterns that landed in each category.
		if (sm->site_category != NULL) {
			int* counts = ivector(sm->cat_count);
			for (int i = 0; i < sm->sp->count; i++){
				counts[sm->get_site_category(sm, i)]++;
			}
			fprintf(out, "patterns rate\n");
			for (int i = 0; i < sm->cat_count; i++){
				fprintf(out, "%d %f\n", counts[i], sm->cat_rates[i]);
			}
			free(counts);
			return;
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
		fprintf(out, "proportion rate %lu\n", Parameter_size(sm->proportions));
		for (int i = 0; i < Parameter_size(sm->proportions); i++){
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

// Loggable interface: the per-category "rates" and "proportions", each yielding
// cat_count columns. Column headers are 0-based (rates.0 … rates.{K-1}), matching
// the multi-element Parameter convention.
static size_t _SiteModel_log_count(Model* model, const char* quantity){
	SiteModel* sm = model->obj;
	if(strcasecmp(quantity, "rates") == 0 || strcasecmp(quantity, "proportions") == 0){
		return sm->cat_count;
	}
	return 0;
}

static void _SiteModel_log_name(Model* model, const char* quantity, size_t i,
                                struct StringBuffer* out){
	// Prefix with the model id so columns stay unique across several site models,
	// e.g. "sitemodel.rates.0".
	StringBuffer_append_format(out, "%s.%s.%zu", model->name, quantity, i);
}

static void _SiteModel_log_value(Model* model, const char* quantity, size_t i,
                                 const char* format, struct StringBuffer* out){
	SiteModel* sm = model->obj;
	sm->get_proportions(sm);  // refresh cat_rates/cat_proportions
	double value = strcasecmp(quantity, "rates") == 0 ? sm->get_rate(sm, i)
	                                                   : sm->get_proportion(sm, i);
	StringBuffer_append_format(out, format != NULL ? format : "%f", value);
}

int _get_site_category(SiteModel* sm, const int pattern){
	return 0;
}

void _no_gradient( SiteModel *sm, const double* ingrad, double* grad ){

}

double _no_derivative( SiteModel *sm, const double* ingrad, Parameter* p ){
	return 0;
}

double _gamma_shape_derivative( SiteModel *sm, const double* ingrad ){
	double derivsum = 0;
	double sum = 0;
	size_t catCount = sm->cat_count;
	size_t variableCat = sm->proportions != NULL ? catCount - 1 : catCount;
	double* proportions = sm->get_proportions(sm);
	double shape = Parameters_value(sm->rates, 0);
	size_t j = (variableCat == catCount ? 0 : 1); // ignore category 0 with r_0=0
	double* temp = dvector(variableCat*2);
	const double eps = pow(DBL_EPSILON, 1.0/3.0);
#ifndef GSL_DISABLED
	gsl_error_handler_t* handler = gsl_set_error_handler_off();
#endif
	size_t i = 0;
	for ( ; i < variableCat; i++, j++) {
		double prob = (2.0*i + 1.0)/(2.0*variableCat);
		double xp = shape*(1.0 + eps);
		double xm = shape*(1.0 - eps);
#ifndef GSL_DISABLED
		double plus = gsl_cdf_gamma_Qinv( prob, xp, 1.0/xp );
		double minus = gsl_cdf_gamma_Qinv( prob, xm, 1.0/xm );
		temp[i+variableCat] = gsl_cdf_gamma_Qinv( prob, shape, 1.0/shape);
#else
		double plus = qgamma( prob, xp, xp );
		double minus = qgamma( prob, xm, xm );
		temp[i+variableCat] = qgamma( prob, shape, shape);
#endif
		temp[i] = (plus - minus)/(xp-xm);
		if(!isfinite(plus) || !isfinite(minus) || !isfinite(temp[i+variableCat])){
			break;
		}
		
		derivsum += temp[i] * proportions[j];
		sum += temp[i+variableCat] * proportions[j];
	}
#ifndef GSL_DISABLED
	gsl_set_error_handler(handler);
#endif
	double shape_gradient = NAN;
	if(i == variableCat){
		shape_gradient = 0.0;
		j = (variableCat == catCount ? 0 : 1); // ignore category 0 with r_0=0
		for (size_t i = 0; i < variableCat; i++, j++) {
			double prob = (2.0*i + 1.0)/(2.0*variableCat);
			double deriv_rate = temp[i]/sum - (temp[i+variableCat]*derivsum)/sum/sum;
			shape_gradient += ingrad[j] * deriv_rate * proportions[j];
		}
	}
	free(temp);
	return shape_gradient;
}

double _gamma_inv_derivative( SiteModel *sm, const double* ingrad ){
	double* proportions = sm->get_proportions(sm);
	// assert(proportions[0] == Parameters_value(sm->proportions->parameters, 0));
	
	// with Gamma
	if(Parameters_count(sm->rates) > 0){
		size_t variableCatCount = sm->cat_count - 1;
		double shape = Parameters_value(sm->rates, 0);
		double sum_rate = 0;
		for (size_t i = 0; i < variableCatCount; i++) {
			double prob = (2.0*i + 1.0)/(2.0*variableCatCount);
#ifndef GSL_DISABLED
			sum_rate += gsl_cdf_gamma_Qinv( prob, shape, 1.0/shape );
#else
			sum_rate += qgamma( prob, shape, shape );
#endif
		}
		double pinv_gradient = 0;
		for (size_t i = 0; i < variableCatCount; i++) {
			double prob = (2.0*i + 1.0)/(2.0*variableCatCount);
#ifndef GSL_DISABLED
			pinv_gradient += ingrad[i+1] * gsl_cdf_gamma_Qinv( prob, shape, 1.0/shape );
#else
			pinv_gradient += ingrad[i+1] * qgamma( prob, shape, shape );
#endif
		}
		return ingrad[0] + pinv_gradient/(sum_rate*(1.0 - proportions[0]));
	}
	// +I only
	else{
		return ingrad[0] + ingrad[1]/proportions[1];
	}
}

double _gamma_derivative( SiteModel *sm, const double* ingrad, Parameter* p ){
	// pinv (with no rate categories, p can only be the proportions parameter)
	if(sm->proportions != NULL &&
	   (Parameters_count(sm->rates) == 0 || Parameters_at(sm->rates, 0) != p)){
		return _gamma_inv_derivative(sm, ingrad);
	}
	else{
		return _gamma_shape_derivative(sm, ingrad);
	}
	return 0;
}

void _gamma_gradient( SiteModel *sm, const double* ingrad, double* grad ){
	size_t offset = 0;
	if(Parameters_count(sm->rates) > 0){
		grad[offset++] = _gamma_shape_derivative(sm, ingrad);
	}
	if(sm->proportions != NULL){
		grad[offset] = _gamma_inv_derivative(sm, ingrad);
	}
}


double _weibull_shape_derivative( SiteModel *sm, const double* ingrad ){
	double derivsum = 0;
	double sum = 0;
	size_t catCount = sm->cat_count;
	size_t variableCat = sm->proportions != NULL ? catCount - 1 : catCount;
	double* proportions = sm->get_proportions(sm);
	double shape = Parameters_value(sm->rates, 0);
	double shape2 = shape*shape;
	size_t j = (variableCat == catCount ? 0 : 1); // ignore category 0 with r_0=0
	for (size_t i = 0; i < variableCat; i++, j++) {
		double prob = (2.0*i + 1.0)/(2.0*variableCat);
		derivsum += -pow(-log(1.0 - prob), 1.0/shape) * log(-log(1.0 - prob))/shape2*proportions[j];
		sum += pow(-log(1.0 - prob), 1.0/shape)*proportions[j];
	}

	double shape_gradient = 0;
	j = (variableCat == catCount ? 0 : 1); // ignore category 0 with r_0=0
	for (size_t i = 0; i < variableCat; i++, j++) {
		double prob = (2.0*i + 1.0)/(2.0*variableCat);
		double deriv = -pow(-log(1.0 - prob), 1.0/shape) * log(-log(1.0 - prob))/shape2;
		double deriv_rate = deriv/sum - (pow(-log(1.0 - prob), 1.0/shape)*derivsum)/sum/sum;
		shape_gradient += ingrad[j] * deriv_rate * proportions[j];
	}
	return shape_gradient;
}

double _weibull_inv_derivative( SiteModel *sm, const double* ingrad ){
	double* proportions = sm->get_proportions(sm);
	// assert(proportions[0] == Parameters_value(sm->proportions->parameters, 0));
	
	// with Weibull
	if(Parameters_count(sm->rates) > 0){
		size_t variableCatCount = sm->cat_count - 1;
		double shape = Parameters_value(sm->rates, 0);
		double sum_rate = 0;
		for (size_t i = 0; i < variableCatCount; i++) {
			double prob = (2.0*i + 1.0)/(2.0*variableCatCount);
			sum_rate += icdf_weibull_1(prob, shape);
		}
		double pinv_gradient = 0;
		for (size_t i = 0; i < variableCatCount; i++) {
			double prob = (2.0*i + 1.0)/(2.0*variableCatCount);
			pinv_gradient += ingrad[i+1] * icdf_weibull_1(prob, shape);
		}
		return ingrad[0] + pinv_gradient/(sum_rate*(1.0 - proportions[0]));
	}
	// +I only
	else{
		return ingrad[0] + ingrad[1]/proportions[1];
	}
}

double _weibull_derivative( SiteModel *sm, const double* ingrad, Parameter* p ){
	// pinv (with no rate categories, p can only be the proportions parameter)
	if(sm->proportions != NULL &&
	   (Parameters_count(sm->rates) == 0 || Parameters_at(sm->rates, 0) != p)){
		return _weibull_inv_derivative(sm, ingrad);
	}
	else{
		return _weibull_shape_derivative(sm, ingrad);
	}
	return 0;
}
void _weibull_gradient( SiteModel *sm, const double* ingrad, double* grad ){
	size_t offset = 0;
	if(Parameters_count(sm->rates) > 0){
		grad[offset++] = _weibull_shape_derivative(sm, ingrad);
	}
	if(sm->proportions != NULL){
		grad[offset] = _weibull_inv_derivative(sm, ingrad);
	}
}

// Mean quadrature places category k at the conditional mean of its bin,
// r_k = K (c_k - c_{k-1})/propVariable, where the c_i are the cumulative partial
// expectations of the rate distribution. Neither family gives an elementary
// derivative of c_i with respect to the shape -- the gamma boundary is an inverse
// CDF, and the Weibull shape enters the order of an incomplete gamma -- so the
// masses are differentiated by a central difference, the same treatment the median
// quadrature gives the gamma quantile. Unlike the median rates these need no
// renormalization, so there is no quotient rule here.
double _mean_quadrature_shape_derivative( SiteModel *sm, const double* ingrad ){
	size_t catCount = sm->cat_count;
	size_t variableCat = sm->proportions != NULL ? catCount - 1 : catCount;
	double* proportions = sm->get_proportions(sm);
	double propVariable = sm->proportions != NULL ? 1.0 - proportions[0] : 1.0;
	double shape = Parameters_value(sm->rates, 0);
	const double eps = pow(DBL_EPSILON, 1.0/3.0);
	double xp = shape*(1.0 + eps);
	double xm = shape*(1.0 - eps);

	double* plus = dvector(variableCat + 1);
	double* minus = dvector(variableCat + 1);
	double shape_gradient = NAN;
	if(_cumulative_partial_means(sm->distribution, xp, variableCat, plus) &&
	   _cumulative_partial_means(sm->distribution, xm, variableCat, minus)){
		shape_gradient = 0.0;
		size_t j = (variableCat == catCount ? 0 : 1); // ignore category 0 with r_0=0
		for (size_t i = 0; i < variableCat; i++, j++) {
			double deriv_rate = variableCat*((plus[i+1] - plus[i]) - (minus[i+1] - minus[i]))
			                    /((xp - xm)*propVariable);
			shape_gradient += ingrad[j] * deriv_rate * proportions[j];
		}
	}
	free(plus);
	free(minus);
	return shape_gradient;
}

double _mean_quadrature_inv_derivative( SiteModel *sm, const double* ingrad ){
	double* proportions = sm->get_proportions(sm);
	// The invariant proportion only scales the conditional means, r_k = K (c_k -
	// c_{k-1})/(1 - p0), so dr_k/dp0 = r_k/(1 - p0).
	double propVariable = 1.0 - proportions[0];
	double pinv_gradient = 0;
	for (size_t i = 1; i < sm->cat_count; i++) {
		pinv_gradient += ingrad[i] * sm->cat_rates[i] * proportions[i]/propVariable;
	}
	return ingrad[0] + pinv_gradient;
}

double _mean_quadrature_derivative( SiteModel *sm, const double* ingrad, Parameter* p ){
	if(sm->proportions != NULL && Parameters_at(sm->rates, 0) != p){
		return _mean_quadrature_inv_derivative(sm, ingrad);
	}
	return _mean_quadrature_shape_derivative(sm, ingrad);
}

void _mean_quadrature_gradient( SiteModel *sm, const double* ingrad, double* grad ){
	size_t offset = 0;
	grad[offset++] = _mean_quadrature_shape_derivative(sm, ingrad);
	if(sm->proportions != NULL){
		grad[offset] = _mean_quadrature_inv_derivative(sm, ingrad);
	}
}

bool _update_gamma_approx_quantile(SiteModel *sm){
	if ( sm->need_update ) {
        return _gamma_approx_quantile(sm);
    }
	return true;
}

bool _update_nothing(SiteModel *sm){
	return true;
}

// (Gamma or/and Invariant) or one rate
// should not be used directly
SiteModel * new_SiteModel_with_parameters( const Parameters *params, Parameter* proportions, const size_t cat_count, distribution_t distribution, bool invariant, quadrature_t quad, rate_parameterization_t rate_parameterization){
	SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	sm->site_category = NULL;
	sm->sp = NULL;
	sm->get_site_category = _get_site_category;

	sm->distribution = distribution;
	sm->invariant = invariant;
	sm->quadrature = quad;
	sm->rate_parameterization = rate_parameterization;

	sm->cat_count = cat_count;
	if (invariant) sm->cat_count++;
	
	sm->cat_rates       = NULL;
	sm->cat_proportions = NULL;
	sm->proportions = proportions;
	if(proportions != NULL){
		Parameter_set_model(proportions, MODEL_SITEMODEL);
	}

	// sm->epsilon = 1.e-6;
    
	sm->rates = NULL;
	
	if(Parameters_count(params) > 0){
		sm->rates = new_Parameters(Parameters_count(params));
		Parameters_add_parameters(sm->rates, params);
		if (Parameters_name2(params) != NULL) {
			Parameters_set_name2(sm->rates, Parameters_name2(params));
		}
	}
	for(size_t i = 0; i < Parameters_count(params); i++){
		Parameter_set_model(Parameters_at(sm->rates, i), MODEL_SITEMODEL);
	}
	sm->mu    = NULL;
	sm->class_rate = NULL;
	sm->class_rate_parameterization = FREE_CLASS_PARAMETERIZATION_RATE;

	sm->cat_rates       = dvector(sm->cat_count);
	sm->cat_proportions = dvector(sm->cat_count);
	sm->gradient = _no_gradient;
	sm->derivative = _no_derivative;
	sm->update = _update_nothing;

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
		sm->update = _update_gamma_approx_quantile;
	}
    
	// Discrete gamma or prop of invariant or both.
	// The DISTRIBUTION_DISCRETE case is the pure +I model, whose only parameter is
	// the invariant proportion (_weibull_inv_derivative then takes its "+I only"
	// branch). A discrete model with free rates only happens to share the shape of
	// that model when it has two categories; none of the free-rate
	// parameterizations has an analytic derivative here, so it must keep
	// _no_gradient like every other category count.
	// Mean quadrature builds its categories from partial expectations rather than
	// quantiles, so both families share a derivative of their own; the median
	// formulas below would silently return the wrong number for it.
	if (quad == QUADRATURE_QUANTILE_MEAN &&
		(distribution == DISTRIBUTION_GAMMA || distribution == DISTRIBUTION_WEIBULL)) {
		sm->gradient = _mean_quadrature_gradient;
		sm->derivative = _mean_quadrature_derivative;
	}
	else if (distribution == DISTRIBUTION_WEIBULL ||
		(distribution == DISTRIBUTION_DISCRETE && Parameters_count(sm->rates) == 0 &&
		 Parameter_size(sm->proportions) == 2)) {
		sm->gradient = _weibull_gradient;
		sm->derivative = _weibull_derivative;
	}
	else if (distribution == DISTRIBUTION_GAMMA){
		sm->gradient = _gamma_gradient;
		sm->derivative = _gamma_derivative;
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

// Cumulative partial expectations c_i = \int_0^{b_i} r f(r) dr of a unit-mean rate
// distribution at the equal-probability bin boundaries b_i = F^{-1}(i/K), for
// i = 0..K, so that c_0 = 0, c_K = 1 and the conditional mean of bin k is
// K (c_k - c_{k-1}). Those rates satisfy the unit-mean constraint by construction:
// the sum telescopes. Returns false if any mass came out non-finite, which the
// caller reports rather than propagating into the likelihood.
bool _cumulative_partial_means(distribution_t distribution, double alpha,
                               size_t cat_count, double* masses){
	masses[0] = 0.0;
	masses[cat_count] = 1.0;
	if(distribution == DISTRIBUTION_GAMMA){
		// For Gamma(alpha, 1/alpha) the partial expectation is the regularized
		// incomplete gamma one shape up, P(alpha+1, alpha b), and the boundary b
		// itself has to be found by inverting the CDF at every alpha.
		for (size_t i = 1; i < cat_count; i++) {
			double boundary = qgamma((double)i/cat_count, alpha, alpha);
			if(!isfinite(boundary)) return false;
			masses[i] = gammp(alpha + 1.0, boundary*alpha);
		}
	}
	else{
		// The substitution u = (r/lambda)^alpha maps the Weibull to a standard
		// exponential, under which the boundaries u_i = -log(1 - i/K) are constants
		// and the unit-mean scale lambda = 1/Gamma(1 + 1/alpha) cancels. The shape
		// only enters through the order 1 + 1/alpha of the incomplete gamma, so no
		// quantile inversion is needed -- and, unlike the Weibull median rates, no
		// quantile can overflow on the way.
		double order = 1.0 + 1.0/alpha;
		for (size_t i = 1; i < cat_count; i++) {
			masses[i] = gammp(order, -log1p(-(double)i/cat_count));
		}
	}
	for (size_t i = 1; i < cat_count; i++) {
		if(!isfinite(masses[i])) return false;
	}
	return true;
}

// Category 0 of an invariant-shaped model is a point mass pinned at rate 0. A
// free class (+F) unpins it: the mass sits at r_f and contributes c = p_f r_f to
// the unit mean, so the variable categories -- which the quadrature has just
// normalised to carry all of it, sum_{k>0} p_k r_k = 1 -- are scaled down to
// carry 1 - c instead. That single factor is the whole difference between +I and
// +F, whatever quadrature produced the rates, because every branch above
// imposes the same constraint (docs/models/sitemodel.md; the derivation is in
// sitemodel.pdf, "Gamma plus a free rate class").
//
// The increment chart asks for r_f = A*g + d instead, where g is the largest of
// the rates the quadrature produced and A = 1 - c is the very factor being
// solved for -- the class is placed relative to a bulk that its own weight moves.
// Substituting c = p_f r_f and solving the resulting linear equation gives
//
//     r_f = (g + d)/(1 + p_f g),   A = (1 - p_f d)/(1 + p_f g),
//
// so the domain is p_f d < 1 rather than p_f r_f < 1, and d = 0 puts the class
// exactly on the fastest variable category (r_f = A*g) rather than at 0.
//
// Returns false, exactly as a non-finite quantile does, when the parameters
// leave the domain: c >= 1 would make the variable rates negative or zero.
static bool _apply_free_class( SiteModel *sm ) {
	// Without an invariant class there is no category 0 to unpin: cat_count does
	// not include one, so scaling from index 1 and writing index 0 would rewrite a
	// perfectly good +G. SiteModel_set_class_rate asserts this, but the assert is
	// compiled out of a release build.
	if (sm->class_rate == NULL || !sm->invariant) return true;
	const double proportion = sm->cat_proportions[0];
	const double x = Parameter_value(sm->class_rate);
	double rate, contribution;
	switch (sm->class_rate_parameterization) {
		case FREE_CLASS_PARAMETERIZATION_CONTRIBUTION:
			contribution = x;
			rate = contribution/proportion;
			break;
		case FREE_CLASS_PARAMETERIZATION_INCREMENT: {
			// The fastest variable rate. Taken as a maximum rather than as the last
			// category because only the quantile rules are guaranteed to come out
			// ordered, and a caller's simplex could reorder the "discrete" bins.
			double g = 0.0;
			for (size_t i = 1; i < (size_t)sm->cat_count; i++) {
				if (sm->cat_rates[i] > g) g = sm->cat_rates[i];
			}
			rate = (g + x)/(1.0 + proportion*g);
			contribution = proportion*rate;
			break;
		}
		case FREE_CLASS_PARAMETERIZATION_RATE:
		default:
			rate = x;
			contribution = proportion*rate;
			break;
	}
	// A weight that has underflowed to 0 makes r_f = c/p_f infinite; the
	// finiteness test catches that as well as a non-finite parameter.
	if (!(contribution < 1.0) || !isfinite(rate) || rate < 0.0) return false;
	const double scale = 1.0 - contribution;
	for (size_t i = 1; i < (size_t)sm->cat_count; i++) {
		sm->cat_rates[i] *= scale;
	}
	sm->cat_rates[0] = rate;
	return true;
}

#ifndef GSL_DISABLED
// Fill quantiles[offset .. offset+count-1] with the Beta inverse-CDF evaluated at
// i/denominator. GSL's betainv is iterative and gives up ("inverse failed to
// converge") for shape pairs that are small or wildly asymmetric -- alpha ~ 0.02 is
// already enough, and a line search on alpha or beta walks straight into that region.
// Under the default GSL handler that failure aborts the process, so the handler is
// disabled here and the failure reported the same way the gamma inverse-CDF below
// reports one: return false, which _update_gamma_approx_quantile propagates to the
// tree likelihood as a NaN so the optimizer simply rejects the point.
//
// A non-converged call can also return a finite but meaningless value, so the
// quantiles are checked for finiteness and monotonicity: the caller turns successive
// differences into category proportions, and a non-monotonic sequence would silently
// produce negative proportions.
static bool _beta_quantiles(double* quantiles, int offset, int count, int denominator,
                            double alpha, double beta) {
	gsl_error_handler_t* handler = gsl_set_error_handler_off();
	bool ok = true;
	double previous = 0.0;
	for (int i = 0; i < count; i++) {
		const double q = gsl_cdf_beta_Pinv((double)i / denominator, alpha, beta);
		if (!isfinite(q) || q < previous || q > 1.0) {
			ok = false;
			break;
		}
		quantiles[i + offset] = q;
		previous = q;
	}
	gsl_set_error_handler(handler);
	return ok;
}
#endif

bool _gamma_approx_quantile( SiteModel *sm ) {
	double propVariable = 1.0;
	int cat = (sm->invariant ? 1 : 0);
	const int nCat = sm->cat_count - cat;
	double* quantiles = dvector(sm->cat_count); // cat_count includes invariant sites if specified

#ifndef GSL_DISABLED
	// proportions are estimated from beta distribution
	if (sm->quadrature == QUADRATURE_BETA || sm->quadrature == QUADRATURE_KUMARASWAMY){
		double shape_alpha = Parameters_value(sm->rates, Parameters_count(sm->rates)-2);
		double shape_beta = Parameters_value(sm->rates, Parameters_count(sm->rates)-1);
		
		// proportion of invariant is estimated from beta
		// cat can be equal to 0 or 1
		if(sm->proportions == NULL){
			if(sm->quadrature == QUADRATURE_BETA){
				if(!_beta_quantiles(quantiles, 0, sm->cat_count, sm->cat_count,
				                    shape_alpha, shape_beta)){
					free(quantiles);
					sm->need_update = false;
					return false;
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
			if(!_beta_quantiles(quantiles, cat, sm->cat_count - cat,
			                    sm->cat_count - cat, shape_alpha, shape_beta)){
				free(quantiles);
				sm->need_update = false;
				return false;
			}
			
			const double* proportions = Parameter_values(sm->proportions);
			propVariable = proportions[1];
			sm->cat_proportions[0] = proportions[0];
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
#else
	if(sm->quadrature == QUADRATURE_DISCRETE){
#endif
		const double* values = Parameter_values(sm->proportions);
		propVariable = 0;
		for (int i = cat; i < sm->cat_count; i++) {
			propVariable += values[i];
		}
		double sum = 0;
		for (int i = cat; i < sm->cat_count; i++) {
			quantiles[i] = (sum + values[i]/2.0)/propVariable; // pick midpoint
			sum += values[i];
		}
		memcpy(sm->cat_proportions, values, sizeof(double)*Parameter_size(sm->proportions));
	}
	// if the dimension of the simplex is then only pinv is estimated
	// and the remaining categories split 1.0-pinv equally
	// That's the traditional +G+I or +I
	else if (sm->proportions != NULL && Parameter_size(sm->proportions) == 2) {
		const double* values = Parameter_values(sm->proportions);
		// +I
		if(sm->distribution == DISTRIBUTION_DISCRETE){
			memcpy(sm->cat_proportions, values, sizeof(double)*2);
			sm->cat_rates[0] = 0;
			sm->cat_rates[1] = 1.0 / sm->cat_proportions[1];
			sm->need_update = false;
			free(quantiles);
			// A free class here makes this a two-point mixture rather than +I;
			// the JSON parser does not offer it (that model is +R with two
			// categories), but the C API can set one and the scaling holds.
			return _apply_free_class(sm);
		}
		// +G+I
		else{
			cat = 1;
			sm->cat_proportions[0] = values[0];
			propVariable = values[1];
			int gammaCat = sm->cat_count - 1;
			for (int i = 0; i < gammaCat; i++) {
				quantiles[i+cat] = (2.0 * i + 1.0) / (2.0 * gammaCat);
				sm->cat_proportions[i+cat] = propVariable/gammaCat;
			}
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
		size_t i = 0;
		// Clamp quantile probabilities away from {0,1} before feeding them to the
		// inverse-CDFs below. When the Beta/Kumaraswamy quadrature collapses (e.g.
		// the Kumaraswamy "a" driven towards 0 during optimization) the interior
		// quantiles underflow to exactly 0.0, and gsl_cdf_gamma_Qinv(0,...) returns
		// +Inf. That Inf later multiplies a 0 category proportion in the mean,
		// producing a NaN that silently poisons the likelihood. Clamping keeps the
		// rate large but finite; categories with ~0 proportion then contribute ~0.
		const double quantile_eps = 1e-12;
		for (size_t j = cat; j < (size_t)sm->cat_count; j++) {
			if (quantiles[j] < quantile_eps) quantiles[j] = quantile_eps;
			else if (quantiles[j] > 1.0 - quantile_eps) quantiles[j] = 1.0 - quantile_eps;
		}
		if(sm->distribution == DISTRIBUTION_GAMMA){
#ifndef GSL_DISABLED
			gsl_error_handler_t* handler = gsl_set_error_handler_off();
			for ( ; i < sm->cat_count - cat; i++) {
				sm->cat_rates[i + cat] = gsl_cdf_gamma_Qinv( quantiles[i+cat], alpha, 1.0/alpha );
				// printf("quantiles[%d] = %f, alpha = %f, sm->cat_rates[%d] = %f, sm->cat_proportions[%d] = %f\n", i+cat, quantiles[i+cat], alpha, i+cat, sm->cat_rates[i + cat], i+cat, sm->cat_proportions[i + cat]);
				// A failing quantile shows up as +Inf at least as often as NaN, and an
				// Inf survives to poison the mean below (Inf/Inf = NaN) instead of
				// being caught here, so test for any non-finite value.
				if(!isfinite(sm->cat_rates[i + cat])){
					break;
				}
			}
			gsl_set_error_handler(handler);
#else
			for ( ; i < sm->cat_count - cat; i++) {
				sm->cat_rates[i + cat] = qgamma( quantiles[i+cat], alpha, alpha );
				if(!isfinite(sm->cat_rates[i + cat])){
					break;
				}
			}
#endif
		}
		else if(sm->distribution == DISTRIBUTION_WEIBULL){
			for ( ; i < sm->cat_count - cat; i++) {
				// Unit mean Weibull
				// sm->cat_rates[i + cat] = gsl_cdf_weibull_Qinv( sm->cat_proportions[i + cat], alpha, 1.0/exp(gammln(1.0 + 1.0/alpha)) );
				// Fix lambda:=1
				sm->cat_rates[i + cat] = icdf_weibull_1( quantiles[i+cat], alpha);
				if(!isfinite(sm->cat_rates[i + cat])){
					break;
				}

			}
		}
#ifndef GSL_DISABLED		
		else if(sm->distribution == DISTRIBUTION_LOGNORMAL){
			for ( ; i < sm->cat_count - cat; i++) {
				sm->cat_rates[i + cat] = gsl_cdf_lognormal_Qinv( quantiles[i+cat], -alpha*alpha/2, alpha );
			}
		}
		else if(sm->distribution == DISTRIBUTION_BETA){
			for ( ; i < sm->cat_count - cat; i++) {
				sm->cat_rates[i + cat] = gsl_cdf_beta_Qinv( quantiles[i+cat], alpha, Parameters_value(sm->rates, 1) );
			}
		}
#endif

		if(i != sm->cat_count - cat){
			free(quantiles);
			sm->need_update = false;
			return false;
		}
		
		if (sm->quadrature == QUADRATURE_BETA || sm->quadrature == QUADRATURE_KUMARASWAMY) {
			for (int i = 0; i < sm->cat_count - cat; i++) {
				mean += sm->cat_rates[i + cat]*sm->cat_proportions[i + cat];
			}
		}
		// QUADRATURE_DISCRETE
		else if ( (sm->proportions != NULL && Parameter_size(sm->proportions) == sm->cat_count)){
			for (int i = 0; i < nCat; i++) {
				mean += sm->cat_rates[i + cat]*sm->cat_proportions[i + cat];
			}
		}
		// +G or +G+I
		else if(sm->proportions == NULL || (sm->proportions != NULL && Parameter_size(sm->proportions) != sm->cat_count)){
			for (int i = 0; i < sm->cat_count - cat; i++) {
				mean += sm->cat_rates[i + cat];
			}
			mean = (propVariable * mean) / nCat;
		}
		
		for (int i = 0; i < sm->cat_count - cat; i++) {
			sm->cat_rates[i + cat] /= mean;
		}
	}
	// mean: category k is represented by the conditional mean of the distribution over
	// its bin instead of by a quantile. The unit-mean constraint then holds by
	// construction, so there is no renormalization by the sample mean.
	else{
		double* masses = dvector(nCat + 1);
		bool ok = _cumulative_partial_means(sm->distribution, alpha, nCat, masses);
		if(ok){
			sm->cat_rates[0] = 0;
			for (int i = 0; i < nCat; i++) {
				// The invariant class contributes nothing to the weighted mean, so the
				// constraint reads sum_{k>0} p_k r_k = 1 with p_k = propVariable/K and
				// the conditional means are scaled up by the variable proportion.
				sm->cat_rates[i + cat] = (masses[i + 1] - masses[i])*nCat/propVariable;
				sm->cat_proportions[i + cat] = propVariable/nCat;
			}
		}
		free(masses);
		free(quantiles);
		sm->need_update = false;
		return ok && _apply_free_class(sm);
	}
	free(quantiles);
	sm->need_update = false;
	return _apply_free_class(sm);
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

void SiteModel_set_mu(SiteModel *sm, Parameter* mu){
	sm->mu = mu;
	Parameter_set_model(sm->mu, MODEL_SITEMODEL);
	mu->refCount++;
}

void SiteModel_set_class_rate(SiteModel* sm, Parameter* class_rate,
                              free_class_parameterization_t parameterization){
	assert(sm->proportions != NULL);
	sm->class_rate = class_rate;
	sm->class_rate_parameterization = parameterization;
	Parameter_set_model(sm->class_rate, MODEL_SITEMODEL);
	class_rate->refCount++;
	// The analytic derivatives installed by new_SiteModel_with_parameters are the
	// +I ones: they assume category 0 sits at rate 0 and that the variable rates
	// carry the whole unit mean. Both are false here, and a wrong gradient is
	// worse than none, so hand the model back to the derivative-free optimizers
	// until the +F derivatives are implemented.
	sm->gradient = _no_gradient;
	sm->derivative = _no_derivative;
}

#pragma mark -
// MARK: Discrete SiteModel

static void _update_rates_discrete(SiteModel *sm) {
    if (!sm->need_update) return;
    switch (sm->rate_parameterization) {
        case RATE_PARAMETERIZATION_RATE_SHAPE:
            _calculate_rates_discrete_rate_shape(sm);
            return;
        case RATE_PARAMETERIZATION_MEAN_CONTRIBUTION:
            _calculate_rates_discrete_mean_contribution(sm);
            return;
        case RATE_PARAMETERIZATION_RATE_INCREMENTS:
            _calculate_rates_discrete_increments(sm);
            return;
        case RATE_PARAMETERIZATION_RATE_RATIOS:
            _calculate_rates_discrete_ratios(sm);
            return;
        case RATE_PARAMETERIZATION_AUTO:
            break;
    }
    // A C API caller may leave the parameterization unset, in which case it is
    // inferred from the shape of the rate parameters: a simplex is a rate shape, a
    // plain vector a sequence of increments, and a pair of parameters the ratios and
    // the top rate they hang off. The mean-contribution simplex is not reachable this
    // way -- its parameter looks exactly like a rate shape -- which is why the JSON
    // parser names the parameterization instead of inferring it.
    if (Parameters_count(sm->rates) == 1) {
        if (Parameters_at(sm->rates, 0)->simplex) _calculate_rates_discrete_rate_shape(sm);
        else _calculate_rates_discrete_increments(sm);
    } else {
        _calculate_rates_discrete_ratios(sm);
    }
}

double _get_rate_discrete( SiteModel *sm, const int index ){
	_update_rates_discrete(sm);
	return sm->cat_rates[index];
}

double _get_proportion_discrete( SiteModel *sm, const int index ){
	_update_rates_discrete(sm);
	return sm->cat_proportions[index];
}

double *_get_proportions_discrete( SiteModel *sm ){
	_update_rates_discrete(sm);
	return sm->cat_proportions;
}

void _calculate_rates_discrete( SiteModel *sm ) {
	int j = 0;
	int cat = (sm->invariant ? 1 : 0);
	int cat_count = Parameter_size(sm->proportions);
	memset(sm->cat_rates, 0, sizeof(double)*cat_count);
	const double* cat_proportions = Parameter_values(sm->proportions);
	memcpy(sm->cat_proportions, cat_proportions, sizeof(double)*cat_count);

	sm->cat_rates[cat] = Parameters_value(sm->rates, j++);
	sm->cat_rates[cat+1] = 1;
	double sum = sm->cat_rates[cat]*sm->cat_proportions[cat] + sm->cat_proportions[cat+1];
	if(sm->invariant) sum += sm->cat_proportions[0];
	double prod = 1;
	const double* rates_multipliers = Parameter_values(Parameters_at(sm->rates, 1));
	for (int i = cat+2; i < cat_count; i++, j++ ) {
		prod *= rates_multipliers[j];
		sm->cat_rates[i] = prod;
		sum += prod*sm->cat_proportions[i];
	}
	
	for (int i = 0; i < cat_count; i++ ) {
		sm->cat_rates[i] /= sum;
		// printf("cat_rates[%d] = %f\n", i, sm->cat_rates[i]);
	}
	sm->need_update = false;
}

// Rate-shape simplex ("rate_shape"; see docs/models/sitemodel.md). The free
// parameter is the simplex x of relative rates; the unit mean is imposed
// afterwards by dividing by the weighted mean,
//
//     r_k = x_k / sum_j p_j x_j .
//
// Compare _calculate_rates_discrete_mean_contribution, its dual.
void _calculate_rates_discrete_rate_shape( SiteModel *sm ) {
	int cat_count = Parameter_size(sm->proportions);
	memset(sm->cat_rates, 0, sizeof(double)*cat_count);
	const double* cat_proportions = Parameter_values(sm->proportions);
	memcpy(sm->cat_proportions, cat_proportions, sizeof(double)*cat_count);

	Parameter* parameter_simplex = Parameters_at(sm->rates, 0);
	const double* cat_rates = Parameter_values(parameter_simplex);
	memcpy(sm->cat_rates, cat_rates, sizeof(double)*cat_count);
	
	double norm = 0;
	for (size_t i = 0; i < cat_count; i++ ) {
		norm += sm->cat_rates[i]*sm->cat_proportions[i];
	}
	
	for (size_t i = 0; i < cat_count; i++ ) {
		sm->cat_rates[i] /= norm;
		// printf("cat_rates[%zu] = %f\n", i, sm->cat_rates[i]);
	}
	sm->need_update = false;
}

// Mean-contribution simplex ("mean_contribution"; see docs/models/sitemodel.md).
// The free parameter is the simplex s of per-category contributions to the mean,
// s_k = p_k r_k. Because sum_k s_k = 1 by construction, so is the unit mean
// sum_k p_k r_k, and the rates are recovered by dividing by the individual
// proportions:
//
//     r_k = s_k / p_k .
//
// This is the dual of _calculate_rates_discrete_rate_shape, where the rate *shape*
// is the simplex and the unit mean is imposed afterwards by dividing by the
// weighted mean sum_j p_j x_j. Here there is no such denominator, so r_k depends
// on s_k and p_k alone rather than on every other category, and the small-
// denominator blow-up of that parameterization cannot occur.
//
// With an invariant class category 0 is held at r_0 = 0 and contributes nothing
// to the mean, so the constraint reads sum_{k>0} p_k r_k = 1 and s covers only
// the cat_count-1 variable categories.
void _calculate_rates_discrete_mean_contribution( SiteModel *sm ) {
	size_t cat_count = sm->cat_count;
	const double* cat_proportions = Parameter_values(sm->proportions);
	memcpy(sm->cat_proportions, cat_proportions, sizeof(double)*cat_count);
	memset(sm->cat_rates, 0, sizeof(double)*cat_count);

	size_t cat = (sm->invariant ? 1 : 0);
	const double* contributions = Parameter_values(Parameters_at(sm->rates, 0));
	for (size_t i = cat; i < cat_count; i++ ) {
		sm->cat_rates[i] = contributions[i - cat]/sm->cat_proportions[i];
	}
	sm->need_update = false;
}

// Ordered increments ("rate_increments"; see docs/models/sitemodel.md). The free
// parameter is the vector of gaps between consecutive categories, so the raw
// rates are its running sum and the unit mean is imposed afterwards as in
// _calculate_rates_discrete_rate_shape:
//
//     rt_1 = theta_1,  rt_k = rt_{k-1} + theta_k,  r_k = rt_k / sum_j p_j rt_j .
//
// Positive gaps make the rates strictly increasing, which breaks the K!
// label-switching symmetry of the unordered parameterizations. The Jacobian of
// the raw map is the lower-triangular matrix of ones whatever the parameter
// values, so nothing compounds -- better conditioned than the equivalent chain
// of multiplicative factors in _calculate_rates_discrete_ratios.
//
// Note that theta has K elements for K-1 degrees of freedom: the normalisation
// is scale-free, so the likelihood is exactly flat along theta -> c*theta. Pin
// the scale (fix one element, or put theta on a simplex) before reading anything
// curvature-based off a fit.
//
// An invariant class can be prepended, as for the mean-contribution simplex: the
// running sum then starts at category 1 and category 0 is left pinned at rate 0,
// so theta carries one element per *variable* category -- one fewer than the
// proportions simplex. The invariant class contributes nothing to the weighted
// mean, so the constraint reads sum_{k>0} p_k r_k = 1 and the normalisation is
// unchanged.
void _calculate_rates_discrete_increments( SiteModel *sm ) {
	int cat_count = Parameter_size(sm->proportions);
	memset(sm->cat_rates, 0, sizeof(double)*cat_count);
	const double* cat_proportions = Parameter_values(sm->proportions);
	memcpy(sm->cat_proportions, cat_proportions, sizeof(double)*cat_count);

	Parameter* rates = Parameters_at(sm->rates, 0);
	const double* cat_rates = Parameter_values(rates);
	size_t first = (sm->invariant ? 1 : 0);
	sm->cat_rates[first] = cat_rates[0];

	double norm = sm->cat_rates[first]*sm->cat_proportions[first];
	for (size_t i = first + 1; i < cat_count; i++ ) {
		sm->cat_rates[i] = cat_rates[i - first] + sm->cat_rates[i-1];
		norm += sm->cat_rates[i]*sm->cat_proportions[i];
	}

	for (size_t i = first; i < cat_count; i++ ) {
		sm->cat_rates[i] /= norm;
	}
	sm->need_update = false;
}

// Ordered ratios ("rate_ratios" plus "top_rate"). The free parameters are a vector
// theta of K-1 ratios of each category to the next one up, and a scalar for the top
// raw rate, from which the sequence is built downwards and then normalised:
//
//     rt_{K-1} = top,  rt_{k-1} = theta_{k-1} * rt_k,
//     r_k = rt_k / sum_j p_j rt_j .
//
// With theta in (0,1) the rates increase, so this is ordered like
// _calculate_rates_discrete_increments; it is the multiplicative counterpart,
// and its parameters are the scale-free ratios between adjacent categories
// rather than absolute gaps. That also makes it the descending form of the
// cumulative-product construction of _calculate_rates_discrete.
//
// The trade-off against the increments is conditioning: the raw
// derivatives are products of the other ratios, so they grow or shrink
// geometrically in K. The top rate is redundant for the same reason the
// increments carry a spare degree of freedom -- normalisation cancels it -- so
// the same caveat applies. No invariant class either: the chain runs down to
// category 0.
void _calculate_rates_discrete_ratios( SiteModel *sm ) {
	int cat_count = Parameter_size(sm->proportions);
	memset(sm->cat_rates, 0, sizeof(double)*cat_count);
	const double* cat_proportions = Parameter_values(sm->proportions);
	memcpy(sm->cat_proportions, cat_proportions, sizeof(double)*cat_count);

	Parameter* prop_parameter = Parameters_at(sm->rates, 0);
	Parameter* last_parameter = Parameters_at(sm->rates, 1);
	const double* prop_rates = Parameter_values(prop_parameter);
	sm->cat_rates[cat_count-1] = Parameter_value(last_parameter);
	
	double norm = sm->cat_rates[cat_count-1] * sm->cat_proportions[cat_count-1];
	for (size_t i = cat_count-1; i >= 1; i-- ) {
		sm->cat_rates[i-1] = prop_rates[i-1] * sm->cat_rates[i];
		norm += sm->cat_rates[i-1]*sm->cat_proportions[i-1];
	}
	
	for (size_t i = 0; i < cat_count; i++ ) {
		sm->cat_rates[i] /= norm;
	}
	sm->need_update = false;
}

#pragma region CAT

// The category rates are empirical, so they are normalised to a mean of one over
// the *alignment* rather than over a distribution: every site contributes the rate
// of the category it was assigned to, weighted by the number of sites its pattern
// stands for. Before fasttree_cat has run, every pattern sits in category 0 and
// the mean is simply that category's rate.
bool _cat_update(SiteModel* sm){
	const double* rates = Parameter_values(Parameters_at(sm->rates, 0));
	double avg = 0;
	for (int i = 0; i < sm->sp->count; i++ ) {
		avg += rates[sm->get_site_category(sm, i)] * sm->sp->weights[i];
	}
	avg /= sm->sp->nsites;
	// A zero or non-finite mean would send every category rate to inf/nan and take
	// the likelihood with it. Report a failed update instead, which the tree
	// likelihood already turns into a NAN log-likelihood.
	if (!isfinite(avg) || avg <= 0.0) {
		return false;
	}
	for (int i = 0; i < sm->cat_count; i++ ) {
		sm->cat_rates[i] = rates[i]/avg;
	}
	sm->need_update = false;
	return true;
}

bool _update_cat(SiteModel *sm){
	if ( sm->need_update ) {
		return _cat_update(sm);
	}
	return true;
}

double _get_rate_cat( SiteModel *sm, const int index ){
	if ( sm->need_update ) {
		_cat_update(sm);
	}
	return sm->cat_rates[index] * (sm->mu == NULL ? 1.0 : Parameter_value(sm->mu));
}

SiteModel * new_CATSiteModel_with_parameters( const Parameters *params,  const size_t cat_count, SitePattern* sp){
	assert(sp != NULL);
	// One rate per category, held as a single vector parameter the way every other
	// rate parameterization holds its own: _cat_update indexes it by category, and
	// a pattern may be assigned to any of them.
	assert(Parameters_count(params) == 1);
	assert(Parameter_size(Parameters_at(params, 0)) == cat_count);
	SiteModel *sm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(sm);
	// The pattern weights drive the normalisation and the assignment is one entry
	// per pattern, so the site model keeps the pattern set alive for as long as it
	// needs it: it takes its own reference and drops it in free_SiteModel.
	sm->sp = sp;
	sp->ref_count++;
	sm->site_category = ivector(sp->count);
	sm->get_site_category = _get_site_category_CAT;

	sm->distribution = -1;
	sm->invariant = false;
	sm->quadrature = -1;
	sm->rate_parameterization = RATE_PARAMETERIZATION_AUTO;

	sm->cat_count = cat_count;

	sm->cat_rates = dvector(sm->cat_count);
	// Not mixture weights: a site sits in exactly one category rather than being
	// averaged over all of them, so the single block the likelihood integrates
	// carries a weight of 1. Kept as a full vector rather than NULL because
	// get_proportions hands it straight to integrate_partials and to the loggers,
	// which index it by category.
	sm->cat_proportions = dvector(sm->cat_count);
	for (size_t i = 0; i < cat_count; i++ ) {
		sm->cat_proportions[i] = 1.0;
	}
	sm->proportions = NULL;

	sm->rates = NULL;

	if(Parameters_count(params) > 0){
		sm->rates = new_Parameters(Parameters_count(params));
		Parameters_add_parameters(sm->rates, params);
		if (Parameters_name2(params) != NULL) {
			Parameters_set_name2(sm->rates, Parameters_name2(params));
		}
	}
	for(size_t i = 0; i < Parameters_count(params); i++){
		Parameter_set_model(Parameters_at(sm->rates, i), MODEL_SITEMODEL);
	}
	sm->mu    = NULL;
	sm->class_rate = NULL;
	sm->class_rate_parameterization = FREE_CLASS_PARAMETERIZATION_RATE;

	sm->get_rate        = _get_rate_cat;
	sm->get_proportion  = _get_proportion;
	sm->get_proportions = _get_proportions;
	sm->update = _update_cat;
	// The category rates are read off the alignment, not from a differentiable
	// discretization, so there is nothing to differentiate here.
	sm->gradient = _no_gradient;
	sm->derivative = _no_derivative;

	sm->integrate   = false;
	sm->need_update = true;

	return sm;
}

#pragma endregion

SiteModel * clone_SiteModel( const SiteModel *sm ){
	return clone_SiteModel_with(sm);
}

// The pattern set and the per-pattern category assignment, shared by both clone
// entry points. Both are NULL on everything but the empirical CAT model, but they
// are read unconditionally (get_site_category is copied along with the rest of the
// vtable), so they must be set either way.
static void _clone_SiteModel_categories( const SiteModel *sm, SiteModel *newsm ){
	newsm->sp = sm->sp;
	if ( sm->sp != NULL ){
		sm->sp->ref_count++;
	}
	newsm->site_category = NULL;
	if ( sm->site_category != NULL ){
		newsm->site_category = clone_ivector(sm->site_category, sm->sp->count);
	}
}

SiteModel * clone_SiteModel_with( const SiteModel *sm ){
	SiteModel *newsm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(newsm);

	newsm->rates = NULL;

	newsm->cat_count = sm->cat_count;

	if ( sm->rates != NULL ){
		newsm->rates = clone_Parameters(sm->rates);
	}

	// A CAT site model reads the pattern weights and carries one category per
	// pattern; the clone needs both or it dereferences whatever malloc left behind.
	// The pattern set is shared rather than deep-copied -- it is immutable here and
	// reference counted.
	_clone_SiteModel_categories(sm, newsm);

	newsm->proportions = NULL;
	if ( sm->proportions != NULL ){
		newsm->proportions = clone_Parameter(sm->proportions);
	}

	newsm->mu = NULL;
	if ( sm->mu != NULL ){
		newsm->mu = clone_Parameter(sm->mu);
	}

	newsm->class_rate = NULL;
	newsm->class_rate_parameterization = sm->class_rate_parameterization;
	if ( sm->class_rate != NULL ){
		newsm->class_rate = clone_Parameter(sm->class_rate);
	}

	newsm->cat_rates = clone_dvector(sm->cat_rates, sm->cat_count);
	
	newsm->cat_proportions = NULL;
	if ( sm->cat_proportions != NULL ){
		newsm->cat_proportions = clone_dvector(sm->cat_proportions, sm->cat_count);
	}
	
	newsm->integrate = sm->integrate;
	newsm->need_update = false;

	newsm->get_rate        = sm->get_rate;
	newsm->get_proportion  = sm->get_proportion;
	newsm->get_proportions = sm->get_proportions;
    newsm->get_site_category = sm->get_site_category;
	newsm->update = sm->update;
	
	newsm->distribution = sm->distribution;
	newsm->invariant = sm->invariant;
	newsm->quadrature = sm->quadrature;
	newsm->rate_parameterization = sm->rate_parameterization;

	newsm->gradient = sm->gradient;
	newsm->derivative = sm->derivative;

	// newsm->epsilon = sm->epsilon;
	return newsm;
}

SiteModel * clone_SiteModel_with_parameters( const SiteModel *sm, Parameter* props, const Parameters* params, Parameter* mu ){
	SiteModel *newsm = (SiteModel *)malloc(sizeof(SiteModel));
	assert(newsm);
	newsm->proportions = props;
	newsm->distribution = sm->distribution;
	newsm->invariant = sm->invariant;
	newsm->quadrature = sm->quadrature;
	newsm->rate_parameterization = sm->rate_parameterization;

	newsm->rates = NULL;

	newsm->cat_count = sm->cat_count;

	if ( sm->rates != NULL ){
		newsm->rates = new_Parameters(Parameters_count(sm->rates));
		for (int i = 0; i < Parameters_count(sm->rates); i++) {
			Parameters_add(newsm->rates, Parameters_at(params, i));
		}
	}

	// See clone_SiteModel_with: CAT's pattern set and category assignment.
	_clone_SiteModel_categories(sm, newsm);

	newsm->mu = NULL;
	if ( mu != NULL ){
		newsm->mu = mu;
		mu->refCount++;
	}

	// The rate of a free class is attached by the caller (_site_model_clone),
	// which resolves it against the hashtable as it does every other parameter.
	newsm->class_rate = NULL;
	newsm->class_rate_parameterization = sm->class_rate_parameterization;

	newsm->cat_rates = clone_dvector(sm->cat_rates, sm->cat_count);
	
	newsm->cat_proportions = NULL;
	if ( sm->cat_proportions != NULL ){
		newsm->cat_proportions = clone_dvector(sm->cat_proportions, sm->cat_count);
	}
	
	newsm->integrate = sm->integrate;
	newsm->need_update = false;

	newsm->get_rate        = sm->get_rate;
	newsm->get_proportion  = sm->get_proportion;
	newsm->get_proportions = sm->get_proportions;
    newsm->get_site_category = sm->get_site_category;
	newsm->update = sm->update;
	
	newsm->gradient = sm->gradient;
	newsm->derivative = sm->derivative;

	// newsm->epsilon = sm->epsilon;
	return newsm;
}


void free_SiteModel( SiteModel *sm ){
	if ( sm->rates != NULL ) free_Parameters(sm->rates);
	if ( sm->mu != NULL ) free_Parameter(sm->mu);
	if ( sm->class_rate != NULL ) free_Parameter(sm->class_rate);
	if ( sm->cat_proportions != NULL ) free(sm->cat_proportions);
	if ( sm->site_category != NULL ) free(sm->site_category);
	// Only the empirical CAT model holds a pattern set, and it took a reference for
	// it in new_CATSiteModel_with_parameters.
	if ( sm->sp != NULL ) free_SitePattern(sm->sp);
	free_Parameter(sm->proportions);
	free(sm->cat_rates);
	free(sm);
}

// Build the internal 2-simplex X = [p, 1-p] the site model consumes, driven by
// the proportion parameter S = p (in (0,1)) through the "proportions" (S -> X)
// transform. X listens on S and is registered in the hashtable.
static Parameter* _build_invariant_simplex(Parameter* proportion,
                                           const char* simplex_id, Hashtable* hash){
	Transform* transform =
		new_SimplexTransform_with_parameter("proportions", proportion);
	double* values = dvector(transform->dim);
	transform->get(transform, values);
	Parameter* simplex =
		new_Parameter2(simplex_id, values, transform->dim, new_Constraint(0.0, 1.0));
	simplex->transform = transform;
	simplex->simplex = true;
	proportion->listeners->add_parameter(proportion->listeners, simplex);
	Hashtable_add(hash, Parameter_name(simplex), simplex);
	free(values);
	return simplex;
}

// Sugar for the weight of the extra point-mass class in category 0, shared by
// "proportion_invariant" (+I, the class is pinned at rate 0) and
// "free_class_proportion" (+F, the class has a rate of its own): <id> always
// names the scalar proportion S = p in (0,1) (so it logs as a single value and
// can carry e.g. a Beta prior). The internal 2-simplex X = [p, 1-p] the site
// model consumes is auto-named "<id>.simplex" and driven by S through the
// "proportions" (S -> X) transform. S is supplied one of two ways (the two
// feeding modes of S -> X):
//
//   - a raw box-constrained (0,1) leaf: a bare number, or {"lower":0,"upper":1,
//     "x":p}. Optimized directly (box-respecting MLE, e.g. brent).
//   - reparameterised through an unconstrained leaf: {..., "x":{...}} (logit) or
//     a "&reference" to such a parameter. Gradient/VI/HMC-ready.
//
// Returns the simplex X (what the site model stores as its proportions).
static Parameter* new_class_proportion_from_json(json_node* node, const char* key,
                                                 Hashtable* hash){
	Parameter* proportion;
	if (node->node_type == MJSON_STRING || node->node_type == MJSON_OBJECT) {
		// new_Parameter_from_json handles both a raw (0,1) leaf ({lower,upper,x:p})
		// and a logit-coupled scalar ({...,"x":{...}}), plus "&references".
		proportion = new_Parameter_from_json(node, hash);
		if (node->node_type == MJSON_OBJECT) {
			Hashtable_add(hash, Parameter_name(proportion), proportion);
		}
	}
	else {
		// Bare number: a raw (0,1) leaf named after the key that supplied it.
		double p = atof((char*)node->value);
		proportion = new_Parameter2(key, &p, 1, new_Constraint(0.0, 1.0));
		Hashtable_add(hash, Parameter_name(proportion), proportion);
	}

	if (Parameter_size(proportion) != 1) {
		fprintf(stderr, "%s must be a scalar in (0,1)\n", key);
		exit(2);
	}
	const double* pv = Parameter_values(proportion);
	if (pv[0] <= 0.0 || pv[0] >= 1.0) {
		fprintf(stderr, "%s must be in (0,1), got %f\n", key, pv[0]);
		exit(2);
	}

	char simplex_id[256];
	snprintf(simplex_id, sizeof(simplex_id), "%s.simplex", Parameter_name(proportion));
	return _build_invariant_simplex(proportion, simplex_id, hash);
}

// The JSON key naming the rate of a free class (+F) *is* its parameterization, so
// the key and the chart come as a pair and exactly one of the three may appear.
static const struct {
	const char* key;
	free_class_parameterization_t parameterization;
} free_class_rate_keys[] = {
	{"free_class_rate", FREE_CLASS_PARAMETERIZATION_RATE},
	{"free_class_contribution", FREE_CLASS_PARAMETERIZATION_CONTRIBUTION},
	{"free_class_increment", FREE_CLASS_PARAMETERIZATION_INCREMENT},
};
#define FREE_CLASS_RATE_KEY_COUNT \
	(sizeof(free_class_rate_keys) / sizeof(free_class_rate_keys[0]))

// The rate of a free class (+F), in whichever chart named it. Accepts the same
// spellings as any other scalar parameter -- a bare number, an object (raw or
// reparameterised through an unconstrained leaf), or a "&reference" -- and checks
// it against the domain of its chart: a contribution is a share of the unit mean
// and lives in (0,1), a rate is any non-negative number (0 being exactly the
// invariant class), an increment is any non-negative number (0 putting the class
// level with the fastest variable category).
//
// Only the contribution can be bounded into its domain by a box, which is the
// reason to prefer it: the rate chart's real constraint is the joint p_f r_f < 1,
// which no box on r_f alone expresses. The increment's is the equally joint
// p_f d < 1, but that one binds only at rates far beyond anything a fit reaches
// (d = 1/p_f puts the whole unit mean in the free class), so in practice the
// lower bound is the only one that does any work. Both are caught by
// _apply_free_class rather than by a constraint.
static Parameter* new_free_class_rate_from_json(
	json_node* node, const char* key, free_class_parameterization_t parameterization,
	Hashtable* hash){
	const bool is_contribution =
		(parameterization == FREE_CLASS_PARAMETERIZATION_CONTRIBUTION);
	Parameter* parameter;
	if (node->node_type == MJSON_STRING || node->node_type == MJSON_OBJECT) {
		parameter = new_Parameter_from_json(node, hash);
		if (node->node_type == MJSON_OBJECT) {
			Hashtable_add(hash, Parameter_name(parameter), parameter);
		}
	}
	else {
		double x = atof((char*)node->value);
		parameter = new_Parameter2(
			key, &x, 1, new_Constraint(0.0, is_contribution ? 1.0 : INFINITY));
		Hashtable_add(hash, Parameter_name(parameter), parameter);
	}

	if (Parameter_size(parameter) != 1) {
		fprintf(stderr, "%s must be a scalar\n", key);
		exit(2);
	}
	const double x = Parameter_values(parameter)[0];
	if (is_contribution) {
		if (x <= 0.0 || x >= 1.0) {
			fprintf(stderr, "%s is the free class's share of the unit mean and must "
			                "be in (0,1), got %f\n", key, x);
			exit(2);
		}
		check_constraint(parameter, 0.0, 1.0, 1.0e-6, 1.0 - 1.0e-6);
	}
	else {
		if (x < 0.0) {
			fprintf(stderr, "%s must be non-negative, got %f\n", key, x);
			exit(2);
		}
		check_constraint(parameter, 0.0, INFINITY, 0.0, 100.0);
	}
	return parameter;
}

// The JSON key holding the free rate parameter of a "discrete" site model *is* the
// parameterization it drives, so a key name and a parameterization name are the same
// string and exactly one of these keys may appear. That leaves nothing to infer from
// the shape of the parameter: the shape is checked against the key
// (_check_rate_parameterization) rather than used to guess what the key meant, so a
// simplex handed to "rate_increments" is an error instead of a silently different
// model. "rate_ratios" additionally needs "top_rate", the rate its chain hangs off.
static const struct {
	const char* key;
	rate_parameterization_t parameterization;
} rate_parameterization_keys[] = {
	{"rate_shape", RATE_PARAMETERIZATION_RATE_SHAPE},
	{"mean_contribution", RATE_PARAMETERIZATION_MEAN_CONTRIBUTION},
	{"rate_increments", RATE_PARAMETERIZATION_RATE_INCREMENTS},
	{"rate_ratios", RATE_PARAMETERIZATION_RATE_RATIOS},
};
#define RATE_PARAMETERIZATION_KEY_COUNT \
	(sizeof(rate_parameterization_keys) / sizeof(rate_parameterization_keys[0]))

static void _check_rate_parameterization(rate_parameterization_t parameterization,
                                         const Parameters* rates,
                                         const Parameter* proportions, size_t cat,
                                         bool invariant){
	char prefix[64];
	snprintf(prefix, sizeof(prefix), "sitemodel: \"%s\"",
	         rate_parameterization_strings[parameterization]);
	if (proportions == NULL) {
		fprintf(stderr, "%s requires a \"proportions\" simplex\n", prefix);
		exit(2);
	}
	// The mean-contribution simplex and the ordered increments start at category 1
	// and leave category 0 pinned at rate 0; the other two write their first rate
	// there, so an invariant class would be overwritten.
	if (invariant && parameterization != RATE_PARAMETERIZATION_MEAN_CONTRIBUTION &&
		parameterization != RATE_PARAMETERIZATION_RATE_INCREMENTS) {
		fprintf(stderr, "%s does not support an invariant category\n", prefix);
		exit(2);
	}
	size_t dim = cat + (invariant ? 1 : 0);
	if (Parameter_size(proportions) != dim) {
		fprintf(stderr, "%s expects a proportions simplex of dimension %zu, got %zu\n",
		        prefix, dim, Parameter_size(proportions));
		exit(2);
	}

	size_t expected = dim;
	switch (parameterization) {
		case RATE_PARAMETERIZATION_MEAN_CONTRIBUTION:
			// The invariant class contributes nothing to the mean, so the simplex
			// covers the variable categories only.
			expected = cat;
			// fall through
		case RATE_PARAMETERIZATION_RATE_SHAPE:
			if (Parameters_count(rates) != 1 || !Parameters_at(rates, 0)->simplex) {
				fprintf(stderr, "%s expects a simplex parameter\n", prefix);
				exit(2);
			}
			if (Parameter_size(Parameters_at(rates, 0)) != expected) {
				fprintf(stderr, "%s expects a simplex of dimension %zu, got %zu\n",
				        prefix, expected, Parameter_size(Parameters_at(rates, 0)));
				exit(2);
			}
			break;
		case RATE_PARAMETERIZATION_RATE_INCREMENTS:
			// As for the mean-contribution simplex, the invariant class is not part of
			// the running sum, so there is one increment per variable category.
			expected = cat;
			if (Parameters_count(rates) != 1 || Parameters_at(rates, 0)->simplex) {
				fprintf(stderr, "%s expects a plain (non-simplex) vector parameter\n",
				        prefix);
				exit(2);
			}
			if (Parameter_size(Parameters_at(rates, 0)) != expected) {
				fprintf(stderr, "%s expects %zu increments, got %zu\n", prefix, expected,
				        Parameter_size(Parameters_at(rates, 0)));
				exit(2);
			}
			break;
		case RATE_PARAMETERIZATION_RATE_RATIOS:
			if (Parameters_count(rates) != 2) {
				fprintf(stderr, "%s requires \"top_rate\", the rate of the last\n"
				                "category its chain hangs off\n", prefix);
				exit(2);
			}
			if (Parameters_at(rates, 0)->simplex) {
				fprintf(stderr, "%s expects a plain (non-simplex) vector parameter\n",
				        prefix);
				exit(2);
			}
			if (Parameter_size(Parameters_at(rates, 0)) != expected - 1) {
				fprintf(stderr, "%s expects %zu ratios, got %zu\n", prefix, expected - 1,
				        Parameter_size(Parameters_at(rates, 0)));
				exit(2);
			}
			if (Parameter_size(Parameters_at(rates, 1)) != 1) {
				fprintf(stderr, "%s expects \"top_rate\" to be a scalar, got %zu "
				                "elements\n", prefix,
				        Parameter_size(Parameters_at(rates, 1)));
				exit(2);
			}
			break;
		case RATE_PARAMETERIZATION_AUTO:
			break;
	}
}

// The key holding the parameter of a parametric distribution or of a Beta/Kumaraswamy
// quadrature cannot be marked JSON_REQUIRED in the schema, because which one is
// required depends on "distribution" and "quadrature". Enforce it here instead, so a
// missing key is a parse error naming the key and what asked for it, rather than a
// NULL dereference inside new_Parameter_from_json.
static json_node* _get_required_json_node(json_node* node, const char* key,
                                          const char* selector,
                                          const char* selector_value){
	json_node* child = get_json_node(node, key);
	if (child == NULL) {
		json_die(node, "\"%s\" is required for \"%s\": \"%s\"", key, selector,
		         selector_value);
	}
	return child;
}

// "categories" counts the *variable* categories, so the weights simplex carries one
// more element when an invariant class is prepended to them. Checked wherever a
// "proportions" is accepted: its elements are copied into cat_proportions, which is
// sized from "categories", so a longer one is a write past the end of that buffer and
// a shorter one leaves the tail categories uninitialised.
static void _check_proportions_dimension(const Parameter* proportions, size_t cat,
                                         bool invariant){
	size_t dim = cat + invariant;
	if (Parameter_size(proportions) != dim) {
		fprintf(stderr,
		        "sitemodel: \"proportions\" expects a simplex of dimension %zu "
		        "(\"categories\": %zu%s), got %zu\n",
		        dim, cat, invariant ? ", plus one for the invariant class" : "",
		        Parameter_size(proportions));
		exit(2);
	}
}

Model* new_SiteModel_from_json(json_node*node, Hashtable*hash){
	static const json_field schema[] = {
		{"a", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},        // Kumaraswamy a
		{"alpha", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},    // Beta alpha
		{"b", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},        // Kumaraswamy b
		{"beta", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},     // Beta beta
		{"categories", JSON_OPTIONAL, JSON_NUMBER},
		{"distribution", JSON_OPTIONAL, JSON_STRING},
		// The free rate class (+F): its weight, and its rate in one of the three
		// charts (see _apply_free_class).
		{"free_class_contribution", JSON_OPTIONAL, JSON_OBJECT | JSON_NUMBER | JSON_STRING},
		{"free_class_increment", JSON_OPTIONAL, JSON_OBJECT | JSON_NUMBER | JSON_STRING},
		{"free_class_proportion", JSON_OPTIONAL, JSON_OBJECT | JSON_NUMBER | JSON_STRING},
		{"free_class_rate", JSON_OPTIONAL, JSON_OBJECT | JSON_NUMBER | JSON_STRING},
		{"invariant", JSON_OPTIONAL, JSON_BOOL},
		{"mean_contribution", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},  // discrete
		{"mu", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},
		// The category rates of the empirical CAT model, selected by "sitepattern".
		{"parameters", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},
		{"proportion_invariant", JSON_OPTIONAL, JSON_OBJECT | JSON_NUMBER | JSON_STRING},
		{"proportions", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},
		{"quadrature", JSON_OPTIONAL, JSON_STRING},
		{"rates", JSON_FORBIDDEN, JSON_ANY},
		{"rate_increments", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},    // discrete
		{"rate_ratios", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},        // discrete
		{"rate_shape", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},         // discrete
		{"scale", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},    // lognormal scale
		{"shape", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},    // gamma/Weibull shape
		{"sitepattern", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING},
		{"top_rate", JSON_OPTIONAL, JSON_OBJECT | JSON_STRING}, // with rate_ratios
	};
	json_validate(node, schema, sizeof(schema) / sizeof(schema[0]));
	// The ratios build their sequence downwards from the rate of the last category,
	// so neither key means anything without the other.
	json_validate_co_required(node, "rate_ratios", "top_rate", NULL);

	char* id = get_json_node_value_string(node, "id");
	json_node* distribution_node = get_json_node(node, "distribution");
	json_node* mu_node = get_json_node(node, "mu");
	char* sp_ref = get_json_node_value_string(node, "sitepattern");
	json_node* proportions_node = get_json_node(node, "proportions");
	json_node* proportion_invariant_node = get_json_node(node, "proportion_invariant");
	json_node* discretization_node = get_json_node(node, "quadrature");

	size_t cat = 1;
	
	Parameters* rates = new_Parameters(1);
	Parameter* proportions = NULL;
	
	distribution_t distribution = DISTRIBUTION_UNIFORM;
	quadrature_t quad = QUADRATURE_QUANTILE_MEDIAN;
	rate_parameterization_t rate_parameterization = RATE_PARAMETERIZATION_AUTO;
	bool invariant = false;


	// Both name the same weights: "proportion_invariant" is the scalar p of the
	// 2-element simplex [p, 1-p] the parser builds from it, "proportions" is the
	// simplex spelled out. Rejected here rather than in the distribution branch so
	// that a model declaring no distribution does not quietly keep one and drop the
	// other.
	if (proportion_invariant_node != NULL && proportions_node != NULL) {
		fprintf(stderr, "sitemodel: specify either \"proportions\" or "
		                "\"proportion_invariant\", not both\n");
		exit(2);
	}

	// The free rate class (+F) is the invariant class with its rate unpinned: it
	// takes the same weight, spelled "free_class_proportion" rather than
	// "proportion_invariant" because nothing about it is invariant, plus a rate in
	// one of three charts. The keys stand or fall together, and the weight cannot
	// also come from "proportion_invariant" or a "proportions" simplex -- each of
	// those names the weight of a *different* category 0.
	json_node* free_class_proportion_node =
		get_json_node(node, "free_class_proportion");
	json_node* free_class_rate_node = NULL;
	free_class_parameterization_t free_class_parameterization =
		FREE_CLASS_PARAMETERIZATION_RATE;
	const char* free_class_rate_key = NULL;
	for (size_t i = 0; i < FREE_CLASS_RATE_KEY_COUNT; i++) {
		json_node* candidate = get_json_node(node, free_class_rate_keys[i].key);
		if (candidate == NULL) continue;
		free_class_rate_node = candidate;
		free_class_rate_key = free_class_rate_keys[i].key;
		free_class_parameterization = free_class_rate_keys[i].parameterization;
	}
	const bool free_class =
		(free_class_proportion_node != NULL || free_class_rate_node != NULL);
	if (free_class) {
		// Exactly one chart for the rate: r_f itself, the class's contribution to
		// the unit mean c = p_f r_f, or its increment above the fastest variable
		// rate. Any two of them would over-determine it.
		json_validate_xor(node, "free_class_rate", "free_class_contribution",
		                  "free_class_increment", NULL);
		if (free_class_proportion_node == NULL) {
			json_die(node, "\"free_class_proportion\" is required for a free rate "
			               "class: a rate with no weight is not a mixture component");
		}
		if (proportion_invariant_node != NULL || proportions_node != NULL) {
			json_die(node, "a free rate class takes its weight from "
			               "\"free_class_proportion\"; \"%s\" weights a different "
			               "category 0",
			         proportion_invariant_node != NULL ? "proportion_invariant"
			                                           : "proportions");
		}
	}

	// Which of the mutually exclusive rate keys was given, if any -- none is the pure
	// +I model. Read before anything is built so that a key meant for "discrete" can
	// be rejected on any other model rather than silently ignored.
	json_node* rates_node = NULL;
	for (size_t i = 0; i < RATE_PARAMETERIZATION_KEY_COUNT; i++) {
		json_node* candidate = get_json_node(node, rate_parameterization_keys[i].key);
		if (candidate == NULL) continue;
		if (rates_node != NULL) {
			json_die(node, "\"%s\" and \"%s\" are alternative parameterizations of the "
			               "same category rates; define at most one",
			         rate_parameterization_strings[rate_parameterization],
			         rate_parameterization_keys[i].key);
		}
		rates_node = candidate;
		rate_parameterization = rate_parameterization_keys[i].parameterization;
	}
	json_node* top_rate_node = get_json_node(node, "top_rate");
	if (rates_node != NULL && distribution_node == NULL) {
		json_die(node, "\"%s\" parameterizes the free rates of \"distribution\": "
		               "\"discrete\", which this site model does not declare",
		         rate_parameterization_strings[rate_parameterization]);
	}

	if (distribution_node != NULL) {
		char* distribution_name = get_json_node_value_string(node, "distribution");
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

		if(get_json_node(node, "categories") == NULL){
			fprintf(stderr, "sitemodel: \"categories\" must be specified for distribution %s\n", distribution_name);
			exit(2);
		}
		else{
			cat = get_json_node_value_size_t(node, "categories", 4);
		}
		
		// +G/W/L+I or +I
		if (proportion_invariant_node != NULL) {
			proportions = new_class_proportion_from_json(proportion_invariant_node,
		                                             "proportion_invariant", hash);
			invariant = true;
			Parameter_set_model(proportions, MODEL_SITEMODEL);
			if(distribution == DISTRIBUTION_DISCRETE){
				cat = 1; // +I
			}
		}
		// +G/W/L+F: the same weights as +I -- category 0 is a point mass of weight
		// p_f and the rest of the distribution shares 1-p_f -- with the rate of that
		// point mass attached below instead of being pinned at 0.
		else if (free_class) {
			// "discrete" has no distribution to put the free class beside: with no
			// rate key it is +I (a point mass and one variable rate, so a free class
			// makes it a two-point mixture) and with one it is +R, whose rates are
			// all free already.
			if (distribution == DISTRIBUTION_DISCRETE) {
				json_die(node, "a free rate class needs a rate distribution to sit "
				               "beside; \"discrete\" already has free rates (+R), or "
				               "is the invariant model (+I)");
			}
			proportions = new_class_proportion_from_json(free_class_proportion_node,
			                                             "free_class_proportion", hash);
			invariant = true;
			Parameter_set_model(proportions, MODEL_SITEMODEL);
		}
		// +R or +R+I or +G+I+D
		// We need to specify invariant=true for +R+I and +G+I+D
		else if (proportions_node != NULL) {
			invariant = get_json_node_value_bool(node, "invariant", false);

			if (proportions_node->node_type == MJSON_STRING) {
				char* ref = (char*)proportions_node->value;
				proportions = safe_get_reference_parameter(ref, hash, id);
				proportions->refCount++;
			}
			else{
				proportions = new_Parameter_from_json(proportions_node, hash);
				Hashtable_add(hash, Parameter_name(proportions), proportions);
			}
			Parameter_set_model(proportions, MODEL_SITEMODEL);

			_check_proportions_dimension(proportions, cat, invariant);
		}
		
		if(discretization_node != NULL && distribution != DISTRIBUTION_DISCRETE){
			char* method = discretization_node->value;
			if(strcasecmp("gausslaguerre", method) == 0){
				quad = QUADRATURE_GAUSS_LAGUERRE;
				if (proportions_node != NULL) {
					fprintf(stderr, "Gauss-Laguerre quadrature does not need proportions to be specified (%s)\n", proportions_node->key);
					exit(13);
				}
				// Its nodes and weights come from the quadrature rule, so there is no
				// separately weighted extra class to add, invariant or free.
				if (free_class) {
					fprintf(stderr, "Gauss-Laguerre quadrature does not support a free rate class\n");
					exit(13);
				}
				if (invariant) {
					fprintf(stderr, "Gauss-Laguerre quadrature does not support invariant sites\n");
					exit(13);
				}
			}
			else if(strcasecmp("median", method) == 0){
				quad = QUADRATURE_QUANTILE_MEDIAN;
			}
			else if(strcasecmp("mean", method) == 0){
				quad = QUADRATURE_QUANTILE_MEAN;
				// The conditional mean of a bin is a partial expectation, not a
				// quantile: it is derived per family and only the gamma and the
				// Weibull have one here. Without this the other families would
				// silently be discretized as a gamma.
				if(distribution != DISTRIBUTION_GAMMA && distribution != DISTRIBUTION_WEIBULL){
					json_die(node, "\"quadrature\": \"mean\" places each category at the "
					               "conditional mean of its bin, which is only "
					               "implemented for \"gamma\" and \"weibull\", not "
					               "\"%s\"", distribution_name);
				}
			}
			else if(strcasecmp("discrete", method) == 0){
				quad = QUADRATURE_DISCRETE;
			}
			else if(strcasecmp("beta", method) == 0){
				quad = QUADRATURE_BETA;
			}
			else if(strcasecmp("kumaraswamy", method) == 0){
				quad = QUADRATURE_KUMARASWAMY;
			}
			else{
				fprintf(stderr, "Cannot not recognize quadrature method %s\n", method);
				exit(13);
			}
		}

		if (proportions_node != NULL && distribution != DISTRIBUTION_DISCRETE &&
			quad != QUADRATURE_DISCRETE) {
			json_die(node, "\"proportions\" weights each category individually, which "
			               "needs \"quadrature\": \"discrete\"; \"%s\" takes the "
			               "weights from the discretisation itself (use "
			               "\"proportion_invariant\" for an invariant class)",
			         discretization_node == NULL ? "median"
			                                     : (char*)discretization_node->value);
		}
		
		if (distribution == DISTRIBUTION_DISCRETE) {
			if(rates_node != NULL){
				Parameters_move(rates, new_Parameter_from_json(rates_node, hash));
				// The ratios alone do not pin the sequence: it is built downwards from
				// the rate of the last category, which lives under its own key rather
				// than in a positional slot alongside them.
				if (rate_parameterization == RATE_PARAMETERIZATION_RATE_RATIOS &&
					top_rate_node != NULL) {
					Parameters_move(rates,
					                new_Parameter_from_json(top_rate_node, hash));
				}
			}
		}
		else{
			if (rates_node != NULL) {
				fprintf(stderr, "sitemodel: \"%s\" parameterizes the free rates of "
				                "\"distribution\": \"discrete\"; the categories of a\n"
				                "\"%s\" model come from its quantile function "
				                "instead\n",
				        rate_parameterization_strings[rate_parameterization],
				        distribution_name);
				exit(2);
			}
			if(distribution == DISTRIBUTION_GAMMA || distribution == DISTRIBUTION_WEIBULL){
				json_node* shape_node = _get_required_json_node(node, "shape",
				                                                "distribution",
				                                                distribution_name);
				Parameter* shape_parameter = new_Parameter_from_json(shape_node, hash);
				Parameters_move(rates, shape_parameter);
			}
			else if(distribution == DISTRIBUTION_LOGNORMAL){
				json_node* scale_node = _get_required_json_node(node, "scale",
				                                                "distribution",
				                                                distribution_name);
				Parameter* scale_parameter = new_Parameter_from_json(scale_node, hash);
				Parameters_move(rates, scale_parameter);
			}

			// quad is only set when "quadrature" is present, so the node naming the
			// method the parameters are missing from is there to be quoted back.
			const char* method =
				(discretization_node == NULL ? "" : (char*)discretization_node->value);

			if(quad == QUADRATURE_BETA){
				json_node* alpha_node = _get_required_json_node(node, "alpha",
				                                                "quadrature", method);
				Parameter* alpha_parameter = new_Parameter_from_json(alpha_node, hash);

				json_node* beta_node = _get_required_json_node(node, "beta",
				                                               "quadrature", method);
				Parameter* beta_parameter = new_Parameter_from_json(beta_node, hash);

				// GSL's iterative Beta inverse-CDF stops converging once either shape
				// gets small -- 0.02 already fails -- and a line search on alpha or
				// beta bounded only by 0 walks right into that region. Unless the user
				// pinned an explicit positive lower bound, keep both shapes above the
				// awkward zone. This is a usability bound, not a guarantee: the inverse
				// also fails for extreme *ratios* (alpha=1000 with beta=0.5), so
				// _beta_quantiles remains the ultimate safety net.
				Parameter* shapes[2] = {alpha_parameter, beta_parameter};
				for (int i = 0; i < 2; i++) {
					if (Parameter_lower(shapes[i]) <= 0.0) {
						Parameter_set_lower(shapes[i], BETA_QUADRATURE_SHAPE_LOWER);
						if (Parameter_value(shapes[i]) < BETA_QUADRATURE_SHAPE_LOWER) {
							Parameter_set_value(shapes[i], BETA_QUADRATURE_SHAPE_LOWER);
						}
					}
				}

				Parameters_move(rates, alpha_parameter);
				Parameters_move(rates, beta_parameter);
			}
			else if(quad == QUADRATURE_KUMARASWAMY){
				json_node* a_node = _get_required_json_node(node, "a", "quadrature",
				                                            method);
				Parameter* a_parameter = new_Parameter_from_json(a_node, hash);

				json_node* b_node = _get_required_json_node(node, "b", "quadrature",
				                                            method);
				Parameter* b_parameter = new_Parameter_from_json(b_node, hash);

				// As "a" -> 0 the interior quantile B(1/K)^(1/a) of the Kumaraswamy
				// inverse-CDF underflows to exactly 0, and the downstream inverse-CDF
				// (e.g. gsl_cdf_gamma_Qinv(0,...)) returns +Inf, poisoning the
				// likelihood with a NaN. Unless the user pinned an explicit positive
				// lower bound, set one a few times above that underflow cliff. The
				// cliff depends on the category count K and on b: the smallest CDF
				// base is B(1/K) = 1 - (1 - 1/K)^(1/b), and B^(1/a) underflows once
				// (1/a)*ln(B) < ln(DBL_MIN) ~ -700. The clamp in
				// _gamma_approx_quantile remains the ultimate safety net.
				if (Parameter_lower(a_parameter) <= 0.0) {
					const double K = (double)cat;
					const double b_init = Parameter_value(b_parameter);
					double base_min = 1.0 - pow(1.0 - 1.0 / K, 1.0 / b_init);
					base_min = fmax(base_min, 1e-300);
					const double a_cliff = -log(base_min) / 700.0;
					double a_lower = fmin(fmax(5.0 * a_cliff, 1e-3), 0.1);
					Parameter_set_lower(a_parameter, a_lower);
					if (Parameter_value(a_parameter) < a_lower) {
						Parameter_set_value(a_parameter, a_lower);
					}
				}

				Parameters_move(rates, a_parameter);
				Parameters_move(rates, b_parameter);
			}
		}

		// The shape of the parameter is checked against the key that named it, so a
		// mis-sized or mis-typed one is a parse error rather than a read off the end
		// of the parameter or a silently different model. The pure +I model has no
		// rate parameter to parameterise.
		if (rates_node != NULL) {
			_check_rate_parameterization(rate_parameterization, rates, proportions, cat,
			                             invariant);
		}

		for (int i = 0; i < Parameters_count(rates); i++) {
			Parameter* p = Parameters_at(rates, i);
			Parameter_set_model(p, MODEL_SITEMODEL);
			Hashtable_add(hash, Parameters_name(rates, i), p);
		}
	}
	// allow +I only without distribution=discrete
	else if (proportion_invariant_node != NULL) {
		proportions = new_class_proportion_from_json(proportion_invariant_node,
		                                             "proportion_invariant", hash);
		invariant = true;
		Parameter_set_model(proportions, MODEL_SITEMODEL);
		distribution = DISTRIBUTION_DISCRETE;
	}
	// The weights belong to the categories of a rate distribution, and "categories"
	// is only read when one is declared. Without it there is a single category of
	// weight 1 and the simplex would be dropped on the floor.
	else if (proportions_node != NULL) {
		json_die(node, "\"proportions\" gives the weight of each category of a "
		               "\"distribution\", which this site model does not declare");
	}
	// Unlike +I, a free class is not a model on its own: on its own it is a
	// two-point mixture, which is what "distribution": "discrete" already is.
	else if (free_class) {
		json_die(node, "a free rate class is an extra component beside a rate "
		               "\"distribution\", which this site model does not declare");
	}

	if(distribution == DISTRIBUTION_GAMMA || distribution == DISTRIBUTION_WEIBULL){
		Parameter* p = Parameters_at(rates, 0);
		// double lower = fmax(Parameter_lower(p), SITEMODEL_ALPHA_MIN);
		// double upper = fmin(Parameter_upper(p), SITEMODEL_ALPHA_MAX);
		Parameter_set_flower(p, SITEMODEL_ALPHA_MIN);
		Parameter_set_fupper(p, SITEMODEL_ALPHA_MAX);
	}
	
	SiteModel* sm = NULL;
	
	if (sp_ref != NULL) {
		// "sitepattern" selects the empirical CAT model, which reads its categories
		// off the alignment. Everything describing a rate *distribution* would be
		// silently dropped on the floor here, so refuse it rather than quietly
		// fitting a different model than the one that was asked for.
		if (distribution_node != NULL || proportions_node != NULL ||
		    proportion_invariant_node != NULL || discretization_node != NULL ||
		    free_class) {
			json_die(node, "\"sitepattern\" selects the empirical CAT site model, "
			               "whose categories come from the alignment; it cannot also "
			               "take a rate \"distribution\" or its weights");
		}
		int cat = get_json_node_value_int(node, "categories", 4);
		SitePattern* sp = Hashtable_get(hash, sp_ref+1);
		if (sp == NULL) {
			json_die(node, "\"sitepattern\" refers to \"%s\", which is not defined",
			         sp_ref + 1);
		}
		json_node* params_node = get_json_node(node, "parameters");
		if (params_node == NULL) {
			json_die(node, "the CAT site model reads its category rates from "
			               "\"parameters\"");
		}
		// One vector parameter of "categories" elements, like every other rate
		// parameterization. The rates are indexed by category, so a mis-sized one is
		// a read off the end of the parameter rather than a different model.
		Parameters_move(rates, new_Parameter_from_json(params_node, hash));
		Parameter* cat_rates = Parameters_at(rates, 0);
		if (Parameter_size(cat_rates) != (size_t)cat) {
			json_die(node, "\"parameters\" holds %zu rate(s) but \"categories\" is %d",
			         Parameter_size(cat_rates), cat);
		}
		Hashtable_add(hash, Parameter_name(cat_rates), cat_rates);

		sm = new_CATSiteModel_with_parameters(rates, cat, sp);
	}
	else {
		sm = new_SiteModel_with_parameters(rates, proportions, cat, distribution, invariant, quad, rate_parameterization);
	}
	
	
	if (sm->rates != NULL && Parameters_name2(sm->rates) != NULL) {
		Hashtable_add(hash, Parameters_name2(sm->rates), sm->rates);
	}
	
	if (mu_node != NULL) {
		sm->mu = new_Parameter_from_json(mu_node, hash);
		Parameter_set_model(sm->mu, MODEL_SITEMODEL);
		check_constraint(sm->mu, 0, INFINITY, 0.001, 100);
		Hashtable_add(hash, Parameter_name(sm->mu), sm->mu);
	}

	// Unpin category 0. Before new_SiteModel2, which registers the model as a
	// listener of every parameter the site model holds.
	if (free_class) {
		Parameter* class_rate = new_free_class_rate_from_json(
			free_class_rate_node, free_class_rate_key, free_class_parameterization,
			hash);
		SiteModel_set_class_rate(sm, class_rate, free_class_parameterization);
		free_Parameter(class_rate);  // the site model took its own reference
	}

	Model* msm = new_SiteModel2(id, sm);
	
	msm->print = _SiteModel_print;
	free_Parameters(rates);
	return msm;
}
