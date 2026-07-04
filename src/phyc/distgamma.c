// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#include "distgamma.h"

#include <math.h>
#include <strings.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>

#include "distmodel.h"
#include "matrix.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"

// NOTE: GSL uses shape and scale parameterization

// Derivative of the regularized lower incomplete gamma P(a, x) wrt the shape a.
//
// P/Q are evaluated by the same series (x < a+1) and continued fraction
// (x >= a+1) as Numerical Recipes' gser/gcf, and the derivative wrt a is
// propagated analytically through each. Both branches share
//   L = -x + a*log(x) - lgamma(a),  dL/da = log(x) - psi(a)
// with P = S*exp(L) (series sum S) or Q = h*exp(L) (continued fraction h).
double gamma_p_grad_a(double a, double x) {
    if (x <= 0.0) return 0.0;
    const double psi = gsl_sf_psi(a);
    const double expL = exp(-x + a * log(x) - lgamma(a));
    const double dLda = log(x) - psi;

    if (x < a + 1.0) {
        // Series: S = sum_n t_n, t_0 = 1/a, t_n = t_{n-1} * x/(a+n).
        // dS/da = -sum_n t_n * H_n,  H_n = sum_{k=0}^{n} 1/(a+k).
        double ap = a;
        double del = 1.0 / a;   // t_0
        double sum = del;       // S
        double H = 1.0 / a;     // H_0
        double dsum = -del * H; // dS/da
        for (int n = 0; n < 1000; n++) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            H += 1.0 / ap;
            dsum -= del * H;
            if (fabs(del) < fabs(sum) * 1e-15) break;
        }
        // P = sum*expL  =>  dP/da = expL*(dS/da + S*dL/da)
        return expL * (dsum + sum * dLda);
    }

    // Continued fraction (Lentz) for Q, with derivatives wrt a propagated
    // alongside each quantity (d* denotes d/da).
    const double FPMIN = 1e-300;
    const double EPS = 1e-15;
    double b = x + 1.0 - a, db = -1.0;
    double c = 1.0 / FPMIN, dc = 0.0;
    double d = 1.0 / b, dd = -d * d * db;  // d = 1/b
    double h = d, dh = dd;
    for (int i = 1; i <= 1000; i++) {
        double an = -i * (i - a), dan = (double)i;
        b += 2.0;  // db unchanged (-1)
        double E = an * d + b;             // next denominator (d is prev 1/D)
        double dE = dan * d + an * dd + db;
        if (fabs(E) < FPMIN) { E = FPMIN; dE = 0.0; }
        double cprev = c, dcprev = dc;
        c = b + an / cprev;
        dc = db + (dan * cprev - an * dcprev) / (cprev * cprev);
        if (fabs(c) < FPMIN) { c = FPMIN; dc = 0.0; }
        d = 1.0 / E;
        dd = -d * d * dE;
        double del = d * c;
        double ddel = dd * c + d * dc;
        double hprev = h, dhprev = dh;
        h = hprev * del;
        dh = dhprev * del + hprev * ddel;
        // Both the value (del -> 1) and its derivative (ddel -> 0) must converge;
        // breaking on the value alone truncates the slower derivative recurrence.
        if (fabs(del - 1.0) < EPS && fabs(ddel) < EPS) break;
    }
    // Q = h*expL  =>  dQ/da = expL*(dh + h*dL/da);  dP/da = -dQ/da
    return -expL * (dh + h * dLda);
}

double DistributionModel_log_gamma(DistributionModel* dm){
    if(!dm->need_update) return dm->lp;
    dm->lp = 0.0;
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* beta = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double betaValue = beta[0];
		if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
			betaValue = 1.0/betaValue;
		}

        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->lp += log(gsl_ran_gamma_pdf(values[i] - dm->shift, *alpha, betaValue));
            }
        }
    }
    // multiple parameter and len(x) >= 1 and \sum_i len(x_i) == len(parameter)
    else{
        size_t index = 0;
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                double betaValue = beta[index];
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    betaValue = 1.0/betaValue;
                }
                dm->lp += log(gsl_ran_gamma_pdf(values[i] - dm->shift, alpha[index], betaValue));
                index++;
            }
        }
    }
    dm->need_update = false;
    return dm->lp;
}

void DistributionModel_gamma_gradient(DistributionModel* dm, Parameters* parameters){
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* beta = Parameters_at(dm->parameters, 1);
    
    const double* alphaValues = Parameter_values(alpha);
    const double* betaValues = Parameter_values(beta);
    size_t mask = -(Parameter_size(alpha) != 1);

    Parameter* alphax = Parameters_depends(parameters, alpha);
    Parameter* betax = Parameters_depends(parameters, beta);
    if(alphax != NULL){
        // A shared scalar alpha (mask == 0) folds every element onto tempp[0] so
        // backward sees the full sum, not just the last element.
        memset(dm->tempp, 0, Parameter_size(alpha) * sizeof(double));
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                size_t idx = index & mask;
                double betaValue = betaValues[idx];
                double xs = xValues[j] - dm->shift;
                double g;
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    g = log(xs) + log(betaValue) - gsl_sf_psi(alphaValues[idx]);
                }
                else{
                    g = log(xs) - log(betaValue) - gsl_sf_psi(alphaValues[idx]);
                }
                dm->tempp[idx] += g;
                alpha->grad[idx] += g;
                index++;
            }
        }
        if(alpha != alphax){
            alpha->transform->backward(alpha->transform, dm->tempp);
        }
    }

    if(betax != NULL){
        memset(dm->tempp, 0, Parameter_size(beta) * sizeof(double));
        size_t index = 0;
        for(size_t k = 0; k < Parameters_count(dm->x); k++){
            Parameter* x = Parameters_at(dm->x, k);
            size_t dim = Parameter_size(x);
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < dim; j++){
                size_t idx = index & mask;
                double alphaValue = alphaValues[idx];
                double betaValue = betaValues[idx];
                double xs = xValues[j] - dm->shift;
                double g;
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    g = alphaValue/betaValue - xs;
                }
                else{
                    g = (xs - alphaValue*betaValue)/(betaValue*betaValue);
                }
                dm->tempp[idx] += g;
                beta->grad[idx] += g;
                index++;
            }
        }

        if(beta != betax){
            beta->transform->backward(beta->transform, dm->tempp);
        }
    }

    size_t index = 0;
    for(size_t k = 0; k < Parameters_count(dm->x); k++){
        Parameter* x = Parameters_at(dm->x, k);
        Parameter* xx = Parameters_depends(parameters, x);
        size_t sizeX = Parameter_size(x);
        if (xx != NULL) {
            const double* xValues = Parameter_values(x);
            for(size_t j = 0; j < sizeX; j++){
                size_t idx = index & mask;
                double alphaValue = alphaValues[idx];
                double betaValue = betaValues[idx];
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_SCALE) {
                    betaValue = 1.0/betaValue;
                }
                // x_k owns its own grad/transform: index into them locally (j)
                // while index & mask tracks the global position for the params.
                double xs = xValues[j] - dm->shift;
                dm->tempp[j] = (alphaValue-1.0)/xs - betaValue;
                x->grad[j] += dm->tempp[j];
                index++;
            }
            if(x != xx){
                x->transform->backward(x->transform, dm->tempp);
            }
        }
        else{
            index += sizeX;
        }
    }
}

static void DistributionModel_gamma_sample(DistributionModel* dm){
    const double* alpha = Parameter_values(Parameters_at(dm->parameters, 0));
    const double* beta = Parameter_values(Parameters_at(dm->parameters, 1));
    size_t dimX = Parameters_count(dm->x);

    // single parameter and len(x_i) >= 1 (e.g. prior) 
    if(Parameter_size(Parameters_at(dm->parameters, 0)) == 1){
        double betaValue = beta[0];
        if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
            betaValue = 1.0/betaValue;
        }
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                dm->tempx[i] = gsl_ran_gamma(dm->rng, *alpha, betaValue);
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
    // multiple parameter and len(x) >= 1 and \sum_i len(x_i) == len(parameter)
    else{
        size_t index = 0;
        for(size_t j = 0; j < dimX; j++){
            Parameter* x = Parameters_at(dm->x, j);
            size_t dim = Parameter_size(x);
            const double* values = Parameter_values(x);
            for (size_t i = 0; i < dim; i++) {
                double betaValue = beta[index];
                if (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE) {
                    betaValue = 1.0/betaValue;
                }
                dm->tempx[i] = gsl_ran_gamma(dm->rng, alpha[index], betaValue);
                index++;
            }
            Parameter_set_values(x, dm->tempx);
        }
    }
}

// Reparameterized sampling for variational use. The Gamma has no explicit
// reparameterization, so we draw an ordinary sample (rejection sampling); the
// implicit reparameterization gradient (see below) differentiates the sample x
// itself, so no auxiliary noise needs to be stored.
static void DistributionModel_gamma_rsample(DistributionModel* dm){
    DistributionModel_gamma_sample(dm);
}

// Implicit reparameterization gradient (Figurnov, Mohamed & Mnih, NeurIPS 2018)
// of the log objective wrt the Gamma parameters, given the downstream gradient
// x->grad[i] = dL/dx_i. Writing x = shift + theta*z with z ~ Gamma(a, 1):
//   scale theta:  dx/dtheta = z                         (explicit, z fixed)
//   shape a:      dx/da = theta * dz/da,
//                 dz/da = -dP(a, z)/da / pdf_std(a, z)   (implicit, theta fixed)
// where pdf_std(a, z) = z^{a-1} e^{-z} / Gamma(a). For the rate parameterization
// (beta = 1/theta), dx/dbeta = -z/beta^2.
static void DistributionModel_gamma_rgradient(DistributionModel* dm){
    Parameter* alpha = Parameters_at(dm->parameters, 0);
    Parameter* beta = Parameters_at(dm->parameters, 1);
    const double* alphaValues = Parameter_values(alpha);
    const double* betaValues = Parameter_values(beta);
    bool rate = (dm->parameterization == DISTRIBUTION_GAMMA_SHAPE_RATE);
    size_t dimX = Parameters_count(dm->x);
    size_t index = 0;
    for (size_t j = 0; j < dimX; j++) {
        Parameter* x = Parameters_at(dm->x, j);
        size_t dim = Parameter_size(x);
        const double* xValues = Parameter_values(x);
        size_t base = index;
        // shape (alpha): implicit reparameterization term.
        for (size_t i = 0; i < dim; i++) {
            double a = alphaValues[index];
            double theta = rate ? 1.0 / betaValues[index] : betaValues[index];
            double z = (xValues[i] - dm->shift) / theta;
            double dzda = -gamma_p_grad_a(a, z) / gsl_ran_gamma_pdf(z, a, 1.0);
            dm->tempp[i] = x->grad[i] * theta * dzda;
            alpha->grad[index] += dm->tempp[i];
            index++;
        }
        if (alpha->transform != NULL) {
            alpha->transform->backward(alpha->transform, dm->tempp);
        }
        // scale / rate (beta): explicit term, z held fixed.
        index = base;
        for (size_t i = 0; i < dim; i++) {
            double theta = rate ? 1.0 / betaValues[index] : betaValues[index];
            double z = (xValues[i] - dm->shift) / theta;
            double dxdbeta = rate ? -z / (betaValues[index] * betaValues[index]) : z;
            dm->tempp[i] = x->grad[i] * dxdbeta;
            beta->grad[index] += dm->tempp[i];
            index++;
        }
        if (beta->transform != NULL) {
            beta->transform->backward(beta->transform, dm->tempp);
        }
    }
}

DistributionModel* new_GammaDistributionModel_with_parameters(Parameters* parameters, Parameters* x, distribution_parameterization parameterization){
	DistributionModel* dm = new_DistributionModel(parameters, x);
	dm->type = DISTRIBUTION_GAMMA;
	dm->logP = DistributionModel_log_gamma;
    dm->gradient = DistributionModel_gamma_gradient;
	dm->sample = DistributionModel_gamma_sample;
	dm->rsample = DistributionModel_gamma_rsample;
	dm->rgradient = DistributionModel_gamma_rgradient;
	dm->parameterization = parameterization;
    dm->shift = 0;
    dm->support[0] = 0;
    dm->support[1] = INFINITY;
	return dm;
}

Model* new_GammaDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
    size_t paramCount = Parameter_size(Parameters_at(x, 0));
    
    char* file = get_json_node_value_string(node, "file");
    Parameters* parameters = new_Parameters(2);
    Parameter* alpha = NULL;
    Parameter* beta = NULL;
    
    distribution_parameterization parameterization = DISTRIBUTION_GAMMA_SHAPE_RATE;
    bool scale = false;
    
    char* parameterization_string = get_json_node_value_string(node, "parameterization");
    if(parameterization_string!= NULL && strcasecmp(parameterization_string, "scale") == 0){
        scale = true;
    }

    // empirical
    if (file != NULL) {
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, x);
        double* alphaValues = malloc(sizeof(double)*paramCount);
        double* betaValues = malloc(sizeof(double)*paramCount);
        
        for (int i = 0; i < paramCount; i++) {
            const double* vec = Vector_data(samples[i]);
            double m = mean(vec, Vector_length(samples[i]));
            double v = variance(vec, Vector_length(samples[i]), m);
            alphaValues[i] = m*m/v;
            if (scale) {
                betaValues[i] = v/m;
            }
            else{
				betaValues[i] = m/v;
            }
            free_Vector(samples[i]);
        }
		
		alpha = new_Parameter2("alpha", alphaValues, paramCount, new_Constraint(0, INFINITY));
        beta = new_Parameter2("beta", betaValues, paramCount, new_Constraint(0, INFINITY));

        free(samples);
    }
    else if(get_json_node(node, "parameters") == NULL){
        double* alphaValues = malloc(sizeof(double)*paramCount);
        double* betaValues = malloc(sizeof(double)*paramCount);
        for (size_t i = 0; i < paramCount; i++) {
            alphaValues[i] = 0;
            betaValues[i] = 1;
        }
        alpha = new_Parameter2("alpha", alphaValues, paramCount, new_Constraint(0, INFINITY));
        beta = new_Parameter2("beta", betaValues, paramCount, new_Constraint(0, INFINITY));
        
        free(alphaValues);
        free(betaValues);
    }
    else{
        json_node* parametersNode = get_json_node(node, "parameters");
        json_node* alphaNode = get_json_node(parametersNode, "shape");
        json_node* betaNode = get_json_node(parametersNode, "rate");
        if(betaNode == NULL){
            betaNode = get_json_node(parametersNode, "scale");
            parameterization = DISTRIBUTION_GAMMA_SHAPE_SCALE;
        }

        if(betaNode == NULL){
            fprintf(stderr, "Gamma distribution should be parametrized with shape and (scale or rate)\n");
            exit(13);
        }
        
        alpha = distmodel_parse_parameter(alphaNode, hash, "", 0, INFINITY);
        beta = distmodel_parse_parameter(betaNode, hash, "", 0, INFINITY);
    }

    Parameters_move(parameters, alpha);
    Parameters_move(parameters, beta);

    DistributionModel* dm = new_GammaDistributionModel_with_parameters(parameters, x, parameterization);
    
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel2(id, dm);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");

    free_Parameters(x);
    free_Parameters(parameters);
    
    return model;
}
