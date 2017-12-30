//
//  gamvi.c
//  viphy
//
//  Created by Mathieu Fourment on 12/04/2017.
//  Copyright © 2017 University of Technology Sydney. All rights reserved.
//

#include "gamvi.h"

#include <math.h>
#include <string.h>

#include "phyc/matrix.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

// Uses the Generalized Reparameterization Gradient

void init_meanfield_gamma(variational_t* var){
	// for reuse
	for (int i = 0; i < Parameters_count(var->parameters); i++) {
		Parameter_store(Parameters_at(var->parameters, i));
	}
}

void grad_elbo_gamma_meanfield(variational_t* var, double* grads){
	if (var->initialized == false) {
		init_meanfield_gamma(var);
		var->initialized = true;
	}

	int dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;

    memset(grads, 0, sizeof(double)*dim*2);
    
    double* z = dvector(dim);
    
    for (int i = 0; i < var->grad_samples; i++) {
		for (int j = 0; j < dim; j++) {
			double alpha = exp(Parameters_value(var->var_parameters, j));
			double beta = exp(Parameters_value(var->var_parameters, j+dim));
			
			// gsl uses shape and scale to parametrize gamma distribution
			z[j] = gsl_ran_gamma(var->rng, alpha, 1.0/beta);
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter_set_value(p, z[j]);
        }
		
		double logP = posterior->logP(posterior);
        
		for (int k = 0; k < dim; k++){
			Parameter* p = Parameters_at(var->parameters, k);
			
            double dlogP =  posterior->dlogP(posterior, p);
			
			double alpha = exp(Parameters_value(var->var_parameters, k));
			double beta = exp(Parameters_value(var->var_parameters, k+dim));
			
			double psi_alpha = gsl_sf_psi(alpha);
			double psi1_alpha = gsl_sf_psi_1(alpha);
			double psi2_alpha = gsl_sf_psi_n(2, alpha);
			
			double epsilon = (log(z[k]) - psi_alpha + log(beta))/sqrt(psi1_alpha);
			
			double h_alpha = z[k] * (epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha);
			double u_alpha = epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha + (psi2_alpha/(2.0*psi1_alpha));
			
			double dlogqdz = (alpha-1.0)/z[k] - beta;
			
			// alpha
			double g_rep = dlogP * h_alpha;
			double dlogqda = log(beta) - psi_alpha + log(z[k]);
			double g_corr = logP*(dlogqdz*h_alpha + dlogqda + u_alpha);
			grads[k] += g_rep + g_corr;
			
			// beta
			// g_corr == 0
			// g_rep
			grads[k+dim] += -dlogP*z[k]/beta;
        }
    }
    
    for (int k = 0; k < dim*2; k++) {
        grads[k] /= var->grad_samples;
    }
    

	// Gradient entropy
	for (int k = 0; k < dim; k++) {
		double alpha = exp(Parameters_value(var->var_parameters, k));
		double beta = exp(Parameters_value(var->var_parameters, k+dim));
		grads[k] += 1.0 + (1.0-alpha) * gsl_sf_psi_1(alpha);
		grads[k+dim] += -1.0/beta;
	}

    free(z);
}

double elbo_gamma_meanfield(variational_t* var){
	if (var->initialized == false) {
		init_meanfield_gamma(var);
		var->initialized = true;
	}
    size_t dim = Parameters_count(var->parameters);
    double elbo = 0;
	Model* posterior = var->posterior;
    
    for (int i = 0; i < var->elbo_samples; i++) {
		for (int j = 0; j < dim; j++) {
			double alpha = exp(Parameters_value(var->var_parameters, j));
			double beta = exp(Parameters_value(var->var_parameters, j+dim));
			
			double z = gsl_ran_gamma(var->rng, alpha, 1.0/beta);
			
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter_set_value(p, z);
        }
        
        elbo += posterior->logP(posterior);
    }
    elbo /= var->elbo_samples;
    
    // Entropy of gamma: H(X) = alpha - ln(beta) lngamma(alpha) + (1-a)psi(alpha)
	for (int j = 0; j < dim; j++) {
		double alpha = exp(Parameters_value(var->var_parameters, j));
		double lbeta = Parameters_value(var->var_parameters, j+dim);
        elbo += alpha - lbeta + gsl_sf_lngamma(alpha) + (1.0-alpha) * gsl_sf_psi(alpha);
    }
    return elbo;
}

