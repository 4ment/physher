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

#include "matrix.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

// Uses the Generalized Reparameterization Gradient

void init_meanfield_gamma(variational_t* var){
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
	
	// for reuse
	for (int i = 0; i < dim; i++) {
		Parameter_store(Parameters_at(var->parameters, i));
	}

//	Parameters_set_all_value(var->var_parameters, 2);
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		
//		if(Parameter_value(p) < 1.0e-6){
//			Parameter* var_p_alpha = Parameters_at(var->var_parameters, i);
//			Parameter* var_p_beta = Parameters_at(var->var_parameters, i+dim);
//			Parameter_set_estimate(var_p_alpha, false);
//			Parameter_set_value(var_p_alpha, 0); // alpha == 1
//			Parameter_set_value(var_p_beta, log(20));
//			// This is an exponential distribution with mean 1/20
//			return;
//		}

		Parameter* var_p_mu = Parameters_at(var->var_parameters, i);
		double dlogP;
		double d2logP = Model_second_derivative(posterior, Parameters_at(var->parameters, i), &dlogP, 0.001);
		double mu = Parameters_value(var->parameters, i); // mean = mode of normal
		double v = -1.0/d2logP; // variance of normal
//		printf("%f %s\n", mu, Parameter_name(var_p_mu));
		double q_beta = log(0.5*(mu + sqrt(mu*mu + 4.0*v))/(2.0*v));
		double q_alpha = log((mu * q_beta + 1.0)*2);
		
		if (isnan(q_alpha)) {
			q_alpha = 0;
		}
		if (isnan(q_beta)) {
			q_beta = log(1000);
		}
		
		Parameters_set_value(var->var_parameters, i, q_alpha);
		Parameters_set_value(var->var_parameters, i+dim, q_beta);
//					printf("%f %f  %f %f %s\n", mu, d2logP, q_alpha, q_beta, p->name);
		//			Parameter_set_estimate(var_p_mu, false);
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
//			printf("%d a: %f b: %f z: %f\n", j, alpha, beta, z[j]);
        }
		
		double logP = posterior->logP(posterior);
        
		for (int k = 0; k < dim; k++){
			Parameter* p = Parameters_at(var->parameters, k);

			double dlogP = Model_first_derivative(posterior, p, 0.001);
//            double dlogP =  posterior->dlogP(posterior, p);
			
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
    
    // Entropy of gamma: H(X) = alpha - ln(beta) + lngamma(alpha) + (1-alpha)psi(alpha)
	for (int j = 0; j < dim; j++) {
		double alpha = exp(Parameters_value(var->var_parameters, j));
		double lbeta = Parameters_value(var->var_parameters, j+dim);
        elbo += alpha - lbeta + gsl_sf_lngamma(alpha) + (1.0-alpha) * gsl_sf_psi(alpha);
    }
    return elbo;
}

