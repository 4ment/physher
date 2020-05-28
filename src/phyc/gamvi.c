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

		// Laplus initialization
		double d2logP = posterior->d2logP(posterior, Parameters_at(var->parameters, i));
		double map = Parameters_value(var->parameters, i);
        
        double rate = map * -d2logP;
        double shape = rate*map + 1;
		
        double q_alpha = log(shape);
        double q_beta = log(rate);
        
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

//MARK: block gamma meanfield

void klqp_block_meanfield_gamma_sample1(variational_block_t* var, double* jacobian){
    size_t dim = Parameters_count(var->parameters);
    double* z = Vector_data(var->etas);
    
    for (int j = 0; j < dim; j++) {
        Parameter* p = Parameters_at(var->parameters, j);
        double shape = Parameters_value(var->var_parameters[0], j);
        double rate = Parameters_value(var->var_parameters[1], j);
        
		// gsl uses shape and scale to parametrize gamma distribution
		z[j] = gsl_ran_gamma(var->rng, shape, 1.0/rate);
		Parameter_set_value(p, z[j]);
    }
    
    if (jacobian != NULL) {
        for(int s = 0; s < var->simplex_count; s++){
            Simplex* simplex = var->simplices[s]->obj;
            const double* constrained_values = simplex->get_values(simplex);
            double stick = 1;
            for(int k = 0; k < simplex->K-1; k++){
                *jacobian += log(stick);
                stick -= constrained_values[k];
            }
        }
    }
}

// Entropy of gamma: H(X) = alpha - ln(beta) + lngamma(alpha) + (1-alpha)psi(alpha)
double klqp_block_meanfield_gamma_entropy(variational_block_t* var){
    double entropy = 0;
    size_t dim = Parameters_count(var->parameters);
    for (size_t i = 0; i < dim; i++) {
        double alpha = Parameters_value(var->var_parameters[0], i);
		double beta = Parameters_value(var->var_parameters[1], i);
        entropy += alpha - log(beta) + gsl_sf_lngamma(alpha) + (1.0-alpha) * gsl_sf_psi(alpha);
    }
    return entropy;
}

// parameters are the variational parameters
void klqp_block_meanfield_gamma_sample2(variational_block_t* var, const Parameters* parameters){
    size_t dim = Parameters_count(var->parameters);
    size_t opt_param_dim = Parameters_count(parameters);
    
    // find first shape
    size_t shape_idx = 0;
    for ( ; shape_idx < opt_param_dim; shape_idx++) {
        if(Parameters_at(parameters, shape_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (shape_idx == opt_param_dim) return;
    
    double* z = Vector_data(var->etas);
    
    for (int idx = 0; idx < dim; idx++) {
        Parameter* p = Parameters_at(var->parameters, idx);
        double shape = Parameters_value(var->var_parameters[0], idx);
        double rate = Parameters_value(var->var_parameters[1], idx);
        
        z[idx] = gsl_ran_gamma(var->rng, shape, 1.0/rate);
        Parameter_set_value(p, z[idx]);
    }
}

void klqp_block_meanfield_gamma_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads){
    
    // find first shape and rate
    size_t shape_idx = 0;
    size_t rate_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
    for ( ; shape_idx < opt_param_dim; shape_idx++) {
        if(Parameters_at(parameters, shape_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (shape_idx == opt_param_dim) return;
    
    for ( ; rate_idx < opt_param_dim; rate_idx++) {
        if(Parameters_at(parameters, rate_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
    
    Model* posterior = var->posterior;
    size_t dim = Parameters_count(var->parameters);
    size_t simplex_parameter_count = var->simplex_parameter_count;
    double* z = Vector_data(var->etas);
//    int idx = 0;
//    if (simplex_parameter_count > 0) {
//        for(int s = 0; s < var->simplex_count; s++){
//            Simplex* simplex = var->simplices[s]->obj;
//            for(int k = 0; k < simplex->K-1; k++){
//                Parameter* p = Parameters_at(var->parameters, idx);
//                double mu = Parameters_value(var->var_parameters[0], idx);
//                double sigma = Parameters_value(var->var_parameters[1], idx);
//                double dlogP = posterior->dlogP(posterior, p);
//                double zeta = etas[idx] * sigma + mu;
//                double gldits = 1.0/zeta + 1.0/(zeta - 1.0); // grad log det transform of stick
//                const double gldit = grad_log_det_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p));
//                double grad_mu = dlogP * grad_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p)) + gldit + gldits;
//                grads[shape_idx + idx] += grad_mu;
//                grads[rate_idx + idx] += grad_mu * etas[idx] * sigma;
//                idx++;
//            }
//        }
//    }
    
    double logP = posterior->logP(posterior);
			
	for (int k = 0; k < dim; k++){
		Parameter* p = Parameters_at(var->parameters, k);
		double dlogP =  posterior->dlogP(posterior, p);
		
		double alpha = Parameters_value(var->var_parameters[0], k);
		double beta = Parameters_value(var->var_parameters[1], k);
		double lnbeta = log(beta);
		
		double psi_alpha = gsl_sf_psi(alpha);
		double psi1_alpha = gsl_sf_psi_1(alpha);
		double psi2_alpha = gsl_sf_psi_n(2, alpha);
		
		double epsilon = (log(z[k]) - psi_alpha + lnbeta)/sqrt(psi1_alpha);
		
		double h_alpha = z[k] * (epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha);
		double u_alpha = (epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha) + (psi2_alpha/(2.0*psi1_alpha));
		
		double h_beta = -z[k]/beta;
		double dlogqdz = (alpha-1.0)/z[k] - beta;
		
		// alpha
		double g_rep = dlogP * h_alpha;
		double dlogqda = lnbeta - psi_alpha + log(z[k]);
		double g_corr = logP*(dlogqdz*h_alpha + dlogqda + u_alpha);
		grads[k+shape_idx] += g_rep + g_corr;
		
		// beta
		// g_corr == 0
		// g_rep
		grads[k+rate_idx] += dlogP * h_beta;
		
		if(!var->use_entropy){
			double logQ = log(gsl_ran_gamma_pdf(z[k], alpha, 1.0/beta));
			double g_rep_q_a = dlogqda * h_alpha;
			double g_corr_q_a = logQ*(dlogqdz*h_alpha + dlogqda + u_alpha);
			
			double u_beta = -1.0/beta;
			double dlogqdb = alpha/beta - z[k];
			double g_rep_q_b = dlogqdb * h_beta;
			
			grads[k+shape_idx] -= g_rep_q_a + g_corr_q_a;
			grads[k+rate_idx]  -= g_rep_q_b;
		}
	}
}

// Entropy: H(X) = alpha - ln(beta) + lngamma(alpha) + (1-alpha)psi(alpha)
void klqp_block_meanfield_gamma_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads){
    // find first shape
    size_t rate_idx = 0;
	size_t shape_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
	for ( ; shape_idx < opt_param_dim; shape_idx++) {
        if(Parameters_at(parameters, shape_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (shape_idx == opt_param_dim) return;
	
	for ( ; rate_idx < opt_param_dim; rate_idx++) {
        if(Parameters_at(parameters, rate_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
	
    
    size_t dim = Parameters_count(var->parameters);
    for (int j = 0; j < dim; j++) {
        double shape = Parameters_value(var->var_parameters[0], j);
		double rate = Parameters_value(var->var_parameters[1], j);
		grads[j+shape_idx] += 1.0 + (1.0-shape) * gsl_sf_psi_1(shape);
		grads[j+rate_idx] += -1.0/rate;
    }
}

double klqp_block_meanfield_gamma_logP(variational_block_t* var, double* values){
    size_t dim = Parameters_count(var->parameters);
    double logP = 0;
    for (size_t i = 0; i < dim; i++) {
        double shape = Parameters_value(var->var_parameters[0], i);
        double rate = Parameters_value(var->var_parameters[1], i);
        logP += log(gsl_ran_gamma_pdf(values[i], shape, 1.0/rate));
    }
    return logP;
}


double klqp_block_meanfield_gamma_logQ(variational_block_t* var, double* values){
	return klqp_block_meanfield_gamma_logP(var, values);
}


void klqp_block_meanfield_gamma_sample(variational_block_t* var, double* values){
    size_t dim = Parameters_count(var->parameters);
    for (size_t i = 0; i < dim; i++) {
        double shape = Parameters_value(var->var_parameters[0], i);
        double rate = Parameters_value(var->var_parameters[1], i);
        values[i] = gsl_ran_gamma(var->rng, shape, 1.0/rate);
    }
}

bool klqp_block_meanfield_gamma_sample_some(variational_block_t* var, const Parameters* parameters, double* values){
    size_t dim = Parameters_count(var->parameters);
    size_t dim2 = Parameters_count(parameters);
    for (size_t i = 0; i < dim2; i++) {
        Parameter* p = Parameters_at(parameters, i);
        for (int j = 0; j < dim; j++) {
            Parameter* p2 = Parameters_at(var->parameters, j);
            if (p == p2) {
				double shape = Parameters_value(var->var_parameters[0], i);
				double rate = Parameters_value(var->var_parameters[1], i);
				values[i] = gsl_ran_gamma(var->rng, shape, 1.0/rate);
                break;
            }
        }
    }
    return true;
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

//			double dlogP = Model_first_derivative(posterior, p, 0.001);
            double dlogP =  posterior->dlogP(posterior, p);
			
			double alpha = exp(Parameters_value(var->var_parameters, k));
            double lbeta = Parameters_value(var->var_parameters, k+dim);
            double beta = exp(lbeta);
			
			double psi_alpha = gsl_sf_psi(alpha);
			double psi1_alpha = gsl_sf_psi_1(alpha);
			double psi2_alpha = gsl_sf_psi_n(2, alpha);
			
			double epsilon = (log(z[k]) - psi_alpha + lbeta)/sqrt(psi1_alpha);
			
			double h_alpha = z[k] * (epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha);
			double u_alpha = epsilon*psi2_alpha/(2.0*sqrt(psi1_alpha)) + psi1_alpha + (psi2_alpha/(2.0*psi1_alpha));
			
			double dlogqdz = (alpha-1.0)/z[k] - beta;
			
			// alpha
			double g_rep = dlogP * h_alpha;
			double dlogqda = lbeta - psi_alpha + log(z[k]);
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

// check success
bool variational_sample_gamma_meanfield(variational_t* var, double* values){
	if(!var->ready_to_sample){
		var->finalize(var);
		var->ready_to_sample = true;
	}
	int dim = Parameters_count(var->parameters);
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		if (!p->estimate) {
			values[i] = 0;
			continue;
		}
		double alpha = exp(Parameters_value(var->var_parameters, i));
		double beta = exp(Parameters_value(var->var_parameters, i+dim));
		values[i] = gsl_ran_gamma(var->rng, alpha, 1.0/beta);
	}
	return true;
}

// assumes that there is transformation
double variational_gamma_meanfield_parameters_logP(variational_t* var, const Parameters* parameters){
	size_t dim = Parameters_count(var->parameters);
	size_t dim2 = Parameters_count(parameters);
	double logP = 0;
	for (int i = 0; i < dim2; i++) {
		Parameter* p = Parameters_at(parameters, i);
		for (int j = 0; j < dim; j++) {
			Parameter* p2 = Parameters_at(var->parameters, j);
			if (p == p2) {
				double alpha = exp(Parameters_value(var->var_parameters, i));
				double beta = exp(Parameters_value(var->var_parameters, i+dim));
				double zeta = Parameter_value(p);
				logP += log(gsl_ran_gamma_pdf(zeta, alpha, 1.0/beta));
				break;
			}
		}
	}
	return logP;
}

// assumes that there is transformation
double variational_gamma_meanfield_logP(variational_t* var, double* values){
	int dim = Parameters_count(var->parameters);
	double logP = 0;
	for (int i = 0; i < dim; i++) {
		double alpha = exp(Parameters_value(var->var_parameters, i));
		double beta = exp(Parameters_value(var->var_parameters, i+dim));
		logP += log(gsl_ran_gamma_pdf(values[i], alpha, 1.0/beta));
	}
	return logP;
}

// assumes that there is transformation
bool variational_sample_some_gamma_meanfield(variational_t* var, const Parameters* parameters, double* values){
	if(!var->ready_to_sample){
		var->finalize(var);
		var->ready_to_sample = true;
	}
	size_t dim = Parameters_count(var->parameters);
	size_t dim2 = Parameters_count(parameters);
	for (int i = 0; i < dim2; i++) {
		Parameter* p = Parameters_at(parameters, i);
		for (int j = 0; j < dim; j++) {
			Parameter* p2 = Parameters_at(var->parameters, j);
			if (p == p2) {
				if (!p->estimate) {
					values[i] = 0;
				}
				else{
					double alpha = exp(Parameters_value(var->var_parameters, i));
					double beta = exp(Parameters_value(var->var_parameters, i+dim));
					values[i] = gsl_ran_gamma(var->rng, alpha, 1.0/beta);
				}
				break;
			}
		}
	}
	return true;
}

void meanfield_gamma_log_samples(variational_t* var, FILE* file){
    size_t dim = Parameters_count(var->parameters);
    Model* posterior = var->posterior;
    double* samples = dvector(dim);
    fprintf(file, "sample,logP,logQ");
    
    int offset = 0; // don't log the simplex unconstrained parameters
    if (var->simplex_count > 0) {
        for(int s = 0; s < var->simplex_count; s++){
            Model* msimplex = var->simplices[s];
            Simplex* simplex = msimplex->obj;
            for(int k = 0; k < simplex->K; k++){
                fprintf(file, ",%s.%d", msimplex->name, k);
            }
            offset += simplex->K-1;
        }
    }
    
    for (int j = offset; j < dim; j++) {
        fprintf(file, ",%s", Parameters_name(var->parameters, j));
    }
    fprintf(file, "\n");
    
    for (int i = 0; i < var->log_samples; i++) {
        
        double logQ = 0;
        double jacobian = 0.0;
        for (int j = 0; j < dim; j++) {
            Parameter* p = Parameters_at(var->parameters, j);
            double mu = exp(Parameters_value(var->var_parameters, j));
            double sigma = exp(Parameters_value(var->var_parameters, j+dim));
            
            
            samples[j] = gsl_ran_gamma(var->rng, mu, 1.0/sigma);
            logQ += log(gsl_ran_gamma_pdf(samples[j], mu, 1.0/sigma));
            Parameter_set_value(p, samples[j]);
        }
        
        double logP = posterior->logP(posterior) + jacobian;
        fprintf(file, "%d,%f,%f", i, logP, logQ);
        
        if (var->simplex_count > 0) {
            for(int s = 0; s < var->simplex_count; s++){
                Simplex* simplex = var->simplices[s]->obj;
                const double* constrained_values = simplex->get_values(simplex);
                for(int k = 0; k < simplex->K; k++){
                    fprintf(file, ",%e",constrained_values[k]);
                }
            }
        }
        
        for (int j = offset; j < dim; j++) {
            fprintf(file, ",%e", samples[j]);
        }
        fprintf(file, "\n");
    }
    free(samples);
}
