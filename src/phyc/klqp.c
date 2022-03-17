//
//  klqp.c
//  physher
//
//  Created by Mathieu Fourment on 28/3/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "klqp.h"

#include <string.h>
#include <strings.h>
#include <math.h>

#include "matrix.h"
#include "gaussian.h"
#include "transforms.h"
#include "solve.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#define MY_PI acos(-1.0)
#define LOG_TWO_PI (log(2.0)+log(MY_PI))

//MARK: Meanfield

//MARK: normal meanfield

void klqp_meanfield_normal_init(variational_t* var){
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
	
	// for reuse
	for (int i = 0; i < Parameters_count(var->parameters); i++) {
		Parameter_store(Parameters_at(var->parameters, i));
	}
	
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		// lower(p) == 0 and upper(p) == inf
		if (!isinf(Parameter_lower(p))) {
			// fix mu when value is small
			// mu will be set such that mu-sigma*sigma = -10
			if(Parameter_value(p) < 1.0e-6){
//				Parameters_set_value(var->var_parameters, i, -5);
//				Parameters_set_value(var->var_parameters, i+dim, log(0.5));
				Parameter* var_p_mu = Parameters_at(var->var_parameters, i);
				Parameter_set_estimate(var_p_mu, false);
			}
			else{
				// Laplus: matching mode and second derivative of the loglikelihood and lognormal distributions
				//	sigma = sqrt(-1/(f''(m)*m^2))
				//	mu    = log(m)+sigma^2
				double d2logP = posterior->d2logP(posterior, Parameters_at(var->parameters, i));
				double map = Parameters_value(var->parameters, i);
				double q_sigma = sqrt(-1.0/(d2logP*map*map));
				double q_mu = log(map) + q_sigma*q_sigma;

				// d2logp could be non negative are create a nan
				if (isnan(q_mu)) {
					q_mu = 0;
				}
				if (isnan(q_sigma)) {
					q_sigma = 0.1;
				}
				// check mean is not infinity
				// For short branch b=0.000032, mu=1614.813030 sigma=40.313458 => infinity
				double mean = exp(q_mu + q_sigma*q_sigma/2);
				if(isinf(mean) || mean > 10.0*map){
					q_mu = -6;
					q_sigma = 0.5;
				}
				Parameters_set_value(var->var_parameters, i, q_mu);
				Parameters_set_value(var->var_parameters, i+dim, log(q_sigma));
			}
		}
	}
}

void klqp_meanfield_normal_finalize(variational_t* var){
    size_t dim = Parameters_count(var->parameters);
    
    for (int i = 0; i < dim; i++) {
        Parameter_restore(Parameters_at(var->parameters, i));
    }
    
    for (int i = 0; i < dim; i++) {
        double sigma = exp(Parameters_value(var->var_parameters, i+dim));
        if(!Parameters_at(var->var_parameters, i)->estimate){
            Parameters_set_value(var->var_parameters, i, -10 + sigma);
        }
        //        double mu = Parameters_value(var->var_parameters, i);
        //        if(!Parameters_at(var->var_parameters, i)->estimate){
        //            mu = -10 + sigma;
        //        }
        ////        double mean = exp(mu + sigma*sigma/2.0);
        ////        double mode = exp(mu - sigma*sigma);
        ////        printf("%f %f %f %f %s\n", Parameters_value(var->parameters, i), mode,
        ////               mu, sigma, Parameters_name(var->parameters, i));
    }
}

//MARK: block KL(Q||P)

double variational_klqp_elbo(variational_t* var){
    double elbo = 0;
    Model* posterior = var->posterior;
    
    double entropy = 0;
    for(int i = 0; i < var->block_count; i++){
        variational_block_t* block = var->blocks[i];
        if (block->use_entropy) {
            entropy += block->entropy(block);
        }
        else{
            // estimate entropy
        }
    }
    int inf_count = 0;
    
    for (int i = 0; i < var->elbo_samples; i++) {
        double jacobian = 0.0;
        double logQ = 0.0;
        for(int j = 0; j < var->block_count; j++){
            variational_block_t* block = var->blocks[j];
            block->sample1(block, &jacobian);
            if(!block->use_entropy){
                logQ += block->logQ(block, Vector_data(block->etas));
            }
        }
        double logP = posterior->logP(posterior);
        if(!isinf(logP) && !isnan(logP))elbo += logP -logQ + jacobian;
		else{
			inf_count++;
			i--;
		}
		if (inf_count == 10) {
			return elbo;
		}
    }
    return elbo/(var->elbo_samples-inf_count) + entropy;
}

double variational_klqp_elbo_multi(variational_t* var){
    double elbo = 0;
    Model* posterior = var->posterior;
    
    double* elbos = dvector(var->elbo_multi);
    
    for (int s = 0; s < var->elbo_samples; s++) {
        double mean = 0;
        int inf_count = 0;
        for (int k = 0; k < var->elbo_multi; k++) {
            double jacobian = 0.0;
            double logQ = 0;
            for(int i = 0; i < var->block_count; i++){
                variational_block_t* block = var->blocks[i];
                block->sample1(block, &jacobian);
                logQ += block->logQ(block, Vector_data(block->etas));
            }
            double logP = posterior->logP(posterior);
            double ratio = logP - logQ + jacobian;
            if(!isinf(ratio) && !isnan(ratio)) elbos[k-inf_count] = ratio;
            else inf_count++;
        }
        if(inf_count == var->elbo_multi){
            free(elbos);
            return elbo;
        }
        double max = elbos[0];
        for (int k = 1; k < var->elbo_multi - inf_count; k++) {
            if(elbos[k] > max) max = elbos[k];
        }
        for (int k = 0; k < var->elbo_multi - inf_count; k++) {
            mean += exp(elbos[k] - max);
        }

        mean /= var->elbo_multi - inf_count;
        elbo += log(mean) + max;
    }
    free(elbos);
    return elbo/var->elbo_samples;
}


void variational_klqp_grad_elbo(variational_t* var, const Parameters* parameters, double* grads){
    size_t dim = Parameters_count(parameters);
    memset(grads, 0, sizeof(double)*Parameters_count(parameters));
	size_t failures = 0;
    for (int i = 0; i < var->grad_samples; i++) {
        // sample within each block. etas should be stored for jacobians
        for(int j = 0; j < var->block_count; j++){
            variational_block_t* block = var->blocks[j];
            block->sample2(block, parameters);
        }
        
		size_t j = 0;
        for( ; j < var->block_count; j++){
            variational_block_t* block = var->blocks[j];
            block->grad_elbo(block, parameters, grads);
			if (isnan(grads[0])) {
				failures++;
				break;
			}
        }
		if (j != var->block_count) {
			i--;
			if (failures == 10) {
				break;
			}
		}
    }
    
    for (int j = 0; j < dim; j++) {
        grads[j] /= var->grad_samples;
    }
    
    // Add gradient entropy if needed
    for(int j = 0; j < var->block_count; j++){
        variational_block_t* block = var->blocks[j];
        if (block->use_entropy) {
            block->grad_entropy(block, parameters, grads);
        }
    }
}

//MARK: block normal meanfield

void klqp_block_meanfield_normal_sample1(variational_block_t* var, double* jacobian){
    size_t dim = Parameters_count(var->parameters);
    double* etas = Vector_mutable_data(var->etas);
    
    for (int j = 0; j < dim; j++) {
        Parameter* p = Parameters_at(var->parameters, j);
        Parameter* var_p_mu = Parameters_at(var->var_parameters[0], j);
        Parameter* var_p_sigma = Parameters_at(var->var_parameters[1], j);
        
        if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
        
        double mu;
        double sigma = Parameter_value(var_p_sigma);
        
        if(Parameter_estimate(var_p_mu)){
            mu = Parameter_value(var_p_mu);
        }
        else{
            //                if(Parameter_value(p) < 1.0e-6){
            mu = -10 + sigma * sigma;
            //                }
            //                else{
            //                    mu = log(Parameter_value(p)) + sigma*sigma;
            //                }
        }
        etas[j] = rnorm();
        double zeta = etas[j] * sigma + mu;
        double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), jacobian);
        Parameter_set_value(p, theta);
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

double klqp_block_meanfield_normal_entropy(variational_block_t* var){
    // Entropy: 0.5(1 + log(2 \pi s^2))
    // 0.5(1+log(2\pi) + sum(log(sigma))
    size_t count = 0;
    double entropy = 0;
    size_t dim = Parameters_count(var->parameters);
    for (size_t i = 0; i < dim; i++) {
        if(Parameters_estimate(var->var_parameters[1], i)){
            entropy += log(Parameters_value(var->var_parameters[1], i));
            count++;
        }
    }
    entropy += 0.5 * count * (1.0 + LOG_TWO_PI);
    return entropy;
}

// parameters are the variational parameters
void klqp_block_meanfield_normal_sample2(variational_block_t* var, const Parameters* parameters){
    size_t dim = Parameters_count(var->parameters);
    size_t opt_param_dim = Parameters_count(parameters);
    
    // find first mu
    size_t mu_idx = 0;
    for ( ; mu_idx < opt_param_dim; mu_idx++) {
        if(Parameters_at(parameters, mu_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (mu_idx == opt_param_dim) return;
    
    double* etas = Vector_mutable_data(var->etas);
    
    for (int idx = 0; idx < dim; idx++) {
        Parameter* p = Parameters_at(var->parameters, idx);
        double mu = Parameters_value(var->var_parameters[0], idx);
        double sigma = Parameters_value(var->var_parameters[1], idx);
        
        etas[idx] = rnorm();
        double zeta = etas[idx] * sigma + mu;
        double theta = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
        Parameter_set_value(p, theta);
    }
}

void klqp_block_meanfield_normal_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads){
    
    // find first mu and sigma
    size_t mu_idx = 0;
    size_t sigma_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
    for ( ; mu_idx < opt_param_dim; mu_idx++) {
        if(Parameters_at(parameters, mu_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (mu_idx == opt_param_dim) return;
    
    for ( ; sigma_idx < opt_param_dim; sigma_idx++) {
        if(Parameters_at(parameters, sigma_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
    
    Model* posterior = var->posterior;
    size_t dim = Parameters_count(var->parameters);
    size_t simplex_parameter_count = var->simplex_parameter_count;
    const double* etas = Vector_data(var->etas);
    int idx = 0;
    if (simplex_parameter_count > 0) {
        for(int s = 0; s < var->simplex_count; s++){
            Simplex* simplex = var->simplices[s]->obj;
            for(int k = 0; k < simplex->K-1; k++){
                Parameter* p = Parameters_at(var->parameters, idx);
                double mu = Parameters_value(var->var_parameters[0], idx);
                double sigma = Parameters_value(var->var_parameters[1], idx);
				double dlogP;
				if(var->derivative_type == DERIVATIVE_ANALYTICAL){
                	dlogP = posterior->dlogP(posterior, p);
				}
				else{
					dlogP = Model_first_derivative(posterior, p, var->numerical_eps);
				}
                double zeta = etas[idx] * sigma + mu;
                double gldits = 1.0/zeta + 1.0/(zeta - 1.0); // grad log det transform of stick
                const double gldit = grad_log_det_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p));
                double grad_mu = dlogP * grad_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p)) + gldit + gldits;
                grads[mu_idx + idx] += grad_mu;
                grads[sigma_idx + idx] += grad_mu * etas[idx] * sigma;
				
				if(!var->use_entropy){
					double sigma2 = sigma*sigma;
					grads[mu_idx + idx] -= -(mu - zeta)/sigma2;
					grads[sigma_idx + idx] -= (mu*mu - sigma2 -2.0*mu*zeta + zeta*zeta)/(sigma2*sigma);
				}
				
                idx++;
            }
        }
    }
    
    for ( ; idx < dim; idx++) {
        Parameter* p = Parameters_at(var->parameters, idx);
        double mu = Parameters_value(var->var_parameters[0], idx);
        double sigma = Parameters_value(var->var_parameters[1], idx);
        
        double dlogP;
		if(var->derivative_type == DERIVATIVE_ANALYTICAL){
			dlogP = posterior->dlogP(posterior, p);
		}
		else{
			dlogP = Model_first_derivative(posterior, p, var->numerical_eps);
		}
        double zeta = etas[idx]*sigma + mu;
        double gldit = grad_log_det_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p));
        double grad_mu = dlogP * grad_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p)) + gldit;
        grads[mu_idx + idx] += grad_mu;
        grads[sigma_idx + idx] += grad_mu * etas[idx] * sigma;
		
		if(!var->use_entropy){
			double sigma2 = sigma*sigma;
			grads[mu_idx + idx] -= -(mu - zeta)/sigma2;
			grads[sigma_idx + idx] -= (mu*mu - sigma2 -2.0*mu*zeta + zeta*zeta)/(sigma2*sigma);
		}
    }
}

void klqp_block_meanfield_normal_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads){
    // Entropy: 0.5(1 + log(2 \pi s^2))
    // 0.5(1+log(2\pi) + sum(log(sigma))
    
    // find first sigma
    size_t sigma_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
    for ( ; sigma_idx < opt_param_dim; sigma_idx++) {
        if(Parameters_at(parameters, sigma_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
    // not beeing optimized
    if (sigma_idx == opt_param_dim) return;
    
    size_t dim = Parameters_count(var->parameters);
    for (int j = 0; j < dim; j++) {
        grads[j + sigma_idx] += 1.0;
    }
}

double klqp_block_meanfield_normal_logP(variational_block_t* var, const double* values){
    size_t dim = Parameters_count(var->parameters);
    double logP = 0;
    for (size_t i = 0; i < dim; i++) {
        Parameter* p = Parameters_at(var->parameters, i);
        double mu = Parameters_value(var->var_parameters[0], i);
        double sd = Parameters_value(var->var_parameters[1], i);
        double zeta = transform2(values[i], Parameter_lower(p), Parameter_upper(p), &logP);
        logP += dnorml(zeta, mu, sd);
    }
    return logP;
}


double klqp_block_meanfield_normal_logQ(variational_block_t* var, const double* values){
    size_t dim = Parameters_count(var->parameters);
    double logP = 0;
    for (size_t i = 0; i < dim; i++) {
        double sd = Parameters_value(var->var_parameters[1], i);
		logP += log(gsl_ran_gaussian_pdf(values[i]*sd, sd));
    }
    return logP;
}


void klqp_block_meanfield_normal_sample(variational_block_t* var, double* values){
    size_t dim = Parameters_count(var->parameters);
    
    for (int j = 0; j < dim; j++) {
        Parameter* p = Parameters_at(var->parameters, j);
        Parameter* var_p_mu = Parameters_at(var->var_parameters[0], j);
        Parameter* var_p_sigma = Parameters_at(var->var_parameters[1], j);
        
        if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
        
        double mu;
        double sigma = Parameter_value(var_p_sigma);
        
        if(Parameter_estimate(var_p_mu)){
            mu = Parameter_value(var_p_mu);
        }
        else{
            //                if(Parameter_value(p) < 1.0e-6){
            mu = -10 + sigma * sigma;
            //                }
            //                else{
            //                    mu = log(Parameter_value(p)) + sigma*sigma;
            //                }
        }
        double zeta = rnorm() * sigma + mu;
        values[j] = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
    }
}

bool klqp_block_meanfield_normal_sample_some(variational_block_t* var, const Parameters* parameters, double* values){
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
                    double mu = Parameters_value(var->var_parameters[0], j);
                    double sd = Parameters_value(var->var_parameters[1], j);
                    const double zeta = rnorm() * exp(sd) + mu;
                    values[i] = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
                }
                break;
            }
        }
    }
    return true;
}

//MARK: block normal fullrank

void klqp_block_fullrank_normal_sample1(variational_block_t* var, double* jacobian){
    size_t dim = Parameters_count(var->parameters);
    double* etas = Vector_mutable_data(var->etas);
    for (int j = 0; j < dim; j++) {
        etas[j] = rnorm();
    }
    
    size_t row = 0;
    for (int j = 0; j < dim; j++) {
        double temp = 0;
        // multiply L_j and eta
        for (int k = 0; k < j+1; k++) {
            if(Parameters_estimate(var->var_parameters[1], row)){
                temp += Parameters_value(var->var_parameters[1], row)*etas[k];
            }
            row++;
        }
        double zeta = temp + Parameters_value(var->var_parameters[0], j); // add mu
        Parameter* p = Parameters_at(var->parameters, j);
        double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), jacobian);
        Parameter_set_value(p, theta);
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

double klqp_block_fullrank_normal_entropy(variational_block_t* var){
    size_t dim = Parameters_count(var->parameters);
    double entropy = 0;
    int r = -1;
    int count = 0;
    for (int d = 0; d < dim; ++d) {
        r += d+1;
//        if(Parameters_estimate(var->var_parameters, r)){
            double tmp = fabs(Parameters_value(var->var_parameters[1], r));
            if (tmp != 0.0) entropy += log(tmp);
            count++;
//        }
    }
    entropy +=  0.5 * (1.0 + LOG_TWO_PI) * count;
    return entropy;
}

// parameters are the variational parameters
void klqp_block_fullrank_normal_sample2(variational_block_t* var, const Parameters* parameters){
    size_t dim = Parameters_count(var->parameters);
    double* etas = Vector_mutable_data(var->etas);
    for (int j = 0; j < dim; j++) {
        etas[j] = rnorm();
    }
    
    size_t row = 0;
    for (int j = 0; j < dim; j++) {
        double temp = 0;
        // multiply L_j and eta
        for (int k = 0; k < j+1; k++) {
//            if(Parameters_estimate(var->var_parameters[1], row)){
                temp += Parameters_value(var->var_parameters[1], row)*etas[k];
//            }
            row++;
        }
        double zeta = temp + Parameters_value(var->var_parameters[0], j); // add mu
        etas[j+dim] = zeta;
        Parameter* p = Parameters_at(var->parameters, j);
        double theta = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
        Parameter_set_value(p, theta);
    }
}

void klqp_block_fullrank_normal_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads){
    
    // find first mu and sigma
    size_t mu_idx = 0;
    size_t sigma_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
    for ( ; mu_idx < opt_param_dim; mu_idx++) {
        if(Parameters_at(parameters, mu_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (mu_idx == opt_param_dim) return;
    
    for ( ; sigma_idx < opt_param_dim; sigma_idx++) {
        if(Parameters_at(parameters, sigma_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
    
    const double* etas = Vector_data(var->etas);
    Model* posterior = var->posterior;
    size_t simplex_parameter_count = var->simplex_parameter_count;
    size_t  dim = Parameters_count(var->parameters);
    size_t row = 0;
    for (int k = 0; k < dim; k++){
        Parameter* p = Parameters_at(var->parameters, k);
        double dlogP;
  		if(var->derivative_type == DERIVATIVE_ANALYTICAL){
			dlogP = posterior->dlogP(posterior, p);
		}
		else{
			dlogP = Model_first_derivative(posterior, p, var->numerical_eps);
		}
        double zeta = etas[dim+k];
        const double gldit = grad_log_det_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p));
        double grad_mu = dlogP * grad_inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p)) + gldit;
        if (k < simplex_parameter_count) {
            grad_mu += 1.0/zeta + 1.0/(zeta - 1.0); // grad log det transform of stick
        }
        grads[mu_idx + k] += grad_mu;
        for (int j = 0; j < k+1; j++) {
            grads[sigma_idx + row] += grad_mu*etas[j];
            row++;
        }
    }
}

void klqp_block_fullrank_normal_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads){
    // find first sigma
    size_t sigma_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
    for ( ; sigma_idx < opt_param_dim; sigma_idx++) {
        if(Parameters_at(parameters, sigma_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
    // not beeing optimized
    if (sigma_idx == opt_param_dim) return;
    
    size_t dim = Parameters_count(var->parameters);
    int diag = -1;
    for (int i = 0; i < dim; i++) {
        diag += i+1;
//        if(Parameters_estimate(var->var_parameters, i)){
            grads[diag + sigma_idx] += 1.0/Parameters_value(var->var_parameters[1], diag);
//        }
    }
}

double klqp_block_fullrank_normal_logP(variational_block_t* var, const double* values){
    size_t dim = Parameters_count(var->parameters);
    gsl_vector * mu = gsl_vector_calloc(dim);
    gsl_matrix * L = gsl_matrix_calloc(dim, dim);
    gsl_vector * x = gsl_vector_calloc(dim);
    gsl_vector * work = gsl_vector_calloc(dim);
    
    double logJac = 0;
    size_t row = 0;
    for (int i = 0; i < dim; i++) {
        Parameter* p = Parameters_at(var->parameters, i);
        double zeta = transform(values[i], Parameter_lower(p), Parameter_upper(p));
        gsl_vector_set(x, i, zeta);
        logJac += zeta;
        gsl_vector_set(mu, i, Parameters_value(var->var_parameters[0], i));
        for (int j = 0; j <= i; j++) {
            gsl_matrix_set(L, i, j, Parameters_value(var->var_parameters[1], row));
            row++;
        }
    }
    double logP = 0;
    gsl_ran_multivariate_gaussian_log_pdf (x,  mu, L, &logP, work);
    
    gsl_vector_free(mu);
    gsl_matrix_free(L);
    gsl_vector_free(x);
    gsl_vector_free(work);
    
    return logP + logJac;
}

double klqp_block_fullrank_normal_logQ(variational_block_t* var, const double* values){
    size_t dim = Parameters_count(var->parameters);
    gsl_vector * mu = gsl_vector_calloc(dim);
    gsl_matrix * L = gsl_matrix_calloc(dim, dim);
    gsl_vector * x = gsl_vector_calloc(dim);
    gsl_vector * work = gsl_vector_calloc(dim);
    
    size_t row = 0;
    for (int i = 0; i < dim; i++) {
        gsl_vector_set(x, i, values[i]);
        gsl_vector_set(mu, i, Parameters_value(var->var_parameters[0], i));
        for (int j = 0; j <= i; j++) {
            gsl_matrix_set(L, i, j, Parameters_value(var->var_parameters[1], row));
            row++;
        }
    }
    double logP = 0;
    gsl_ran_multivariate_gaussian_log_pdf (x,  mu, L, &logP, work);
    
    gsl_vector_free(mu);
    gsl_matrix_free(L);
    gsl_vector_free(x);
    gsl_vector_free(work);
    
    return logP;
}


void klqp_block_fullrank_normal_sample(variational_block_t* var, double* values){
    size_t dim = Parameters_count(var->parameters);
    double* etas = Vector_mutable_data(var->etas);
    for (int j = 0; j < dim; j++) {
        etas[j] = rnorm();
    }
    
    size_t row = 0;
    for (int j = 0; j < dim; j++) {
        double temp = 0;
        // multiply L_j and eta
        for (int k = 0; k < j+1; k++) {
            if(Parameters_estimate(var->var_parameters[1], row)){
                temp += Parameters_value(var->var_parameters[1], row)*etas[k];
            }
            row++;
        }
        double zeta = temp + Parameters_value(var->var_parameters[0], j); // add mu
        Parameter* p = Parameters_at(var->parameters, j);
        double theta = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
        values[j] = theta;
    }
}

bool klqp_block_fullrank_normal_sample_some(variational_block_t* var, const Parameters* parameters, double* values){
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
                    double mu = Parameters_value(var->var_parameters[0], j);
                    double sd = Parameters_value(var->var_parameters[1], j);
                    const double zeta = rnorm() * exp(sd) + mu;
                    values[i] = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
                }
                break;
            }
        }
    }
    return true;
}


#pragma mark - normal meanfield

double klqp_meanfield_normal_elbo(variational_t* var){
	if (var->initialized == false) {
 		klqp_meanfield_normal_init(var);
		var->initialized = true;
	}
	//	printf("elbo\n");
	double elbo = 0;
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
	
	// Entropy: 0.5(1 + log(2 \pi s^2))
	// 0.5(1+log(2\pi) + sum(log(sigma))
	size_t count = 0;
	double entropy = 0;
	for (size_t i = 0; i < dim; i++) {
		if(Parameters_estimate(var->var_parameters, dim+i)){
			entropy += Parameters_value(var->var_parameters, dim+i);
			count++;
		}
	}
	entropy += 0.5 * count * (1.0 + LOG_TWO_PI);
	int inf_count = 0;
	
	for (int i = 0; i < var->elbo_samples; i++) {
		double jacobian = 0.0;
		
		for (int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
			
			if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
			
			double mu;
			double sigma = exp(Parameter_value(var_p_sigma));
			
			if(Parameter_estimate(var_p_mu)){
				mu = Parameter_value(var_p_mu);
			}
			else{
				//				if(Parameter_value(p) < 1.0e-6){
				mu = -10 + sigma * sigma;
				//				}
				//				else{
				//					mu = log(Parameter_value(p)) + sigma*sigma;
				//				}
			}
			double zeta = rnorm() * sigma + mu;
			double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), &jacobian);
			Parameter_set_value(p, theta);
		}
		
		if (var->simplex_count > 0) {
			for(int s = 0; s < var->simplex_count; s++){
				Simplex* simplex = var->simplices[s]->obj;
				const double* constrained_values = simplex->get_values(simplex);
				double stick = 1;
				for(int k = 0; k < simplex->K-1; k++){
					jacobian += log(stick);
					stick -= constrained_values[k];
				}
			}
		}
		double logP = posterior->logP(posterior);
		if(!isinf(logP))elbo += logP + jacobian;
		else inf_count++;
		
		if(isinf(elbo) || isnan(elbo)){
			return elbo;
		}
	}
	return elbo/(var->elbo_samples-inf_count) + entropy;
}

double klqp_meanfield_normal_elbo_multi(variational_t* var){
	if (var->initialized == false) {
		klqp_meanfield_normal_init(var);
		var->initialized = true;
	}
	double elbo = 0;
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
	
	double* elbos = dvector(var->elbo_multi);
	
	for (int i = 0; i < var->elbo_samples; i++) {
		double mean = 0;
		int inf_count = 0;
		for (int k = 0; k < var->elbo_multi; k++) {
			double jacobian = 0.0;
			double logQ = 0;
			for (int j = 0; j < dim; j++) {
				Parameter* p = Parameters_at(var->parameters, j);
				Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
				Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
				
				if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
				
				double mu;
				double sigma = exp(Parameter_value(var_p_sigma));
				
				if(Parameter_estimate(var_p_mu)){
					mu = Parameter_value(var_p_mu);
				}
				else{
					//				if(Parameter_value(p) < 1.0e-6){
					mu = -10 + sigma * sigma;
					//				}
					//				else{
					//					mu = log(Parameter_value(p)) + sigma*sigma;
					//				}
				}
				double eta = rnorm();
				double zeta = eta * sigma + mu;
				logQ += dnorml(zeta, mu, sigma);
				double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), &jacobian);
				Parameter_set_value(p, theta);
			}
			
			if (var->simplex_count > 0) {
				for(int s = 0; s < var->simplex_count; s++){
					Simplex* simplex = var->simplices[s]->obj;
					const double* constrained_values = simplex->get_values(simplex);
					double stick = 1;
					for(int k = 0; k < simplex->K-1; k++){
						jacobian += log(stick);
						stick -= constrained_values[k];
					}
				}
			}
			
			double logP = posterior->logP(posterior);
			double ratio = logP + jacobian - logQ;
			if(!isinf(ratio)) elbos[k-inf_count] = ratio;
			else inf_count++;
		}
		if(inf_count == var->elbo_multi){
			free(elbos);
			return elbo;
		}
		double max = elbos[0];
		for (int k = 1; k < var->elbo_multi - inf_count; k++) {
			if(elbos[k] > max) max = elbos[k];
		}
		for (int k = 0; k < var->elbo_multi - inf_count; k++) {
			mean += exp(elbos[k] - max);
		}

		mean /= var->elbo_multi - inf_count;
		elbo += log(mean) + max;
	}
	free(elbos);
	return elbo/var->elbo_samples;
}

void klqp_meanfield_normal_grad_elbo(variational_t* var, const Parameters* parameters, double* grads){
	if (var->initialized == false) {
 		klqp_meanfield_normal_init(var);
		// save for later
		//		for (int i = 0; i < Parameters_count(var->parameters); i++) {
		//			Parameter_store(Parameters_at(var->parameters, i));
		//		}
		var->initialized = true;
	}
	
	Model* posterior = var->posterior;
	size_t dim = Parameters_count(var->parameters);
	double* eta = dvector(dim);
	double* zeta = dvector(dim);
	memset(grads, 0, sizeof(double)*dim*2);
	
	int simplex_parameter_count = 0;
	if (var->simplex_count > 0) {
		for(int s = 0; s < var->simplex_count; s++){
			Simplex* simplex = var->simplices[s]->obj;
			simplex_parameter_count += simplex->K - 1;
		}
	}
	
	// process simplices first
	for (int i = 0; i < var->grad_samples; i++) {
		if (var->simplex_count > 0) {
			int idx = 0;
			for(int s = 0; s < var->simplex_count; s++){
				Simplex* simplex = var->simplices[s]->obj;
				for(int k = 0; k < simplex->K-1; k++){
					Parameter* p = Parameters_at(var->parameters, idx);
					double mu = Parameters_value(var->var_parameters, idx);
					double sigma = exp(Parameters_value(var->var_parameters, idx+dim));
					
					eta[idx] = rnorm();
					zeta[idx] = eta[idx] * sigma + mu;
					double theta = inverse_transform2(zeta[idx], Parameter_lower(p), Parameter_upper(p));
					Parameter_set_value(p, theta);
					idx++;
				}
			}
			idx = 0;
			for(int s = 0; s < var->simplex_count; s++){
				Simplex* simplex = var->simplices[s]->obj;
				for(int k = 0; k < simplex->K-1; k++){
					Parameter* p = Parameters_at(var->parameters, idx);
					double dlogP = posterior->dlogP(posterior, p);
					double gldits = 1.0/zeta[idx] + 1.0/(zeta[idx]-1.0); // grad log det transform of stick
					const double gldit = grad_log_det_inverse_transform(zeta[idx], Parameter_lower(p), Parameter_upper(p));
					double grad_mu = dlogP * grad_inverse_transform(zeta[idx], Parameter_lower(p), Parameter_upper(p)) + gldit + gldits;
					grads[idx] += grad_mu;
					grads[dim+idx] += grad_mu * eta[idx] * exp(Parameters_value(var->var_parameters, dim+idx));
					idx++;
				}
			}
		}
		
		for ( int j = simplex_parameter_count; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
			
			if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
			
			double mu;
			double sigma = exp(Parameter_value(var_p_sigma));
			
			if(Parameter_estimate(var_p_mu)){
				mu = Parameter_value(var_p_mu);
			}
			else{
				//				if(Parameter_value(p) < 1.0e-6){
				mu = -10 + sigma * sigma;
				//				}
				//				else{
				//					mu = log(Parameter_value(p)) + sigma*sigma;
				//				}
			}
			
			eta[j] = rnorm();
			zeta[j] = eta[j] * sigma + mu;
			
			double theta = inverse_transform2(zeta[j], Parameter_lower(p), Parameter_upper(p));
			//						printf("%f %f %f %d\n", mu, sigma, theta, Parameters_at(var->var_parameters, j)->estimate);
			Parameter_set_value(p, theta);
		}
		
		for ( int j = simplex_parameter_count; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
			
			if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
			
			//			double dlogP = Model_first_derivative(posterior, p, 0.001);
			double dlogP = posterior->dlogP(posterior, p);
			//			printf("%f %f %s\n", dlogP, Model_first_derivative(posterior, p, 0.001), Parameter_name(p));
			//			printf("dlogP %f %e\n",dlogP, Parameter_value(p));
			const double gldit = grad_log_det_inverse_transform(zeta[j], Parameter_lower(p), Parameter_upper(p));
			double grad_mu = dlogP * grad_inverse_transform(zeta[j], Parameter_lower(p), Parameter_upper(p)) + gldit;
			grads[j] += grad_mu;
			grads[dim+j] += grad_mu * eta[j] * exp(Parameters_value(var->var_parameters, dim+j));
			//			printf("b %d %e %f %f %f\n", j, Parameter_value(p), dlogP, grads[j], grads[j+dim]);
		}
	}
	
	for (int j = 0; j < dim; j++) {
		
		Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
		Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
		
		if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
		
		grads[j] /= var->grad_samples;
		grads[dim+j] = grads[dim+j]/var->grad_samples + 1.0;
		//printf("gradient %f %f\n", grads[i] ,grads[dim+i] );
	}
	//printf("\n\n");
	free(eta);
	free(zeta);
	
	if (var->file != NULL) {
		fprintf(var->file, "%zu", var->iter++);
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%e", Parameters_value(var->var_parameters, i));
		}
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%e", grads[i]);
		}
		fprintf(var->file, "\n");
	}
}


// TODO: should use a jacobian fucntion instead transform().
// In this function it works for transform == log but not anything else
double klqp_meanfield_normal_logP(variational_t* var, const double* values){
	size_t dim = Parameters_count(var->parameters);
	double logP = 0;
	for (size_t i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		double mu = Parameters_value(var->var_parameters, i);
		double sd = exp(Parameters_value(var->var_parameters, i+dim));
		double zeta = transform(values[i], Parameter_lower(p), Parameter_upper(p));
		logP += dnorml(zeta, mu, sd) - zeta;
	}
	return logP;
}

// check success
bool klqp_meanfield_normal_sample(variational_t* var, double* values){
	if(!var->ready_to_sample){
		var->finalize(var);
		var->ready_to_sample = true;
	}
	size_t dim = Parameters_count(var->parameters);
	for (size_t i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		if (!p->estimate) {
			values[i] = 0;
			continue;
		}
		double mu = Parameters_value(var->var_parameters, i);
		double sd = exp(Parameters_value(var->var_parameters, i+dim));
		const double zeta = rnorm() * sd + mu;
		values[i] = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
	}
	return true;
}

bool klqp_meanfield_normal_sample_some(variational_t* var, const Parameters* parameters, double* values){
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
					double mu = Parameters_value(var->var_parameters, j);
					double sd = Parameters_value(var->var_parameters, j+dim);
					const double zeta = rnorm() * exp(sd) + mu;
					values[i] = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
				}
				break;
			}
		}
	}
	return true;
}

void klqp_meanfield_normal_log_samples(variational_t* var, FILE* file){
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
			Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
			
			if(!Parameter_estimate(var_p_mu) && !Parameter_estimate(var_p_sigma)) continue;
			
			double mu;
			double sigma = exp(Parameter_value(var_p_sigma));
			
			if(Parameter_estimate(var_p_mu)){
				mu = Parameter_value(var_p_mu);
			}
            else{
				mu = -10 + sigma * sigma;
			}
            double zeta = rnorm() * sigma + mu;
            logQ += dnorml(zeta, mu, sigma);
            samples[j] = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), &jacobian);
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

//MARK: lognormal meanfield

// elbo is calculated using the entropy of the variational distribution
double klqp_meanfield_lognormal_elbo(variational_t* var){
    if (var->initialized == false) {
        klqp_meanfield_normal_init(var);
        var->initialized = true;
    }

    double elbo = 0;
    size_t dim = Parameters_count(var->parameters);
    Model* posterior = var->posterior;
    
    int inf_count = 0;
    // entropy = log_2(sigma exp(mu + 1/2) sqrt(2 pi))
    // entropy (nats) = mu + 1/2(log(2 pi e sigma^2))
    double entropy = log(2.0*M_PI*M_E)/2.0*dim;
    for (int j = 0; j < dim; j++) {
        double mu = Parameters_value(var->var_parameters, j);
        double sigma = exp(Parameters_value(var->var_parameters, j+dim));
        entropy += mu + log(sigma*sigma)/2.0;
    }
    
    for (int i = 0; i < var->elbo_samples; i++) {
//        double logQ = 0;
        for (int j = 0; j < dim; j++) {
            Parameter* p = Parameters_at(var->parameters, j);
            double mu = Parameters_value(var->var_parameters, j);
            double sigma = exp(Parameters_value(var->var_parameters, j+dim));
            double theta = gsl_ran_lognormal(var->rng, mu, sigma);
//            logQ += log(gsl_ran_lognormal_pdf(theta, mu, sigma));
            Parameter_set_value(p, theta);
        }
        
        double logP = posterior->logP(posterior);
        if(!isinf(logP)) elbo += logP;// - logQ;
        else inf_count++;
        
        if(isinf(elbo) || isnan(elbo)){
            return elbo;
        }
    }
    return elbo/(var->elbo_samples-inf_count) + entropy;
}

void klqp_meanfield_lognormal_grad_elbo(variational_t* var, const Parameters* parameters, double* grads){
    if (var->initialized == false) {
        //klqp_meanfield_normal_init(var);
        // save for later
        //        for (int i = 0; i < Parameters_count(var->parameters); i++) {
        //            Parameter_store(Parameters_at(var->parameters, i));
        //        }
        var->initialized = true;
    }
    
    Model* posterior = var->posterior;
    size_t dim = Parameters_count(var->parameters);
    double* thetas = dvector(dim);
    memset(grads, 0, sizeof(double)*dim*2);
    // entropy = mu + 1/2(log(2 pi e sigma^2))
    // dE/d\sigma = 1/sigma
    // dE/d\mu = 1
    
    for (int i = 0; i < var->grad_samples; i++) {
        for ( int j = 0; j < dim; j++) {
            Parameter* p = Parameters_at(var->parameters, j);
            double mu = Parameters_value(var->var_parameters, j);
            double sigma = exp(Parameters_value(var->var_parameters, j+dim));
            
            thetas[j] = gsl_ran_lognormal(var->rng, mu, sigma);
            Parameter_set_value(p, thetas[j]);
        }
        
        for ( int j = 0; j < dim; j++) {
            Parameter* p = Parameters_at(var->parameters, j);
            double mu = Parameters_value(var->var_parameters, j);
            double sigma = exp(Parameters_value(var->var_parameters, j+dim));
            
            double dlogP = posterior->dlogP(posterior, p);
            grads[j] += dlogP * thetas[j];// + 1.0;
            grads[dim+j] += dlogP * thetas[j] * (log(thetas[j]) - mu)/sigma;// + (log(thetas[j]) - mu)/sigma + 1.0/sigma;
        }
    }
    
    for (int j = 0; j < dim; j++) {
        double sigma = exp(Parameters_value(var->var_parameters, j+dim));
        // last terms are entropy derivatives
        grads[j] = grads[j]/var->grad_samples + 1.0;
        grads[dim+j] = grads[dim+j]/var->grad_samples + 1.0/(sigma);
    }
    free(thetas);
    
    if (var->file != NULL) {
        fprintf(var->file, "%zu", var->iter++);
        for(int i = 0; i < Parameters_count(var->var_parameters); i++){
            fprintf(var->file, ",%e", Parameters_value(var->var_parameters, i));
        }
        for(int i = 0; i < Parameters_count(var->var_parameters); i++){
            fprintf(var->file, ",%e", grads[i]);
        }
        fprintf(var->file, "\n");
    }
}

// assumes that values are in R+
double klqp_meanfield_lognormal_logP(variational_t* var, const double* values){
    size_t dim = Parameters_count(var->parameters);
    double logP = 0;
    for (size_t i = 0; i < dim; i++) {
        double mu = Parameters_value(var->var_parameters, i);
        double sigma = exp(Parameters_value(var->var_parameters, i+dim));
        logP += log(gsl_ran_lognormal_pdf(values[i], mu, sigma));
    }
    return logP;
}

// check success
bool klqp_meanfield_lognormal_sample(variational_t* var, double* values){
    if(!var->ready_to_sample){
        var->finalize(var);
        var->ready_to_sample = true;
    }
    size_t dim = Parameters_count(var->parameters);
    for (size_t i = 0; i < dim; i++) {
        Parameter* p = Parameters_at(var->parameters, i);
        if (!p->estimate) {
            values[i] = 0;
            continue;
        }
        double mu = Parameters_value(var->var_parameters, i);
        double sigma = exp(Parameters_value(var->var_parameters, i+dim));
        values[i] = gsl_ran_lognormal(var->rng, mu, sigma);
    }
    return true;
}

bool klqp_meanfield_lognormal_sample_some(variational_t* var, const Parameters* parameters, double* values){
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
                    double mu = Parameters_value(var->var_parameters, j);
                    double sigma = Parameters_value(var->var_parameters, j+dim);
                    values[i] = gsl_ran_lognormal(var->rng, mu, sigma);
                }
                break;
            }
        }
    }
    return true;
}

void klqp_meanfield_lognormal_log_samples(variational_t* var, FILE* file){
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
            Parameter* var_p_mu = Parameters_at(var->var_parameters, j);
            Parameter* var_p_sigma = Parameters_at(var->var_parameters, j+dim);
            
            double mu = Parameter_value(var_p_mu);
            double sigma = exp(Parameter_value(var_p_sigma));
            
            if(i==0)printf("%f %f\n", mu, sigma);
            samples[j] = gsl_ran_lognormal(var->rng, mu, sigma);
            logQ += dnorml(samples[j], mu, sigma);
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


//MARK: Fullrank

void klqp_fullrank_normal_init(variational_t* var){
	
	for (int i = 0; i < Parameters_count(var->parameters); i++) {
		Parameter_store(Parameters_at(var->parameters, i));
	}
	return;
	
	Model* posterior = var->posterior;
	size_t dim = Parameters_count(var->parameters);
	size_t n = 2*dim+(dim*(dim-1))/2;
	
	// mu vector
	for (int i = 0; i < dim; i++) {
		double mu = log(Parameters_value(var->parameters, i));
		if (isnan(mu) || isinf(mu)) {
			mu = 1;
		}
		Parameters_set_value(var->var_parameters, i, mu); // use median
	}
	
	for (int i = dim; i < n; i++) {
		Parameters_set_value(var->var_parameters, i, 1);
	}
	return;
	double* d2lnl2 = dvector(dim*dim);
	double* hessian = dvector(dim*dim);
	double* diag = dvector(dim);
	
	for (int i = 0; i < dim; i++) {
		d2lnl2[i*dim+i] = Model_second_derivative(posterior, Parameters_at(var->parameters, i), NULL, 0.0001);
		for (int j = i+1; j < dim; j++) {
			d2lnl2[i*dim+j] = d2lnl2[j*dim+i] = Model_mixed_derivative(posterior, Parameters_at(var->parameters, i), Parameters_at(var->parameters, j));
		}
	}
	
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			hessian[i*dim+j] = -d2lnl2[i*dim+j];
		}
	}
	
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			printf("%f ", d2lnl2[i*dim+j]);
		}
		printf("\n");
	}
	
	inverse2(hessian, dim);
	
	printf("hessian\n");
	printf("\n");
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			printf("%e ", hessian[i*dim+j]);
		}
		printf("\n");
	}
	
	printf("\n");
	for (int i = 0; i < dim; i++) printf("%s ", Parameters_name(var->parameters, i));
	printf("\n");
	for (int i = 0; i < dim; i++) {
		printf("%s ", Parameters_name(var->parameters, i));
		for (int j = 0; j < dim; j++) {
			printf("%e ", hessian[i*dim+j]/(sqrt(hessian[i*dim+i])*sqrt(hessian[j*dim+j])));
		}
		printf("\n");
	}
	
	cholesky(dim, hessian, diag);
	
	//  Transform elements in L matrix
	int d = dim;
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j <= i; j++) {
			Parameters_set_value(var->var_parameters, d++, exp(hessian[dim*j+i]));
		}
	}
	
	// Set elements between distant nodes to INFINITY
	//	make_sparse(nodes, dim, params);
	free(hessian);
	free(d2lnl2);
	free(diag);
}

double klqp_fullrank_normal_elbo(variational_t* var){
	if (var->initialized == false) {
		klqp_fullrank_normal_init(var);
		// save for later
		for (int i = 0; i < Parameters_count(var->parameters); i++) {
			Parameter_store(Parameters_at(var->parameters, i));
		}
		var->initialized = true;
	}
	
	Model* posterior = var->posterior;
	double elbo = 0;
	size_t dim = Parameters_count(var->parameters);
	double* eta = dvector(dim);
	
	// Entropy
	double entropy = 0;
	int r = dim-1;
	int count = 0;
	for (int d = 0; d < dim; ++d) {
		r += d+1;
//		if(Parameters_estimate(var->var_parameters, r)){
			double tmp = fabs(Parameters_value(var->var_parameters, r));
			if (tmp != 0.0) entropy += log(tmp);
			count++;
//		}
	}
	entropy +=  0.5 * (1.0 + LOG_TWO_PI) * count;
	
	for (int i = 0; i < var->elbo_samples; i++) {
		double jacobian = 0.0;
		for (int j = 0; j < dim; j++) {
//			if(Parameters_estimate(var->var_parameters, j)){
				eta[j] = rnorm();
//			}
		}
		
		size_t row = dim;
		for (int j = 0; j < dim; j++) {
			double temp = 0;
			// multiply L_j and eta
			for (int k = 0; k < j+1; k++) {
				if(Parameters_estimate(var->var_parameters, row)){
					temp += Parameters_value(var->var_parameters, row)*eta[k];
				}
				row++;
			}
			double zeta = temp + Parameters_value(var->var_parameters, j); // add mu
			Parameter* p = Parameters_at(var->parameters, j);
			double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), &jacobian);
			Parameter_set_value(p, theta);
		}
        
        for(int s = 0; s < var->simplex_count; s++){
            Simplex* simplex = var->simplices[s]->obj;
            const double* constrained_values = simplex->get_values(simplex);
            double stick = 1;
            for(int k = 0; k < simplex->K-1; k++){
                jacobian += log(stick);
                stick -= constrained_values[k];
            }
        }
		
		double logP = posterior->logP(posterior);
		elbo += logP + jacobian;
		
		if(isinf(elbo)){
			free(eta);
			fprintf(stderr, "elbo_fullrank elbo %f logP %f jacobian %f\n", elbo, logP, jacobian);
			return elbo;
		}
		
	}
	free(eta);
	return elbo/var->elbo_samples + entropy;
}

void klqp_fullrank_normal_grad_elbo(variational_t* var, const Parameters* parameters, double* grads){
	if (var->initialized == false) {
		klqp_fullrank_normal_init(var);
		var->initialized = true;
	}
	
	Model* posterior = var->posterior;
	size_t dim = Parameters_count(var->parameters);
	size_t grad_dim = 2*dim + (dim*dim-dim)/2;
	double* eta = dvector(dim);
	double* zeta = dvector(dim);
	memset(grads, 0, sizeof(double)*grad_dim);
	
    size_t simplex_parameter_count = 0;
    for(int s = 0; s < var->simplex_count; s++){
        Simplex* simplex = var->simplices[s]->obj;
        simplex_parameter_count += simplex->K - 1;
    }
    
	for (int i = 0; i < var->grad_samples; i++) {
		for (int j = 0; j < dim; j++) {
//			if(Parameters_estimate(var->var_parameters, j)){
				eta[j] = rnorm();
//			}
		}
		size_t row = dim;
		for (int j = 0; j < dim; j++) {
			double temp = 0;
			// multiply L_j and eta
			for (int k = 0; k < j+1; k++) {
//				if(Parameters_estimate(var->var_parameters, row)){
					temp += Parameters_value(var->var_parameters, row)*eta[k];
//				}
				row++;
			}
			zeta[j] = temp + Parameters_value(var->var_parameters, j); // add mu
			
			Parameter* p = Parameters_at(var->parameters, j);
			double theta = inverse_transform2(zeta[j], Parameter_lower(p), Parameter_upper(p));
			Parameter_set_value(p, theta);
		}
		
		row = dim;
		for (int k = 0; k < dim; k++){
			Parameter* p = Parameters_at(var->parameters, k);
			//			double dlogP = Model_first_derivative(posterior, p, 0.001);
			double dlogP = posterior->dlogP(posterior, p);
			
			const double gldit = grad_log_det_inverse_transform(zeta[k], Parameter_lower(p), Parameter_upper(p));
			double grad_mu = dlogP * grad_inverse_transform(zeta[k], Parameter_lower(p), Parameter_upper(p)) + gldit;
            if (k < simplex_parameter_count) {
                grad_mu += 1.0/zeta[k] + 1.0/(zeta[k]-1.0); // grad log det transform of stick
            }
			grads[k] += grad_mu;
			for (int j = 0; j < k+1; j++) {
//				if(Parameters_estimate(var->var_parameters, row)){
					grads[row] += grad_mu*eta[j];
//				}
				row++;
			}
		}
	}
	
	for (int i = 0; i < grad_dim; i++) {
		grads[i] /= var->grad_samples;
	}
	
	size_t diag = dim-1;
	for (int i = 0; i < dim; i++) {
		diag += i+1;
//		if(Parameters_estimate(var->var_parameters, i)){
			grads[diag] += 1.0/Parameters_value(var->var_parameters, diag);
//		}
	}
	free(eta);
	free(zeta);
	
	if (var->file != NULL) {
		fprintf(var->file, "%zu", var->iter++);
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%e", Parameters_value(var->var_parameters, i));
		}
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%e", grads[i]);
		}
		fprintf(var->file, "\n");
	}
}

double klqp_fullrank_normal_logP(variational_t* var, const double* values){
	size_t dim = Parameters_count(var->parameters);
	gsl_vector * mu = gsl_vector_calloc(dim);
	gsl_matrix * L = gsl_matrix_calloc(dim, dim);
	gsl_vector * x = gsl_vector_calloc(dim);
	gsl_vector * work = gsl_vector_calloc(dim);
	
	double logJac = 0;
	size_t row = dim;
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		double zeta = transform(values[i], Parameter_lower(p), Parameter_upper(p));
		gsl_vector_set(x, i, zeta);
		logJac += zeta;
		gsl_vector_set(mu, i, Parameters_value(var->var_parameters, i));
		for (int j = 0; j <= i; j++) {
			gsl_matrix_set(L, i, j, Parameters_value(var->var_parameters, row));
			row++;
		}
	}
	double logP = 0;
	gsl_ran_multivariate_gaussian_log_pdf (x,  mu, L, &logP, work);
	
	gsl_vector_free(mu);
	gsl_matrix_free(L);
	gsl_vector_free(x);
	gsl_vector_free(work);
	
	return logP + logJac;
}

bool klqp_fullrank_normal_sample(variational_t* var, double* values){
	if(!var->ready_to_sample){
		var->finalize(var);
		var->ready_to_sample = true;
	}
	size_t dim = Parameters_count(var->parameters);
	
	gsl_vector * mu = gsl_vector_calloc(dim);
	gsl_matrix * L = gsl_matrix_calloc(dim, dim);
	gsl_vector * samples = gsl_vector_calloc(dim);
	
	size_t row = dim;
	for (int i = 0; i < dim; i++) {
		gsl_vector_set(mu, i, Parameters_value(var->var_parameters, i));
		for (int j = 0; j <= i; j++) {
			gsl_matrix_set(L, i, j, Parameters_value(var->var_parameters, row));
			row++;
		}
	}
	
	gsl_ran_multivariate_gaussian(var->rng, mu, L, samples);
	
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		values[i] = inverse_transform2(gsl_vector_get(samples, i), Parameter_lower(p), Parameter_upper(p));
	}
	
	gsl_vector_free(mu);
	gsl_vector_free(samples);
	gsl_matrix_free(L);
	return true;
}

void klqp_fullrank_log_samples(variational_t* var, FILE* file){
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
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
	
	gsl_vector * mu = gsl_vector_calloc(dim);
	gsl_matrix * L = gsl_matrix_calloc(dim, dim);
	gsl_vector * gsl_samples = gsl_vector_calloc(dim);
	gsl_vector * work = gsl_vector_calloc(dim);
	double* samples = dvector(dim);
	
	size_t row = dim;
	for (int i = 0; i < dim; i++) {
		gsl_vector_set(mu, i, Parameters_value(var->var_parameters, i));
		for (int j = 0; j <= i; j++) {
			gsl_matrix_set(L, i, j, Parameters_value(var->var_parameters, row));
			row++;
		}
	}
	
	for (int k = 0; k < var->log_samples; k++) {
		gsl_ran_multivariate_gaussian(var->rng, mu, L, gsl_samples);
		
		double logQ = 0;
		gsl_ran_multivariate_gaussian_log_pdf (gsl_samples,  mu, L, &logQ, work);
		
		double jacobian = 0.0;
		for (int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			samples[j] = inverse_transform(gsl_vector_get(gsl_samples, j), Parameter_lower(p), Parameter_upper(p), &jacobian);
			Parameter_set_value(p, samples[j]);
		}
		
		double logP = posterior->logP(posterior) + jacobian;
		fprintf(file, "%d,%f,%f", k, logP, logQ);
		
		if (var->simplex_count > 0) {
			for(int s = 0; s < var->simplex_count; s++){
				Simplex* simplex = var->simplices[s]->obj;
				const double* constrained_values = simplex->get_values(simplex);
				for(int k = 0; k < simplex->K; k++){
					fprintf(file, ",%f",constrained_values[k]);
				}
			}
		}
		
		for (int j = offset; j < dim; j++) {
			fprintf(file, ",%f", samples[j]);
		}
		fprintf(file, "\n");
	}
	free(samples);
	gsl_vector_free(mu);
	gsl_vector_free(gsl_samples);
	gsl_matrix_free(L);
	gsl_vector_free(work);
}

