//
//  weibullvi.c
//  physher
//
//  Created by mathieu on 21/5/20.
//  Copyright © 2020 Mathieu Fourment. All rights reserved.
//

#include "weibullvi.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

//MARK: block weibull meanfield

double qweibull(double p, double scale, double shape){
	return scale*pow(-log(1.0 - p), 1.0/shape);
}

void klqp_block_meanfield_weibull_sample1(variational_block_t* var, double* jacobian){
    size_t dim = Parameters_count(var->parameters);
    double* etas = Vector_mutable_data(var->etas);
    
    for (int j = 0; j < dim; j++) {
        Parameter* p = Parameters_at(var->parameters, j);
        double scale = Parameters_value(var->var_parameters[0], j);
        double shape = Parameters_value(var->var_parameters[1], j);
        
        etas[j] = gsl_ran_flat(var->rng, 0, 1);
		double z = qweibull(etas[j], scale, shape);
		Parameter_set_value(p, z);
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

// Entropy: H(X) = M_EULER (1-1/shape) + ln(scale/shape) + 1
double klqp_block_meanfield_weibull_entropy(variational_block_t* var){
    double entropy = 0;
    size_t dim = Parameters_count(var->parameters);
    for (size_t i = 0; i < dim; i++) {
        double scale = Parameters_value(var->var_parameters[0], i);
		double shape = Parameters_value(var->var_parameters[1], i);
		entropy += M_EULER * (1.0 - 1.0/shape) + log(scale/shape) + 1.0;
    }
    return entropy;
}

// parameters are the variational parameters
void klqp_block_meanfield_weibull_sample2(variational_block_t* var, const Parameters* parameters){
    size_t dim = Parameters_count(var->parameters);
    size_t opt_param_dim = Parameters_count(parameters);
    
    // find first shape
    size_t scale_idx = 0;
    for ( ; scale_idx < opt_param_dim; scale_idx++) {
        if(Parameters_at(parameters, scale_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (scale_idx == opt_param_dim) return;
    
    double* etas = Vector_mutable_data(var->etas);
    
    for (int idx = 0; idx < dim; idx++) {
        Parameter* p = Parameters_at(var->parameters, idx);
        double scale = Parameters_value(var->var_parameters[0], idx);
        double shape = Parameters_value(var->var_parameters[1], idx);
        
        etas[idx] = gsl_ran_flat(var->rng, 0, 1);
		double z = qweibull(etas[idx], scale, shape);
        Parameter_set_value(p, z);
    }
}

void klqp_block_meanfield_weibull_grad_elbo(variational_block_t* var, const Parameters* parameters, double* grads){
    
    // find first shape and rate
    size_t scale_idx = 0;
    size_t shape_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
    for ( ; scale_idx < opt_param_dim; scale_idx++) {
        if(Parameters_at(parameters, scale_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (scale_idx == opt_param_dim) return;
    
    for ( ; shape_idx < opt_param_dim; shape_idx++) {
        if(Parameters_at(parameters, shape_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
    
    Model* posterior = var->posterior;
    size_t dim = Parameters_count(var->parameters);
    size_t simplex_parameter_count = var->simplex_parameter_count;
    const double* etas = Vector_data(var->etas);
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
    
			
	for (int k = 0; k < dim; k++){
		Parameter* p = Parameters_at(var->parameters, k);
		double dlogP =  posterior->dlogP(posterior, p);
		
		double scale = Parameters_value(var->var_parameters[0], k);
		double shape = Parameters_value(var->var_parameters[1], k);
		double z = Parameter_value(p);
		
		double dtransdscale = pow(-log(1.0 - etas[k]), 1.0/shape);
		double dtransdshape = -log(-log(1.0 - etas[k]))*z/(shape*shape);
		grads[k+scale_idx] += dlogP * dtransdscale;
		grads[k+shape_idx] += dlogP * dtransdshape;
		
		if(!var->use_entropy){
//			f(x) = shape/scale * (x/scale)**(shape-1) * exp(-x/scale)**shape
//			log(f) = log(shape) - log(scale) + (shape-1)*log(x/scale) - shape*x/scale
			grads[k+scale_idx] -= (shape*(z - scale)/(scale*scale)) * dtransdscale;
			grads[k+shape_idx] -= (-z/scale + log(z/scale) + 1.0/shape) * dtransdshape;
			
		}
	}
}


// Entropy: H(X) = M_EULER (1-1/shape) + ln(scale/shape) + 1
void klqp_block_meanfield_weibull_grad_entropy(variational_block_t* var, const Parameters* parameters, double* grads){
    // find first shape
    size_t scale_idx = 0;
	size_t shape_idx = 0;
    size_t opt_param_dim = Parameters_count(parameters);
	for ( ; scale_idx < opt_param_dim; scale_idx++) {
        if(Parameters_at(parameters, scale_idx) == Parameters_at(var->var_parameters[0], 0)){
            break;
        }
    }
    // not beeing optimized
    if (scale_idx == opt_param_dim) return;
	
	for ( ; shape_idx < opt_param_dim; shape_idx++) {
        if(Parameters_at(parameters, shape_idx) == Parameters_at(var->var_parameters[1], 0)){
            break;
        }
    }
	
    
    size_t dim = Parameters_count(var->parameters);
    for (int j = 0; j < dim; j++) {
        double scale = Parameters_value(var->var_parameters[0], j);
		double shape = Parameters_value(var->var_parameters[1], j);
		grads[j+scale_idx] += 1.0/scale;
		grads[j+shape_idx] += M_EULER/(shape*shape) - 1.0/shape;
    }
}

double klqp_block_meanfield_weibull_logP(variational_block_t* var, const double* values){
    size_t dim = Parameters_count(var->parameters);
    double logP = 0;
    for (size_t i = 0; i < dim; i++) {
        double scale = Parameters_value(var->var_parameters[0], i);
        double shape = Parameters_value(var->var_parameters[1], i);
		double z = scale*pow(-log(1.0 - values[i]), 1.0/shape);
        logP += log(gsl_ran_weibull_pdf(z, scale, shape));
    }
    return logP;
}


double klqp_block_meanfield_weibull_logQ(variational_block_t* var, const double* values){
	return klqp_block_meanfield_weibull_logP(var, values);
}


void klqp_block_meanfield_weibull_sample(variational_block_t* var, double* values){
    size_t dim = Parameters_count(var->parameters);
    for (size_t i = 0; i < dim; i++) {
        double shape = Parameters_value(var->var_parameters[0], i);
        double rate = Parameters_value(var->var_parameters[1], i);
        values[i] = gsl_ran_gamma(var->rng, shape, 1.0/rate);
    }
}

bool klqp_block_meanfield_weibull_sample_some(variational_block_t* var, const Parameters* parameters, double* values){
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
