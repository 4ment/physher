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
#include <tgmath.h>

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

void klqp_meanfield_normal_init(variational_t* var){
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
	
	// for reuse
	for (int i = 0; i < Parameters_count(var->parameters); i++) {
		Parameter_store(Parameters_at(var->parameters, i));
	}
	
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		// fix mu when value is small
		// mu will be set such that mu-sigma*sigma = -10
		if(Parameter_value(p) < 1.0e-6){
			//			printf("%e %s\n", Parameter_value(p), Parameter_name(p));
			Parameter* var_p_mu = Parameters_at(var->var_parameters, i);
			Parameter* var_p_sigma = Parameters_at(var->var_parameters, i+dim);
			Parameter_set_estimate(var_p_mu, false);
		}
		else{
			Parameter* var_p_mu = Parameters_at(var->var_parameters, i);
			double dlogP;
			double d2logP = posterior->d2logP(posterior, Parameters_at(var->parameters, i)); //Model_second_derivative(posterior, Parameters_at(var->parameters, i), &dlogP, 0.001);
			double mu = Parameters_value(var->parameters, i); // mean = mode of normal
			double v = -1.0/d2logP; // variance of normal
			
			double q_var = log(sqrt(log(exp(log(v)-2.0*log(mu))+1.0)));
			double q_mu = log(mu) - q_var*0.5;
			
			// d2logp could be non negative are create a nan
			if (isnan(q_mu)) {
				q_mu = 1;
			}
			if (isnan(q_var)) {
				q_var = 0;
			}
			Parameters_set_value(var->var_parameters, i, q_mu);
			Parameters_set_value(var->var_parameters, i+dim, q_var);
			//			printf("%f %f  %f %f %s\n", mu, d2logP, q_mu, q_var, p->name);
			//			Parameter_set_estimate(var_p_mu, false);
		}
	}
	
	return;
	
	for (int i = 0; i < dim; i++) {
		//		printf("%s %e", Parameters_name(var->parameters, i), Parameters_value(var->parameters, i));
		//		printf("bl: %e", Parameters_value(var->parameters, i));
		
		double dlogP;
		double d2logP = Model_second_derivative(posterior, Parameters_at(var->parameters, i), &dlogP, 0.001);
		double mu = Parameters_value(var->parameters, i); // mean = mode of normal
		double v = -1.0/d2logP; // variance of normal
		
		double q_var = log(sqrt(log(exp(log(v)-2.0*log(mu))+1.0)));
		double q_mu = log(mu) - q_var*0.5;
		
		// e.g. small branch: derivative vanishes at 0
		
		if(isnan(q_var)){
			q_var = -5;
		}
		if(mu < 0.0001 && dlogP < -100){
			q_mu = log(mu) + 1;
			q_var = 1;
		}
		else{
			q_mu = log(mu) - q_var*0.5;
		}
		//		if(d2logP >= 0 || isnan(q_var)){
		//			q_var = 1;
		//			q_mu = log(mu)+1;
		//		}
		//		if(isinf(q_mu) || isnan(q_mu)){
		//			q_mu = 1;
		//		}
		//		if(isinf(q_var) || isnan(q_var)){
		//			q_var = 1;
		//		}
		
		Parameters_set_value(var->var_parameters, i, q_mu);
		Parameters_set_value(var->var_parameters, i+dim, q_var);
		//		printf(" m: %f d: %f qm: %f qv: %f v: %f d2: %f\n", mu, dlogP, q_mu, q_var, v, d2logP);
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
		//		double mu = Parameters_value(var->var_parameters, i);
		//		if(!Parameters_at(var->var_parameters, i)->estimate){
		//			mu = -10 + sigma;
		//		}
		////		double mean = exp(mu + sigma*sigma/2.0);
		////		double mode = exp(mu - sigma*sigma);
		////		printf("%f %f %f %f %s\n", Parameters_value(var->parameters, i), mode,
		////			   mu, sigma, Parameters_name(var->parameters, i));
	}
}


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
		
		double logP = posterior->logP(posterior);
		if(!isinf(logP))elbo += logP + jacobian;
		else inf_count++;
		
		if(isinf(elbo) || isnan(elbo)){
			return elbo;
		}
	}
	return elbo/(var->elbo_samples-inf_count) + entropy;
}

void klqp_meanfield_normal_grad_elbo(variational_t* var, double* grads){
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
	
	for (int i = 0; i < var->grad_samples; i++) {
		for ( int j = 0; j < dim; j++) {
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
		
		for ( int j = 0; j < dim; j++) {
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
double klqp_meanfield_normal_logP(variational_t* var, double* values){
	int dim = Parameters_count(var->parameters);
	double logP = 0;
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		double mu = Parameters_value(var->var_parameters, i);
		double sd = exp(Parameters_value(var->var_parameters, i+dim));
		double zeta = transform(values[i], Parameter_lower(p), Parameter_upper(p));
		logP += dnorml(zeta, mu, sd) - zeta;
	}
	return logP;
}

double klqp_meanfield_normal_logP_parameters(variational_t* var, const Parameters* parameters){
	size_t dim = Parameters_count(var->parameters);
	size_t dim2 = Parameters_count(parameters);
	double logP = 0;
	for (int i = 0; i < dim2; i++) {
		Parameter* p = Parameters_at(parameters, i);
		for (int j = 0; j < dim; j++) {
			Parameter* p2 = Parameters_at(var->parameters, j);
			if (p == p2) {
				double mu = Parameters_value(var->var_parameters, i);
				double sd = exp(Parameters_value(var->var_parameters, i+dim));
				double zeta = transform(Parameter_value(p), Parameter_lower(p), Parameter_upper(p));
				logP += dnorml(zeta, mu, sd) - zeta;
				break;
			}
		}
	}
	return logP;
}

// check success
bool klqp_meanfield_normal_sample(variational_t* var, double* values){
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
	size_t grad_dim = Parameters_count(var->var_parameters);
	
	double* eta = dvector(grad_dim);
	
	// Entropy
	double entropy = 0;
	int r = dim-1;
	int count = 0;
	for (int d = 0; d < dim; ++d) {
		r += d+1;
		if(Parameters_estimate(var->var_parameters, r)){
			double tmp = fabs(Parameters_value(var->var_parameters, r));
			if (tmp != 0.0) entropy += log(tmp);
			count++;
		}
	}
	entropy +=  0.5 * (1.0 + LOG_TWO_PI) * count;
	
	for (int i = 0; i < var->elbo_samples; i++) {
		double jacobian = 0.0;
		for (int j = 0; j < grad_dim; j++) {
			if(Parameters_estimate(var->var_parameters, j)){
				eta[j] = rnorm();
			}
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

void klqp_fullrank_normal_grad_elbo(variational_t* var, double* grads){
	if (var->initialized == false) {
		klqp_fullrank_normal_init(var);
		var->initialized = true;
	}
	
	Model* posterior = var->posterior;
	size_t dim = Parameters_count(var->parameters);
	size_t grad_dim = 2*dim + (dim*dim-dim)/2;
	double* eta = dvector(grad_dim);
	double* zeta = dvector(grad_dim);
	memset(grads, 0, sizeof(double)*grad_dim);
	
	for (int i = 0; i < var->grad_samples; i++) {
		for (int j = 0; j < grad_dim; j++) {
			if(Parameters_estimate(var->var_parameters, j)){
				eta[j] = rnorm();
			}
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
			
			grads[k] += grad_mu;
			for (int j = 0; j < k+1; j++) {
				if(Parameters_estimate(var->var_parameters, row)){
					grads[row] += grad_mu*eta[j];
				}
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
		if(Parameters_estimate(var->var_parameters, i)){
			grads[diag] += 1.0/Parameters_value(var->var_parameters, diag);
		}
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

double klqp_fullrank_normal_logP(variational_t* var, double* values){
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

double klqp_fullrank_normal_logP_parameters(variational_t* var, const Parameters* parameters){
	//TODO: need to rethink this function. Does not make sense in multivariate case
	error("variational_fullrank_parameters_logP\n");
	return 0;
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
		values[i] = gsl_vector_get(samples, i);
	}
	
	gsl_vector_free(mu);
	gsl_vector_free(samples);
	gsl_matrix_free(L);
	return true;
}

