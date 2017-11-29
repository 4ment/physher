//
//  transforms.c
//  viphy
//
//  Created by Mathieu Fourment on 4/05/2017.
//  Copyright © 2017 University of Technology Sydney. All rights reserved.
//

#include "transforms.h"

#include <math.h>
#include "matrix.h"
#include "hessenberg.h"

double log1p_exp(double a) {
	// like log_sum_exp below with b=0.0; prevents underflow
	if (a > 0.0)
		return a + log1p(exp(-a));
	return log1p(exp(a));
}

double softplus(double x){
    return log(1.0+exp(x));
}

// logit
double logit(const double x){
    return log(x/(1.0-x));
}

// inverse logit
double inverse_logit(const double y){
    return 1.0/(1+exp(-y));
}

double grad_inverse_logit(const double y){
    const double exp_y = exp(y);
    return exp_y/pow(exp_y+1.0, 2);
}

double grad_log_det_inverse_logit(double y){
    const double exp_y = exp(y);
    return (1.0-exp_y)/(exp_y+1.0);
}

double transform(double x, double lb, double ub){
    // lower bound
    if(!isinf(lb) && isinf(ub)){
        return log(x-lb);
    }
    // upper bound
    else if(!isinf(ub) && isinf(lb)){
        return log(ub-x);
    }
	else if(!isinf(ub) && !isinf(lb)){
		return logit(x-lb)/(ub-lb);
	}
    return x;
}

double inverse_transform(double z, double lb, double ub, double* lp){
    // lower bound
    if(!isinf(lb) && isinf(ub)){
        *lp += z;
        return exp(z) + lb;
    }
    // upper bound
    else if(!isinf(ub) && isinf(lb)){
        *lp += z;
        return ub - exp(z);
    }
    // no transform
    else if(isinf(ub) && isinf(lb)){
        return z;
    }
    
    double inv_logit_x;
    if (z > 0) {
        double exp_minus_x = exp(-z);
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        //printf("%f %f\n", inv_logit_x, x);
        *lp += log(ub - lb) - z - 2 * log1p(exp_minus_x);
        // Prevent x from reaching one unless it really really should.
        if (!isinf(z) && inv_logit_x == 1)
            inv_logit_x = 1 - 1e-15;
    } else {
        double exp_x = exp(z);
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        *lp += log(ub - lb) + z - 2 * log1p(exp_x);
        // Prevent x from reaching zero unless it really really should.
        if (!isinf(z) && inv_logit_x== 0)
            inv_logit_x = 1e-15;
    }
    return lb + (ub - lb) * inv_logit_x;
    
}

double inverse_transform2(double z, double lb, double ub){
    // lower bound
    if(!isinf(lb) && isinf(ub)){
        return exp(z) + lb;
    }
    // upper bound
    else if(!isinf(ub) && isinf(lb)){
        return ub - exp(z);
    }
    // no transform
    else if(isinf(ub) && isinf(lb)){
        return z;
    }
    
    double inv_logit_x;
    if (z > 0) {
        double exp_minus_x = exp(-z);
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        // Prevent x from reaching one unless it really really should.
        if (!isinf(z) && inv_logit_x == 1)
            inv_logit_x = 1 - 1e-15;
    } else {
        double exp_x = exp(z);
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        //printf("%f %f\n", inv_logit_x, x);
        // Prevent x from reaching zero unless it really really should.
        if (!isinf(z) && inv_logit_x== 0)
            inv_logit_x = 1e-15;
    }
    return lb + (ub - lb) * inv_logit_x;
    
}

double grad_inverse_transform(double z, double lb, double ub){
    // lower bound
    if(!isinf(lb) && isinf(ub)){
        return exp(z);
    }
    // upper bound
    else if(!isinf(ub) && isinf(lb)){
        return -exp(z);
    }
    // no transform
    else if(isinf(ub) && isinf(lb)){
        return 1.0;
    }
    
    double inv_logit_x;
    if (z > 0) {
        double exp_minus_x = exp(-z);
        inv_logit_x = exp_minus_x / pow((1.0 + exp_minus_x), 2);
    } else {
        double exp_x = exp(z);
        inv_logit_x = exp_x / pow((1.0 + exp_x), 2);
    }
    return (ub - lb) * inv_logit_x;
    
}

// Calculate gradient of log determinant of Jacobian
// the det should always be positive because all the terms are positive so no sign(det) term
double grad_log_det_inverse_transform(double z, double lb, double ub){
    if(lb == 0.0 && isinf(ub)){
        return 1.0;
    }
    // no transform
    else if(isinf(ub) && isinf(lb)){
        return 0;
    }
    
    double exp_x = exp(z);
    return (1.0-exp_x)/(1.0+exp_x);
}

// Simplex

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_simplex2(double* x, const double* y, int K){
	double stick = 1.0;
	int k = 0;
	for(; k < K-1; k++){
		double z = inverse_logit(y[k] - log(K-k-1));
		x[k] = stick * z;
		stick -= x[k];
	}
	x[k] = stick;
}

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_simplex(double* x, const double* y, int K, double* lp){
	double stick = 1.0;
	int k = 0;
	for(; k < K-1; k++){
		double adj_y_k = y[k] - log(K-k-1);
		double z = inverse_logit(adj_y_k);
		x[k] = stick * z;
		*lp += log(stick);
		*lp -= log1p_exp(-adj_y_k);
		*lp -= log1p_exp(adj_y_k);
		stick -= x[k];
	}
	x[k] = stick;
}

void grad_inverse_transform_simplex(double* dz, const double* z, int K){
    
}

/*void grad_log_det_inverse_transform_simplex(const double* x, double* dz, int K){
	//dx[0] = 1.0/(x[0]-1.0) + 1.0/(x[0]+x[1]-1.0);
	//dx[1] = 1.0/(x[0]+x[1]-1.0);
	dz[K-1] = 0;
	for(int k = K-2; k >=0; --k){
		dz[k] = -1.0;
		for(int i = 0; i <= k; --i){
			dz[k] += x[k];
		}
		dz[k] = 1.0/dz[k] + dz[k+1];
	}
}*/

void grad_log_det_inverse_transform_simplex(const double* z, double* dz, int K){
	for(int k = 0; k < K-1; k++){
		dz[k] = 1.0/z[k] + 1.0/(z[k]-1.0);
	}
}

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_msimplex(double* x, const double* y, int K){
	double sum = 0;
	int k = 0;
	for(; k < K-1; k++){
		x[k] = exp(y[k]);
		sum += x[k];
	}
	x[k] = 1.0/(1.0+sum);
	for(int i = 0; i < K-1; i++){
		x[i] /= x[k];
	}
}

// x: array of containing K elements
// y: array containing K-1 elements with support (-Inf, Inf)
// K: dimension of K-simplex
void inverse_transform_msimplex2(double* x, const double* y, int K, double* lp){
	double sum = 0;
	double sum_one = 1.0;
	int dim = K-1;
	int k = 0;
	for(; k < dim; k++){
		x[k] = exp(y[k]);
		sum_one += x[k];
	}
	x[k] = 1.0/(sum_one);
	for(int i = 0; i < dim; i++){
		x[i] /= x[k];
	}
	
	double inv_sum_one_squared = 1.0/(sum_one * sum_one);
	double* m = dvector(dim*dim);
	k = 0;
	for(int i = 0; i < dim; i++){
		double temp = 1;
		int j = 0;
		for(; j < i; j++) temp += x[j];
		j++;
		for(; j < dim; j++)	temp += x[j];
		m[dim*i+i] = x[i] * temp/inv_sum_one_squared;
		
		for(int j = i+1; j < dim; j++){
			m[dim*i+j] = -exp(y[i]+y[j])/ inv_sum_one_squared;
			m[dim*j+i] = m[dim*i+j];
		}
	}
	*lp = log(LUDecompose_det2(m, dim));
	free(m);
}

