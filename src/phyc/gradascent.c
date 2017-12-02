//
//  gradascent.c
//  physher
//
//  Created by Mathieu Fourment on 30/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "gradascent.h"


int compare (const void * a, const void * b){
	return ( *(double*)a - *(double*)b );
}

opt_result optimize_stochastic_gradient(Parameters* parameters, opt_func f, opt_grad_func grad_f, double eta, void *data, OptStopCriterion *stop, double *fmin){
	size_t dim = Parameters_count(parameters);
	double tau = 1;
	double pre_factor  = 0.9;
	double post_factor = 0.1;
	double elbo = 0;
	double elbo_prev = -INFINITY;
	double elbo_best = -INFINITY;
	double tol_rel_obj = 0.0001;
	int eval_elbo = 100;
	int max_conv = 3;
	int conv = 0;
	stop->iter = 0;
	double *elbos = calloc(stop->iter_max/eval_elbo, sizeof(double));
	double* grads = calloc(dim, sizeof(double));
	double *history_grad_squared = calloc(dim, sizeof(double));
	
	opt_result result = OPT_SUCCESS;

	while(stop->iter++ < stop->iter_max){
		grad_f(parameters, grads, data);
		
		double eta_scaled = eta / sqrt(stop->iter);
		
		// Update step-size
		if (stop->iter == 1) {
			for (int i = 0; i < dim; i++) {
				history_grad_squared[i] = grads[i]*grads[i];
			}
		} else {
			for (int i = 0; i < dim; i++) {
				history_grad_squared[i] = pre_factor * history_grad_squared[i] + post_factor * grads[i]*grads[i];
			}
		}
		// ascent
		for (int i = 0; i < dim; i++) {
			if (Parameters_estimate(parameters, i)) {
				double v = Parameters_value(parameters, i) + eta_scaled * grads[i] / (tau + sqrt(history_grad_squared[i]));
				Parameters_set_value(parameters, i, v);
			}
		}
		
		if (stop->iter % eval_elbo == 0) {
			//            for (int j = 0; j < dim; j++) {
			//                Parameter* p = Parameters_at(var->parameters, j);
			//                //Parameter_set_value(p, exp(parameters[i]));
			//                Parameter_set_value(p, exp(var->var_parameters[j]-var->var_parameters[j+dim]*var->var_parameters[j+dim])); // mode
			//            }
			
			elbo_prev = elbo;
			elbo = f(parameters, NULL, data);
			printf("%zu ELBO: %f (%f)\n", stop->iter, elbo, elbo_prev);
			if (isnan(elbo)) {
				result = OPT_FAIL;
				break;
			}
			//            printf("%d ELBO: %f (%f) logL: %f\n",iter, elbo, elbo_prev, var->posterior->logP(var->posterior));
			//            double cubo = cubo_meanfield(tlk, nodes, nodeCount, elbo_samples, parameters,r);
			//            printf("%d ELBO: %f (%f) cubo: %f\n",iter, elbo, elbo_prev, cubo);
			
			if (elbo > elbo_best){
				elbo_best = elbo;
			}
			double delta_elbo = fabs((elbo_prev - elbo) / elbo);
			size_t eval_count = stop->iter/eval_elbo;
			elbos[eval_count-1] = delta_elbo;
			qsort (elbos, eval_count, sizeof(double), compare);
			size_t median = eval_count/2;
			if(elbos[median] < tol_rel_obj && conv == max_conv){
				//                for(int i = 0; i < parameterCount/2; i++){
				//                    printf("%d grad: %f %f\n",i, grads[i], grads[dim+i]);
				//                }
				printf("ELBO converged: %f < %f  %f  %zu %zu\n",elbos[median], tol_rel_obj, delta_elbo, median, stop->iter);
				break;
			}
			else if(elbos[median] < tol_rel_obj){
				conv++;
			}
		}
	}
	free(grads);
	free(elbos);
	free(history_grad_squared);
	
	return result;
}