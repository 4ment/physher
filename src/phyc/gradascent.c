//
//  gradascent.c
//  physher
//
//  Created by Mathieu Fourment on 30/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "gradascent.h"

#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <signal.h>

#include "matrix.h"
#include "checkpoint.h"


#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif

bool sgd_interrupted = false;

void sgd_signal_callback_handler( int signum ) {
	printf("Caught signal %d\n",signum);
	if ( sgd_interrupted == true ) {
		exit(SIGINT);
	}
	sgd_interrupted = true;
}

#if defined (PTHREAD_ENABLED)

typedef struct threadpool_sg_t{
	pthread_t *threads;
	pthread_mutex_t lock;
	Model** models; // variational models
	Parameters** parameters; // associated parameters
	size_t model_count; // number of models available
	double* etas;
	double* elbos;
	OptStopCriterion *stop;
	opt_func f;
	opt_grad_func grad_f;
	size_t index_model;
	size_t count; // read-write
	size_t total; // read-only
	void(*reset)(void*);
}threadpool_sg_t;

static void * _threadpool_sg_worker( void *threadpool ){
	threadpool_sg_t* pool = (threadpool_sg_t *)threadpool;
	
	size_t index_model = 0;
	pthread_mutex_lock(&(pool->lock));
	index_model = pool->index_model++;
	pthread_mutex_unlock(&(pool->lock));
	
	Model* model = pool->models[index_model];
	Parameters* parameters = pool->parameters[index_model];
	
	while ( 1 ) {
		size_t count;
		pthread_mutex_lock(&(pool->lock));
		if( pool->count == pool->total ){
			pthread_mutex_unlock(&(pool->lock));
			break;
		}
		count = pool->count++;
		pthread_mutex_unlock(&(pool->lock));
		
		// should model reset here
		OptStopCriterion stop = *pool->stop; // copy
		stop.iter_max = 50;
		stop.frequency_check = 1;
		
		for (int j = 0; j < Parameters_count(parameters); j++) {
			Parameter_restore(Parameters_at(parameters, j));
		}
		pool->reset(model);
		
		optimize_stochastic_gradient(parameters, pool->f, pool->grad_f, pool->etas[count], model, &stop, 0, pool->elbos+count, NULL);
	}
	pthread_exit(NULL);
	return NULL;
}
#endif

opt_result optimize_stochastic_gradient_adapt(Parameters* parameters, opt_func f, opt_grad_func grad_f, void(*reset)(void*),
											  double* etas, size_t eta_count, void *data,
											  OptStopCriterion *stop, int verbose, double *best_eta, size_t nthreads){
	opt_result result = OPT_SUCCESS;
	
	for (size_t i = 0; i < Parameters_count(parameters); i++) {
		Parameter_store(Parameters_at(parameters, i));
	}
	nthreads = dmin(nthreads, eta_count);
	
	if (nthreads == 1) {
		
		double* elbos = dvector(eta_count);
		for (size_t i = 0; i < eta_count; i++) {
			for (int j = 0; j < Parameters_count(parameters); j++) {
				Parameter_restore(Parameters_at(parameters, j));
			}
			reset(data);
			
			OptStopCriterion stop_local = *stop; // copy
			stop_local.iter_max = 50;
			stop_local.frequency_check = 0;
			fprintf(stdout, "Trying eta: %f ", etas[i]);
			optimize_stochastic_gradient(parameters, f, grad_f, etas[i], data, &stop_local, 0, elbos+i, NULL);
			fprintf(stdout, "ELBO: %f\n", elbos[i]);
		}
		
		for (size_t i = 0; i < Parameters_count(parameters); i++) {
			Parameter_restore(Parameters_at(parameters, i));
		}
		
		double elbo_max = -INFINITY;
		int elbo_max_index = -1;
		for (int i = 0; i < eta_count; i++) {
			if(!isnan(elbos[i]) && !isnan(elbos[i]) && elbos[i] > elbo_max){
				elbo_max = elbos[i];
				elbo_max_index = i;
			}
		}
		free(elbos);
	
		if (isinf(elbo_max)) {
			result = OPT_FAIL;
		}
		else{
			*best_eta = etas[elbo_max_index];
		}
		return result;
	}

#if defined (PTHREAD_ENABLED)
	Hashtable* hash = new_Hashtable_string(10);
	hashtable_set_key_ownership( hash, false );
	hashtable_set_value_ownership( hash, false );
	
	threadpool_sg_t* threadpool = malloc(sizeof(threadpool_sg_t));
	assert(threadpool);
	threadpool->threads = malloc(nthreads*sizeof(pthread_t));
	threadpool->index_model = 0;
	threadpool->count = 0;
	threadpool->total = eta_count;
	threadpool->models = malloc(nthreads*sizeof(Model*));
	threadpool->parameters = malloc(nthreads*sizeof(Model*));
	threadpool->model_count = nthreads;
	threadpool->etas = etas;
	threadpool->elbos = dvector(eta_count);
	threadpool->stop = stop;
	threadpool->f = f;
	threadpool->grad_f = grad_f;
	threadpool->reset = reset;
	
	Model* model = data;
	threadpool->models[0] = model;
	threadpool->parameters[0] = parameters;
	
	for (size_t i = 1; i < threadpool->model_count; i++) {
		threadpool->models[i] = model->clone(model, hash);
		threadpool->parameters[i] = Hashtable_get(hash, Parameters_name2(parameters)); // fish parameters out
		Hashtable_empty(hash);
	}
	free_Hashtable(hash);
	
	pthread_mutex_init(&(threadpool->lock), NULL);
	
	for ( size_t i = 0; i < nthreads; i++ ) {
		pthread_create( &(threadpool->threads[i]), NULL, _threadpool_sg_worker, threadpool );
	}
	for ( size_t i = 0; i < nthreads; i++ ) {
		pthread_join(threadpool->threads[i], NULL);
	}
	
	for (size_t i = 0; i < Parameters_count(parameters); i++) {
		Parameter_restore(Parameters_at(parameters, i));
	}
	
	double elbo_max = -INFINITY;
	int elbo_max_index = -1;
	for (int i = 0; i < eta_count; i++) {
		if(!isnan(threadpool->elbos[i]) && !isnan(threadpool->elbos[i]) && threadpool->elbos[i] > elbo_max){
			elbo_max = threadpool->elbos[i];
			elbo_max_index = i;
		}
	}
	
	if (isinf(elbo_max)) {
		result = OPT_FAIL;
	}
	else{
		*best_eta = etas[elbo_max_index];
	}
	
	free(threadpool->threads);
	pthread_mutex_destroy(&(threadpool->lock));
	
	for (size_t i = 1; i < threadpool->model_count; i++) {
		threadpool->models[i]->free(threadpool->models[i]);
	}
	printf("-------\n");
	printf("eta\telbo\n");
	for (int i = 0; i < eta_count; i++) {
		printf("%f\t%f\n", etas[i], threadpool->elbos[i]);
	}
	printf("-------\n");
	free(threadpool->models);
	free(threadpool->parameters);
	free(threadpool->elbos);
	free(threadpool);
#endif
	
	return OPT_SUCCESS;
}


int compare (const void * a, const void * b){
	return ( *(double*)a - *(double*)b );
}

opt_result optimize_stochastic_gradient_adam(Parameters* parameters, opt_func f, opt_grad_func grad_f, double eta, void *data, OptStopCriterion *stop, int verbose, double *fmin, OptimizerCheckpoint* checkpointer){
	size_t dim = Parameters_count(parameters);
	
	double beta1 = 0.9;
	double beta2 = 0.999;
	double* grads = dvector(dim);
	double* mean_grad = dvector(dim);
	double* var_grad = dvector(dim);
	double eps = 1.e-8;
	
	double elbo = 0;
	double elbo_prev = -INFINITY;
	double elbo_best = -INFINITY;
	double tol_rel_obj = stop->tolfx;
	size_t eval_elbo = stop->frequency_check;
	int max_conv = 3;
	int conv = 0;
	stop->iter = 0;
	double *elbos = calloc(stop->iter_max/eval_elbo, sizeof(double));
	
	opt_result result = OPT_SUCCESS;
	size_t max_failures = 10;
	size_t failure_count = 0;
	double elbo0 = f(parameters, NULL, data);
	if(verbose > 0)
		printf("%zu ELBO: %f\n", stop->iter, elbo0);
	
	while(stop->iter++ < stop->iter_max){
		grad_f(parameters, grads, data);
		double eta_scaled = eta / sqrt(stop->iter);
		double beta1p = pow(beta1, stop->iter);
		double beta2p = pow(beta2, stop->iter);

		int i = 0;
		for (; i < dim; i++) {
			if(isnan(grads[i])){
				break;
			}
		}
		if(i == dim){
			failure_count = 0;
		}
		else if(failure_count < max_failures){
			stop->iter--;
			failure_count++;
			continue;
		}

		for (int i = 0; i < dim; i++) {
			if (Parameters_estimate(parameters, i)) {
				double grad = grads[i];
				mean_grad[i] = beta1*mean_grad[i] + (1.0 - beta1)*grad;
				var_grad[i] = beta2*var_grad[i] + (1.0 - beta2)*grad*grad;
				double hat_mean_grad = mean_grad[i]/(1.0 - beta1p);
				double hat_var_grad = var_grad[i]/(1.0 - beta2p);
				double v = Parameters_value(parameters, i) + eta_scaled * hat_mean_grad/(sqrt(hat_var_grad) + eps);
				Parameters_set_value(parameters, i, v);
			}
		}
		
		if (stop->iter % eval_elbo == 0) {
			elbo_prev = elbo;
			elbo = f(parameters, NULL, data);
			if(verbose > 0)  printf("%zu ELBO: %f (%f)\n", stop->iter, elbo, elbo_prev);
			if (isnan(elbo) || isinf(elbo)) {
				result = OPT_FAIL;
				*fmin = elbo;
				break;
			}
			
			if (elbo > elbo_best){
				elbo_best = elbo;
				*fmin = elbo_best;
			}

			if (stop->iter % checkpointer->frequency == 0){
				checkpoint_save(checkpointer->file, parameters);
			}

			double delta_elbo = fabs((elbo_prev - elbo) / elbo);
			size_t eval_count = stop->iter/eval_elbo;
			elbos[eval_count-1] = delta_elbo;
			qsort (elbos, eval_count, sizeof(double), compare);
			size_t median = eval_count/2;
			if(elbos[median] < tol_rel_obj && conv == max_conv){
				if(verbose > 0)  printf("ELBO converged: %f < %f  %f  %zu %zu\n",elbos[median], tol_rel_obj, delta_elbo, median, stop->iter);
				*fmin = elbo;
				break;
			}
			else if(elbos[median] < tol_rel_obj){
				conv++;
			}
		}
	}
	free(grads);
	free(var_grad);
	free(var_grad);
	free(elbos);
	
	return result;
}

opt_result optimize_stochastic_gradient(Parameters* parameters, opt_func f, opt_grad_func grad_f, double eta, void *data, OptStopCriterion *stop, int verbose, double *fmin, OptimizerCheckpoint* checkpointer){
	size_t dim = Parameters_count(parameters);
	double tau = 1;
	double pre_factor  = 0.9;
	double post_factor = 0.1;
	double elbo = 0;
	double elbo_prev = -INFINITY;
	double elbo_best = -INFINITY;
	double tol_rel_obj = stop->tolfx;
	size_t eval_elbo = stop->frequency_check;
	int max_conv = 3;
	int conv = 0;
	stop->iter = 0;
	double *elbos = NULL;
	if(eval_elbo > 0){
		elbos = calloc(stop->iter_max/eval_elbo, sizeof(double));
	}
	double* grads = calloc(dim, sizeof(double));
	double *history_grad_squared = calloc(dim, sizeof(double));
	
	opt_result result = OPT_SUCCESS;
	
	if(eval_elbo > 0){
		double elbo0 = f(parameters, NULL, data);
		if(verbose > 0)  printf("%zu ELBO: %f\n", stop->iter, elbo0);
	}
	
	signal(SIGINT, sgd_signal_callback_handler);
	sgd_interrupted = false;
	
	
	while(stop->iter++ < stop->iter_max){
		grad_f(parameters, grads, data);
		if(isnan(grads[0])){
			result = OPT_FAIL;
			break;
		}
		
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
                double v = Parameters_value(parameters, i);
                if (Parameters_constraint(parameters, i) != NULL && Parameters_lower(parameters, i) == 0) {
                    v = exp(log(v) + eta_scaled * grads[i] / (tau + sqrt(history_grad_squared[i])));
                }
                else{
                    v += eta_scaled * grads[i] / (tau + sqrt(history_grad_squared[i]));
                }
				Parameters_set_value(parameters, i, v);
			}
		}
		
		if (eval_elbo != 0 && stop->iter % eval_elbo == 0) {
			//            for (int j = 0; j < dim; j++) {
			//                Parameter* p = Parameters_at(var->parameters, j);
			//                //Parameter_set_value(p, exp(parameters[i]));
			//                Parameter_set_value(p, exp(var->var_parameters[j]-var->var_parameters[j+dim]*var->var_parameters[j+dim])); // mode
			//            }
			
			elbo_prev = elbo;
			elbo = f(parameters, NULL, data);
			if(verbose > 0)  printf("%zu ELBO: %f (%f)\n", stop->iter, elbo, elbo_prev);
			if (isnan(elbo) || isinf(elbo)) {
				result = OPT_FAIL;
				*fmin = elbo;
				break;
			}
			//            printf("%d ELBO: %f (%f) logL: %f\n",iter, elbo, elbo_prev, var->posterior->logP(var->posterior));
			//            double cubo = cubo_meanfield(tlk, nodes, nodeCount, elbo_samples, parameters,r);
			//            printf("%d ELBO: %f (%f) cubo: %f\n",iter, elbo, elbo_prev, cubo);
			
			if (elbo > elbo_best){
				elbo_best = elbo;
				*fmin = elbo_best;
			}

			if (stop->iter % checkpointer->frequency == 0){
				checkpoint_save(checkpointer->file, parameters);
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
				if(verbose > 0)  printf("ELBO converged: %f < %f  %f  %zu %zu\n",elbos[median], tol_rel_obj, delta_elbo, median, stop->iter);
				*fmin = elbo;
				break;
			}
			else if(elbos[median] < tol_rel_obj){
				conv++;
			}
			if(sgd_interrupted){
				*fmin = elbo;
				break;
			}
		}
	}
	free(grads);
	if(elbos != NULL) free(elbos);
	free(history_grad_squared);
	if(eval_elbo == 0){
		*fmin = f(parameters, NULL, data);
	}
	
	return result;
}
