/*
 *  optimizer.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/4/11.
 *  Copyright (C) 2016 Mathieu Fourment. All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with this program; if not,
 *  write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "optimizer.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <strings.h>

#include "parameters.h"
#include "brent.h"
#include "matrix.h"
#include "powell.h"
#include "bfgs.h"
#include "frpmrn.h"
#include "gradascent.h"
#include "topologyopt.h"
#include "tree.h"
#include "treelikelihood.h"

#include "vb.h"

static double _logP( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	double logP = model->logP(model);
	if (grad != NULL) {
		for (int i = 0; i < Parameters_count(params); i++) {
			grad[i] = model->dlogP(model, Parameters_at(params, i));
		}
	}
	//	printf("%f\n", logP);
	return logP;
}

double model_negative_logP( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
	double logP = model->logP(model);
//		printf("%f\n", model->logP(model));
	if (grad != NULL) {
		for (int i = 0; i < Parameters_count(params); i++) {
			grad[i] = -model->dlogP(model, Parameters_at(params, i));
		}
	}
	return -logP;
}

static void _gradient( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;
//	printf("%f\n", model->logP(model));
	model->gradient(model, grad);
}

static void _reset(void* data){
	Model* model = (Model*)data;
	model->reset(model);
}

static bool dummy_update_data( void *data, Parameters *p){return false;}

static bool xStop( const Parameters *x,  double *xold, const double tolx);
static bool fxStop(double fx, double *fxold, const double tolfx);

struct _Optimizer{
	opt_algorithm algorithm;
	unsigned int dimension;
	
	opt_func f;
    opt_grad_func grad_f;
	void (*reset)(void*);
	void *data;
//    Model* model;
	Parameters* parameters;
	Model* treelikelihood;

	OptStopCriterion stop;
	
	opt_update_data update;
	int verbosity;
	OptimizerSchedule* schedule;
    
    // for stochastic gradient
    double* etas;
	size_t eta_count;
    bool ascent;
	size_t threads;
};

opt_result topology_optimize(TopologyOptimizer* topopt, double *fmin){
	*fmin = -topopt->optimize(topopt);
	return OPT_SUCCESS;
}

opt_result serial_brent_optimize_tree( Model* mtlk, opt_func f, void *data, OptStopCriterion *stop, double *fmin ){
	Parameters* temp = new_Parameters(1);
	SingleTreeLikelihood* tlk = (SingleTreeLikelihood*)mtlk->obj;
	Tree* tree = tlk->tree;
	Node** nodes = Tree_get_nodes(tree, POSTORDER);
	tlk->node_upper = NULL;
	tlk->use_upper = true;
	tlk->update_upper = true;
	if(Node_distance(Tree_root(tree)->right) != 0){
		double tot = Node_distance(Tree_root(tree)->right) + Node_distance(Tree_root(tree)->left);
		Node_set_distance(Tree_root(tree)->right, 0);
		Node_set_distance(Tree_root(tree)->left, tot);
	}
	// initialize lower and upper
	SingleTreeLikelihood_update_uppers(tlk);
	//for(int j = 0; j < stop->iter_min; j++)
	for(int i = 0; i < Tree_node_count(tree); i++){
		Node* node = nodes[i];//Tree_node(tree, i);
		if(Node_isroot(node) || (Node_isroot(Node_parent(node)) && Node_right(Node_parent(node)) == node)) continue;
		if(tlk->node_upper == NULL) tlk->node_upper = node;

		Parameters_add(temp, node->distance);
		stop->iter = 0;
		stop->f_eval_current = 0;
		stop->count = 0;
		if ( stop->time_max != 0 ) {
			time( &stop->time_start );
		}
#ifdef UPPER_PARTIALS
		printf("brent\n");
#endif
		opt_result status = brent_optimize(temp, f, data, stop, fmin);
		Parameters_pop(temp);
#ifdef UPPER_PARTIALS
		printf("%f %s\n", -*fmin, nodes[i]->name);
#endif
		//tlk->node_upper = nodes[i];
	}
	free_Parameters(temp);
	tlk->use_upper = false;
	return OPT_SUCCESS;
}

opt_result meta_optimize( opt_func f, void *data, OptStopCriterion *stop, double *fmin, OptimizerSchedule* schedule ){
	double lnl = f(NULL, NULL, data);
	double fret = lnl;
	for (stop->iter = 0; stop->iter < stop->iter_max; stop->iter++) {
		double lnl_current = lnl;
		for (int i = 0; i < schedule->count; i++) {
			Optimizer* opt = schedule->optimizers[i];
			double local_fret;
			for (int k = 0; k < schedule->rounds[i]; k++){
				local_fret = fret;
				opt_result status;
				if(opt->treelikelihood != NULL){
					status = serial_brent_optimize_tree(opt->treelikelihood, opt->f, opt->data, &opt->stop, &fret);
				}
				else{
					status = opt_optimize( opt, opt->parameters, &fret);
				}
				bool stopit = schedule->post[i](schedule,local_fret, fret);
				//				printf("%s %f %f -> %f (%f)\n", Parameters_name(parameters, 0), Parameters_value(parameters, 0), lnl, fret, local_fret);
				if(stopit) break;
			}
			if(stop->iter == 2) schedule->rounds[i] = 1;
			// Optimizing branches efficiently
			if(opt->treelikelihood != NULL){
//				SingleTreeLikelihood_update_all_nodes(opt->treelikelihood->obj);
//				fret = _logP(NULL, NULL, data);
				//printf("-== %f\n", _logP(NULL, NULL, data));
				printf("branches %f %f\n", -fret, -lnl);
			}
			// Optimizing any parameters
			else if(opt->parameters != NULL){
				printf("%s %f %f\n", Parameters_name(opt->parameters, 0), -fret, -lnl);
			}
			// Optimizing topology
			else{
				TopologyOptimizer* topopt = opt->data;
				printf("topology %f %f (%d)\n", -fret, -lnl, topopt->moves);
			}
			lnl = fret;
		}
		if(  lnl_current-lnl < stop->tolfx ){
			printf("\n%f %f\n\n", lnl, lnl_current);
			*fmin = lnl;
			return OPT_SUCCESS;
		}
		printf("\n");
	}
	return OPT_MAXITER;
}



Optimizer * new_Optimizer( opt_algorithm algorithm ) {
	Optimizer * opt;
	opt = (Optimizer*) malloc(sizeof(struct _Optimizer));
	assert(opt);
	opt->algorithm = algorithm;
	opt->f = NULL;
    opt->grad_f = NULL;
	opt->data = NULL;
	opt->parameters = NULL;
	opt->treelikelihood = NULL;
	opt->dimension = 0;
	
	opt->stop.iter_min = 1;
	opt->stop.iter_max = 200;
	opt->stop.iter = 0;

	opt->stop.time_max = 0;
	opt->stop.time_start = 0;
	opt->stop.time_end = 0;
	opt->stop.time_current = 0;
	
	opt->stop.f_eval_max = 0;
	opt->stop.f_eval_current = 0;
	
	opt->stop.tolfx = 0;
	opt->stop.tolx = OPT_XTOL;
	opt->stop.tolg = 1.e-5;
	
	opt->stop.oldx = NULL;
	opt->stop.oldfx = 0;
	
	opt->stop.count = 0;
	opt->verbosity = 1;
	opt->update = dummy_update_data;
	opt->schedule = NULL;
    
	opt->stop.frequency_check = 100;
    opt->etas = NULL;
	opt->eta_count = 0;
    opt->ascent = true;
	opt->reset = NULL;
	return opt;
}

// need to change:
// data treelikelihood
// f: maybe different log likleihood for tree with upper
Optimizer* clone_Optimizer(Optimizer *opt, void* data, Parameters* parameters){
	Optimizer* clone = (Optimizer*) malloc(sizeof(struct _Optimizer));
	assert(clone);
	clone->algorithm = opt->algorithm;
	clone->f = opt->f;
	clone->data = data;
	clone->dimension = opt->dimension;
	
	clone->stop.iter_min = opt->stop.iter_min;
	clone->stop.iter_max = opt->stop.iter_max;
	clone->stop.iter = opt->stop.iter;
	clone->stop.frequency_check = opt->stop.frequency_check;
	
	clone->stop.time_max = opt->stop.time_max;
	clone->stop.time_start = opt->stop.time_start;
	clone->stop.time_end = opt->stop.time_end;
	clone->stop.time_current = opt->stop.time_current;
	
	clone->stop.f_eval_max = opt->stop.f_eval_max;
	clone->stop.f_eval_current = opt->stop.f_eval_current;
	
	clone->stop.tolfx = opt->stop.tolfx;
	clone->stop.tolx = opt->stop.tolx;
	clone->stop.tolg = opt->stop.tolg;
	
	clone->stop.oldx = opt->stop.oldx;
	clone->stop.oldfx = opt->stop.oldfx;
	
	clone->stop.count = opt->stop.count;
	clone->verbosity = opt->verbosity;
	clone->update = opt->update;
	clone->schedule = NULL;
	clone->parameters = NULL;
	clone->treelikelihood = NULL;
    
    clone->ascent = opt->ascent;
	clone->etas = clone_dvector(opt->etas, opt->eta_count);
	clone->eta_count = opt->eta_count;
    clone->grad_f = opt->grad_f;
	
	if(opt->schedule != NULL){
		clone->schedule = (OptimizerSchedule*)malloc(sizeof(OptimizerSchedule));
		clone->schedule->capacity = opt->schedule->capacity;
		clone->schedule->count = opt->schedule->count;
		clone->schedule->optimizers = (Optimizer**)calloc(opt->schedule->capacity, sizeof(Optimizer*));
		clone->schedule->post = (OptimizerSchedule_post*)calloc(opt->schedule->capacity, sizeof(OptimizerSchedule_post));
		clone->schedule->rounds = ivector(opt->schedule->capacity);
		memcpy(clone->schedule->rounds, opt->schedule->rounds, sizeof(int)*opt->schedule->count);
		
		for (int i = 0; i < clone->schedule->count; i++) {
			clone->schedule->optimizers[i] = clone_Optimizer(opt->schedule->optimizers[i], data, parameters);
			clone->schedule->post[i] = opt->schedule->post[i];
			
			if(Parameters_count(opt->parameters) > 0){
				opt->parameters = new_Parameters(Parameters_count(opt->parameters));
			}
			
			// Find matching parameters
			for (int j = 0; j < Parameters_count(opt->parameters); j++) {
				int k = 0;
				for (k = 0; k < Parameters_count(parameters); k++) {
					if(strcmp(Parameters_name(opt->parameters, j), Parameters_name(parameters, k)) == 0){
						Parameters_add(clone->parameters, Parameters_at(parameters, k));
						break;
					}
				}
				if(k == Parameters_count(parameters)){
					printf("not found %s\n", Parameters_name(parameters, k));
					for (k = 0; k < Parameters_count(parameters); k++) {
						printf("%s\n", Parameters_name(parameters, k));
					}
					exit(1);
				}
			}
		}
	}
	clone->threads = opt->threads;
	return clone;
}

void free_Optimizer( Optimizer *opt ){
	if(opt->algorithm == OPT_TOPOLOGY){
		TopologyOptimizer* topopt = opt->data;
		free_TopologyOptimizer(topopt);
	}
	opt->data = NULL;
	free_Parameters(opt->parameters);
	if(opt->etas != NULL){
		free(opt->etas);
	}
	//if(opt->treelikelihood != NULL) opt->treelikelihood->free(opt->treelikelihood);
	if(opt->schedule != NULL){
		for (int i = 0; i < opt->schedule->count; i++) {
			free_Optimizer(opt->schedule->optimizers[i]);
		}
		free(opt->schedule->optimizers);
		free(opt->schedule->rounds);
		free(opt->schedule->post);
		free(opt->schedule);
	}
	free(opt);
}

void opt_set_objective_function( Optimizer *opt, opt_func f ){
	if( opt != NULL ){
		opt->f = f;
	}
}

void opt_set_data( Optimizer *opt, void *data ){
	if( opt != NULL ){
		opt->data = data;
	}
}

void opt_set_parameters( Optimizer *opt, const Parameters *parameters ){
	if( opt != NULL ){
		if(opt->parameters == NULL){
			opt->parameters = new_Parameters(Parameters_count(parameters));
			if(Parameters_name2(parameters) != NULL){
				Parameters_set_name2(opt->parameters, Parameters_name2(parameters));
			}
		}
		Parameters_add_parameters(opt->parameters, parameters);
	}
}

void opt_set_treelikelihood( Optimizer *opt, Model* treelikelihood){
	opt->treelikelihood = treelikelihood;
}

void opt_set_max_evaluation( Optimizer *opt, const size_t maxeval ){
	if( opt != NULL ){
		opt->stop.f_eval_max = maxeval;
	}
}

void opt_set_max_iteration( Optimizer *opt, const size_t maxiter ){
	if( opt != NULL ){
		opt->stop.iter_max = maxiter;
	}
}

void opt_set_min_iteration( Optimizer *opt, const size_t miniter ){
	if( opt != NULL ){
		opt->stop.iter_min = miniter;
	}
}

// In seconds by default
void opt_set_time_max( Optimizer *opt, const double maxtime ){
	if( opt != NULL ){
		opt->stop.time_max = maxtime;
	}
}

void opt_set_time_max_minutes( Optimizer *opt, const double maxtime ){
	if( opt != NULL ){
		opt->stop.time_max = maxtime*60;
	}
}

void opt_set_time_max_hours( Optimizer *opt, const double maxtime ){
	if( opt != NULL ){
		opt->stop.time_max = maxtime*3600;
	}
}

void opt_set_time_max_relative( Optimizer *opt, const time_t time, const double factor ){
	if( opt != NULL ){
		opt->stop.time_max = time*factor;
	}
}

void opt_set_tolfx( Optimizer *opt, const double tolfx ){
	if( opt != NULL ){
		opt->stop.tolfx = tolfx;
	}
}

void opt_set_tolx( Optimizer *opt, const double tolx ){
	if( opt != NULL ){
		opt->stop.tolx = tolx;
	}
}

double opt_tolx( Optimizer *opt ){
    return opt->stop.tolx;
}
void opt_set_tolg( Optimizer *opt, const double tolg ){
	if( opt != NULL ){
		opt->stop.tolg = tolg;
	}
}

size_t opt_frequency_check( Optimizer *opt ){
	return opt->stop.frequency_check;
}
void opt_set_frequency_check( Optimizer *opt, const size_t frequency_check ){
	if( opt != NULL ){
		opt->stop.frequency_check = frequency_check;
	}
}

int opt_iterations( Optimizer *opt ){
    return opt->stop.iter;
}

int opt_f_evaluations( Optimizer *opt ){
    return opt->stop.f_eval_current;
}

// Meta optimizer

bool opt_post(OptimizerSchedule* schedule, double before, double after){
	return false;
}

void opt_add_optimizer(Optimizer *opt_meta, Optimizer *opt){
	if(opt_meta->schedule == NULL){
		opt_get_schedule(opt_meta);
	}
	else if(opt_meta->schedule->capacity == opt_meta->schedule->count){
		opt_meta->schedule->capacity++;
		opt_meta->schedule->optimizers = (Optimizer**)realloc(opt_meta->schedule->optimizers, sizeof(Optimizer*)*opt_meta->schedule->capacity);
		opt_meta->schedule->rounds = (int*)realloc(opt_meta->schedule->rounds, sizeof(int)*opt_meta->schedule->capacity);
		opt_meta->schedule->post = (OptimizerSchedule_post*)realloc(opt_meta->schedule->post, sizeof(OptimizerSchedule_post)*opt_meta->schedule->capacity);
	}
	opt_meta->schedule->optimizers[opt_meta->schedule->count] = opt;
	opt_meta->schedule->rounds[opt_meta->schedule->count] = 1;
	opt_meta->schedule->post[opt_meta->schedule->count] = opt_post;
	opt_meta->schedule->count++;
}

OptimizerSchedule* opt_get_schedule(Optimizer *opt_meta){
	if(opt_meta->schedule == NULL){
		opt_meta->schedule = (OptimizerSchedule*)malloc(sizeof(OptimizerSchedule));
		opt_meta->schedule->optimizers = (Optimizer**)malloc(sizeof(Optimizer*));
		opt_meta->schedule->post = (OptimizerSchedule_post*)malloc(sizeof(OptimizerSchedule_post));
		opt_meta->schedule->rounds = ivector(1);
		opt_meta->schedule->capacity = 1;
		opt_meta->schedule->count = 0;
	}
	return opt_meta->schedule;
}


// At least one of these conditions is sufficient to stop the optimization
opt_result opt_check_stop( OptStopCriterion *stop, Parameters *x, double fx ){
	opt_result stopflag = OPT_KEEP_GOING;
	
	if ( stop->count == 0 ) {
		for (int i = 0; i < Parameters_count(x); i++) {
			stop->oldx[i] = Parameters_value(x,i);
		}
		stop->oldfx = fx;
		stop->count++;
		return stopflag;
	}
	
	if( stop->time_max != 0 ){
		time( &stop->time_current );
		bool stop_time = ( difftime( stop->time_current, stop->time_start ) > stop->time_max );
		if( stop_time ){
			stopflag = OPT_MAXTIME;
			//fprintf(stderr, "time max\n");
		}
		//if(stop_time) fprintf(stderr, "Max time reached %f (max=%f)\n", difftime( stop->time_current, stop->time_start), stop->time_max );
		
	}
	
	if ( stop->iter_max != 0 ) {
		bool stop_eval = ( stop->iter > stop->iter_max );
		if( stop_eval ){//fprintf(stderr, "iter max\n");
			stopflag = OPT_MAXITER;
			//fprintf(stderr, "Max iteration reached %d (max=%d)\n", stop->iter, stop->iter_max );
		}
		
	}
	if ( stop->f_eval_max != 0 ) {
		bool stop_eval = ( stop->f_eval_current+1 > stop->f_eval_max );
		if( stop_eval ){//fprintf(stderr, "f eval max\n");
			stopflag = OPT_MAXEVAL;
			//fprintf(stderr, "Max eval reached %d (max=%d)\n", stop->f_eval_current, stop->f_eval_max );
		}
		
	}
	
//	if( xStop( x, stop->oldx, stop->tolx) ){
//		fprintf(stderr, "tolx criterion: %f\n",stop->tolx);
//		return OPT_SUCCESS;
//	}
    
	if( fxStop(fx, &stop->oldfx, stop->tolfx) ){
		return OPT_SUCCESS;
	}
	
	return stopflag;
}

opt_result opt_optimize( Optimizer *opt, Parameters *ps, double *fmin ){
	if ( opt->stop.time_max != 0 ) {
		time( &opt->stop.time_start );
	}
	opt_result result = OPT_SUCCESS;
	
	if(ps != NULL){
		opt->dimension = Parameters_count(ps);
	}
	
	// should probably moved somewhere else
	opt->stop.oldx = dvector(opt->dimension);
    opt->stop.iter = 0;
    opt->stop.f_eval_current = 0;
    opt->stop.count = 0;
	
	switch (opt->algorithm ) {
		case OPT_META:{
			result = meta_optimize( opt->f, opt->data, &opt->stop, fmin, opt->schedule );
			break;
		}
		case OPT_POWELL:{
			result = powell_optimize( ps, opt->f, opt->data, opt->stop, fmin, opt->update );
			break;
		}
		case OPT_BRENT:{
			result = brent_optimize( ps, opt->f, opt->data, &opt->stop, fmin );
			break;
		}
		case OPT_SERIAL_BRENT:{
			if(opt->treelikelihood == NULL){
				result = serial_brent_optimize( ps, opt->f, opt->data, &opt->stop, fmin );
			}
			else{
				result = serial_brent_optimize_tree(opt->treelikelihood, opt->f, opt->data, &opt->stop, fmin);
			}
			break;
		}
		case OPT_BFGS:{
			result = dfpmin_optimize( ps, opt->f, opt->data, opt->stop, fmin );
			break;
		}
		case OPT_CG_FR:{
			result = frprmn_optimize( ps, opt->f, opt->data, opt->stop, fmin, OPT_CG_FR );
			break;
        }
        case OPT_CG_PR:{
            result = frprmn_optimize( ps, opt->f, opt->data, opt->stop, fmin, OPT_CG_PR );
            break;
        }
        case OPT_SG:{
			double eta = *opt->etas;
			opt_result eta_adapt_result = OPT_SUCCESS;
			if(opt->eta_count > 1){
				eta_adapt_result = optimize_stochastic_gradient_adapt(opt->parameters, opt->f, opt->grad_f, opt->reset, opt->etas, opt->eta_count, opt->data, &opt->stop, opt->verbosity, &eta, opt->threads);
			}
			opt->reset(opt->data);
			if(eta_adapt_result == OPT_SUCCESS){
				if(opt->verbosity > 0) printf("Stochastic gradient ascent using eta: %f\n", eta);
				result = optimize_stochastic_gradient(opt->parameters, opt->f, opt->grad_f, eta, opt->data, &opt->stop, opt->verbosity, fmin);
				opt->reset(opt->data);
			}
			else{
				result = OPT_FAIL;
			}
            break;
        }
		case OPT_TOPOLOGY:{
			result = topology_optimize(opt->data, fmin);
			break;
		}
		default:
			result = -100;
			break;
	}
	
	free(opt->stop.oldx);
	return result;
}

opt_result opt_optimize_univariate( Optimizer *opt, Parameter *p, double *fmin ){
	if ( opt->stop.time_max != 0 ) {
		time( &opt->stop.time_start );
	}
	opt_result result = OPT_SUCCESS;
	
	opt->dimension = 1;
	
	// should probably moved somewhere else
	opt->stop.oldx = dvector(opt->dimension);
    opt->stop.iter = 0;
    opt->stop.f_eval_current = 0;
    opt->stop.count = 0;
    
    if( opt->algorithm != OPT_BRENT ){
        error("optimize_univariate only works with Brent algorithm\n");
    }
    Parameters *ps = new_Parameters(1);
    Parameters_add(ps, p);
    
	result = brent_optimize( ps, opt->f, opt->data, &opt->stop, fmin );
	
    free_Parameters(ps);
	free(opt->stop.oldx);
	return result;
}

opt_result opt_maximize( Optimizer *opt, Parameters *ps, double *fmin ){
	opt_result result = opt_optimize(opt, ps, fmin);
    *fmin = - *fmin;
	return result;
}

opt_result opt_maximize_univariate( Optimizer *opt, Parameter *p, double *fmin ){
	opt_result result = opt_optimize_univariate(opt, p, fmin);
    *fmin = - *fmin;
	return result;
}

void opt_set_update_data_function( Optimizer *opt, opt_update_data uf ){
	if( opt != NULL ){
		opt->update = uf;
	}
}



bool xStop( const Parameters *x, double *xold, const double tolx){
	bool stop = true;
	
	for (int i = 0; i < Parameters_count(x); i++){
		//fprintf(stderr, "x=%f xold=%f tolx=%f (%d)\n", Parameters_value(x,i), xold[i], tolx,(fabs(Parameters_value(x,i) - xold[i] ) > tolx));
		if ( fabs(Parameters_value(x,i) - xold[i] ) > tolx ){
			stop = false;
		}
		xold[i] = Parameters_value(x,i);
	}
	
	return stop;
}

bool fxStop( double fx, double *fxold, const double tolfx){
	if ( fabs(fx - *fxold) > tolfx ){
		*fxold = fx;
		return false;
	}
	else{
		*fxold = fx;
		return true;
	}
}


Optimizer* new_Optimizer_from_json(json_node* node, Hashtable* hash){

	const char* algorithm_string = get_json_node_value_string(node, "algorithm");
	if (algorithm_string == NULL) {
		fprintf(stderr, "The `algorithm' key is not specified for object %s\n", get_json_node_value_string(node, "id"));
		exit(13);
	}
	
	if(strcasecmp(algorithm_string, "topology") == 0){
		Optimizer* opt = new_Optimizer(OPT_TOPOLOGY);
		opt->data = new_TopologyOptimizer_from_json(node, hash);
		return opt;
	}
	
	char* allowed[] = {
		"algorithm",
		"eta",
		"frequency_check",
		"list",
		"max",
		"min",
		"model",
		"parameters",
		"precision",
		"rounds",
		"threads",
		"tol",
		"treelikelihood",
		"verbosity"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	
	size_t max = get_json_node_value_size_t(node, "max", 1000);
	size_t min = get_json_node_value_size_t(node, "min", 1);
	double precision = get_json_node_value_double(node, "precision", 0.001);
	Parameters* parameters = new_Parameters(1);
	Optimizer* opt = NULL;
	
	if (strcasecmp(algorithm_string, "meta") == 0) {
		opt = new_Optimizer(OPT_META);
		OptimizerSchedule* schedule = opt_get_schedule(opt);
		json_node* list_node = get_json_node(node, "list");
		for (int i = 0; i < list_node->child_count; i++) {
			json_node* child = list_node->children[i];
			Optimizer* opt_child = new_Optimizer_from_json(child, hash);
			opt_add_optimizer(opt, opt_child);
			int child_rounds = get_json_node_value_int(child, "rounds", 1);
			opt->schedule->rounds[opt->schedule->count-1] = child_rounds;
		}
		opt_set_max_iteration(opt, max);
		opt_set_min_iteration(opt, min);
		opt_set_tolfx(opt, precision);
	}
	// Conjugate gradient
	else if (strcasecmp(algorithm_string, "cg") == 0) {
		opt = new_Optimizer(OPT_CG_PR);
		get_parameters_references(node, hash, parameters);
		opt_set_parameters(opt, parameters);
		opt_set_max_iteration(opt, max);
		opt_set_min_iteration(opt, min);
		opt_set_tolfx(opt, precision);
	}
	else if (strcasecmp(algorithm_string, "brent") == 0 || strcasecmp(algorithm_string, "serial") == 0) {
		if (get_json_node(node, "treelikelihood") != NULL) {
			const char* ref = get_json_node_value_string(node, "treelikelihood");
			opt = new_Optimizer(OPT_SERIAL_BRENT);
			opt->treelikelihood = Hashtable_get(hash, ref+1);
			//opt->treelikelihood->ref_count++;
		}
		else{
			get_parameters_references(node, hash, parameters);
			if (Parameters_count(parameters) == 1) {
				opt = new_Optimizer(OPT_BRENT);
			}
			else{
				opt = new_Optimizer(OPT_SERIAL_BRENT);
			}
			opt_set_parameters(opt, parameters);
		}
		opt_set_max_iteration(opt, max);
		opt_set_min_iteration(opt, min);
		opt_set_tolx(opt, precision);
	}
    // stochastic gradient
    else if(strcasecmp(algorithm_string, "sg") == 0){
        opt = new_Optimizer(OPT_SG);
		opt->stop.frequency_check = get_json_node_value_size_t(node, "frequency_check", 100);
        opt->stop.tolfx = get_json_node_value_double(node, "tol", 0.001);
        json_node* etas = get_json_node(node, "eta");
        if (etas != NULL && etas->node_type == MJSON_ARRAY) {
			opt->etas = dvector(etas->child_count);
			for (int i = 0; i < etas->child_count; i++) {
				json_node* child = etas->children[i];
				opt->etas[i] = atof((char*)child->value);
			}
			opt->eta_count = etas->child_count;
		}
		else{
			opt->etas = dvector(1);
			opt->etas[0] = get_json_node_value_double(node, "eta", 1.0);
			opt->eta_count = 1;
		}
		
		opt_set_max_iteration(opt, max);
		opt_set_min_iteration(opt, min);
    }
	json_node* model_node = get_json_node(node, "model");
    
	if(model_node != NULL){
        const char* ref = (char*)model_node->value;
        Model* model = Hashtable_get(hash, ref+1);
		opt_set_data(opt, model);
		opt_set_objective_function(opt, model_negative_logP);
		opt->grad_f = _gradient;
		
		if(strcasecmp(algorithm_string, "sg") == 0){
			get_parameters_references(node, hash, parameters);
			opt_set_parameters(opt, parameters);
			opt_set_objective_function(opt, _logP);
		}
		opt->reset = _reset;
    }
	opt->threads = get_json_node_value_size_t(node, "threads", 1);
	opt->verbosity = get_json_node_value_int(node, "verbosity", 1);
	free_Parameters(parameters);
	return opt;
}

