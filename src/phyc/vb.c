//
//  vb.c
//  physher
//
//  Created by Mathieu Fourment on 24/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "vb.h"

#include <string.h>

#include "optimizer.h"
#include "matrix.h"
#include "gaussian.h"
#include "transforms.h"
#include "solve.h"
#include "tree.h"

#include <gsl/gsl_randist.h>

#define MY_PI acos(-1.0)
#define LOG_TWO_PI (log(2.0)+log(MY_PI))

void init_fullrank_normal(variational_t* var){return;
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

void init_meanfield_normal(variational_t* var){
	size_t dim = Parameters_count(var->parameters);
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		if(Parameter_value(p) < 1.0e-4){
			p->estimate = false;
		}
	}
	
	return;
	
	Model* posterior = var->posterior;
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

// Mean-field
double elbo_meanfield(variational_t* var){
	if (var->initialized == false) {
		init_meanfield_normal(var);
		var->initialized = true;
	}
//	printf("elbo\n");
	double elbo = 0;
	size_t dim = Parameters_count(var->parameters);
	Model* posterior = var->posterior;
	
	// Entropy: 0.5(1 + log(2 \pi s^2))
	size_t count = 0;
	double entropy = 0;
	for (size_t i = 0; i < dim; i++) {
		if(Parameters_estimate(var->var_parameters, dim+i)){
			entropy += Parameters_value(var->var_parameters, dim+i);
			count++;
		}
	}
	entropy = 0.5 * count * (1.0 + LOG_TWO_PI);
	int inf_count = 0;
	
	for (int i = 0; i < var->elbo_samples; i++) {
		double jacobian = 0.0;
		
		for (int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			if(!p->estimate) continue;
			double var_mean = Parameters_value(var->var_parameters, j);
			double var_sd = Parameters_value(var->var_parameters, dim+j);
			double zeta = rnorm() * exp(var_sd) + var_mean;
			double theta = inverse_transform(zeta, Parameter_lower(p), Parameter_upper(p), &jacobian);
			Parameter_set_value(p, theta);
		}
		
		double logP = posterior->logP(posterior);
		if(!isinf(logP))elbo += posterior->logP(posterior) + jacobian;
		
		if(isinf(elbo)){
			return elbo;
		}
	}
	return elbo/(var->elbo_samples-inf_count) + entropy;
}

void grad_elbo_meanfield(variational_t* var, double* grads){
	if (var->initialized == false) {
		init_meanfield_normal(var);
		// save for later
		for (int i = 0; i < Parameters_count(var->parameters); i++) {
			Parameter_store(Parameters_at(var->parameters, i));
		}
		var->initialized = true;
	}

	Model* posterior = var->posterior;
	size_t dim = Parameters_count(var->parameters);
	double* eta = dvector(dim);
	double* zeta = dvector(dim);
	memset(grads, 0, sizeof(double)*dim*2);
	
	for (int i = 0; i < var->grad_samples; i++) {
		for ( int j = 0; j < dim; j++) {
			if(!Parameters_estimate(var->parameters, j)) continue;
			double var_mean = Parameters_value(var->var_parameters, j);
			double var_sd = Parameters_value(var->var_parameters, dim+j);
			eta[j] = rnorm();
			zeta[j] = eta[j] * exp(var_sd) + var_mean;
			
			Parameter* p = Parameters_at(var->parameters, j);
			double theta = inverse_transform2(zeta[j], Parameter_lower(p), Parameter_upper(p));
			Parameter_set_value(p, theta);
		}
		
		for ( int j = 0; j < dim; j++) {
			Parameter* p = Parameters_at(var->parameters, j);
			if(!p->estimate) continue;
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


double elbo_fullrank(variational_t* var){
	if (var->initialized == false) {
		init_fullrank_normal(var);
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

void grad_elbo_fullrank(variational_t* var, double* grads){
	if (var->initialized == false) {
		init_fullrank_normal(var);
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

double elbo( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;exit(12);
	return model->logP(model);
}

void grad_elbo( Parameters *params, double *grad, void *data ){
	Model* model = (Model*)data;exit(13);
	model->gradient(model, grad);
}

// check success
bool variational_sample_meanfield(variational_t* var, double* values){
	int dim = Parameters_count(var->parameters);
	for (int i = 0; i < dim; i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		if (!p->estimate) {
			values[i] = 0;
			continue;
		}
		double mu = Parameters_value(var->var_parameters, i);
		double sd = Parameters_value(var->var_parameters, i+dim);
		const double zeta = rnorm() * exp(sd) + mu;
		values[i] = inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
	}
	return true;
}

bool variational_sample_some_meanfield(variational_t* var, const Parameters* parameters, double* values){
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

bool variational_sample_fullrank(variational_t* var, double* values){
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

//double variational_sample_meanfield(variational_t* var, const Parameter* p){
//	int index = 0;
//	for (; index < Parameters_count(var->parameters); index++) {
//		if(strcmp(p->name, Parameters_at(var->parameters, index)) == 0) break
//			}
//	const double zeta = rnorm() * var->var_parameters[index+Parameters_count(var->parameters)] + var->var_parameters[index]
//	return inverse_transform2(zeta, Parameter_lower(p), Parameter_upper(p));
//}

Model* new_Variational_from_json(json_node* node, Hashtable* hash){
	const char* posterior_string = get_json_node_value_string(node, "posterior");
	const char* var_string = get_json_node_value_string(node, "var");
	const char* dist_string = get_json_node_value_string(node, "distribution");
	const char* id = get_json_node_value_string(node, "id");
	const char* filename = get_json_node_value_string(node, "log");
	json_node* parameters_node = get_json_node(node, "parameters");
	json_node* tree_node = get_json_node(node, "tree");
	
	variational_t* var = malloc(sizeof(variational_t));
	var->file = NULL;
	var->sample = NULL;
	var->sample_some = NULL;
	var->parameters = new_Parameters(1);
	if(parameters_node){
		get_parameters_references(node, hash, var->parameters);
	}
	if (tree_node != NULL) {
		char* tree_string = get_json_node_value_string(node, "tree");
		Model* mtree = Hashtable_get(hash, tree_string+1);
		Tree* tree = mtree->obj;
		Node** nodes = Tree_get_nodes(tree, POSTORDER);
	
		for(int i = 0; i < Tree_node_count(tree)-2; i++){
			Parameters_add(var->parameters, nodes[i]->distance);
		}
	}
	
	size_t dim = Parameters_count(var->parameters);
	
	var->posterior = Hashtable_get(hash, posterior_string+1);
	var->posterior->ref_count++;
	
	var->elbo_samples = get_json_node_value_size_t(node, "elbosamples", 100);
	var->grad_samples = get_json_node_value_size_t(node, "gradsamples", 1);
	
	StringBuffer* buffer = new_StringBuffer(10);
	if(strcasecmp(var_string, "meanfield") == 0){
		printf("Creating meanfield variational model\n");
		var->var_parameters = new_Parameters(dim*2); // mean + sd
		for(int i = 0; i < dim; i++){
			StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
			StringBuffer_append_string(buffer, ".mean");
			Parameter* p = new_Parameter(buffer->c, 0, NULL);
			p->id = i;
			Parameters_move(var->var_parameters, p);
		}
		for(int i = 0; i < dim; i++){
			StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
			StringBuffer_append_string(buffer, ".sd");
			Parameter* p = new_Parameter(buffer->c, 0, NULL);
			p->id = i+dim;
			Parameters_move(var->var_parameters, p);
		}
		var->elbofn = elbo_meanfield;
		var->grad_elbofn = grad_elbo_meanfield;
		var->f = elbo;
		var->grad_f = grad_elbo;
		var->sample = variational_sample_meanfield;
	}
	else if(strcasecmp(var_string, "fullrank") == 0){
		printf("Creating fullrank variational model\n");
		size_t n = 2*dim+(dim*dim-dim)/2; // mus, vars, covs
		var->var_parameters = new_Parameters(n);
		json_node* init_node = get_json_node(node, "init");
		
		// means
		for(int i = 0; i < dim; i++){
			StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
			StringBuffer_append_string(buffer, ".mean");
			Parameter* p = new_Parameter(buffer->c, 1, NULL);
			Parameters_move(var->var_parameters, p);
		}
	 
		// sparse covariance matrix
		if(tree_node != NULL){
			char* tree_string = get_json_node_value_string(node, "tree");
			Model* mtree = Hashtable_get(hash, tree_string+1);
			Tree* tree = mtree->obj;
			Node** nodes = Tree_get_nodes(tree, POSTORDER);
			for (int i = 0; i < Tree_node_count(tree)-2; i++) {
				Node* nodei = nodes[i];
				for (int j = 0; j < i; j++) {
					Node* nodej = nodes[j];
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_strings(buffer, 2, ".", Parameters_name(var->parameters, j));
					StringBuffer_append_string(buffer, ".cov");
					Parameter* p = new_Parameter(buffer->c, 1, NULL);
					if (Node_parent(nodei) != nodej && Node_parent(nodej) != nodei && Node_sibling(nodei) != nodej) {
						p->estimate = false;
						p->value = 0;
					}
//					p->estimate = false;
//					p->value = 0;
					Parameters_move(var->var_parameters, p);
				}
				StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
				StringBuffer_append_string(buffer, ".var");
				Parameter* p = new_Parameter(buffer->c, 1, NULL);
				Parameters_move(var->var_parameters, p);
			}
		}
		// full covariance matrix
		else{
			size_t row = dim;
			for(int i = 0; i < dim; i++){
				for(int j = 0; j < i; j++){
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_strings(buffer, 2, ".", Parameters_name(var->parameters, j));
					StringBuffer_append_string(buffer, ".cov");
					Parameter* p = new_Parameter(buffer->c, 1, NULL);
					Parameters_move(var->var_parameters, p);
					row++;
				}
				StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
				StringBuffer_append_string(buffer, ".var");
				Parameter* p = new_Parameter(buffer->c, 1, NULL);
				Parameters_move(var->var_parameters, p);
	//			p->estimate = false;
				row++;
			}
		}
		
		// this is run at startup so not useful here
//		if(init_node != NULL){
//			char* model = get_json_node_value_string(init_node, "model");
//			Model* mvar = Hashtable_get(hash, model+1);
//			variational_t* var2 = mvar->obj;
//			for(int i = 0; i < dim; i++){
//				if(strcmp(Parameters_name(var->var_parameters, i), Parameters_name(var2->var_parameters, i)) != 0){
//					printf("%s %s\n", Parameters_name(var->var_parameters, i), Parameters_name(var2->var_parameters, i));
//					exit(81);
//				}
//				Parameters_set_value(var->var_parameters, i, Parameters_value(var2->var_parameters, i));
//			}
//			int r = 0;
//			for(int i = 0; i < dim; i++){
//				r += i;
//				//printf("%s %s\n", Parameters_name(var->var_parameters, r+dim), Parameters_name(var2->var_parameters, i+dim));
//				Parameters_set_value(var->var_parameters, r+dim, Parameters_value(var2->var_parameters, i+dim));
//				r++;
//			}
//		}
		var->elbofn = elbo_fullrank;
		var->grad_elbofn = grad_elbo_fullrank;
		var->f = elbo;
		var->grad_f = grad_elbo;
		var->sample = variational_sample_fullrank;
	}
	
	var->initialized = false;
	var->iter = 0;
	if (filename != NULL) {
		var->file = fopen(filename, "w");
		fprintf(var->file, "iteration");
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%s", Parameters_name(var->var_parameters, i));
		}
		for(int i = 0; i < Parameters_count(var->var_parameters); i++){
			fprintf(var->file, ",%s.grad", Parameters_name(var->var_parameters, i));
		}
		fprintf(var->file, "\n");
	}
	
	free_StringBuffer(buffer);
	var->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	return new_VariationalModel(id, var);
}

static void _variational_model_free( Model *self ){
	if(self->ref_count == 1){
		printf("Free variational model %s\n", self->name);
		variational_t* var = (variational_t*)self->obj;
//		for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
//			printf("%s %f\n", Parameters_name(var->var_parameters, i), Parameters_value(var->var_parameters, i));
//		}
		var->posterior->free(var->posterior);
		free_Parameters(var->var_parameters);
		free_Parameters(var->parameters);
		if (var->file != NULL) fclose(var->file);
		free(var);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static void _variational_model_handle_change( Model *self, Model *model, int index ){
	variational_t* var = (variational_t*)self->obj;
	size_t i = index;
	if(index >= Parameters_count(var->parameters)){
		index -= Parameters_count(var->parameters);
	}
//	self->listeners->fire( self->listeners, self, index );
}

static Model* _variational_model_clone(Model* self, Hashtable *hash){
	//TODO: _variational_model_clone
	error(_variational_model_clone);
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	return NULL;
}

static double _variational_model_logP(Model *self){
	variational_t* var = (variational_t*)self->obj;
	return var->elbofn(var);
}

static void _variational_model_gradient(Model *self, double* grad){
	variational_t* var = (variational_t*)self->obj;
	return var->grad_elbofn(var, grad);
}

static void _variational_model_get_free_parameters(Model* self, Parameters* parameters){
	variational_t* var = (variational_t*)self->obj;
//	Parameters_add_free_parameters(parameters, var->var_parameters);
	Parameters_add_parameters(parameters, var->var_parameters);

}

Model* new_VariationalModel(const char* name, variational_t* var){
	Model *model = new_Model("variational", name, var);
	model->free = _variational_model_free;
	model->clone = _variational_model_clone;
	model->logP = _variational_model_logP;
	model->gradient = _variational_model_gradient;
	model->get_free_parameters = _variational_model_get_free_parameters;
	
	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
		Parameters_at(var->var_parameters, i)->listeners->add( Parameters_at(var->var_parameters, i)->listeners, model );
	}
	return model;
}

void free_Variational(variational_t* var){
//	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
//		printf("%s %f\n", Parameters_name(var->var_parameters, i), Parameters_value(var->var_parameters, i));
//	}
	free_Parameters(var->var_parameters);
    free_Model(var->posterior);
    free_Parameters(var->parameters);
	free(var);
}
