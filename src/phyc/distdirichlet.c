//
//  distdirichlet.c
//  physher
//
//  Created by Mathieu Fourment on 2/04/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distdirichlet.h"

#include <assert.h>

#include <gsl/gsl_randist.h>

#include "matrix.h"
#include "dirichlet.h"
#include "parametersio.h"
#include "descriptivestats.h"
#include "statistics.h"

// Flat dirichlet
double DistributionModel_log_flat_dirichlet(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	dm->lp = log(ddirchlet_flat(dm->simplex->K));
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_flat_dirichlet_with_values(DistributionModel* dm, const double* values){
	if(!dm->need_update) return dm->lp;
	dm->lp = log(ddirchlet_flat(dm->simplex->K));
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_dlog_flat_dirichlet(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_d2log_flat_dirichlet(DistributionModel* dm, const Parameter* p){
	return 0.0;
}

double DistributionModel_log_dirichlet(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	if(dm->x != NULL){
		for (int i = 0; i < Parameters_count(dm->x); i++) {
			dm->tempx[i] = Parameters_value(dm->x, i);
			dm->lp = ddirchletln(dm->tempx, Parameters_count(dm->x), dm->tempp);
		}
	}
	else if(dm->simplex != NULL){
		dm->lp = gsl_ran_dirichlet_lnpdf(dm->simplex->K, dm->tempp, dm->simplex->get_values(dm->simplex));
	}
	else{
		assert(0);
	}
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_dirichlet_with_values(DistributionModel* dm, const double* values){
	if(!dm->need_update) return dm->lp;
	dm->lp = gsl_ran_dirichlet_lnpdf(dm->simplex->K, dm->tempp, values);
	dm->need_update = false;
	return dm->lp;
}

static void DistributionModel_dirichlet_sample(DistributionModel* dm, double* samples){
	gsl_ran_dirichlet(dm->rng, dm->simplex->K, dm->tempp, samples);
}

static double DistributionModel_dirichlet_sample_evaluate(DistributionModel* dm){
	double* samples = dvector(dm->simplex->K);
	gsl_ran_dirichlet(dm->rng, dm->simplex->K, dm->tempp, samples);
	double logP = DistributionModel_log_dirichlet_with_values(dm, samples);
	dm->simplex->set_values(dm->simplex, samples);
	free(samples);
	return logP;
}

//TODO: implement
double DistributionModel_dlog_dirichlet(DistributionModel* dm, const Parameter* p){
	// find corresponding alpha
	// return (alpha-1.0)/Parameter_value(p);
	exit(1);
	return 0;
}

//TODO: implement
double DistributionModel_d2log_dirichlet(DistributionModel* dm, const Parameter* p){
	// find corresponding alpha
	// return -(alpha-1.0)/(Parameter_value(p)*Parameter_value(p));
	exit(1);
	return 0;
}

DistributionModel* new_FlatDirichletDistributionModel(Simplex* simplex){
	DistributionModel* dm = new_DistributionModelSimplex(NULL, simplex);
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_flat_dirichlet;
	dm->logP_with_values = DistributionModel_log_flat_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_flat_dirichlet;
	dm->d2logP = DistributionModel_d2log_flat_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i <simplex->K; i++) {
		dm->tempp[i] = 1;
	}
	return dm;
}

DistributionModel* new_DirichletDistributionModel(const double* alpha, Simplex* simplex){
	Parameters* ps = new_Parameters(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		Parameters_move(ps, new_Parameter("dirichlet.", alpha[i], NULL));
	}
	DistributionModel* dm = new_DistributionModelSimplex(ps, simplex);// len(x)==len(alpha)
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_dirichlet;
	dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
	dm->tempx = dvector(simplex->K);
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		dm->tempp[i] = Parameters_value(dm->parameters, i);
	}
	free_Parameters(ps);
	return dm;
}

DistributionModel* new_DirichletDistributionModel_with_parameters(const Parameters* parameters, Simplex* simplex){
	Parameters* ps = new_Parameters(simplex->K);
	Parameters_add_parameters(ps, parameters);
	DistributionModel* dm = new_DistributionModelSimplex(ps, simplex);// len(x)==len(alpha)
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_dirichlet;
	dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
	dm->tempx = dvector(simplex->K);
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i < simplex->K; i++) {
		dm->tempp[i] = Parameters_value(dm->parameters, i);
	}
	return dm;
}


Model* new_DirichletDistributionModel_from_json(json_node* node, Hashtable* hash){
	json_node* x_node = get_json_node(node, "x");
	char* ref = (char*)x_node->value;
	Model* msimplex = Hashtable_get(hash, ref+1);
	Simplex* simplex = msimplex->obj;
	
	char* id = get_json_node_value_string(node, "id");
	
	Parameters* parameters = new_Parameters(1);
	
	DistributionModel* dm = NULL;
	
	char* file = get_json_node_value_string(node, "file");
	
	// empirical
	if (file != NULL) {
		char* name = msimplex->name;
		char** names = malloc(sizeof(char*)*simplex->K);
		StringBuffer* buffer = new_StringBuffer(10);
		for(int i = 0; i < simplex->K; i++){
			StringBuffer_empty(buffer);
			StringBuffer_append_format(buffer, "%s%d", name, i+1);
			names[i] = StringBuffer_tochar(buffer);
		}
		size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
		Vector** samples = read_log_for_names_t(file, burnin, names, simplex->K);
		
		size_t paramCount = simplex->K;
		double* means = dvector(paramCount);
		double* variances = dvector(paramCount);
		for (int i = 0; i < paramCount; i++) {
			means[i] = dmean(Vector_data(samples[i]), Vector_length(samples[i]));
			variances[i] = variance(Vector_data(samples[i]), Vector_length(samples[i]), means[i]);
		}
		double num = 0;
		double denom = 0;
		for (int i = 0; i < paramCount; i++) {
			num += means[i]*pow(1.0 - means[i], 2);
			denom += means[i]*variances[i]*(1.0 - means[i]);
		}
		double mhat = num/denom - 1.0;
		for (int i = 0; i < paramCount; i++) {
			Parameters_move(parameters, new_Parameter("alpha", mhat*means[i], NULL));
			free(names[i]);
			free_Vector(samples[i]);
		}
		free(samples);
		
		dm = new_DirichletDistributionModel_with_parameters(parameters, simplex);
		
		free(means);
		free(variances);
		free(names);
		free_StringBuffer(buffer);
	}
	else{
		get_parameters_references(node, hash, parameters);
		int i = 0;
		for ( ; i < Parameters_count(parameters); i++) {
			if(Parameters_value(parameters, i) != 1.0) break;
		}
		if(i == Parameters_count(parameters)){
			dm = new_FlatDirichletDistributionModel(simplex);
			dm->lp = log(ddirchlet_flat(dm->simplex->K));
			dm->need_update = false;
		}
		else{
			dm = new_DirichletDistributionModel_with_parameters(parameters, simplex);
		}
		for (int i = 0; i < Parameters_count(parameters); i++) {
			Hashtable_add(hash, Parameters_name(parameters, i), Parameters_at(parameters, i));
		}
	}
	
	Model* model = new_DistributionModel3(id, dm, msimplex);
	
	model->samplable = true;
	dm->shift = get_json_node_value_double(node, "shift", 0);
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	free_Parameters(parameters);
	
	return model;

}