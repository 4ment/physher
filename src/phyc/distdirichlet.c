//
//  distdirichlet.c
//  physher
//
//  Created by Mathieu Fourment on 2/04/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "distdirichlet.h"

#include <assert.h>
#include <strings.h>

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
    for (int i = 0; i < Parameters_count(dm->parameters[0]); i++) {
        dm->tempp[i] = Parameters_value(dm->parameters[0], i);
    }
    dm->lp = gsl_ran_dirichlet_lnpdf(dm->simplex->K, dm->tempp, dm->simplex->get_values(dm->simplex));
	dm->need_update = false;
	return dm->lp;
}

double DistributionModel_log_dirichlet_with_values(DistributionModel* dm, const double* values){
	if(!dm->need_update) return dm->lp;
    for (int i = 0; i < Parameters_count(dm->x); i++) {
        dm->tempp[i] = Parameters_value(dm->parameters[0], i);
    }
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

double DistributionModel_dlog_dirichlet(DistributionModel* dm, const Parameter* p){
	for (int i = 0; i < Parameters_count(dm->x); i++) {
        if( p == Parameters_at(dm->x, i)){
            return (Parameters_value(dm->parameters[0], i)-1.0)/Parameter_value(p);
        }
    }
	return 0;
}

//TODO: implement
double DistributionModel_d2log_dirichlet(DistributionModel* dm, const Parameter* p){
	// find corresponding alpha
	// return -(alpha-1.0)/(Parameter_value(p)*Parameter_value(p));
	exit(1);
	return 0;
}

static void _DistributionModel_error_sample_dirichlet(DistributionModel* dm, double* samples){
    fprintf(stderr, "_DistributionModel_error_sample_dirichlet not implemented\n");
    exit(1);
}

static double _DistributionModel_error_sample_evaluate_dirichlet(DistributionModel* dm){
    fprintf(stderr, "_DistributionModel_error_sample_evaluate_dirichlet not implemented\n");
    exit(1);
}

DistributionModel* new_FlatDirichletDistributionModel(Simplex* simplex){
    Parameters* x = new_Parameters_with_name(Parameters_name2(simplex->parameters), simplex->K-1);
    Parameters_add_parameters(x, simplex->parameters);
	DistributionModel* dm = new_DistributionModel(NULL, 0, x);
    dm->simplex = simplex;
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_flat_dirichlet;
	dm->logP_with_values = DistributionModel_log_flat_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_flat_dirichlet;
	dm->d2logP = DistributionModel_d2log_flat_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
    dm->tempx = dvector(simplex->K);
	dm->tempp = dvector(simplex->K);
	for (int i = 0; i <simplex->K; i++) {
		dm->tempp[i] = 1;
	}
    dm->shift = 0;
	return dm;
}

DistributionModel* new_DirichletDistributionModel_with_parameters(Parameters** parameters, Simplex* simplex){
    Parameters* x = new_Parameters_with_name(Parameters_name2(simplex->parameters), simplex->K-1);
    Parameters_add_parameters(x, simplex->parameters);
    DistributionModel* dm = new_DistributionModel(parameters, 1, x);
    dm->simplex = simplex;
	dm->type = DISTRIBUTION_DIRICHLET;
	dm->logP = DistributionModel_log_dirichlet;
	dm->logP_with_values = DistributionModel_log_dirichlet_with_values;
	dm->dlogP = DistributionModel_dlog_dirichlet;
	dm->d2logP = DistributionModel_d2log_dirichlet;
	dm->ddlogP = DistributionModel_ddlog_0;
	dm->sample = DistributionModel_dirichlet_sample;
	dm->sample_evaluate = DistributionModel_dirichlet_sample_evaluate;
	dm->tempx = dvector(simplex->K);
    if(parameters != 0){
        dm->tempp = dvector(simplex->K);
        for (int i = 0; i < simplex->K; i++) {
            dm->tempp[i] = Parameters_value(dm->parameters[0], i);
        }
    }
    dm->shift = 0;
	return dm;
}


Model* new_DirichletDistributionModel_from_json(json_node* node, Hashtable* hash){
    char* id = get_json_node_value_string(node, "id");
    json_node* x_node = get_json_node(node, "x");
    
//    Parameters* x = new_Parameters(1);
    
    Model* msimplex = NULL;

    if(x_node->node_type == MJSON_STRING){
        char* ref = (char*)x_node->value;
        msimplex = Hashtable_get(hash, ref+1);
        msimplex->ref_count++;
    }
    else if(x_node->node_type == MJSON_OBJECT){
        msimplex = new_SimplexModel_from_json(x_node, hash);
        Hashtable_add(hash, msimplex->name, msimplex);
    }
    
    char* file = get_json_node_value_string(node, "file");
    Parameters** parameters = NULL;
    size_t parameters_dim = 0;
	DistributionModel* dm = NULL;
	
    // empirical
    if (file != NULL) {
        Simplex* simplex = msimplex->obj;
        size_t burnin = get_json_node_value_size_t(node, "burnin", 0);
        Vector** samples = read_log_for_parameters_t(file, burnin, simplex->parameters);
        size_t paramCount = Parameters_count(simplex->parameters);
        parameters_dim = 1;
        parameters = malloc(sizeof(Parameters*));
        parameters[0] = new_Parameters(paramCount);
        
        char* name = msimplex->name;
        char** names = malloc(sizeof(char*)*simplex->K);
        StringBuffer* buffer = new_StringBuffer(10);
        for(int i = 0; i < simplex->K; i++){
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%s%d", name, i+1);
            names[i] = StringBuffer_tochar(buffer);
        }
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
            Parameters_move(parameters[0], new_Parameter("alpha", mhat*means[i], new_Constraint(0, INFINITY)));
            free(names[i]);
            free_Vector(samples[i]);
        }
        free(samples);
        free(means);
        free(variances);
        free(names);
        free_StringBuffer(buffer);
		dm = new_DirichletDistributionModel_with_parameters(parameters, msimplex->obj);
    }
    // Flat dirichlet
    else if(get_json_node(node, "parameters") == NULL){
        dm = new_FlatDirichletDistributionModel(msimplex->obj);
    }
    else{
        json_node* parameters_node = get_json_node(node, "parameters");
        
        if (parameters_node->child_count != 1 && strcasecmp(parameters_node->children[0]->key, "concentration") != 0 && strcasecmp(parameters_node->children[0]->key, "alpha") != 0) {
            fprintf(stderr, "Dirichlet distribution should be parametrized with concentration parameter\n");
            exit(13);
        }

        parameters = distmodel_get_parameters(id, parameters_node, hash, &parameters_dim);
		dm = new_DirichletDistributionModel_with_parameters(parameters, msimplex->obj);
    }
    
    dm->parameterization = 0;
    dm->shift = get_json_node_value_double(node, "shift", dm->shift);
    
    
    Model* model = new_DistributionModel3(id, dm, msimplex);
    msimplex->free(msimplex);
    
    model->samplable = true;
    dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    if(parameters != NULL){
        free_Parameters(parameters[0]);
        free(parameters);
    }
    
    return model;
}
