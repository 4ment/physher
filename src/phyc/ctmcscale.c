//
//  ctmcscale.c
//  physher
//
//  Created by mathieu on 31/3/21.
//  Copyright © 2021 Mathieu Fourment. All rights reserved.
//

#include "ctmcscale.h"

#include <gsl/gsl_sf_gamma.h>

#include "distmodel.h"
#include "mathconstant.h"


double DistributionModel_log_ctmc_scale(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	Tree* tree = dm->tree;
	size_t rateCount = Parameters_count(dm->x);
	Node** nodes = Tree_nodes(tree);
	double totalTreeTime = 0.0;
	double shape = 0.5;
	Tree_update_heights(tree); // make sure node heights are updated
	for(size_t i = 0; i < Tree_node_count(tree); i++){
		if(!Node_isroot(nodes[i])){
			totalTreeTime += Node_time_elapsed(nodes[i]);
		}
	}
	double logNormalization = 0.5 * log(totalTreeTime) - gsl_sf_lngamma(0.5);
	dm->lp = 0.0;
	for(size_t i = 0; i < rateCount; i++){
		dm->lp += logNormalization - shape * log(Parameters_value(dm->x, i)) - Parameters_value(dm->x, i) * totalTreeTime;
	}
	dm->need_update = false;
	
	return dm->lp;
}


double DistributionModel_log_ctmc_scale_with_values(DistributionModel* dm, const double* values){
	Tree* tree = dm->tree;
	Node** nodes = Tree_nodes(tree);
	double totalTreeTime = 0.0;
	double shape = 0.5;
	Tree_update_heights(tree); // make sure node heights are updated
	for(size_t i = 0; i < Tree_node_count(tree); i++){
		if(!Node_isroot(nodes[i])){
			totalTreeTime += Node_time_elapsed(nodes[i]);
		}
	}
	size_t rateCount = Parameters_count(dm->x);
	double logNormalization = shape * log(totalTreeTime) - gsl_sf_lngamma(shape);
	dm->lp = 0.0;
	for(size_t i = 0; i < rateCount; i++){
		dm->lp += logNormalization - shape * log(values[i]) - values[i] * totalTreeTime;
	}
	dm->need_update = false;
	return dm->lp;
}


double DistributionModel_dlog_ctmc_scale(DistributionModel* dm, const Parameter* p){
	fprintf(stderr, "DistributionModel_dlog_ctmc_scale not implemented\n");
	exit(2);
}

double DistributionModel_d2log_ctmc_scale(DistributionModel* dm, const Parameter* p){
	fprintf(stderr, "DistributionModel_d2log_ctmc_scale not implemented\n");
	exit(2);
}

double DistributionModel_ddlog_ctmc_scale(DistributionModel* dm, const Parameter* p1, const Parameter* p2){
	fprintf(stderr, "DistributionModel_ddlog_ctmc_scale not implemented\n");
	exit(2);
}


static void DistributionModel_ctmc_scale_sample(DistributionModel* dm, double* samples){
	fprintf(stderr, "DistributionModel_ctmc_scale_sample not implemented\n");
	exit(2);
}

static double DistributionModel_ctmc_scale_sample_evaluate(DistributionModel* dm){
	fprintf(stderr, "DistributionModel_ctmc_scale_sample_evaluate not implemented\n");
	exit(2);
}

DistributionModel* new_CTMCScale_with_parameters(Parameters* x, Tree* tree){
	DistributionModel* dm = new_DistributionModel(NULL, 0, x);
	dm->type = DISTRIBUTION_CTMC_SCALE;
	dm->parameterization = 0;
	dm->logP = DistributionModel_log_ctmc_scale;
	dm->logP_with_values = DistributionModel_log_ctmc_scale_with_values;
	dm->dlogP = DistributionModel_dlog_ctmc_scale;
	dm->sample = DistributionModel_ctmc_scale_sample;
	dm->sample_evaluate = DistributionModel_ctmc_scale_sample_evaluate;
	dm->d2logP = DistributionModel_d2log_ctmc_scale;
	dm->ddlogP = DistributionModel_ddlog_ctmc_scale;
	dm->tree = tree;
    dm->shift = 0;
	return dm;
}

Model* new_CTMCScaleModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	char* model_key = get_json_node_value_string(node, "tree");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
	
	Model *mtree = Hashtable_get(hash, model_key+1);
	DistributionModel* dm = new_CTMCScale_with_parameters(x, mtree->obj);
	Model* model = new_DistributionModel3(id, dm, mtree);
	mtree->listeners->add(mtree->listeners, model);
	
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
	
	return model;
}
