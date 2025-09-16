//
//  ctmcscale.c
//  physher
//
//  Created by mathieu on 31/3/21.
//  Copyright © 2021 Mathieu Fourment. All rights reserved.
//

#include "ctmcscale.h"

#include "mathconstant.h"
#include "matrix.h"
#include "gradient.h"


double DistributionModel_log_ctmc_scale(DistributionModel* dm){
	if(!dm->need_update) return dm->lp;
	Tree* tree = dm->tree;
	Parameter* x = Parameters_at(dm->x, 0);
	size_t rateCount = Parameter_size(x);
	Node** nodes = Tree_nodes(tree);
	double totalTreeTime = 0.0;
	double shape = 0.5;
	Tree_update_heights(tree); // make sure node heights are updated
	for(size_t i = 0; i < Tree_node_count(tree); i++){
		if(!Node_isroot(nodes[i])){
			totalTreeTime += Node_time_elapsed(nodes[i]);
		}
	}
	double logGammaHalf = 0.57236494292470041501;
	double logNormalization = shape * log(totalTreeTime) - logGammaHalf;
	dm->lp = 0.0;
	for(size_t i = 0; i < rateCount; i++){
		dm->lp += logNormalization - shape * log(Parameter_value_at(x, i)) - Parameter_value_at(x, i) * totalTreeTime;
	}
	dm->need_update = false;
	
	return dm->lp;
}


// double DistributionModel_log_ctmc_scale_with_values(DistributionModel* dm, const double* values){
// 	Tree* tree = dm->tree;
// 	Node** nodes = Tree_nodes(tree);
// 	double totalTreeTime = 0.0;
// 	double shape = 0.5;
// 	Tree_update_heights(tree); // make sure node heights are updated
// 	for(size_t i = 0; i < Tree_node_count(tree); i++){
// 		if(!Node_isroot(nodes[i])){
// 			totalTreeTime += Node_time_elapsed(nodes[i]);
// 		}
// 	}
// 	size_t rateCount = Parameter_size(dm->x);
// 	double logGammaHalf = 0.57236494292470041501;
// 	double logNormalization = shape * log(totalTreeTime) - logGammaHalf;
// 	dm->lp = 0.0;
// 	for(size_t i = 0; i < rateCount; i++){
// 		dm->lp += logNormalization - shape * log(values[i]) - values[i] * totalTreeTime;
// 	}
// 	dm->need_update = false;
// 	return dm->lp;
// }

// Set flags for gradient calculation
void _ctmcscale_model_prepare_gradient(Model* self, const Parameters* ps){
	DistributionModel* dm = self->obj;
	Parameter* x = Parameters_at(dm->x, 0);
	size_t paramCount = Parameters_count(ps);
	bool prepare_tree = false;
	bool prepare_theta = false;
	bool prepare_clock = false;
	dm->prepared_gradient = 0;
	size_t gradient_length = 0;
	for (size_t i = 0; i < paramCount; i++) {
		Parameter* p = Parameters_at(ps, i);
		if (p->model == MODEL_TREE && !prepare_tree) {
			prepare_tree = true;
			dm->prepared_gradient |= GRADIENT_FLAG_TREE_HEIGHTS;
			gradient_length += Tree_tip_count(dm->tree) - 1;
		}
		else if (p->model == MODEL_TREE_TRANSFORM && !prepare_tree) {
			prepare_tree = true;
			dm->prepared_gradient |= GRADIENT_FLAG_TREE_RATIOS;
			gradient_length += Tree_tip_count(dm->tree) - 1;
		}
		else if(p->model == MODEL_BRANCHMODEL && !prepare_clock){
			prepare_clock = true;
			dm->prepared_gradient |= GRADIENT_FLAG_CLOCK_RATE;
			gradient_length += Parameter_size(x);
		}
	}
	if(dm->gradient == NULL){
		dm->gradient = calloc(gradient_length, sizeof(double));
		dm->gradient_length = gradient_length;
	}
	else if (dm->gradient_length < gradient_length) {
		dm->gradient = realloc(dm->gradient, sizeof(double)* gradient_length);
		dm->gradient_length = gradient_length;
	}
}

double _ctmcscale_model_dlogP_prepared(Model *self, const Parameter* p){
	// DistributionModel* dm = (DistributionModel*)self->obj;
	
	// if(dm->gradient_length == 0) return 0.0;

	// if(dm->need_update_gradient){
	// 	DistributionModel_gradient(self);
	// 	dm->need_update_gradient = false;
	// }
	
	// int prepare_tree_height = dm->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS;
	// int prepare_tree_ratio = dm->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS;
	// int prepare_rate = dm->prepared_gradient & GRADIENT_FLAG_CLOCK_RATE;
	// size_t offset = 0;
	
	// if(prepare_tree_height && p->model == MODEL_TREE){
	// 	offset = Tree_tip_count(dm->tree)-1;
	// 	return dm->gradient[Tree_node(dm->tree, p->id)->class_id];
	// }
	// else if(prepare_tree_ratio && p->model == MODEL_TREE_TRANSFORM){
	// 	offset = Tree_tip_count(dm->tree)-1;
	// 	return dm->gradient[p->id];
	// }
	// if (prepare_rate) {
	// 	if (p->model == MODEL_BRANCHMODEL) {
	// 		return dm->gradient[offset + p->id];
	// 	}
	// }
	return 0;
}

// by default we get the derivative wrt reparameterized node heights (ratios)
size_t DistributionModel_initialize_gradient(Model *self, int flags){
	DistributionModel* dm = self->obj;
	Parameter* x = Parameters_at(dm->x, 0);
	dm->prepared_gradient = flags;
	size_t gradient_length = 0;
	
	int prepare_tree_height = dm->prepared_gradient & GRADIENT_FLAG_TREE_HEIGHTS;
	int prepare_tree = dm->prepared_gradient & GRADIENT_FLAG_TREE_RATIOS;
	int prepare_rate = dm->prepared_gradient & GRADIENT_FLAG_CLOCK_RATE;
	
	if(flags == 0){
		prepare_tree = true;
		prepare_tree_height = false;
		prepare_rate = true;
		dm->prepared_gradient = GRADIENT_FLAG_TREE_RATIOS | GRADIENT_FLAG_CLOCK_RATE;
	}
	if(prepare_tree || prepare_tree_height){
		gradient_length += Tree_tip_count(dm->tree) - 1;
	}
	if (prepare_rate) {
		gradient_length += Parameter_size(x);
		
	}

	if(dm->gradient == NULL){
		dm->gradient = calloc(gradient_length, sizeof(double));
		dm->gradient_length = gradient_length;
	}
	else if (dm->gradient_length < gradient_length) {
		dm->gradient = realloc(dm->gradient, sizeof(double)* gradient_length);
		dm->gradient_length = gradient_length;
	}
	return gradient_length;
}

static void _calculate_height_gradient(Tree* tree, double rate, double shape, double totalTreeTime, double* gradient){
	size_t tipCount = Tree_tip_count(tree);
	double temp = shape/totalTreeTime - rate;
	for(size_t i = 0; i < tipCount-1; i++){
		gradient[i] = temp;
	}
	gradient[Tree_root(tree)->class_id] = 2.0*temp;
}

void CTMCModel_gradient(Model *self, int flags, double* gradient){
	DistributionModel* dm = self->obj;
	Tree* tree = dm->tree;
	Parameters* parameters = new_Parameters(1);
	if(flags & GRADIENT_FLAG_TREE_RATIOS){
		Parameters_add_parameters(parameters, get_reparams(tree));
	}
	else if(flags & GRADIENT_FLAG_TREE_HEIGHTS){
		Node** nodes = Tree_nodes(tree);
		size_t nodeCount = Tree_node_count(tree);
		for(size_t i = 0; i < nodeCount; i++){
			if(!Node_isleaf(nodes[i])){
				Parameters_add(parameters, nodes[i]->height);
			}
		}
	}
	
	if(flags & GRADIENT_FLAG_CLOCK_RATE){
		Parameters_add(parameters, Parameters_at(dm->x, 0));
	}

	Parameters_zero_grad(parameters);
	dm->gradient2(dm, parameters);

	size_t offset = 0;
	for(size_t i = 0; i < Parameters_count(parameters); i++){
		Parameter* p = Parameters_at(parameters, i);
		memcpy(gradient + offset, p->grad, sizeof(double)*Parameter_size(p));
		offset += Parameter_size(p);
	}

	free_Parameters(parameters);
}

double DistributionModel_ctmc_gradient(DistributionModel *dm, const Parameters* parameters){
	// only works for strict clock models
	Parameter* rate = Parameters_at(dm->x, 0);
	double xValue = Parameter_value(rate);
	Tree* tree = dm->tree;
	size_t tipCount = Tree_tip_count(tree);
	Node** nodes = Tree_nodes(tree);
	double shape = 0.5;
	Tree_update_heights(tree); // make sure node heights are updated
	size_t offset = 0;
	size_t nodeCount = Tree_node_count(tree);
	double totalTreeTime = 0.0;
	for(size_t i = 0; i < nodeCount; i++){
		if(!Node_isroot(nodes[i])){
			totalTreeTime += Node_time_elapsed(nodes[i]);
		}
	}

	Parameter* xx = Parameters_depends(parameters, rate);
	if(xx != NULL){
		double dLogP = -shape/xValue - totalTreeTime;
		rate->grad[0] += dLogP;
		if(rate != xx){
			// apply chain rule
			rate->transform->backward(rate->transform, &dLogP);
		}
	}

	Parameters* treeModelParameters = new_Parameters(1);
	Parameters* reparam = get_reparams(tree);
	if(reparam != NULL){
		for(size_t i = 0; i < Parameters_count(reparam); i++){
			Parameter* xx = Parameters_depends(parameters, Parameters_at(reparam, i));
			if(xx != NULL) {
				Parameters_add(treeModelParameters, xx);
			}
		}
	}
	else{
		for(size_t i = 0; i < Tree_node_count(tree); i++){
			Node* node = Tree_node(tree, i);
			if(!Node_isleaf(node)){
				Parameter* xx = Parameters_depends(parameters, node->height);
				if(xx != NULL){
					Parameters_add(treeModelParameters, xx);
				}
			}
		}
	}

	if(Parameters_count(treeModelParameters) > 0){
		double* heightGradient = dvector(nodeCount);
		_calculate_height_gradient(tree, xValue, shape, totalTreeTime, heightGradient);
		Tree_height_backward(tree, treeModelParameters, heightGradient);
		free(heightGradient);
	}
	free_Parameters(treeModelParameters);

	return 0;
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


static void DistributionModel_ctmc_scale_sample(DistributionModel* dm){
	fprintf(stderr, "DistributionModel_ctmc_scale_sample not implemented\n");
	exit(2);
}

// static double DistributionModel_ctmc_scale_sample_evaluate(DistributionModel* dm){
// 	fprintf(stderr, "DistributionModel_ctmc_scale_sample_evaluate not implemented\n");
// 	exit(2);
// }

DistributionModel* new_CTMCScale_with_parameters(Parameters* x, Tree* tree){
	DistributionModel* dm = new_DistributionModel(NULL, x);
	dm->type = DISTRIBUTION_CTMC_SCALE;
	dm->parameterization = 0;
	dm->logP = DistributionModel_log_ctmc_scale;
	dm->gradient2 = DistributionModel_ctmc_gradient;
	// dm->logP_with_values = DistributionModel_log_ctmc_scale_with_values;
	dm->dlogP = DistributionModel_dlog_ctmc_scale;
	dm->sample = DistributionModel_ctmc_scale_sample;
	// dm->sample_evaluate = DistributionModel_ctmc_scale_sample_evaluate;
	dm->d2logP = DistributionModel_d2log_ctmc_scale;
	dm->ddlogP = DistributionModel_ddlog_ctmc_scale;
	dm->tree = tree;
    dm->shift = 0;
	return dm;
}
Model* new_CTMCScaleModel(const char* name, DistributionModel* dm, Model* tree){
	Model* model = new_DistributionModel3(name, dm, tree);
	tree->listeners->add(tree->listeners, model);
	model->prepare_gradient = _ctmcscale_model_prepare_gradient;
	model->dlogP = _ctmcscale_model_dlogP_prepared;
	return model;
}

Model* new_CTMCScaleModel_from_json(json_node* node, Hashtable* hash){
	char* id = get_json_node_value_string(node, "id");
	char* model_key = get_json_node_value_string(node, "tree");
    
    json_node* x_node = get_json_node(node, "x");
    Parameters* x = distmodel_get_x(id, x_node, hash);
	
	Model *mtree = Hashtable_get(hash, model_key+1);
	DistributionModel* dm = new_CTMCScale_with_parameters(x, mtree->obj);
	Model* model = new_CTMCScaleModel(id, dm, mtree);
#ifndef GSL_DISABLED
	dm->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
#endif

	free_Parameters(x);

	return model;
}
