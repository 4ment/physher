//
//  vb.c
//  physher
//
//  Created by Mathieu Fourment on 24/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "vb.h"

#include <string.h>
#include <strings.h>

#include "optimizer.h"
#include "tree.h"
#include "utils.h"
#include "transforms.h"
#include "matrix.h"
#include "filereader.h"
#include "distmodel.h"

#include "gamvi.h"
#include "weibullvi.h"

#include "klpq.h"
#include "klqp.h"

variational_block_t* clone_VariationalBlock(variational_block_t* block, Hashtable* hash){
    variational_block_t* clone = malloc(sizeof(variational_block_t));
    
    clone->posterior = NULL;
    clone->simplex_count = 0;
    clone->simplex_parameter_count = 0;
    clone->derivative_type = block->derivative_type;
    clone->numerical_eps = block->numerical_eps;
    
    if(block->simplices != NULL){
        clone->simplices = malloc(block->simplex_count*sizeof(Model*));
        clone->simplex_count = block->simplex_count;
        clone->simplex_parameter_count = block->simplex_parameter_count;
        
        for (int i = 0; i < block->simplex_count; i++) {
            if (Hashtable_exists(hash, block->simplices[i]->name)) {
                clone->simplices[i] = Hashtable_get(hash, block->simplices[i]->name);
                clone->simplices[i]->ref_count++; // it is decremented at the end using free
            }
            else{
                clone->simplices[i] = block->simplices[i]->clone(block->simplices[i], hash);
                Hashtable_add(hash, clone->simplices[i]->name, clone->simplices[i]);
            }
            Parameters_add_parameters(clone->parameters, ((Simplex*)clone->simplices[i])->parameters);
        }
    }
    
    clone->parameters = new_Parameters(Parameters_count(block->parameters));
    for (int i = 0; i < Parameters_count(block->parameters); i++) {
        char* name = Parameters_name(block->parameters, i);
        if (Hashtable_exists(hash, name)) {
            Parameters_add(clone->parameters, Hashtable_get(hash, name));
        }
        else{
            Parameter* p = clone_Parameter(Parameters_at(block->parameters, i));
            Parameters_move(clone->parameters, p);
            Hashtable_add(hash, name, p);
        }
    }
    
    
    clone->var_parameters = malloc(block->var_parameters_count*sizeof(Parameters*));
    clone->var_parameters_count = block->var_parameters_count;
    for(size_t j = 0; j < clone->var_parameters_count; j++){
        clone->var_parameters[j] = new_Parameters(Parameters_count(block->var_parameters[j]));
        Parameters_set_name2(clone->var_parameters[j], Parameters_name2(block->var_parameters[j]));
        Hashtable_add(hash, Parameters_name2(clone->var_parameters[j]), clone->var_parameters[j]);
        
        for (int i = 0; i < Parameters_count(block->var_parameters[j]); i++) {
            char* name = Parameters_name(block->var_parameters[j], i);
            if (Hashtable_exists(hash, name)) {
                Parameters_add(clone->var_parameters[j], Hashtable_get(hash, name));
            }
            else{
                Parameter* p = clone_Parameter(Parameters_at(block->var_parameters[j], i));
                Parameters_move(clone->var_parameters[j], p);
                Hashtable_add(hash, name, p);
            }
        }
    }
    
    clone->sample1 = block->sample1;
    clone->sample2 = block->sample2;
    clone->sample = block->sample;
    clone->entropy = block->entropy;
    clone->grad_elbo = block->grad_elbo;
    clone->grad_entropy = block->grad_entropy;
    clone->logP = block->logP;
    clone->logQ = block->logQ;
    clone->use_entropy = block->use_entropy;
    clone->rng = block->rng;
    clone->initialized = block->initialized;
    clone->etas = clone_Vector(block->etas);
    return clone;
}

variational_block_t* new_VariationalBlock_from_json(json_node* node, Hashtable* hash){
    variational_block_t* var = malloc(sizeof(variational_block_t));
    char* allowed[] = {
        "derivative",
        "distribution",
        "initialize",
        "entropy", // use entropy or not
        "eps",
        "parameters",
        "simplices",
        "tree",
        "x"
    };
    json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
    
    const char* id = get_json_node_value_string(node, "id");
    const char* derivative_string = get_json_node_value_string(node, "derivative");
    const char* dist_string = get_json_node_value_string(node, "distribution");
    json_node* parameters_node = get_json_node(node, "parameters");
    json_node* simplices_node = get_json_node(node, "simplices");
    json_node* tree_node = get_json_node(node, "tree");
    var->use_entropy = get_json_node_value_bool(node, "entropy", true);

    var->derivative_type = DERIVATIVE_ANALYTICAL;
    var->numerical_eps = -INFINITY;
    if(derivative_string != NULL && strcasecmp(derivative_string, "numerical") == 0){
        var->derivative_type = DERIVATIVE_NUMERICAL;
        var->numerical_eps = get_json_node_value_double(node, "eps", 1.e-8);
    }

    var->sample = NULL;
    var->var_parameters = NULL;
	var->parameters = new_Parameters(2);

    var->simplex_count = 0;
    var->simplex_parameter_count = 0;
    var->simplices = NULL;
    
    var->initialized = false;
    var->initialize = NULL;
    var->rng = Hashtable_get(hash, "RANDOM_GENERATOR!@");
    
    if(simplices_node != NULL){
        if (simplices_node->node_type == MJSON_ARRAY) {
            var->simplices = malloc(simplices_node->child_count*sizeof(Model*));
            for (int i = 0; i < simplices_node->child_count; i++) {
                json_node* simplex_node = simplices_node->children[i];
                char* ref = simplex_node->value;
                Model* msimplex = Hashtable_get(hash, ref+1);
                Simplex* simplex = msimplex->obj;
                Parameters_add_parameters(var->parameters, simplex->parameters);
                var->simplices[i] = msimplex;
                msimplex->ref_count++;
                var->simplex_count++;
                var->simplex_parameter_count += simplex->K-1;
            }
        }
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

    var->posterior = NULL; // set by variational_t
    json_node* x_node = get_json_node(node, "x");
	if(x_node != NULL){
		Parameters* parameters = distmodel_get_x(id, x_node, hash);
		Parameters_add_parameters(var->parameters, parameters);
		free_Parameters(parameters);
	}
    
    var->var_parameters = malloc(parameters_node->child_count*sizeof(Parameters*));
    var->var_parameters_count = parameters_node->child_count;
    for (int i = 0; i < parameters_node->child_count; i++) {
        json_node* p_node = parameters_node->children[i];
        json_node* dim_node = get_json_node(p_node, "dimension");
        Parameters* multi_parameter = NULL;
        if(dim_node != NULL){
            multi_parameter = new_MultiParameter_from_json(p_node, hash);
            Hashtable_add(hash, Parameters_name2(multi_parameter), multi_parameter);
        }
        else{
            multi_parameter = new_Parameters(1);
            Parameter* parameter = new_Parameter_from_json(p_node, hash);
            Parameters_move(multi_parameter, parameter);
            Parameters_set_name2(multi_parameter, Parameter_name(parameter));
            Hashtable_add(hash, Parameter_name(parameter), parameter);
        }
        var->var_parameters[i] = multi_parameter;
    }
    
    if (dist_string == NULL || strcasecmp(dist_string, "normal") == 0) {
        var->etas = new_Vector(Parameters_count(var->parameters));
        Vector_resize(var->etas, Parameters_count(var->parameters));
        size_t mu_length = Parameters_count(var->var_parameters[0]);
        size_t sigma_length = Parameters_count(var->var_parameters[1]);
        if (var->var_parameters_count != 2 || mu_length != sigma_length) {
            fprintf(stderr, "Meanfield normal should have 2 parameter vectors of same length: mu (length: %zu) and sigma (length: %zu)\n", mu_length, sigma_length);
            exit(2);
        }
        if (Parameters_count(var->var_parameters[0]) != Parameters_count(var->parameters)) {
            fprintf(stderr, "The length of mu and sigma (%zu) in meanfield model should be equal to the length of input parameters x (%zu)\n", mu_length, Parameters_count(var->parameters));
            exit(2);
        }
        // check order: mu should be first, then sigma
        if (strcasecmp(parameters_node->children[0]->key, "mu") != 0) {
            Parameters* temp = var->var_parameters[0];
            var->var_parameters[0] = var->var_parameters[1];
            var->var_parameters[1] = temp;
        }
        bool backward = get_json_node_value_bool(node, "backward", true); // true is KL(Q||P)
        if(backward){
            var->entropy = klqp_block_meanfield_normal_entropy;
            var->grad_elbo = klqp_block_meanfield_normal_grad_elbo;
            var->grad_entropy = klqp_block_meanfield_normal_grad_entropy;
            var->sample = klqp_block_meanfield_normal_sample;
            var->sample1 = klqp_block_meanfield_normal_sample1;
            var->sample2 = klqp_block_meanfield_normal_sample2;
//            var->sample_some = klqp_block_meanfield_normal_sample_some;
//            var->finalize = klqp_meanfield_normal_finalize;
            var->logP = klqp_block_meanfield_normal_logP;
            var->logQ = klqp_block_meanfield_normal_logQ;
            var->initialize = klqp_block_meanfield_normal_initialize;
//            var->parameters_logP = klqp_block_meanfield_normal_logP_parameters;
//            var->print = klqp_meanfield_normal_log_samples;
            
        }
        else{
//                var->elbofn = klpq_normal_meanfield;
//                var->grad_elbofn = grad_klpq_normal_meanfield;
        }
    }
    else if (strcasecmp(dist_string, "multivariatenormal") == 0) {
        var->etas = new_Vector(Parameters_count(var->parameters)*2); // etas and zetas
        Vector_resize(var->etas, Parameters_count(var->parameters)*2);
        size_t mu_length = Parameters_count(var->var_parameters[0]);
        size_t sigma_length = Parameters_count(var->var_parameters[1]);
        if (var->var_parameters_count != 2 || sigma_length != (mu_length*(mu_length-1))/2+mu_length) {
            fprintf(stderr, "Fullrank normal should have 2 parameter vectors: mu (length: %zu) and sigma (length: %zu)\n", mu_length, sigma_length);
            exit(2);
        }
        // check order: mu should be first, then sigma
        if (strcasecmp(parameters_node->children[0]->key, "mu") != 0) {
            Parameters* temp = var->var_parameters[0];
            var->var_parameters[0] = var->var_parameters[1];
            var->var_parameters[1] = temp;
        }
    
        if (Parameters_count(var->var_parameters[0]) != Parameters_count(var->parameters)) {
            fprintf(stderr, "The length of mu (%zu) in fullrank model should be equal to the length of input parameters x (%zu)\n", mu_length, Parameters_count(var->parameters));
            exit(2);
        }
    
        bool backward = get_json_node_value_bool(node, "backward", true); // true is KL(Q||P)
        if(backward){
            var->entropy = klqp_block_fullrank_normal_entropy;
            var->grad_elbo = klqp_block_fullrank_normal_grad_elbo;
            var->grad_entropy = klqp_block_fullrank_normal_grad_entropy;
            var->sample = klqp_block_fullrank_normal_sample;
            var->sample1 = klqp_block_fullrank_normal_sample1;
            var->sample2 = klqp_block_fullrank_normal_sample2;
//            var->sample_some = klqp_block_meanfield_normal_sample_some;
//            var->finalize = klqp_meanfield_normal_finalize;
            var->logP = klqp_block_fullrank_normal_logP;
            var->logQ = klqp_block_fullrank_normal_logQ;
//            var->parameters_logP = klqp_block_meanfield_normal_logP_parameters;
//            var->print = klqp_meanfield_normal_log_samples;
            
        }
        else{

        }
    }
	else if (strcasecmp(dist_string, "gamma") == 0) {
			var->etas = new_Vector(Parameters_count(var->parameters));
			Vector_resize(var->etas, Parameters_count(var->parameters));
			size_t shape_length = Parameters_count(var->var_parameters[0]);
			size_t rate_length = Parameters_count(var->var_parameters[1]);
			if (var->var_parameters_count != 2 || shape_length != rate_length) {
				fprintf(stderr, "Meanfield gamma should have 2 parameter vectors of same length: shape (length: %zu) and rate (length: %zu)\n", shape_length, rate_length);
				exit(2);
			}
			if (Parameters_count(var->var_parameters[0]) != Parameters_count(var->parameters)) {
				fprintf(stderr, "The length of shape and rate (%zu) in meanfield model should be equal to the length of input parameters x (%zu)\n", shape_length, Parameters_count(var->parameters));
				exit(2);
			}
			// check order: shape should be first, then rate
			if (strcasecmp(parameters_node->children[0]->key, "shape") != 0) {
				Parameters* temp = var->var_parameters[0];
				var->var_parameters[0] = var->var_parameters[1];
				var->var_parameters[1] = temp;
			}
			bool backward = get_json_node_value_bool(node, "backward", true); // true is KL(Q||P)
			if(backward){
				var->entropy = klqp_block_meanfield_gamma_entropy;
				var->grad_elbo = klqp_block_meanfield_gamma_grad_elbo;
				var->grad_entropy = klqp_block_meanfield_gamma_grad_entropy;
				var->sample = klqp_block_meanfield_gamma_sample;
				var->sample1 = klqp_block_meanfield_gamma_sample1;
				var->sample2 = klqp_block_meanfield_gamma_sample2;
	//            var->sample_some = klqp_block_meanfield_normal_sample_some;
	//            var->finalize = klqp_meanfield_normal_finalize;
				var->logP = klqp_block_meanfield_gamma_logP;
				var->logQ = klqp_block_meanfield_gamma_logQ;
	//            var->parameters_logP = klqp_block_meanfield_normal_logP_parameters;
	//            var->print = klqp_meanfield_normal_log_samples;
			}
			else{
	//                var->elbofn = klpq_normal_meanfield;
	//                var->grad_elbofn = grad_klpq_normal_meanfield;
			}
		}
	else if (strcasecmp(dist_string, "weibull") == 0) {
		var->etas = new_Vector(Parameters_count(var->parameters));
		Vector_resize(var->etas, Parameters_count(var->parameters));
		size_t scale_length = Parameters_count(var->var_parameters[0]);
		size_t shape_length = Parameters_count(var->var_parameters[1]);
		if (var->var_parameters_count != 2 || scale_length != shape_length) {
			fprintf(stderr, "Meanfield weibull should have 2 parameter vectors of same length: scale (length: %zu) and shape (length: %zu)\n", scale_length, shape_length);
			exit(2);
		}
		if (Parameters_count(var->var_parameters[0]) != Parameters_count(var->parameters)) {
			fprintf(stderr, "The length of scale and shape (%zu) in meanfield model should be equal to the length of input parameters x (%zu)\n", scale_length, Parameters_count(var->parameters));
			exit(2);
		}
		// check order: scale should be first, then shape
		if (strcasecmp(parameters_node->children[0]->key, "scale") != 0) {
			Parameters* temp = var->var_parameters[0];
			var->var_parameters[0] = var->var_parameters[1];
			var->var_parameters[1] = temp;
		}
		bool backward = get_json_node_value_bool(node, "backward", true); // true is KL(Q||P)
		if(backward){
			var->entropy = klqp_block_meanfield_weibull_entropy;
			var->grad_elbo = klqp_block_meanfield_weibull_grad_elbo;
			var->grad_entropy = klqp_block_meanfield_weibull_grad_entropy;
			var->sample = klqp_block_meanfield_weibull_sample;
			var->sample1 = klqp_block_meanfield_weibull_sample1;
			var->sample2 = klqp_block_meanfield_weibull_sample2;
//            var->sample_some = klqp_block_meanfield_normal_sample_some;
//            var->finalize = klqp_meanfield_normal_finalize;
			var->logP = klqp_block_meanfield_weibull_logP;
			var->logQ = klqp_block_meanfield_weibull_logQ;
//            var->parameters_logP = klqp_block_meanfield_normal_logP_parameters;
//            var->print = klqp_meanfield_normal_log_samples;
		}
		else{
//                var->elbofn = klpq_normal_meanfield;
//                var->grad_elbofn = grad_klpq_normal_meanfield;
		}
	}
    else{
        fprintf(stderr, "distribution %s not recognized\n", dist_string);
        exit(1);
    }
    bool initialize = get_json_node_value_bool(node, "initialize", false);
    if(initialize && var->initialize != NULL){
        var->initialize(var);
        var->initialized = true;
    }

    return var;
}

double variational_logP(variational_t* var, const double* values){
    double logP = 0;
    size_t shift = 0;
    for(int j = 0; j < var->block_count; j++){
        variational_block_t* block = var->blocks[j];
        logP += block->logP(block, values+shift);
        for(int k = 0; k < block->var_parameters_count; k++){
            shift += Parameters_count(block->var_parameters[k]);
        }
    }
    return logP;
}

double variational_logP_parameters(variational_t* var, const Parameters* parameters){
    double logP = 0;
    fprintf(stderr, "variational_logP_parameters no implemented \n");
    exit(3);
    return logP;
}

bool variational_sample(variational_t* var, double* values){
    if(!var->ready_to_sample){
        var->finalize(var);
        var->ready_to_sample = true;
    }
    
    size_t shift = 0;
    for(int j = 0; j < var->block_count; j++){
        variational_block_t* block = var->blocks[j];
        block->sample(block, values+shift);
        for(int k = 0; k < block->var_parameters_count; k++){
            shift += Parameters_count(block->var_parameters[k]);
        }
    }
    return true;
}

void _variational_jsonize(Model* model, json_node* jroot){
    variational_t* var = model->obj;
    for (int i = 0; i < var->block_count; i++) {
        variational_block_t* block = var->blocks[i];
        for (int j = 0; j < block->var_parameters_count; j++) {
            Parameters_to_json(block->var_parameters[j], jroot);
        }
    }
}

void variational_load(variational_t* var, const char* filename){
    char* content = load_file(filename);
    json_node* json = create_json_tree(content);
    free(content);
    for(size_t i = 0; i < var->block_count; i++){
        variational_block_t* block = var->blocks[i];
        for(size_t j = 0; j < block->var_parameters_count; j++){
            const char* var_id = Parameters_name2(block->var_parameters[j]);
            json_node* parameters_node = get_json_node(json, var_id);
            if(parameters_node != NULL){
                json_node* values = get_json_node(parameters_node, "values");
                if(values->child_count != Parameters_count(block->var_parameters[j])){
                    fprintf(stderr, "Could not load values of parameters with ID: %s (%zu %zu)", var_id, values->child_count, Parameters_count(block->var_parameters[j]));
                    continue;
                }
                for(size_t k = 0; k < values->child_count; k++){
                    double value = atof((char*)values->children[k]->value);
                    Parameters_set_value(block->var_parameters[j], k, value);
                }
            }
        }
    }
    json_free_tree(json);
}

Model* new_Variational_from_json2(json_node* node, Hashtable* hash){
    char* allowed[] = {
        "distributions",
        "divergence",
        "elbosamples",
        "elbomulti",
        "gradsamples",
        "init",
        "log",
        "log_samples",
        "parameters",
        "posterior",
        "simplices",
        "tree",
        "var",
        "var_parameters"
    };
    json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
    
    const char* posterior_string = get_json_node_value_string(node, "posterior");
    const char* id = get_json_node_value_string(node, "id");
    const char* filename = get_json_node_value_string(node, "log");
    const char* divergence = get_json_node_value_string(node, "divergence");
    json_node* blocks_node = get_json_node(node, "distributions");
    json_node* init_node = get_json_node(node, "init");
    
    variational_t* var = malloc(sizeof(variational_t));
    
    var->posterior = Hashtable_get(hash, posterior_string+1);
    var->posterior->ref_count++;
    
    var->file = NULL;
    var->sample = NULL;
    var->sample_some = NULL;
    var->parameters = NULL;
    var->var_parameters = NULL;
    
    var->simplex_count = 0;
    var->simplices = NULL;
    
    var->block_count = blocks_node->child_count;
    var->blocks = malloc(blocks_node->child_count*sizeof(variational_block_t*));
    for (int i = 0; i < var->block_count; i++) {
        json_node* block_node = blocks_node->children[i];
        variational_block_t* block = new_VariationalBlock_from_json(block_node, hash);
        block->posterior = var->posterior;
        var->blocks[i] = block;
    }
    
    var->elbo_samples = get_json_node_value_size_t(node, "elbosamples", 100);
    var->grad_samples = get_json_node_value_size_t(node, "gradsamples", 1);
    var->elbo_multi = get_json_node_value_int(node, "elbomulti", 1);
    var->log_samples = get_json_node_value_int(node, "log_samples", 0);
    
    if(divergence == NULL || strcasecmp(divergence, "klqp") == 0){
        if (var->elbo_multi > 1) {
            var->elbofn = variational_klqp_elbo_multi;
        }
        else{
            var->elbofn = variational_klqp_elbo;
        }
        var->grad_elbofn = variational_klqp_grad_elbo;
    }
    var->sample = variational_sample;
    var->logP = variational_logP;
    var->print = klqp_fullrank_log_samples;
    
    var->initialized = false;
    var->ready_to_sample = false;
    var->iter = 0;
    StringBuffer* buffer = new_StringBuffer(10);
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
    
    if (init_node != NULL) {
        variational_load(var, init_node->value);
    }
	
	// Collect every parameter and inform posterior about gradient calculation
	Parameters* x = new_Parameters(10);
	for (size_t i = 0; i < var->block_count; i++) {
		variational_block_t* block = var->blocks[i];
		Parameters_add_parameters(x, block->parameters);
	}
	var->posterior->prepare_gradient(var->posterior, x);
	free_Parameters(x);
    
    return new_VariationalModel(id, var);
}

Model* new_Variational_from_json(json_node* node, Hashtable* hash){
    json_node* blocks_node = get_json_node(node, "distributions");
    if (blocks_node != NULL) {
        return new_Variational_from_json2(node, hash);
    }
    
	char* allowed[] = {
		"backward",
		"distribution",
		"elbosamples",
		"elbomulti",
		"gradsamples",
		"log",
		"log_samples",
		"parameters",
		"posterior",
		"simplices",
		"tree",
		"var",
		"var_parameters"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	const char* posterior_string = get_json_node_value_string(node, "posterior");
	const char* var_string = get_json_node_value_string(node, "var");
	const char* dist_string = get_json_node_value_string(node, "distribution");
	const char* id = get_json_node_value_string(node, "id");
	const char* filename = get_json_node_value_string(node, "log");
	json_node* parameters_node = get_json_node(node, "parameters");
	json_node* tree_node = get_json_node(node, "tree");
	json_node* simplices_node = get_json_node(node, "simplices");
	
	variational_t* var = malloc(sizeof(variational_t));
    var->blocks = NULL;
    var->block_count = 0;
	var->file = NULL;
	var->sample = NULL;
	var->sample_some = NULL;
	var->parameters = new_Parameters(1);
	var->var_parameters = NULL;
	
	var->simplex_count = 0;
	var->simplices = NULL;
	if(simplices_node != NULL){
		if (simplices_node->node_type == MJSON_ARRAY) {
			var->simplices = malloc(simplices_node->child_count*sizeof(Model*));
			for (int i = 0; i < simplices_node->child_count; i++) {
				json_node* simplex_node = simplices_node->children[i];
				char* ref = simplex_node->value;
				Model* msimplex = Hashtable_get(hash, ref+1);
				Simplex* simplex = msimplex->obj;
				Parameters_add_parameters(var->parameters, simplex->parameters);
				var->simplices[i] = msimplex;
				msimplex->ref_count++;
				var->simplex_count++;
			}
		}
	}
	
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
	var->elbo_multi = get_json_node_value_int(node, "elbomulti", 1);
	var->log_samples = get_json_node_value_int(node, "log_samples", 0);
	
	json_node* var_parameters_node = get_json_node(node, "var_parameters");
	// Variational model parameters are provided
	if (var_parameters_node->node_type == MJSON_ARRAY) {
		var->var_parameters = new_Parameters(var_parameters_node->child_count);
		for(size_t i = 0; i < Parameters_count(var->var_parameters); i++){
			json_node* var_p_node = var_parameters_node->children[i];
			Parameter* var_p = new_Parameter_from_json(var_p_node, hash);
			Parameters_move(var->var_parameters, var_p);
		}
	}
	
	StringBuffer* buffer = new_StringBuffer(10);
	if(strcasecmp(var_string, "meanfield") == 0){
//		printf("Creating meanfield variational model\n");
        size_t param_count = 2;
        char** params = malloc(param_count*sizeof(char*));
        
        if (strcasecmp(dist_string, "normal") == 0 || strcasecmp(dist_string, "lognormal") == 0) {
            params[0] = String_clone(".mu");
            params[1] = String_clone(".sigma");
        }
        else if(strcasecmp(dist_string, "gamma") == 0){
            params[0] = String_clone(".alpha");
            params[1] = String_clone(".beta");
        }
        else{
            fprintf(stderr, "meanfield distribution should be: normal, lognormal, or gamma");
            exit(1);
        }
        if(var->var_parameters == NULL){
            var->var_parameters = new_Parameters(dim*2); // mean + sd
            for(size_t i = 0; i < dim; i++){
                StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
                StringBuffer_append_string(buffer, params[0]);
                Parameter* p = new_Parameter(buffer->c, 0, NULL);
                p->id = (int)i;
                Parameters_move(var->var_parameters, p);
            }
            for(size_t i = 0; i < dim; i++){
                StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
                StringBuffer_append_string(buffer, params[1]);
                Parameter* p = new_Parameter(buffer->c, 0, NULL);
                p->id = (int)(i+dim);
                Parameters_move(var->var_parameters, p);
            }
        }
        free_cmatrix(params, param_count);
        
		if (strcasecmp(dist_string, "normal") == 0) {
			bool backward = get_json_node_value_bool(node, "backward", true); // true is KL(Q||P)
			if(backward){
				var->elbofn = klqp_meanfield_normal_elbo;
				if (var->elbo_multi > 1) {
					var->elbofn = klqp_meanfield_normal_elbo_multi;
				}
				var->grad_elbofn = klqp_meanfield_normal_grad_elbo;
			}
			else{
				var->elbofn = klpq_normal_meanfield;
				var->grad_elbofn = grad_klpq_normal_meanfield;
			}
			var->sample = klqp_meanfield_normal_sample;
			var->sample_some = klqp_meanfield_normal_sample_some;
			var->finalize = klqp_meanfield_normal_finalize;
			var->logP = klqp_meanfield_normal_logP;
            var->print = klqp_meanfield_normal_log_samples;
		}
        else if (strcasecmp(dist_string, "gamma") == 0) {
			var->elbofn = elbo_gamma_meanfield;
			var->grad_elbofn = grad_elbo_gamma_meanfield;
			
			var->sample = variational_sample_gamma_meanfield;
			var->sample_some = variational_sample_some_gamma_meanfield;
			var->finalize = NULL; //TODO: implement function
			var->logP = variational_gamma_meanfield_logP;
            var->print = meanfield_gamma_log_samples;
		}
        else if (strcasecmp(dist_string, "lognormal") == 0) {
            var->elbofn = klqp_meanfield_lognormal_elbo;
            var->grad_elbofn = klqp_meanfield_lognormal_grad_elbo;
            var->sample = klqp_meanfield_lognormal_sample;
            var->sample_some = klqp_meanfield_lognormal_sample_some;
            //TODO: implement lognormal finalize
            var->finalize = klqp_meanfield_normal_finalize;
            var->logP = klqp_meanfield_lognormal_logP;
            var->print = klqp_meanfield_lognormal_log_samples;
        }
		else{
			fprintf(stderr, "distribution %s not recognized\n", dist_string);
			exit(1);
		}
	}
	else if(strcasecmp(var_string, "fullrank") == 0){
//		printf("Creating fullrank variational model\n");
		if(var->var_parameters == NULL){
			size_t n = 2*dim+(dim*dim-dim)/2; // mus, vars, covs
			var->var_parameters = new_Parameters(n);
			
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
						Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
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
					Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
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
						Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
						Parameters_move(var->var_parameters, p);
						row++;
					}
					StringBuffer_set_string(buffer, Parameters_name(var->parameters, i));
					StringBuffer_append_string(buffer, ".var");
					Parameter* p = new_Parameter(buffer->c, 0.1, NULL);
					Parameters_move(var->var_parameters, p);
		//			p->estimate = false;
					row++;
				}
			}
		}
		
		var->elbofn = klqp_fullrank_normal_elbo;
		var->grad_elbofn = klqp_fullrank_normal_grad_elbo;
		var->sample = klqp_fullrank_normal_sample;
		var->logP = klqp_fullrank_normal_logP;
        var->print = klqp_fullrank_log_samples;
	}
	
	// variational parameters were not provided 
	if (var_parameters_node->node_type == MJSON_STRING) {
		Parameters_set_name2(var->var_parameters, get_json_node_value_string(node, "var_parameters"));
		Hashtable_add(hash, Parameters_name2(var->var_parameters), var->var_parameters);
	}
	
	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
		Hashtable_add(hash, Parameters_name(var->var_parameters, i), Parameters_at(var->var_parameters, i));
	}
	
	var->initialized = false;
	var->ready_to_sample = false;
	var->iter = 0;
	
	// Inform posterior about gradient calculation
	var->posterior->prepare_gradient(var->posterior, var->parameters);
	
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


void free_Variational(variational_t* var){
	var->posterior->free(var->posterior);
	free_Parameters(var->var_parameters);
	free_Parameters(var->parameters);
	if(var->simplices != NULL){
		for(int i = 0; i < var->simplex_count; i++){
			var->simplices[i]->free(var->simplices[i]);
		}
		free(var->simplices);
	}
    if(var->blocks != NULL){
        for (int i = 0; i < var->block_count; i++) {
            variational_block_t* block = var->blocks[i];
            for (int j = 0; j < block->var_parameters_count; j++) {
                free_Parameters(block->var_parameters[j]);
            }
			if(block->simplices != NULL){
				for (int j = 0; j < block->simplex_count; j++) {
					block->simplices[j]->free(block->simplices[j]);
				}
				free(block->simplices);
			}
            free(block->var_parameters);
            free_Parameters(block->parameters);
            if(block->etas != NULL) free_Vector(block->etas);
            free(block);
        }
//        var->posterior->free(var->posterior);
        free(var->blocks);
    }
	if (var->file != NULL) fclose(var->file);
	free(var);
}

static void _variational_model_free( Model *self ){
	if(self->ref_count == 1){
//		printf("Free variational model %s\n", self->name);
		variational_t* var = (variational_t*)self->obj;
		free_Variational(var);
		free_Model(self);
	}
	else{
		self->ref_count--;
	}
}

static void _variational_model_handle_change( Model *self, Model *model, int index ){
	variational_t* var = (variational_t*)self->obj;
	if(index >= Parameters_count(var->parameters)){
		index -= Parameters_count(var->parameters);
	}
//	self->listeners->fire( self->listeners, self, index );
}

static Model* _variational_model_clone(Model* self, Hashtable *hash){
	if (Hashtable_exists(hash, self->name)) {
		return Hashtable_get(hash, self->name);
	}
	variational_t* var = self->obj;
	
	Model* posterior = var->posterior->clone(var->posterior, hash);
	Hashtable_add(hash, posterior->name, posterior);
    
    
    variational_t* clone = malloc(sizeof(variational_t));
    clone->posterior = posterior;
    clone->parameters = NULL;
    clone->var_parameters = NULL;
    clone->blocks = NULL;
    clone->block_count = 0;

    // old
    if(var->parameters != NULL){
        clone->parameters = new_Parameters(Parameters_count(var->parameters));
        for (int i = 0; i < Parameters_count(var->parameters); i++) {
            char* name = Parameters_name(var->parameters, i);
            if (Hashtable_exists(hash, name)) {
                Parameters_add(clone->parameters, Hashtable_get(hash, name));
            }
            else{
                Parameter* p = clone_Parameter(Parameters_at(var->parameters, i));
                Parameters_move(clone->parameters, p);
                Hashtable_add(hash, name, p);
            }
        }
    }
	
    // old
    if(var->var_parameters != NULL){
        clone->var_parameters = new_Parameters(Parameters_count(var->var_parameters));
        if (Parameters_name2(var->var_parameters) != NULL) {
            Parameters_set_name2(clone->var_parameters, Parameters_name2(var->var_parameters));
        }
        Hashtable_add(hash, Parameters_name2(clone->var_parameters), clone->var_parameters);
        
        for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
            char* name = Parameters_name(var->var_parameters, i);
            if (Hashtable_exists(hash, name)) {
                Parameters_add(clone->var_parameters, Hashtable_get(hash, name));
            }
            else{
                Parameter* p = clone_Parameter(Parameters_at(var->var_parameters, i));
                Parameters_move(clone->var_parameters, p);
                Hashtable_add(hash, name, p);
            }
        }
    }
    
    if(var->blocks != 0){
        clone->block_count = var->block_count;
        clone->blocks = malloc(var->block_count*sizeof(variational_block_t*));
        for(size_t i = 0; i < var->block_count; i++){
            clone->blocks[i] = clone_VariationalBlock(var->blocks[i], hash);
            clone->blocks[i]->posterior = posterior;
        }
    }

	clone->file = NULL;
	clone->sample = var->sample;
	clone->sample_some = var->sample_some;
	clone->elbofn = var->elbofn;
	clone->grad_elbofn = var->grad_elbofn;
	clone->elbo_samples = var->elbo_samples;
	clone->grad_samples = var->grad_samples;
	clone->logP = var->logP;
	clone->initialized = var->initialized;
	clone->finalize = var->finalize;
	clone->iter = var->iter;
	clone->ready_to_sample = var->ready_to_sample;
	clone->rng = var->rng;
    clone->elbo_multi = var->elbo_multi;
    
    clone->simplex_count = var->simplex_count;
    clone->simplices = NULL;
    if(var->simplices != NULL){
        clone->simplices = malloc(var->simplex_count*sizeof(Model*));
        for (int i = 0; i < var->simplex_count; i++) {
            if (Hashtable_exists(hash, var->simplices[i]->name)) {
                clone->simplices[i] = Hashtable_get(hash, var->simplices[i]->name);
                clone->simplices[i]->ref_count++; // it is decremented at the end using free
            }
            else{
                clone->simplices[i] = var->simplices[i]->clone(var->simplices[i], hash);
                Hashtable_add(hash, clone->simplices[i]->name, clone->simplices[i]);
            }
            Parameters_add_parameters(clone->parameters, ((Simplex*)clone->simplices[i])->parameters);
        }
    }
    
	Model* mclone = new_VariationalModel(self->name, clone);

	return mclone;
}

static double _variational_model_logP(Model *self){
	variational_t* var = (variational_t*)self->obj;
	return var->elbofn(var);
}

static void _variational_model_gradient(Model *self, const Parameters* parameters, double* grad){
	variational_t* var = (variational_t*)self->obj;
	var->grad_elbofn(var, parameters, grad);
}

void _variational_model_reset(Model* self){
	variational_t* var = (variational_t*)self->obj;
	if(var->initialized){
		// Restore the parameters of the posterior
		for (int i = 0; i < Parameters_count(var->parameters); i++) {
			Parameter_restore(Parameters_at(var->parameters, i));
		}
		var->ready_to_sample = false;
		var->initialized = false;
//		printf("reset restore\n");
	}
//	else{
//		
//		printf("reset store\n");
//	}

}

void _variational_model_sample(Model* self, double* samples, double* logP){
	variational_t* var = (variational_t*)self->obj;
	var->sample(var, samples);
	if(logP != NULL){
		*logP = var->logP(var, samples);
	}
}

static void _variational_print(Model* self, FILE* file){
	variational_t* var = (variational_t*)self->obj;
	if(var->log_samples > 0){
        var->print(var, file);
		return;
	}

	fprintf(file, "muraw=c(");
	size_t i = 0;
	for (; i < Parameters_count(var->parameters); i++) {
		double value = Parameters_value(var->var_parameters, i);
		fprintf(file,"%f%s", value, i ==Parameters_count(var->parameters)-1 ? "":",");
	}
	
	fprintf(file, ")\n\n");
	fprintf(file, "sigma=c(");
	for (; i < Parameters_count(var->var_parameters); i++) {
		fprintf(file, "%f%s", Parameters_value(var->var_parameters, i), i ==Parameters_count(var->var_parameters)-1 ? "":",");
	}
	fprintf(file, ")\n");
	
	fprintf(file, "mu=c(");
	size_t dim = Parameters_count(var->parameters);
	for (i = 0; i < Parameters_count(var->parameters); i++) {
		Parameter* p = Parameters_at(var->parameters, i);
		double sd = Parameters_value(var->var_parameters, i+dim);
		double value = Parameters_value(var->var_parameters, i);
		if(Parameter_lower(p) == 0){
			value = exp(value + sd*sd/2.0);
		}
		fprintf(file,"%f%s", value, i ==Parameters_count(var->parameters)-1 ? "":",");
	}
	fprintf(file, ")\n");
	
//	size_t row = dim;
//	for (int j = 0; j < dim; j++) {
//		double temp = 0;
//		// multiply L_j and eta
//		for (int k = 0; k < j+1; k++) {
//			if(Parameters_estimate(var->var_parameters, row)){
//				temp += Parameters_value(var->var_parameters, row)*eta[k];
//			}
//			row++;
//		}
//		zeta[j] = temp + Parameters_value(var->var_parameters, j); // add mu
//
//		Parameter* p = Parameters_at(var->parameters, j);
//		double theta = inverse_transform2(zeta[j], Parameter_lower(p), Parameter_upper(p));
//		Parameter_set_value(p, theta);
//	}
}

Model* new_VariationalModel(const char* name, variational_t* var){
	Model *model = new_Model(MODEL_VARIATIONAL, name, var);
	model->free = _variational_model_free;
	model->clone = _variational_model_clone;
	model->logP = _variational_model_logP;
	model->gradient = _variational_model_gradient;
	model->reset = _variational_model_reset;
	model->print = _variational_print;
	model->sample = _variational_model_sample;
	model->samplable = true;
    model->jsonize = _variational_jsonize;
	
	for (int i = 0; i < Parameters_count(var->var_parameters); i++) {
		Parameters_at(var->var_parameters, i)->listeners->add( Parameters_at(var->var_parameters, i)->listeners, model );
	}
	return model;
}
