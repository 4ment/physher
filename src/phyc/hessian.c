//
//  hessian.c
//  physher
//
//  Created by Mathieu Fourment on 5/12/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "hessian.h"

#include "matrix.h"


void calculate_hessian(Hessian* hessian){
	size_t dim = Parameters_count(hessian->parameters);
	
	for (int i = 0; i < dim; i++) {
		Parameter* p1 = Parameters_at(hessian->parameters, i);
		hessian->hessian[i*dim+i] = hessian->hessian[i*dim+i] = Model_second_derivative(hessian->likelihood, p1, NULL, 0.0001);
		for (int j = 0; j < i; j++) {
			Parameter* p2 = Parameters_at(hessian->parameters, j);
			hessian->hessian[i*dim+j] = hessian->hessian[j*dim+i] = Model_mixed_derivative(hessian->likelihood, p1, p2);
		}
	}
}

void print_hessian(Hessian* hessian){
	size_t dim = Parameters_count(hessian->parameters);
	for (int i = 0; i < dim; i++) {
		fprintf(stdout, ",%s", Parameters_name(hessian->parameters, i));
	}
	fprintf(stdout, "\n");
	
	for (int i = 0; i < dim; i++) {
		fprintf(stdout, "%s", Parameters_name(hessian->parameters, i));
		for (int j = 0; j < dim; j++) {
			fprintf(stdout, ",%e", hessian->hessian[i*dim+j]);
		}
		fprintf(stdout, "\n");
	}
}

void _free_Hessian(Hessian* hessian){
	free(hessian->hessian);
	free_Parameters(hessian->parameters);
	hessian->likelihood->free(hessian->likelihood);
	free(hessian);
}

Hessian* new_Hessian_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"model",
		"parameters"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* ref = get_json_node_value_string(node, "model");
	Hessian* hessian = malloc(sizeof(Hessian));
	hessian->likelihood = Hashtable_get(hash, ref+1);
	hessian->likelihood->ref_count++;
	hessian->parameters = new_Parameters(1);
	get_parameters_references(node, hash, hessian->parameters);
	hessian->hessian = dvector(Parameters_count(hessian->parameters));
	hessian->calculate = calculate_hessian;
	hessian->free = _free_Hessian;
	return hessian;
}

