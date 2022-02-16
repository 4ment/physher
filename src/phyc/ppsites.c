//
//  ppsites.c
//  physher
//
//  Created by Mathieu Fourment on 28/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "ppsites.h"

#include "treelikelihood.h"
#include "matrix.h"

// partials:    number of patterns * number of categpries * number of states
// frequencies: number of states
// posteriors:  [number of patterns x number of categories]
double ** SingleTreeLikelihood_posterior_sites( SingleTreeLikelihood *tlk ){
	
	tlk->calculate(tlk);
	const double *partials = tlk->partials[tlk->current_partials_indexes[Node_id(Tree_root(tlk->tree))]][Node_id(Tree_root(tlk->tree))];
	const double *frequencies = tlk->get_root_frequencies(tlk);
	const double *root_lk   = tlk->pattern_lk; // P[D_i]
	const double *rate_prop = tlk->sm->get_proportions(tlk->sm); // P[R_j]
	
	double **posteriors = dmatrix(tlk->sp->count, tlk->sm->cat_count);
	int v = 0;
	int i = 0;
	
	// P[ D_i| r_i=r_j ]
	for ( int j = 0; j < tlk->sm->cat_count; j++) {
		for ( i = 0; i < tlk->sp->count; i++) {
			// sum pi_x P[ D_i| r_i=r_j, D_i^r=x ]
			for ( int x = 0; x < tlk->m->nstate; x++ ) {
				posteriors[i][j] += frequencies[x] * partials[v];
				v++;
			}
			// P[ r_i=r_j| D_i ] = P[ D_i| r_i=r_j ] * P[ r_i=r_j ] / P[D_i]
			posteriors[i][j] = posteriors[i][j] *rate_prop[j]/exp(root_lk[i]);
		}
	}
	
	return posteriors;
}

void posteriors_sites_calculator_from_json(json_node* node, Hashtable* hash){
	char* allowed [] = {
		"file",
		"model",
		"verbosity"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* ref = get_json_node_value_string(node, "model");
	Model* mtlk = Hashtable_get(hash, ref+1);
	const char *file = get_json_node_value_string(node, "file");
	int verbosity = get_json_node_value_int(node, "verbosity", 0);
	SingleTreeLikelihood *tlk = mtlk->data;
	
	int i,s,c,best_index;
	double conditionalMean,best_prob;
	FILE *pfile = fopen(file, "w");
	
	for(int i = 0; i < Parameters_count(tlk->sm->rates); i++){
		if(strcmp(Parameters_name(tlk->sm->rates, i), "sitemodel.alpha") == 0){
			fprintf(pfile, "Alpha %e\n", Parameters_value(tlk->sm->rates, i));
		}
		else if(strcmp(Parameters_name(tlk->sm->rates, i), "sitemodel.pinv") == 0){
			fprintf(pfile, "Pinv %e\n", Parameters_value(tlk->sm->rates, i));
		}
		else{
			fprintf(pfile, "%s %e\n", Parameters_name(tlk->sm->rates, i), Parameters_value(tlk->sm->rates, i));
		}
	}
	
	fprintf(pfile, "Proportions\n");
	for ( i = 0; i < tlk->sm->cat_count; i++ ) {
		fprintf(pfile, "%f ", tlk->sm->get_proportion(tlk->sm, i ));
	}
	fprintf(pfile, "\nRates\n");
	for ( i = 0; i < tlk->sm->cat_count; i++ ) {
		fprintf(pfile, "%e ",tlk->sm->get_rate(tlk->sm, i));
	}
	fprintf(pfile, "\n");
	
	
	double **posteriors = SingleTreeLikelihood_posterior_sites(tlk);
	for ( s = 0; s < tlk->sp->nsites; s++ ) {
		i = tlk->sp->indexes[s];
		fprintf(pfile, "%d", s);
		conditionalMean = 0.0;
		
		best_index = 0;
		best_prob = posteriors[i][0];
		for ( c = 0; c < tlk->sm->cat_count; c++ ) {
			if( posteriors[i][c] > best_prob ){
				best_index = c;
				best_prob = posteriors[i][c];
			}
			fprintf(pfile, ",%e", posteriors[i][c]);
			conditionalMean += posteriors[i][c] * tlk->sm->get_rate(tlk->sm, c);
		}
		fprintf(pfile, ",%e,%d\n", conditionalMean, best_index);
	}
	
	fprintf(pfile, "\n");
	
	free_dmatrix(posteriors, tlk->sp->count);
	fclose(pfile);
}
