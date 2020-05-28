/*
 *  asr.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 16/07/2014.
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

#include "asr.h"

#include <stdio.h>
#include <strings.h>

#include "treelikelihood.h"
#include "matrix.h"
#include "sequenceio.h"

// Marginal reconstruction of ancestral sequence
void _marginal_reconstruction( SingleTreeLikelihood *tlk ){
	Node** nodes = Tree_get_nodes(tlk->tree, POSTORDER);
	size_t node_count = Tree_node_count(tlk->tree);
	SitePattern *sp = tlk->sp;
	int sp_count = sp->count;
	int nstate = tlk->sm->nstate;
	double *backup = dvector(tlk->partials_size);
	double *Lj = dvector(sp_count);

	double* spare_partials = tlk->partials[0][tlk->upper_partial_indexes[Tree_root(tlk->tree)->id]];
	double* spare_root_partials = tlk->root_partials + sp_count*nstate;
	double* spare_pattern_lk = tlk->pattern_lk + sp_count;
	
	for (size_t k = 0; k < node_count; k++) {
		Node* node = nodes[k];
		if( !Node_isleaf(node) ){
			size_t nodeID = Node_id(node);
			size_t index = tlk->mapping[nodeID];
			double* partials = tlk->partials[tlk->current_partials_indexes[nodeID]][nodeID];
			
			memcpy(backup, partials, sizeof(double)*tlk->partials_size);
			
			for ( int j = 0; j < nstate; j++ ) {
				
				for ( int i = 0; i < tlk->partials_size; i++ ){
					if( (i-j) % nstate != 0 ) partials[i] = 0;
				}
				
				if( Node_isroot(node) ){
					spare_partials = partials;
				}
				else {
					tlk->calculate_branch_likelihood(tlk, spare_partials, tlk->upper_partial_indexes[nodeID], nodeID, nodeID );
				}
				
				if(tlk->sm->integrate){
					tlk->integrate_partials(tlk, spare_partials, tlk->sm->get_proportions(tlk->sm), spare_root_partials );
				}
				else{
					memcpy(spare_root_partials, spare_partials, sizeof(double)*sp_count*nstate);
				}
				tlk->node_log_likelihoods( tlk, spare_root_partials, tlk->get_root_frequencies(tlk), spare_pattern_lk);
				
				memcpy(partials, backup, sizeof(double)*tlk->partials_size);
				
				if ( j == 0 ) {
					for ( int i = 0; i < sp_count; i++ ) {
						sp->patterns[i][index] = 0;
					}
					memcpy(Lj, spare_pattern_lk, sp_count*sizeof(double));
				}
				else {
					for ( int i = 0; i < sp_count; i++ ) {
						if( spare_pattern_lk[i] > Lj[i] ){
							Lj[i] = spare_pattern_lk[i];
							sp->patterns[i][index] = j;
						}
					}
				}
			}
			
			if( sp->names[index] != NULL ){
				if (strcmp(sp->names[index], node->name) != 0) {
					free(sp->names[index]);
					sp->names[index] = String_clone( node->name );
				}
			}
			else{
				sp->names[index] = String_clone( node->name );
			}
		}
	}
	free(backup);
	free(Lj);
}

void asr_marginal( SingleTreeLikelihood *tlk ){
	tlk->calculate(tlk);
	tlk->node_upper = NULL;
	tlk->use_upper = true;
	tlk->update_upper = true;
	SingleTreeLikelihood_update_uppers(tlk);
    SitePattern *sp = tlk->sp;
    size_t tip_count = Tree_tip_count(tlk->tree);
	size_t node_count = Tree_node_count(tlk->tree);
	Node** nodes = Tree_get_nodes(tlk->tree, POSTORDER);
	
	// could be called multiple times
    if( tlk->sp->size != node_count ){
        for ( int i = 0; i < tlk->sp->count; i++ ) {
            sp->patterns[i] = realloc(sp->patterns[i], node_count*sizeof(uint8_t));
        }
        sp->size = Tree_node_count(tlk->tree);
        sp->names = realloc(sp->names, node_count*sizeof(char*));
        for ( size_t i = tip_count; i < node_count; i++ ) {
            sp->names[i] = NULL;
        }
		for ( size_t i = 0; i < node_count; i++ ) {
			if(!Node_isleaf(nodes[i]))
				tlk->mapping[Node_id(nodes[i])] = Node_class_id(nodes[i]) + tip_count;
		}
    }
    
    _marginal_reconstruction(tlk);
	tlk->use_upper = false;
}

// Method described in TreeTime paper
void estimate_GTR(SingleTreeLikelihood* tlk){
	asr_marginal(tlk);
	size_t nstate = tlk->sp->nstate;
	Node** nodes = Tree_get_nodes(tlk->tree, POSTORDER);
	double* m = dvector(nstate);
	double* T = dvector(nstate);
	double** N = dmatrix(nstate, nstate);
	size_t rootID = tlk->mapping[Tree_root(tlk->tree)->id];
	for(int i = 0; i < tlk->sp->count; i++){
		m[tlk->sp->patterns[i][rootID]]++;
	}
	
	for (int k = 0; k < Tree_node_count(tlk->tree); k++) {
		if(!Node_isleaf(nodes[k])){
			Node* node = nodes[k];
			size_t nodeID2 = Node_id(nodes[k]);
			size_t nodeID = tlk->mapping[nodeID2];
			size_t nodeLeftID = tlk->mapping[Node_id(nodes[k]->left)];
			size_t nodeRightID = tlk->mapping[Node_id(nodes[k]->right)];
			
			for(size_t p = 0; p < tlk->sp->count; p++){
//				int i = 0;
////				size_t idxi = p*tlk->sp->count;
//				for(; i < nstate; i++){
//					if(tlk->sp->patterns[p][nodeID] == 1){
//						T[i] += Node_distance(nodes[k]->left) + Node_distance(nodes[k]->right);
//						break;
//					}
////					idxi++;
//				}
////				size_t idxj = p*tlk->sp->count;
//				for(int j = 0; j < nstate; j++){
//					if(tlk->sp->patterns[p][nodeLeftID] == 1){
//						N[i][j]++;
////						N[j][i]++;
//						break;
//					}
////					idxj++;
//				}
////				idxj = p*tlk->sp->count;
//				for(int j = 0; j < nstate; j++){
//					if(tlk->sp->patterns[p][nodeRightID] == 1){
//						N[i][j]++;
////						N[j][i]++;
//						break;
//					}
////					idxj++;
//				}
//				printf("%d %d\n", tlk->sp->patterns[p][nodeID], tlk->sp->patterns[p][nodeRightID]);
//				printf("%d %d\n", tlk->sp->patterns[p][nodeID], tlk->sp->patterns[p][nodeLeftID]);
				if(tlk->sp->patterns[p][nodeRightID] < nstate)
				N[tlk->sp->patterns[p][nodeID]][tlk->sp->patterns[p][nodeRightID]]++;
				if(tlk->sp->patterns[p][nodeLeftID] < nstate)
				N[tlk->sp->patterns[p][nodeID]][tlk->sp->patterns[p][nodeLeftID]]++;
			}
		}
	}

	double** Q = dmatrix(nstate, nstate);
	double* pi = dvector(nstate);
	double pc = 0.1; // pseudo count
	for (int i = 0; i < nstate; i++) {
		pi[i] = 1.0/nstate;
	}
	for(int l = 0; l < 3; l++){
		for (int i = 0; i < nstate; i++) {
			for (int j = i+1; j < nstate; j++) {
				Q[i][j] = Q[j][i] = (N[i][j] + N[j][i] + 2.0*pc)/(pi[i]*T[j] + pi[j]*T[i] + 2.0*pc);
			}
		}
		double sum_pi = 0;
		for (int i = 0; i < nstate; i++) {
			double sum = 0;
			for (int j = 0; j < nstate; j++) {
				sum += N[i][j];
			}
			pi[i] = sum + m[i] + pc;
			sum = 0;
			double sum2 = 0;
			for (int j = 0; j < nstate; j++) {
				sum += Q[i][j]*T[j];
				sum2 += m[j] + pc;
			}
			pi[i] /= sum + sum2;
			sum_pi += pi[i];
			Q[i][i] = 0;
		}
		
		// normalize pi
		for (int i = 0; i < nstate; i++) {
			pi[i] /= sum_pi;
		}
		
		double subst = 0.0;
		for (int i = 0; i < nstate; i++) {
			double sum = 0;
			for (int j = 0; j < nstate; j++) {
				Q[i][j] *= pi[j];
				sum += Q[i][j];
			}
			Q[i][i] = -sum;
			subst += sum*pi[i];
		}

		for (int i = 0; i < nstate; i++) {
			for (int j = 0; j < nstate; j++) {
				Q[i][j] /= subst;
			}
		}
		print_dmatrix(stdout, Q, 4, 4, ',');
		print_dvector(pi, 4);
	}
	free(pi);
	free(T);
	free(m);
	free_dmatrix(N, nstate);
	free_dmatrix(Q, nstate);
}

void asr_calculator_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"algorithm",
		"file",
		"model",
		"verbosity"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* algorithm = get_json_node_value_string(node, "algorithm");
	char* ref = get_json_node_value_string(node, "model");
	Model* mtlk = Hashtable_get(hash, ref+1);
	const char *file = get_json_node_value_string(node, "file");
	int verbosity = get_json_node_value_int(node, "verbosity", 0);
	
	SingleTreeLikelihood *tlk = mtlk->obj;
	
	if(verbosity > 0) printf("Inferring ancestral sequences...");
	if(algorithm == NULL || strcasecmp(algorithm, "marginal") == 0){
		asr_marginal(tlk);
	}
	else{
		fprintf(stderr, "Only marginal ancestral reconstruction\n");
		exit(2);
	}
	Sequences *sequences = SitePattern_to_Sequences( tlk->sp );
	Sequences_save_fasta(sequences, file);
	free_Sequences(sequences);
	
	if(verbosity > 0) printf("done\n");
}
