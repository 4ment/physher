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
