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


#include "treelikelihood.h"
#include "matrix.h"
#include "sequenceio.h"


void _traverse2( SingleTreeLikelihood *tlk, Node *node, double *partials, int *index ){
    if( node == NULL ) return;
    _traverse2(tlk, Node_left(node), partials, index);
    _traverse2(tlk, Node_right(node), partials, index);
    double *backup = dvector(tlk->partials_size);
    
    if( !Node_isleaf(node) ){
        int sp_count = tlk->sp->count;
        int nstate = tlk->sm->nstate;
        SitePattern *sp = tlk->sp;
        //tlk->calculate_upper(tlk, node);
        //tlk->node_log_likelihoods_upper( tlk, node );
        //tlk->integrate_partials(tlk, tlk->partials[Node_id(node)], tlk->sm->get_proportions(tlk->sm), partials );
        
        for ( int i = 0; i < sp_count; i++ ) {
            memcpy(backup, tlk->partials[Node_id(node)], sizeof(double)*tlk->partials_size);
            double tot = 0;
            for ( int j = 0; j < nstate; j++ ) {
                for ( int k = 0; k < nstate; k++ ){
                    if( j != k ) tlk->partials[Node_id(node)][i*nstate+k] = 0;
                    //else tlk->partials[Node_id(node)][i*nstate+k] = 1;
                }
                SingleTreeLikelihood_update_all_nodes(tlk);
                tlk->calculate(tlk);
                partials = tlk->pattern_lk + tlk->sp->count;
                tot += partials[i*nstate+j];
                //tot += tlk->node_pattern_lk[i*nstate+j];
                memcpy(tlk->partials[Node_id(node)], backup, sizeof(double)*tlk->partials_size);
            }
            int max_j  = 0;
            double max_Lj = partials[i*nstate]/tot;
            //double max_Lj = tlk->node_pattern_lk[i*nstate]/tot;
            for ( int j = 1; j < nstate; j++ ) {
                double Lj = partials[i*nstate+j]/tot;
                //double Lj = tlk->node_pattern_lk[i*nstate+j]/tot;
                if( Lj > max_Lj ){
                    max_Lj = Lj;
                    max_j  = j;
                }
            }
            sp->patterns[i][*index] = max_j;
        }
        sp->names[*index] = String_clone( node->name );
        ++(*index);
        free(backup);
        SingleTreeLikelihood_update_all_nodes(tlk);
    }
}
bool _calculate_partials( SingleTreeLikelihood *tlk, Node *n  ){
	bool updated = tlk->update_nodes[ Node_id(n) ];

	
	if( !Node_isleaf(n) ){
		bool update_child1 = _calculate_partials( tlk, Node_left(n) );
		bool update_child2 = _calculate_partials( tlk, Node_right(n) );
		
		if( update_child1 || update_child2 ){
			int indx_child1 = Node_id(Node_left(n));
			int indx_child2 = Node_id(Node_right(n));
            
			tlk->update_partials( tlk, Node_id(n), indx_child1, indx_child1, indx_child2, indx_child2);
			
			if( Node_isroot(n) ){
				tlk->integrate_partials(tlk, tlk->partials[Node_id(n)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
				
				tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
			}
			updated = true;
		}
	}
	return updated;
}

void _traverse3( SingleTreeLikelihood *tlk, Node *node, double *partials, int *index ){
    if( node == NULL ) return;
    _traverse3(tlk, Node_left(node), partials, index);
    _traverse3(tlk, Node_right(node), partials, index);
    
    
    if( !Node_isleaf(node) ){
        double *backup = dvector(tlk->partials_size);
        int sp_count = tlk->sp->count;
        int nstate = tlk->sm->nstate;
        uint8_t **patterns = tlk->sp->patterns;
        
        for ( int i = 0; i < sp_count; i++ ) {
            memcpy(backup, tlk->partials[Node_id(node)], sizeof(double)*tlk->partials_size);
            
            int max_j  = 0;
            double max_Lj = -INFINITY;
        
            for ( int j = 0; j < nstate; j++ ) {
                
                for ( int k = 0; k < nstate; k++ ){
                    if( j != k )
                        tlk->partials[Node_id(node)][i*nstate+k] = 0;
                }
                
                if( Node_isroot(node) ){
                    tlk->integrate_partials(tlk, tlk->partials[Node_id(node)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
                    tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
                }
                else {
                    SingleTreeLikelihood_update_one_node(tlk, node);
                    _calculate_partials(tlk, Tree_root(tlk->tree));
                }
                
                double Lj = tlk->pattern_lk[i];
                
                if( Lj > max_Lj ){
                    max_Lj = Lj;
                    max_j  = j;
                }
                memcpy(tlk->partials[Node_id(node)], backup, sizeof(double)*tlk->partials_size);
                SingleTreeLikelihood_update_all_nodes(tlk);
                tlk->calculate(tlk);
            }
            
            patterns[i][*index] = max_j;
        }
        tlk->sp->names[*index] = String_clone( node->name );
        ++(*index);
        free(backup);
    }
}

// Marginal reconstruction of ancestral sequence
void _traverse( SingleTreeLikelihood *tlk, Node *node, int *index ){
    
    if( !Node_isleaf(node) ){
        _traverse(tlk, Node_left(node), index);
        _traverse(tlk, Node_right(node), index);
    
        double *backup = dvector(tlk->partials_size);
        
        int sp_count = tlk->sp->count;
        int nstate = tlk->sm->nstate;
        SitePattern *sp = tlk->sp;
        
        double *Lj = dvector(sp_count);
        
        memcpy(backup, tlk->partials[Node_id(node)], sizeof(double)*tlk->partials_size);
        
        for ( int j = 0; j < nstate; j++ ) {
            
            for ( int i = 0; i < tlk->partials_size; i++ ){
                if( (i-j) % nstate != 0 ) tlk->partials[Node_id(node)][i] = 0;
            }
            
            if( Node_isroot(node) ){
                tlk->integrate_partials(tlk, tlk->partials[Node_id(node)], tlk->sm->get_proportions(tlk->sm), tlk->root_partials );
                tlk->node_log_likelihoods( tlk, tlk->root_partials, tlk->get_root_frequencies(tlk), tlk->pattern_lk);
            }
            else {
                SingleTreeLikelihood_update_one_node(tlk, node);
                _calculate_partials(tlk, Tree_root(tlk->tree));
            }
            
            memcpy(tlk->partials[Node_id(node)], backup, sizeof(double)*tlk->partials_size);
            
            if ( j == 0 ) {
                for ( int i = 0; i < sp_count; i++ ) {
                    sp->patterns[i][*index] = 0;
                }
                memcpy(Lj, tlk->pattern_lk, sp_count*sizeof(double));
            }
            else {
                for ( int i = 0; i < sp_count; i++ ) {
                    if( tlk->pattern_lk[i] > Lj[i] ){
                        Lj[i] = tlk->pattern_lk[i];
                        sp->patterns[i][*index] = j;
                    }
                }
            }
            
        }
        
        SingleTreeLikelihood_update_all_nodes(tlk);
        tlk->calculate(tlk);
        
        if( sp->names[*index] != NULL ){
            free(sp->names[*index]);
            sp->names[*index] = NULL;
        }
        sp->names[*index] = String_clone( node->name );
        ++(*index);
        free(backup);
        free(Lj);
    }
}

void asr_marginal( SingleTreeLikelihood *tlk ){
    tlk->calculate(tlk);
    SitePattern *sp = tlk->sp;
    int index = Tree_tip_count(tlk->tree);
    double *partials = dvector(tlk->sp->count*tlk->sm->nstate);
    
    if( tlk->sp->size != Tree_node_count(tlk->tree) ){
        for ( int i = 0; i < tlk->sp->count; i++ ) {
            sp->patterns[i] = realloc(sp->patterns[i], Tree_node_count(tlk->tree)*sizeof(uint8_t));
        }
        sp->size = Tree_node_count(tlk->tree);
        sp->names = realloc(sp->names, Tree_node_count(tlk->tree)*sizeof(char*));
        for ( int i = Tree_tip_count(tlk->tree); i < Tree_node_count(tlk->tree); i++ ) {
            sp->names[i] = NULL;
        }
    }
    
    _traverse(tlk, Tree_root(tlk->tree), &index);
    free(partials);
}

void asr_marginal_calculator_from_json(json_node* node, Hashtable* hash){
	char* ref = get_json_node_value_string(node, "model");
	Model* mtlk = Hashtable_get(hash, ref+1);
	const char *file = get_json_node_value_string(node, "file");
	int verbosity = get_json_node_value_int(node, "verbosity", 0);
	
	SingleTreeLikelihood *tlk = mtlk->data;
	
	if(verbosity > 0) printf("Inferring ancestral sequences...");
	asr_marginal(tlk);
	Sequences *sequences = SitePattern_to_Sequences( tlk->sp );
	Sequences_save_fasta(sequences, file);
	free_Sequences(sequences);
	
	if(verbosity > 0) printf("done\n");
	
	Node **nodes = Tree_nodes(tlk->tree);
	for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
		tlk->mapping[i] = get_sequence_index(tlk->sp, nodes[i]->name);
	}
	tlk->mapping[Tree_root(tlk->tree)->id] = -1;
	
	tlk->sm->need_update = true;
}
