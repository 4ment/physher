/*
 *  physim.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 9/09/13.
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

#include "physim.h"

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "sequence.h"
#include "branchmodel.h"
#include "sitemodel.h"
#include "tree.h"
#include "random.h"
#include "matrix.h"
#include "geneticcode.h"
#include "sitepattern.h"



// for codons len is the number of codons
Sequences * Sequence_simulate( Tree *tree, SiteModel *sm, BranchModel *bm, DataType *datatype, unsigned len, bool keep_internal ){
    int nstate = sm->nstate;
    double *p = dvector(nstate*nstate);
    int *rates = ivector(len);
    
    double *freqs = sm->m->_freqs;
    
    double accum = 0;
    for ( int i = 0; i < nstate; i++ ) {
        accum += freqs[i];
    }
    
    if(accum != 1.0 && fabs(1.0-accum) < 0.001){
        accum = 0;
        for ( int i = 1; i < nstate; i++ ) {
            accum += freqs[i];
        }
        freqs[0] = 1.0 - accum;
    }
    else if( accum != 1.0 ){
        fprintf(stderr, "Frequencies must add up to 1 (%e)\n", accum);
        for ( int i = 0; i < nstate; i++ ) {
            fprintf(stderr, "%d %e\n", i, freqs[i]);
        }
        exit(1);
    }
    
    sm->m->need_update = true;
    sm->need_update = true;
    
    if( sm->type != SITEMODEL_NONE ){
        double accum = 0;
        for ( int i = 0; i < sm->cat_count; i++ ) {
            accum += sm->get_proportion(sm,i);
        }
        if( accum < 0.000001 ){
            printf("Proportions in SiteModel must add up to 1\n");
            for ( int i = 0; i < sm->cat_count; i++ ) {
                fprintf(stderr, "%d %e\n",i, sm->get_proportion(sm,i));
            }
            exit(1);
        }
        
        for ( int i = 0; i < len; i++ ) {
            rates[i] = roulette_wheel(sm->get_proportions(sm), sm->cat_count);
            //fprintf(stderr,"%d,", rates[i]);
        }
        //fprintf(stderr,"\n");
    }
    
    Sequences *seqs = new_Sequences( Tree_node_count(tree));
    seqs->datatype = datatype;
    
    seqs->length = len;
    if( datatype->type == DATA_TYPE_CODON ){
        seqs->length *= 3;
    }
    
    Node **nodes = Tree_get_nodes(tree, PREORDER);
    // This is the root sequence
    Sequence *current = new_Sequence2(Node_name(Tree_root(tree)), seqs->length);
    seqs->seqs[Tree_root(tree)->postorder_idx] = current;
    
    
    int residue_length = strlen(seqs->datatype->state_string(seqs->datatype,0));
    char *residue = (char*)malloc((residue_length+1)*sizeof(char));
    assert(residue);
    residue[residue_length] = '\0';
    
    for ( int i = 0; i < len; i++ ) {
        const char *temp = seqs->datatype->state_string(seqs->datatype, roulette_wheel(freqs, nstate));
        strcpy(&current->seq[i*residue_length], temp);
    }
    current->length = len;
    
    
    Sequence *parent  = NULL;
    // Iterate over nodes
    for ( int i = 1; i < Tree_node_count(tree); i++ ) {
        
        current = new_Sequence2(Node_name(nodes[i]), seqs->length);
        seqs->seqs[ nodes[i]->postorder_idx ] = current;
        parent = seqs->seqs[ Node_parent(nodes[i])->postorder_idx ];
        
        // Iterate over sites
        for ( int j = 0; j < len; j++ ) {
            double t = sm->get_rate(sm, rates[j]);
            if ( bm == NULL ) {
                t *= Node_distance(nodes[i]);
            }
            else {
                t *= Node_time_elapsed(nodes[i]) * bm->get(bm, nodes[i]);
            }
            sm->m->p_t(sm->m, t, p);

            memcpy(residue, &parent->seq[j], residue_length*sizeof(char));
            int parent_state = seqs->datatype->encoding_string(seqs->datatype, residue);
            const char *temp = seqs->datatype->state_string(seqs->datatype, roulette_wheel(&p[parent_state*nstate], nstate));
            strcpy(&current->seq[j*residue_length], temp);
            current->length = len;
            
        }
    }
    seqs->size = Tree_node_count(tree);
    
    if( !keep_internal ){
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            if ( !Node_isleaf(nodes[i]) ) {
                free_Sequence(seqs->seqs[nodes[i]->postorder_idx]);
                seqs->seqs[nodes[i]->postorder_idx] = NULL;
                //Sequences_delete(seqs, nodes[i]->postorder_idx);
            }
        }
        Sequences_pack2(seqs);
    }
    
    seqs->aligned  = true;
    seqs->datatype->genetic_code  = sm->m->gen_code;
    seqs->datatype->type = sm->m->dtype;
    
    free(residue);
    free(p);
    free(rates);
    return seqs;
}

