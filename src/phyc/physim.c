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
#include <strings.h>

#include "sequence.h"
#include "sequenceio.h"
#include "branchmodel.h"
#include "sitemodel.h"
#include "tree.h"
#include "random.h"
#include "matrix.h"
#include "geneticcode.h"
#include "sitepattern.h"
#include "lognormal.h"
#include "exponential.h"



// for codons len is the number of codons
Sequences * Sequence_simulate( Tree *tree, SiteModel *sm, BranchModel *bm, DataType *datatype, unsigned len, bool keep_internal ){
    int nstate = sm->nstate;
    double *p = dvector(nstate*nstate);
    int *rates = ivector(len);
    
    const double *freqs = sm->m->get_frequencies(sm->m);
    
    double accum = 0;
    for ( int i = 0; i < nstate; i++ ) {
        accum += freqs[i];
    }
    
    if(fabs(1.0-accum) > 0.001){
        fprintf(stderr, "Frequencies must add up to 1 (%e)\n", accum);
        for ( int i = 0; i < nstate; i++ ) {
            fprintf(stderr, "%d %e\n", i, freqs[i]);
        }
        exit(1);
    }
    
    sm->m->need_update = true;
    sm->need_update = true;
    
    if( sm->distribution != DISTRIBUTION_UNIFORM ){
        double accum = 0;
        for ( int i = 0; i < sm->cat_count; i++ ) {
            accum += sm->get_proportion(sm,i);
        }
        if( accum < 0.000001 ){
            fprintf(stderr, "Proportions in SiteModel must add up to 1\n");
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
    seqs->datatype->type = sm->m->datatype->type;
    
    free(residue);
    free(p);
    free(rates);
    return seqs;
}

char * Node_get_string_from_info2( const Node *node, const char *str ){
    char *value = NULL;
    if ( node->info != NULL ) {
        char *ptr = node->info;
        if ( String_contains_str(ptr, str) ) {
            StringBuffer *buffer = new_StringBuffer(20);
            ptr += String_index_of_str(ptr, str) + strlen(str);
            bool range = false;
            if(*ptr=='{') range = true;
            
            if( range ){
                while( *ptr != '}' ){
                    StringBuffer_append_char(buffer, *ptr);
                    ptr++;
                }
                StringBuffer_append_char(buffer, *ptr);
            }
            else {
                while( *ptr != ',' && *ptr != ']' ){
                    
                    StringBuffer_append_char(buffer, *ptr);
                    ptr++;
                }
            }
            if (buffer->length !=0 ) {
                value = String_clone(buffer->c);
            }
            
            
            free_StringBuffer(buffer);
        }
    }
    return value;
}

void SimulateSequences_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"branchmodel",
		"datatype",
		"distribution",
		"format",
		"internal",
		"length",
		"output",
		"scaler",
		"sitemodel",
		"tree",
		"verbosity"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
    char* output = get_json_node_value_string(node, "output");
    char* format = get_json_node_value_string(node, "format");
    json_node* tree_node = get_json_node(node, "tree");
    json_node* sm_node = get_json_node(node, "sitemodel");
    json_node* bm_node = get_json_node(node, "branchmodel");
    json_node* datatype_node = get_json_node(node, "datatype");
    int length = get_json_node_value_int(node, "length", 1000);
    bool internal = get_json_node_value_bool(node, "internal", false);
    double scaler = get_json_node_value_double(node, "scaler", 1.0);
    char* rate_dist = get_json_node_value_string(node, "distribution");
	int verbosity = get_json_node_value_int(node, "verbosity", 1);

	if(format == NULL){
		format = "fasta";
	}
    
    Model* mtree = NULL;
    Model* msm = NULL;
    Model* mbm = NULL;
    
    if (tree_node->node_type == MJSON_STRING) {
        char* ref = (char*)tree_node->value;
        // check it starts with a &
        mtree = Hashtable_get(hash, ref+1);
        mtree->ref_count++;
    }
    else{
        char* id = get_json_node_value_string(tree_node, "id");
        mtree = new_TreeModel_from_json(tree_node, hash);
        Hashtable_add(hash, id, mtree);
    }
    
    if (sm_node->node_type == MJSON_STRING) {
        char* ref = (char*)sm_node->value;
        // check it starts with a &
        msm = Hashtable_get(hash, ref+1);
        msm->ref_count++;
    }
    else{
        char* id = get_json_node_value_string(sm_node, "id");
        msm = new_SiteModel_from_json(sm_node, hash);
        Hashtable_add(hash, id, msm);
    }
    
    BranchModel* bm = NULL;
    if(bm_node != NULL){
        if (bm_node->node_type == MJSON_STRING) {
            char* ref = (char*)bm_node->value;
            // check it starts with a &
            mbm = Hashtable_get(hash, ref+1);
            mbm->ref_count++;
        }
        else{
            //			sm = new_SiteModel_from_json(sm_node, hash);
        }
        bm = mbm->obj;
    }
    
    DataType* datatype = NULL;
    if(datatype_node->node_type == MJSON_STRING && (strcasecmp((char*)datatype_node->value, "nucleotide") == 0 ||
                                                    strcasecmp((char*)datatype_node->value, "codon") == 0 || strcasecmp((char*)datatype_node->value, "aa") == 0)){
        datatype = new_DataType_from_json(datatype_node, hash);
    }
    else if (datatype_node->node_type == MJSON_STRING && Hashtable_exists(hash, (char*)datatype_node->value)) {
        datatype = Hashtable_get(hash, (char*)datatype_node->value);
        datatype->ref_count++;
    }
    else{
        datatype = new_DataType_from_json(datatype_node, hash);
        Hashtable_add(hash, datatype->name, datatype);
    }
    Tree* tree = mtree->obj;
    
    if(scaler != 1.0){
        Tree_scale_distance(tree, scaler);
    }
    
    if ( rate_dist != NULL ) {
        double *rates = dvector(Tree_node_count(tree)-1);
        Node **nodes = Tree_get_nodes(tree, POSTORDER);
        
        // rates are provided in the tree
        if ( strcasecmp(rate_dist, "tree") == 0 ) {
            for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
                rates[i] = Node_get_double_from_info(nodes[i], "rate=");
            }
            
            //StringBuffer_append_format(buffer, "\nRates from tree file\n");
        }
        else {
            // lognormal
            double* means = NULL;
            double* logsigmas = NULL;
            // exponential
            double* lambdas = NULL;
            
            unsigned l = 0;
            if (strcasecmp(rate_dist, "lognormal") == 0) {
                
                char* mean_str = get_json_node_value_string(node, "mean");
                if( mean_str == NULL ) error("Could not read the mean of the lognormal distribution of categories");
                
                means = String_to_double_array(mean_str, ',', &l);
                for (int i = 0; i < l; i++) {
                    if( means[i] <= 0 ) error("Mean cannot be negative");
                }
                free(mean_str);

                char* logsigma_str = get_json_node_value_string(node, "sd");
                if( logsigma_str == NULL ) error("Could not read the standard deviation of the distribution of categories [-S]");
                unsigned l2 = 0;
                logsigmas = String_to_double_array(logsigma_str, ',', &l2);
                for (int i = 0; i < l2; i++) {
                    if( logsigmas[i] <= 0 ) error("logsigma cannot be negative");
                }
                free(logsigma_str);
                
                if(l!=l2){
                    error("The number of means and sds don't match");
                }
                
            }
            else if(strcasecmp(rate_dist, "exponential") == 0){
                char* lambda_str = get_json_node_value_string(node, "lambda");
                if( lambda_str == NULL ) error("Could not read the lambda of the exponential distribution");

                lambdas = String_to_double_array(lambda_str, ',', &l);
                for (int i = 0; i < l; i++) {
                    if( lambdas[i] <= 0 ) error("Lambdas cannot be negative");
                }
                free(lambda_str);
            }
            
            if(l > 1){
                int* counts = calloc(l, sizeof(int));
                for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                    if( !Node_isroot(nodes[i]) ){
                        counts[Node_get_int_from_info(nodes[i], "class=")]++;
                    }
                }
                double* temp_rates = dvector(Tree_node_count(tree)-1);
                
                for (int i = 0; i < l; i++) {
                    if(strcasecmp(rate_dist, "lognormal") == 0){
                        lognormal_discretize(log(means[i]), logsigmas[i], temp_rates, counts[i]);
//                            StringBuffer_append_format(buffer, "\n%d Lognormal distribution mean = %f log(stdev) = %f", i, means[i], logsigmas[i]);
                    }
                    else{
                        exponential_discretize(lambdas[i], temp_rates, counts[i]);
//                            StringBuffer_append_format(buffer, "\n%d Exponential distribution mean = %f", i, means[i]);
                    }
                    randomize_dvector( temp_rates, counts[i]);
                    int k = 0;
                    for ( int j = 0; j < Tree_node_count(tree); j++ ) {
                        int class = Node_get_int_from_info(nodes[j], "class=");
                        if( !Node_isroot(nodes[j]) && class == i){
                            rates[j] = temp_rates[k];
                            k++;
                        }
                    }
                }
//                    StringBuffer_append_string(buffer, "\n");
                free(temp_rates);
                free(counts);
            }
            else{
                if(strcasecmp(rate_dist, "lognormal") == 0){
                    lognormal_discretize(log(means[0]), logsigmas[0], rates, Tree_node_count(tree)-1);
                }
                else{
                    exponential_discretize(means[0], rates, Tree_node_count(tree)-1);
                }
                randomize_dvector( rates, Tree_node_count(tree)-1);
            }
            if(means != NULL)free(means);
            if(logsigmas != NULL) free(logsigmas);
            if(lambdas != NULL) free(lambdas);
            
        
            StringBuffer *info = new_StringBuffer(100);
            for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    StringBuffer_empty(info);
                    StringBuffer_append_format(info, "%e", rates[i]);
                    Node_set_annotation(nodes[i], "rate", info->c);
                    
                    StringBuffer_empty(info);
                    int class = Node_get_int_from_info(nodes[i], "class=");
                    StringBuffer_append_format(info, "%d", class);
                    Node_set_annotation(nodes[i], "class", info->c);
                }
            }
            
            // setup heights
            for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                double height = Node_get_double_from_info(nodes[i], "height=");
                Node_set_height(nodes[i], height);
                
                StringBuffer_empty(info);
                StringBuffer_append_format(info, "%f", height);
                Node_set_annotation(nodes[i], "height", info->c);
                
                char *cal = Node_get_string_from_info2(nodes[i], "cal_height=");
                if( cal != NULL ){
                    Node_set_annotation(nodes[i], "cal_height", cal);
                    free(cal);
                }
            }
            free_StringBuffer(info);
        }
        
        for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
            Node_set_distance(nodes[i], Node_distance(nodes[i])*rates[i]);
        }
        
//        double min = dmin_vector(rates, Tree_node_count(tree)-1);
//        double max = dmax_vector(rates, Tree_node_count(tree)-1);
//        double mean = dmean(rates, Tree_node_count(tree)-1);
//        double median = dmedian(rates, Tree_node_count(tree)-1);
//        StringBuffer_append_format(buffer, "\nRate min = %f max = %f mean = %f median = %f\n", min,max,mean,median);
        
        free(rates);
        free(rate_dist);
    }
    
    Sequences *sim = Sequence_simulate(tree, msm->obj, bm, datatype, length, internal);
    
    if(strcasecmp(format, "fasta") == 0){
        Sequences_save_fasta(sim, output);
    }
    else if(strcasecmp(format, "nexus") == 0){
        Sequences_save_nexus(sim, output);
    }
    else if(strcasecmp(format, "phylip") == 0){
        Sequences_save_phylip(sim, output);
    }
    
    int polymorphisms = 0;
    for ( int i = 0; i < sim->length; i++ ) {
        char c = sim->seqs[0]->seq[i];
        int j = 1;
        for ( ; j < sim->size; j++ ) {
            if( c != sim->seqs[j]->seq[i] ) break;
        }
        if( j != sim->size ){
            polymorphisms++;
        }
    }
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    double tree_length = 0;
    for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
        tree_length += Node_distance(nodes[i]);
    }
	if (verbosity > 0) {
		printf("Tree length %f\n", tree_length);
		fprintf(stdout, "Number of polymorphic sites: %d/%d (%f)\n", polymorphisms, length,((double)polymorphisms/length) );
	}
	
    free_Sequences(sim);
    mtree->free(mtree);
    msm->free(msm);
    if(mbm != NULL) mbm->free(mbm);
}
