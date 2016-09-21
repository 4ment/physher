/*
 *  localclock.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 1/12/11.
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

#include "localclock.h"

#include <math.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif

#include "ga.h"
#include "treelikelihood.h"
#include "optimize.h"
#include "utils.h"
#include "matrix.h"
#include "treeio.h"
#include "random.h"
#include "modelselection.h"
#include "filereader.h"
#include "descriptivestats.h"
#include "combinatorics.h"

#include "heterotachy.h"

//#define PTHREAD_ENABLED 1
//#define OPENMP_TASK_ENABLED 1

static void _localclock_termination_callback( GA *ga );

// reset and return the TreeLikelihood with index index
static SingleTreeLikelihood * localclock_reset_paramters( ClockSearch *searcher, unsigned index ){
	Tree_vector_to_heights( searcher->current_heights, searcher->pool->tlks[index]->tree );
	BranchModel_vector_to_rates(searcher->pool->tlks[index]->bm, searcher->current_rates);
	Tree_constraint_heights(searcher->pool->tlks[index]->tree);
	SingleTreeLikelihood_update_all_nodes(searcher->pool->tlks[index]);
	
	return searcher->pool->tlks[index];
}


#pragma mark -
// MARK: Localclock search



#if defined (PTHREAD_ENABLED)
static void _localclock_greedy_search_pthreads( ClockSearch *search, const int current );
#else

#ifdef OPENMP_TASK_ENABLED
static void _localclock_greedy_search_task( ClockSearch *search, const int current, const int start, int *count );
#else
static void _localclock_greedy_search( ClockSearch *search, const int current );
#endif

#endif

static void _localclock_exhaustive_search( ClockSearch *search, int *count );

static double _localclock_optimize( ClockSearch *search );

ClockSearch * new_LocalClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads ){
	ClockSearch *search = new_ClockSearch(tlk, nThreads);
	
	search->pDerivedObj = search; // point to itself
	
	search->optimize = _localclock_optimize;
	search->free     = free_ClockSearch;
	search->init     = ClockSearch_init;
	
	search->algorithm = HETEROTACHY_LOCALCLOCK_GREEDY;

	search->verbosity = 2;
	
	return search;
}

static void quick_check( ClockSearch *search );

static double _localclock_optimize( ClockSearch *search  ){
	if ( !search->initialized ) {
		search->init(search);
	}
	
	Tree *tree = search->pool->tlks[0]->tree;
	BranchModel *bm = search->pool->tlks[0]->bm;
	int nNodes = Tree_node_count(tree);
    Node **nodes = Tree_nodes(tree);
	
	
	search->best_indexes    = realloc( search->best_indexes,    (search->n_rate-1) * sizeof(unsigned) );
	search->current_indexes = realloc( search->current_indexes, (search->n_rate-1) * sizeof(unsigned) );
    
    LocalClock_get_indexes(bm, search->best_indexes);
    memcpy(search->current_indexes, search->best_indexes, (search->n_rate-1) * sizeof(unsigned));
    
	Tree_heights_to_vector(tree, search->best_heights);
    memcpy(search->current_heights, search->best_heights, nNodes * sizeof(double));
    
	BranchModel_rates_to_vector( bm, search->best_rates );
    memcpy(search->current_rates, search->best_rates, search->n_rate * sizeof(double));
	
    search->default_rate = dmean(search->best_rates, Parameters_count(bm->rates) );
    
	double previous_lk = search->best_lk;
    
    int previous_df = search->starting_df;
	double previous_IC = search->IC(search, previous_lk, search->starting_df);
	double IC;
    int previous_n_rate = SingleTreeLikelihood_df_count(search->pool->tlks[0])-previous_df;
    
	for ( ; search->n_rate <= search->max_n_rate; search->n_rate++) {
		
		if ( search->n_rate != search->starting_n_rate ) {
			for ( int j = 0; j < search->pool->count; j++ ) {
				LocalClock_set_number_of_clocks(search->pool->tlks[j]->bm, search->n_rate-1);
				BranchModel_vector_to_rates(search->pool->tlks[j]->bm, search->current_rates);
                
			}
		}
		
//        for ( int j = 0; j < search->pool->count; j++ ) {
//            // fix lower bounds
//            for (int p = 0; p < Parameters_count(search->pool->tlks[j]->bm->rates); p++) {
//                //Parameters_set_lower(search->pool->tlks[j]->bm->rates, p, Parameters_lower(search->pool->tlks[j]->bm->rates, 0) );
//            }
//        }
		
		time_t start_time;
		time_t end_time;
		time(&start_time);
		
		if ( search->verbosity > 0) {
			fprintf(stdout, "---------------------------------\n");
			fprintf(stdout, "Number of local clocks: %d\n", search->n_rate-1);
            
			if ( search->verbosity > 1) {
                fprintf(stdout, "Iteration 0/%d\r", Tree_node_count(tree)-search->n_rate+1);
            }
		}
        
		
        if ( search->algorithm == HETEROTACHY_LOCALCLOCK_GREEDY ) {
#ifdef PTHREAD_ENABLED
            _localclock_greedy_search_pthreads(search, search->n_rate-2);
#else
#ifdef OPENMP_TASK_ENABLED
            // used with tasks in openmp 3.0
            int count = 0;
            #pragma omp parallel
            #pragma omp single nowait
            {
                _localclock_greedy_search_task(search, search->n_rate-2, search->n_rate-2, &count);
            }
#else
            _localclock_greedy_search(search, search->n_rate-2);
#endif
#endif
        }
        else if( search->algorithm == HETEROTACHY_LOCALCLOCK_EXHAUSTIVE ){
            int count = 0;
            _localclock_exhaustive_search(search, &count);
        }
        else {
            error("Error algorithm in localclock\n");
        }
		
        search->df = SingleTreeLikelihood_df_count(search->pool->tlks[0]);
        
        if( search->logFile != NULL ) fprintf(search->logFile, "\n");
        
		if ( search->algorithm == HETEROTACHY_LOCALCLOCK_EXHAUSTIVE ) {
            IC = search->IC(search, search->best_lk, search->df);
            if ( search->verbosity > 0) {
                time(&end_time);
                printf("\n");
                double t = difftime(end_time, start_time);
                printf("Runtime ");
                print_pretty_time(stdout,t);
                
                printf("\n");
                
                printf("%s_%d = %f LnL = %f df = %d\n", INFORMATION_CRITERION[search->ic], previous_n_rate, previous_IC, previous_lk, previous_df);
                printf("%s_%d = %f LnL = %f df = %d\n", INFORMATION_CRITERION[search->ic], search->n_rate, IC, search->best_lk, search->df);
                printf("==============================================\n\n");
			}
            
			if ( IC > previous_IC ) {
				break;
			}
            previous_df = search->df;
            previous_IC = IC;
            previous_n_rate = search->n_rate;
        }
        else if(true) {
            
            IC = search->IC(search, search->best_lk, search->df);
            
            if ( search->verbosity > 0) {
                time(&end_time);
                printf("\n");
                double t = difftime(end_time, start_time);
                printf("Runtime ");
                print_pretty_time(stdout,t);
                
                if ( search->verbosity > 1 ) {
                    printf("Locations: ");
                    for (int i = 0; i < search->n_rate-1; i++ ) {
                        printf(" %s[%d]", Node_name(nodes[search->best_indexes[i]]), search->best_indexes[i]);
                    }
                    printf("\nRates: ");
                    for (int i = 0; i < search->n_rate; i++ ) {
                        printf(" %e", search->best_rates[i]);
                    }
                    
                    printf("\n");
                    printf("Root height %f\n", search->best_heights[Node_id(Tree_root(tree))] );
                }
                
                printf("\n");
                printf("%s_%d = %f LnL = %f df = %d\n", INFORMATION_CRITERION[search->ic], previous_n_rate, previous_IC, previous_lk, previous_df);
                printf("%s_%d = %f LnL = %f df = %d\n", INFORMATION_CRITERION[search->ic], search->n_rate, IC, search->best_lk, search->df);
            }
            
            if ( IC > previous_IC ) {
				break;
			}
            previous_df = search->df;
            previous_IC = IC;
            previous_n_rate = search->n_rate;
            //quick_check(search);
		}
        else {
            
            double p = LRT(previous_lk, search->best_lk, search->n_rate-1, search->n_rate);
            
            if ( search->verbosity > 0) {
                time(&end_time);
                printf("\n");
                double t = difftime(end_time, start_time);
                printf("Runtime ");
                print_pretty_time(stderr,t);
                
                if ( search->verbosity > 1 ) {
                    printf("Locations: ");
                    for (int i = 0; i < search->n_rate-1; i++ ) {
                        printf(" %s[%d]", Node_name(nodes[search->best_indexes[i]]), search->best_indexes[i]);
                    }
                    
                    printf("\n");
                }
                
                printf("\n");
                printf("LnL0 = %f\n", previous_lk);
                printf("LnL1 = %f\n", search->best_lk);
                printf("p-value = %f %s\n\n", p, (p<search->alpha ? "*" : "") );
            }
            
            if( p > search->alpha ){
                break;
            }
            //quick_check(search);
		}
                
		previous_lk = search->best_lk;
		
		memcpy(search->current_indexes, search->best_indexes,(search->n_rate-1) * sizeof(unsigned) );
		search->current_indexes = realloc( search->current_indexes, search->n_rate * sizeof(unsigned));
		search->best_indexes    = realloc( search->best_indexes,    search->n_rate * sizeof(unsigned));
		search->current_indexes[search->n_rate-1] = search->best_indexes[search->n_rate-1] = 0;
		
		memcpy(search->current_heights, search->best_heights, nNodes * sizeof(double) );
		
		search->best_rates    = realloc( search->best_rates,    (search->n_rate+1) * sizeof(double) );
		search->current_rates = realloc( search->current_rates, (search->n_rate+1) * sizeof(double) );
		// assign the ancestral rate to the new local clock (could inherit the rate from it's parent)
		search->best_rates[search->n_rate] = search->best_rates[0];
		memcpy(search->current_rates, search->best_rates, (search->n_rate+1) * sizeof(double) );
	}
	
	if ( search->n_rate == search->max_n_rate){
		fprintf(stderr, "The potential number of local clocks (%d) is greater than the maximum (%d)\n", search->n_rate, search->max_n_rate );
	}
    // the algorithm stopped after the first iteration and the previous model was a strict clock
    else if ( search->n_rate == 2 ){
		
	}
	else {
        
        // from here I assume that the search has found at least one local clock
        // should do somthing about it
        
        // save the best heights and rates
        search->n_rate--;
        
        memcpy(search->best_heights, search->current_heights, nNodes * sizeof(double) );
        search->best_rates = realloc( search->best_rates, search->n_rate * sizeof(double));
        memcpy(search->best_rates, search->current_rates, search->n_rate * sizeof(double) );
        
        search->best_indexes = realloc( search->best_indexes, (search->n_rate-1) * sizeof(unsigned));
        memcpy(search->best_indexes, search->current_indexes, (search->n_rate-1) * sizeof(unsigned) );
        
        search->best_lk = previous_lk;
    }
    
    // reset BranchModel
    LocalClock_set_number_of_clocks(bm, search->n_rate-1);
    localclock_set_indicators2(bm, search->best_indexes);
    BranchModel_vector_to_rates(bm, search->best_rates);
    // reset tree heights in order to get the root height
    Tree_vector_to_heights( search->best_heights, tree);
    Tree_constraint_heights(tree);
    SingleTreeLikelihood_update_all_nodes(search->pool->tlks[0]);
    
	if( search->logFile != NULL ){
		// we only write the header if it is on write mode.
		// If we are in append mode then the header will be printed by the caller
		if ( search->logFileMode[0] == 'w') {
			fprintf(search->logFile, "\nEnd;\n");
		}
		fclose(search->logFile);
	}
	
	return search->best_lk;
	
}



void quick_check( ClockSearch *search ){
    FileReader *reader = new_FileReader( search->logFileName, 1000);
    StringBuffer *buffer = new_StringBuffer(100);
    
    // keep track of the best trees and their likelihood given a number of clock
    double *lnls = dvector(10);
    char **trees = (char**)malloc(sizeof(char *)*10);
    assert(trees);
    for ( int i = 0; i < 10; i++ ) {
        lnls[i] = -INFINITY;
        trees[i] = NULL;
    }
    
    char *ptr = NULL;
    
    while ( reader->read_line(reader) ){
        ptr = reader->line;
        if ( String_start_with(ptr, "tree", true) ) {
            StringBuffer_empty(buffer);
            
            ptr += String_index_of_str(ptr, "LnL=") + strlen("LnL=");
            while ( isspace(*ptr) ) ptr++;
            while ( (*ptr >= 48 && *ptr <= 57) || *ptr == '.' || *ptr == '-' ) {
                StringBuffer_append_char(buffer, *ptr);
                ptr++;
            }
            double lnl = atof(buffer->c);
            
            int pos = 0;
            int nClock = 1;
            while ( (pos = String_index_of_str(ptr, "local=") ) != -1 ) {
                ptr += pos + strlen("local=");
                if ( *ptr == '1') {
                    nClock++;
                }
            }
            
            if ( nClock == search->n_rate ) {
                ptr = reader->line;
                while ( *ptr != '(' ) ptr++;
                
                int i = 10-1;
                
                // insert lnl and tree if better than the last tree
                if(  lnl > lnls[i] ) {
                    free(trees[i]);
                    trees[i] = NULL;
                    
                    while( i > 0 ){
                        if ( lnl < lnls[i-1] ) {
                            trees[i] = String_clone(ptr);
                            lnls[i]  = lnl;
                            break;
                        }
                        lnls[i] = lnls[i-1];
                        trees[i] = trees[i-1];
                        i--;
                    }
                    // better than any other models
                    if ( i == 0 ) {
                        trees[i] = String_clone(ptr);
                        lnls[i]  = lnl;
                    }
                }
                
            }
        }
    }
    
    free_StringBuffer(buffer);
    free_FileReader(reader);
    
    SingleTreeLikelihood *tlk = search->pool->tlks[0];
    Node **nodes = Tree_nodes(tlk->tree);
    
    unsigned *indicators = uivector(search->n_rate);
    

    for ( int i = 1; i < 10; i++) {
        if ( isinf(lnls[i]) ) break;
        
        memcpy(indicators, search->best_indexes, (search->n_rate-1) * sizeof(unsigned) );
        
        int index = 0;
        int count = 0;
        Tree *t = new_Tree(trees[i], false);
        Node **ns = Tree_nodes(t);
        //printf("count %d\n", Tree_node_count(t));
        for ( ; index < Tree_node_count(t); index++ ) {
            if( Node_isroot(ns[i]) ) continue;
            
            if( String_index_of_str(ns[index]->info, "local=") != -1 ){
                int k = 0;
                for ( ; k < search->n_rate-1; k++ ) {
                    if ( indicators[k] == index ) {
                        break;
                    }
                }
                // we found the index not in common
                if( k == search->n_rate-1){
                    indicators[search->n_rate-1] = index;
                    break;
                }
                count++;
            }
        }

        printf(" - %s [%d] %f\n", Node_name(ns[index]), index, lnls[i]);
        free_Tree(t);
        
//        while ( (pos = String_index_of_str(ptr, "local=") ) != -1 ) {
//            ptr += pos + strlen("local=");
//            if ( *ptr == '1') {
//                int j = 0;
//                for ( ; j < search->nClock; j++ ) {
//                    if ( indicators[j] != index ) {
//                        break;
//                    }
//                }
//                // we found the index not in common
//                if( j == search->nClock){
//                    indicators[search->nClock] = index;
//                    break;
//                }
//            }
//            index++;
//        }
        
//        // test if the node is not a descendent or parent of the other nodes with a local clock
//        int j = 0;
//        for ( ; j < search->n_rate-1; j++ ) {
//            if ( Node_isrelated(nodes[index], nodes[ indicators[j] ]) ) {
//                break;
//            }
//        }
//        
//        if ( j ==  search->n_rate-1 ) {
        
            time_t start_time;
            time_t end_time;
            time(&start_time);
            
            LocalClock_set_number_of_clocks(tlk->bm, search->n_rate);
            for ( int k = 0; k < search->n_rate-1; k++ ) {
                Parameters_set_value(tlk->bm->rates, k, search->best_rates[k]);
            }
            Parameters_set_value(tlk->bm->rates, search->n_rate-1, search->best_rates[0]);
            Tree_vector_to_heights( search->best_heights, tlk->tree );
            localclock_set_indicators2(tlk->bm, indicators);
            
            
            Tree_constraint_heights( tlk->tree );
            SingleTreeLikelihood_update_all_nodes( tlk );
            
            double lk = optimize_singletreelikelihood(tlk);
            
            double new_ic = search->IC(search, lk, search->df+1);
            double old_ic = search->IC(search, search->best_lk, search->df);
            
            //double p = LRT(search->best_lk, lk, search->n_rate-1, search->n_rate);
            
            if ( search->verbosity > 1) {
                time(&end_time);
                fprintf(stdout, "Testing %s\n", Node_name(nodes[index]));
                
                double t = difftime(end_time, start_time);
                fprintf(stderr, "Runtime ");
                print_pretty_time(stderr,t);
                
                if ( search->verbosity > 1 ) {
                    fprintf(stdout, "Locations: ");
                    for (int i = 0; i < search->n_rate; i++ ) {
                        fprintf(stdout, " %s", Node_name(nodes[indicators[i]]));
                    }
                    
                    fprintf(stdout, "\n");
                }
                
                fprintf(stdout, "\n");
                printf("%s_%d = %f LnL = %f df = %d\n", INFORMATION_CRITERION[search->ic], search->n_rate, old_ic, search->best_lk, search->df);
                printf("%s_%d = %f LnL = %f df = %d\n", INFORMATION_CRITERION[search->ic], search->n_rate+1, new_ic, lk, search->df+1);
//                fprintf(stderr, "LnL0 = %f %s = %f\n", search->best_lk);
//                fprintf(stderr, "LnL1 = %f %s = %f\n", lk);
//                fprintf(stderr, "p-value = %f %s\n\n", p, (p<search->alpha ? "*" : "") );
            }
            
//            if( p > search->alpha ){
//                break;
//            }
            if ( new_ic > old_ic ) {
				break;
			}
            
            if( search->logFile != NULL ){
                fprintf(search->logFile, "tree TREE%d [&LnL=%f] = [&R] ", search->tree_count++, lk);
                Node **nodes = Tree_nodes(tlk->tree);
                StringBuffer * buffer = new_StringBuffer(10);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buffer);
                        StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buffer->c);
                        
                        if ( tlk->bm->indicators[i] ) {
                            Node_set_annotation(nodes[i], "local", "1");
                        }
                    }
                }
                Tree_print_nexus_with_annotation(search->logFile, tlk->tree);
                fprintf(search->logFile, "\n");
                free_StringBuffer(buffer);
            }
            
            search->best_lk = lk;
            search->df++;
            
            search->current_indexes = realloc( search->current_indexes, search->n_rate * sizeof(unsigned));
            search->best_indexes    = realloc( search->best_indexes,    search->n_rate * sizeof(unsigned));
            
            memcpy(search->best_indexes,    indicators, search->n_rate * sizeof(unsigned) );
            memcpy(search->current_indexes, indicators, search->n_rate * sizeof(unsigned) );
            
            Tree_heights_to_vector(tlk->tree, search->best_heights);
            memcpy(search->current_heights, search->best_heights, Tree_node_count(tlk->tree) * sizeof(double) );
            
            
            search->n_rate++;
            
            search->best_rates    = realloc( search->best_rates,    search->n_rate * sizeof(double));
            search->current_rates = realloc( search->current_rates, search->n_rate * sizeof(double));
            
            BranchModel_rates_to_vector(tlk->bm, search->best_rates);
            memcpy(search->current_rates, search->best_rates, search->n_rate * sizeof(double) );
            
            indicators = realloc( indicators, search->n_rate * sizeof(unsigned));
        //}
        
    }
    
    free(indicators);
    free(lnls);
    for ( int i = 0; i < 10; i++ ) {
        if( trees[i] != NULL ){
            free(trees[i]);
        }
    }
    free(trees);
}

#pragma mark -
// MARK: Greedy search

#ifdef OPENMP_TASK_ENABLED
//static int busy;

void _localclock_greedy_search_task( ClockSearch *search, const int current, const int start, int *count ) {
	
		if ( current == search->n_rate-1 ) {
			unsigned *indexes_copy = clone_uivector(search->current_indexes, search->n_rate-1);
			
			
			#pragma omp task firstprivate(indexes_copy)
			{
				int tid = 0;
				
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
				SingleTreeLikelihood *tlk = localclock_reset_paramters( search, tid );
				localclock_set_indicators2( tlk->bm, indexes_copy );
                //if( tlk->bm->indicators[Node_id(Node_left(Tree_root(tlk->tree)) )] )
                
				double lk = optimize_singletreelikelihood(tlk);
                //tlk->opt.verbosity = 0;
				
				#pragma omp critical
				{

					(*count)++;
				
					if( search->logFile != NULL ){
                        fprintf(search->logFile, "tree TREE%d [&LnL=%f", search->tree_count++, lk);
                        
                        for(int i = 0; i < search->n_rate-1; i++ ){
                            fprintf(search->logFile,",C%d=%d", i, indexes_copy[i]);
                        }
                        fprintf(search->logFile, "] = [&R] ");
                        
                        Node **nodes = Tree_nodes(tlk->tree);
                        StringBuffer * buffer = new_StringBuffer(10);
                        for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                            Node_empty_annotation(nodes[i]);
                            if( !Node_isroot(nodes[i]) ){
                                StringBuffer_empty(buffer);
                                StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                                Node_set_annotation(nodes[i], "rate", buffer->c);
                                
                                if ( tlk->bm->indicators[i] ) {
                                    Node_set_annotation(nodes[i], "local", "1");
                                }
                            }
                        }
                        Tree_print_nexus_with_annotation(search->logFile, tlk->tree);
                        free_StringBuffer(buffer);
                        
						fprintf(search->logFile, "\n");
					}
					
					if ( lk > search->best_lk ){
						search->best_lk = lk;
						Tree_heights_to_vector(tlk->tree, search->best_heights);
						BranchModel_rates_to_vector(tlk->bm, search->best_rates);
						memcpy(search->best_indexes, indexes_copy, (search->n_rate-1) * sizeof(unsigned) );
					}
					
					if ( search->verbosity > 1) {
						fprintf(stdout, "Iteration %d/%d (lnL: %f)\r", *count, Tree_node_count(search->pool->tlks[0]->tree)-search->n_rate+1, search->best_lk);
					}
	
				}
				free(indexes_copy);
			}
			
		} 
		else {
            // We dont test the root and the right node
            int root_id  = Node_id( Tree_root(search->pool->tlks[0]->tree));
            int right_id = Node_id( Node_right(Tree_root(search->pool->tlks[0]->tree)) );
            int left_id = Node_id( Node_left(Tree_root(search->pool->tlks[0]->tree)) );
            bool child = false;
            
			int i = ( current == start ? 0 : search->current_indexes[current-1]+1);
			for ( ; i < Tree_node_count(search->pool->tlks[0]->tree); i++) {
                if ( i == root_id || (i == right_id && start == 0) ) continue;
                
				// check if local clock has already been assigned
				int j = 0;
				for ( ; j < current; j++) {
                    if ( search->current_indexes[j] == left_id || search->current_indexes[j] == right_id ) child = true; // one the child of the root is alreay assigned
					if ( i == search->current_indexes[j] ) break;
				}
				if( j != current ) continue;
				if ( child && (i == left_id || i == right_id) ) continue; // one child done and the selected one is the other one
                    
				search->current_indexes[current] = i;
				_localclock_greedy_search_task( search, current + 1, start, count );
			}
		}
}

#else

// current is an index (start from 0)
void _localclock_greedy_search( ClockSearch *search, const int current ) {
	int count = 0;
    
    int idx = 0;
    int tot = Tree_node_count(search->pool->tlks[0]->tree)-search->n_rate+1;
    
#pragma omp parallel for shared(idx,tot,count,search)
	for ( int i = 0; i < tot; i++ ) {
		unsigned *indexes = uivector(search->n_rate-1);
        
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        int root_id = Node_id( Tree_root(search->pool->tlks[tid]->tree));
        
#pragma omp critical
		{
            if(current != 0 ){
                memcpy(indexes, search->current_indexes, (search->n_rate-2) * sizeof(unsigned) );
                
                while(1){
                    int j = 0;
                    for ( j = 0; j < current; j++ ) {
                        if ( idx == indexes[j] ) {
                            idx++;
                            break;
                        }
                    }
                    // not done yet
                    if ( j == current ) {
                        break;
                    }
                }
            }
            
            indexes[current] = idx++;
		}
        
        // do not test the root
        if ( indexes[current] == root_id ) {
            free(indexes);
            continue;
        }
        
        SingleTreeLikelihood *tlk = localclock_reset_paramters( search, tid );
        localclock_set_indicators2( tlk->bm, indexes );
        
        SingleTreeLikelihood_update_all_nodes(tlk);
        
        double lk = optimize_singletreelikelihood(tlk);
        
        
#pragma omp critical
        {
            count++;
            if( search->logFile != NULL ){
                fprintf(search->logFile, "tree TREE%d [&LnL=%f] = [&R] ", search->tree_count++, lk);
                
                Node **nodes = Tree_nodes(tlk->tree);
                StringBuffer * buffer = new_StringBuffer(10);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buffer);
                        StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buffer->c);
                        
                        if ( tlk->bm->indicators[i] ) {
                            Node_set_annotation(nodes[i], "local", "1");
                        }
                    }
                }
                Tree_print_nexus_with_annotation(search->logFile, tlk->tree);
                free_StringBuffer(buffer);
                fprintf(search->logFile, "\n");
            }
            
            if ( lk > search->best_lk ){
                search->best_lk = lk;
                Tree_heights_to_vector(tlk->tree, search->best_heights);
                BranchModel_rates_to_vector(tlk->bm, search->best_rates);
                memcpy(search->best_indexes, indexes, (search->n_rate-1) * sizeof(unsigned) );
            }
            
            if ( search->verbosity > 1) {
                fprintf(stdout, "Iteration %d/%d (lnL: %f)\r", count, Tree_node_count(search->pool->tlks[0]->tree)-search->n_rate+1, search->best_lk);
                //fprintf(stderr, "Iteration %d/%d (lnL: %f) %s\n", count, Tree_node_count(search->pool->tlks[0]->tree)-search->n_rate+1, search->best_lk, Tree_node(tlk->tree, indexes[0])->name);
            }
            
            
        }
        
		free(indexes);
	}
	
	
	
}
#endif

#ifdef PTHREAD_ENABLED
typedef struct threadpool_t{
    pthread_t *threads;
    pthread_mutex_t lock;
    ClockSearch *search;      // read-only
    int idx;                  // read-write
    int current;              // read-only in threads
    int index_treelikelihood; // read-write
    int count; // read-write
    int schedule;
}threadpool_t;

// current is the position in indexes that we want to test. If there 3 local clock then current == 2
static void * _threadpool_work( void *threadpool ){
    threadpool_t *pool = (threadpool_t *)threadpool;
    ClockSearch *search = pool->search;
    
    int index_treelikelihood = 0;
    pthread_mutex_lock(&(pool->lock));
    index_treelikelihood = pool->index_treelikelihood++;
    pthread_mutex_unlock(&(pool->lock));
    
    int root_id = Node_id( Tree_root(search->pool->tlks[index_treelikelihood]->tree));
    
    int tot = Tree_node_count(search->pool->tlks[index_treelikelihood]->tree)-search->n_rate+1;
    unsigned *indexes = uivector(search->n_rate-1);
    while( 1 ) {
        
        pthread_mutex_lock(&(pool->lock));
        pool->count++;
        if(pool->current != 0 ){
            memcpy(indexes, search->current_indexes, (search->n_rate-2) * sizeof(unsigned) );
            
            while(1){
                int j = 0;
                for ( j = 0; j < pool->current; j++ ) {
                    if ( pool->idx == indexes[j] ) {
                        pool->idx++;
                        break;
                    }
                }
                // not done yet
                if ( j == pool->current ) {
                    break;
                }
            }
        }
        
        indexes[pool->current] = pool->idx++;
    
        pthread_mutex_unlock(&(pool->lock));
        
        if ( indexes[pool->current] == root_id ) {
            continue;
        }
        
        SingleTreeLikelihood *tlk = localclock_reset_paramters( search, index_treelikelihood );
        localclock_set_indicators2( tlk->bm, indexes );
        
        SingleTreeLikelihood_update_all_nodes(tlk);
        
        double lk = optimize_singletreelikelihood(tlk);
        
        
        pthread_mutex_lock(&(pool->lock));
        {
            if( search->logFile != NULL ){
                fprintf(search->logFile, "tree TREE%d [&LnL=%f] = [&R] ", search->tree_count++, lk);
                
                Node **nodes = Tree_nodes(tlk->tree);
                StringBuffer * buffer = new_StringBuffer(10);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buffer);
                        StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buffer->c);
                        
                        if ( tlk->bm->indicators[i] ) {
                            Node_set_annotation(nodes[i], "local", "1");
                        }
                    }
                }
                Tree_print_nexus_with_annotation(search->logFile, tlk->tree);
                free_StringBuffer(buffer);
                fprintf(search->logFile, "\n");
            }
            
            if ( lk > search->best_lk ){
                search->best_lk = lk;
                Tree_heights_to_vector(tlk->tree, search->best_heights);
                BranchModel_rates_to_vector(tlk->bm, search->best_rates);
                memcpy(search->best_indexes, indexes, (search->n_rate-1) * sizeof(unsigned) );
            }
            
            if ( search->verbosity > 1) {
                fprintf(stdout, "Iteration %d/%d (lnL: %f)\r", pool->count, tot, search->best_lk);
            }
            
            
        }
        if( pool->count == tot ){
            pthread_mutex_unlock(&(pool->lock));
            break;
        }
        pthread_mutex_unlock(&(pool->lock));
        
		
	}
    free(indexes);
    pthread_exit(NULL);
    return NULL;
}

// current is an index (start from 0)
void _localclock_greedy_search_pthreads( ClockSearch *search, const int current ) {
    //printf("Using threads\n");
    threadpool_t *threadpool = malloc(sizeof(threadpool_t));
    assert(threadpool);
    threadpool->count = 0;
    threadpool->current = current;
    threadpool->idx = 0;
    threadpool->search = search;
    threadpool->index_treelikelihood = 0;
    threadpool->schedule = 4;
    
    threadpool->threads = malloc(sizeof(pthread_t)*search->nThreads);
    assert(threadpool->threads);
    
    pthread_mutex_init(&(threadpool->lock), NULL);
                       
    for ( int i = 0; i < search->nThreads; i++ ) {
        pthread_create( &(threadpool->threads[i]), NULL, _threadpool_work, threadpool );
    }
    for ( int i = 0; i < search->nThreads; i++ ) {
        pthread_join(threadpool->threads[i], NULL);
    }
    free(threadpool->threads);
    pthread_mutex_destroy(&(threadpool->lock));
    free(threadpool);
}
#endif


#pragma mark -
// MARK: Exhaustive search

void _localclock_exhaustive_search( ClockSearch *search, int *count ){
    int nnodes = Tree_node_count(search->pool->tlks[0]->tree);
    int k = search->n_rate-1;
    long n = choose(nnodes, k);
    
    int nrow = search->nThreads;
    unsigned **combinations = uimatrix(nrow, k);
    long counter = 0;
    
    int root_id = Node_id( Tree_root(search->pool->tlks[0]->tree));
    
    #pragma omp parallel for
    for ( long i = 0; i < n; i++ ) {
        
        int tid = 0;
        
#ifdef _OPENMP
        tid = omp_get_thread_num();
        //printf("tid %d %d\n",tid,omp_get_num_threads());
#endif
        
        combination_at_index(nnodes, k, i, combinations[tid] );
        
        int j = 0;
        for ( ; j < k; j++ ) {
            if( root_id == combinations[tid][j]) break;
        }
        // combinations[tid] contains at least one root assignment
        if ( j != k ) continue;
        
        SingleTreeLikelihood *tlk = localclock_reset_paramters( search, tid );
        localclock_set_indicators2( tlk->bm, combinations[tid] );
        
        double lk = optimize_singletreelikelihood(tlk);

#pragma omp critical
        {
            
            (*count)++;
            counter++;
            if( search->logFile != NULL ){
                fprintf(search->logFile, "tree TREE%d [&LnL=%f] = [&R] ", search->tree_count++, lk);
                
                Node **nodes = Tree_nodes(tlk->tree);
                StringBuffer * buffer = new_StringBuffer(10);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buffer);
                        StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buffer->c);
                        
                        if ( tlk->bm->indicators[i] ) {
                            Node_set_annotation(nodes[i], "local", "1");
                        }
                    }
                }
                Tree_print_nexus_with_annotation(search->logFile, tlk->tree);
                free_StringBuffer(buffer);
                
                //print_tree_nexus_treelikelihood(search->logFile, tlk);
                fprintf(search->logFile, "\n");
            }
            
            if ( lk > search->best_lk ){
                search->best_lk = lk;
                Tree_heights_to_vector(tlk->tree, search->best_heights);
                BranchModel_rates_to_vector(tlk->bm, search->best_rates);
                memcpy(search->best_indexes, combinations[tid], k * sizeof(unsigned) );
            }
            
            if ( search->verbosity > 1) {
                fprintf(stdout, "Iteration %ld/%ld (lnL: %f)\r", counter, n, search->best_lk);
            }
            
        }
    }
    free_uimatrix(combinations, nrow);
}



#pragma mark -

// MARK: GA
// chromosome_size == number of clocks
/*
 *********************************************************************************************************
 */



static double localclock_fitness( GA *ga, Individual *individual );
static void  localclock_mate_chc( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );
static void localclock_mutate_chc( GA *ga, Individual *individual );

static void _localclock_ga_init( GA *ga, ClockSearch *search );
static void _random_individual( Individual *indiv, int n_nodes );
static bool is_chromosome_ok( const Individual *individual, const int nLocalClock );
static double _localclock_ga_optmize( ClockSearch *search );

static double localclock_fitness2( GA *ga, Individual *individual );
static void  localclock_mate_chc2( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );
static void localclock_mutate_chc2( GA *ga, Individual *individual );
void _localclock_ga_init2( GA *ga, ClockSearch *search );
static double _localclock_ga_optmize2( ClockSearch *search );



ClockSearch * new_GALocalClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize ){
	ClockSearch *search = new_LocalClockSearch(tlk, nThreads);
	GAClockSearch *gaSearch = (GAClockSearch*)malloc( sizeof(GAClockSearch) );
	assert(gaSearch);
    
	search->pDerivedObj = gaSearch;
    gaSearch->pBaseObj = search;
	
	search->optimize = _localclock_ga_optmize;
	search->free     = free_GAClockSearch;
	search->init     = ClockSearch_init;
	
	search->algorithm = HETEROTACHY_LOCALCLOCK_GA;
	
    // GA stuff
   	gaSearch->initGA = _localclock_ga_init;
	gaSearch->popSize = popSize;
	gaSearch->ngen = 200;
	gaSearch->max_no_improvement = 50;
	
    gaSearch->tree_count = 0;
	gaSearch->verbosity = 1;
	gaSearch->feedback_sampling = 1;
    
    gaSearch->pop_individuals = NULL;
    gaSearch->pop_heights     = NULL;
    gaSearch->pop_rates       = NULL;
	
	return search;
}

ClockSearch * new_GALocalClockSearch2( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize ){
    ClockSearch *search = new_GALocalClockSearch(tlk, nThreads, popSize);
    GAClockSearch *gaSearch = (GAClockSearch*)search->pDerivedObj;
    
	
	search->optimize = _localclock_ga_optmize2;
   	gaSearch->initGA = _localclock_ga_init2;
	
	return search;
}


// when the number of nodes is lower than the pop size it might cause the init of the GA go on forever
double _localclock_ga_optmize( ClockSearch *search ){
	
	if ( !search->initialized ) {
		search->init(search);
	}
	
	Tree *tree = search->pool->tlks[0]->tree;
	BranchModel *bm = search->pool->tlks[0]->bm;
	int nNodes = Tree_node_count(tree);
    
	search->best_indexes    = realloc( search->best_indexes,    (search->n_rate-1) * sizeof(unsigned) );
	search->current_indexes = realloc( search->current_indexes, (search->n_rate-1) * sizeof(unsigned) );
    
    //TODO: best_indexes now contains id not postorder_idx
    fprintf(stderr, "Replace postorder_idx with id %s %d\n",__FILE__,__LINE__);
    exit(1);
    
    LocalClock_get_indexes(bm, search->best_indexes);
    memcpy(search->current_indexes, search->best_indexes, (search->n_rate-1) * sizeof(unsigned));
    
	Tree_heights_to_vector(tree, search->best_heights);
    memcpy(search->current_heights, search->best_heights, nNodes * sizeof(double));
    
	BranchModel_rates_to_vector( bm, search->best_rates );
    memcpy(search->current_rates, search->best_rates, search->n_rate * sizeof(double));
	
    search->default_rate = dmean(search->best_rates, Parameters_count(bm->rates) );
	
	GAClockSearch *gasearch = (GAClockSearch*)search->pDerivedObj;
	
	
	// logs all the configurations with the likelihood value and the number of times the GA generated it
	// ex: 0:1:5:6:9: -2555.2 3
	bool logHash = false;
	FILE *hashLogFile = NULL;
	if ( logHash ) {
		StringBuffer *buff = new_StringBuffer(100);
		StringBuffer_append_strings(buff, 2, search->logFileName, ".log");
		char *hashLog = StringBuffer_tochar(buff);
		free_StringBuffer(buff);
		hashLogFile = fopen(hashLog, "w");
		free(hashLog);
		assert(hashLogFile);
	}
	
	int nsites = search->pool->tlks[0]->sp->nsites;
	
	double previous_lk     = search->best_lk;
    int    previous_df     = search->starting_df;
	double previous_IC     = AICc(previous_lk, previous_df, nsites);
    int    previous_n_rate = search->n_rate;
    
	double IC = previous_IC;
	
	for ( ; search->n_rate <= search->max_n_rate; search->n_rate++ ) {
		//ga->mutation_rate = 1.75/(30*sqrt(i));
		fprintf(stdout, "Start GA with %d local clocks (max %d)\n", search->n_rate, search->max_n_rate);
        
        search->df = SingleTreeLikelihood_df_count(search->pool->tlks[0]);
        
		GA *ga = new_GA( GA_CHC, GA_CHROMOSOME_UNSIGNED, gasearch->popSize, search->n_rate );
		
		ga->fitness     = localclock_fitness;
		ga->mate        = localclock_mate_chc;
		ga->mutate      = localclock_mutate_chc;
		ga->diversity   = ga_diversity_fitness; //ga_diversity_atomic;
		ga->termination = ga_default_termination;
        
        ga->mutation_threshold =  0.0001;
		ga->mutation_rate = 0.15;
        
        ga_set_ngeneration(ga,gasearch->ngen);
        
		ga->data = search;
        
		ga->use_max_no_improvement = true;
		ga->max_no_improvement = gasearch->max_no_improvement;
        
		ga_set_n_threads(ga, search->nThreads);
		
        ga->feedback_computeall = true;
		ga->verbosity = gasearch->verbosity;
		ga->feedback_sampling = gasearch->feedback_sampling;
        
		gasearch->initGA(ga, search); // init the Ks but don't compute the fitness
		
		ga_init(ga);
        
		ga_evolve(ga);
		
		search->best_lk = ga->maxfitness;
		
        IC = AICc(search->best_lk, search->df, nsites);
        
		// this is the first round we just want the LnL to be better (should not be absoute)
		if ( search->n_rate == search->starting_n_rate ) {
			if( search->best_lk < previous_lk ){
				free_GA(ga);
				break;
			}
		}
		else {
			fprintf(stdout, "\n");
			fprintf(stdout, "AIC_%d = %f LnL = %f\n", previous_n_rate, previous_IC, previous_lk);
			fprintf(stdout, "AIC_%d = %f LnL = %f\n", search->n_rate, IC, search->best_lk);
			fprintf(stdout, "==============================================\n\n");
			
			if ( IC > previous_IC ) {
				if ( logHash ) {
					fprintf(hashLogFile, "# %d\n", search->n_rate);
					print_ModelTally_Hashtable(hashLogFile, ga->lookup , search->n_rate);
					fprintf(hashLogFile, "\n");
				}
				
				free_GA(ga);
				break;
			}
		}
		
		previous_df     = search->df;
		previous_lk     = search->best_lk;
		previous_IC     = IC;
        previous_n_rate = search->n_rate;
        
        memcpy(search->current_indexes, search->best_indexes,(search->n_rate-1) * sizeof(unsigned) );
		search->current_indexes = realloc( search->current_indexes, search->n_rate * sizeof(unsigned));
		search->best_indexes    = realloc( search->best_indexes,    search->n_rate * sizeof(unsigned));
		search->current_indexes[search->n_rate-1] = search->best_indexes[search->n_rate-1] = 0;
		
		memcpy(search->current_heights, search->best_heights, nNodes * sizeof(double) );
		
		search->best_rates    = realloc( search->best_rates,    (search->n_rate+1) * sizeof(double) );
		search->current_rates = realloc( search->current_rates, (search->n_rate+1) * sizeof(double) );
		// assign the ancestral rate to the new local clock (could inherit the rate from it's parent)
		search->best_rates[search->n_rate] = search->best_rates[0];
		memcpy(search->current_rates, search->best_rates, (search->n_rate+1) * sizeof(double) );
		
		if ( logHash ) {
			fprintf(hashLogFile, "# %d\n", search->n_rate);
			print_ModelTally_Hashtable(hashLogFile, ga->lookup , search->n_rate);
			fprintf(hashLogFile, "\n");
		}
		
		free_GA(ga);
	}
    
	if ( search->n_rate == search->max_n_rate){
		fprintf(stderr, "The potential number of local clocks (%d) is greater than the maximum (%d)\n", search->n_rate, search->max_n_rate );
	}
    // the algorithm stopped after the first iteration and the previous model was another local clock
    else if ( search->n_rate == search->starting_n_rate ){
		
	}
	else {
        
        // from here I assume that the search has found at least one local clock
        // should do somthing about it
        
        // save the best heights and rates
        search->n_rate--;
        
        memcpy(search->best_heights, search->current_heights, nNodes * sizeof(double) );
        search->best_rates = realloc( search->best_rates, search->n_rate * sizeof(double));
        memcpy(search->best_rates, search->current_rates, search->n_rate * sizeof(double) );
        
        search->best_indexes = realloc( search->best_indexes, (search->n_rate-1) * sizeof(unsigned));
        memcpy(search->best_indexes, search->current_indexes, (search->n_rate-1) * sizeof(unsigned) );
        
        search->best_lk = previous_lk;
    }
    
    // reset BranchModel
    LocalClock_set_number_of_clocks(bm, search->n_rate-1);
    localclock_set_indicators2(bm, search->best_indexes);
    BranchModel_vector_to_rates(bm, search->best_rates);
    // reset tree heights in order to get the root height
    Tree_vector_to_heights( search->best_heights, tree);
    Tree_constraint_heights(tree);
    SingleTreeLikelihood_update_all_nodes(search->pool->tlks[0]);
    
	
	if( search->logFile != NULL ){
		// we only write the header if it is on write mode.
		// If we are in append mode then the header will be printed by the caller
		if ( search->logFileMode[0] == 'w') {
			fprintf(search->logFile, "\nEnd;\n");
		}
		
		fclose(search->logFile);
	}
	
	if ( logHash ) {
		fclose(hashLogFile);
	}
    
	return search->best_lk;
}

double _localclock_ga_optmize2( ClockSearch *search ){
	
	if ( !search->initialized ) {
		search->init(search);
	}
	
	Tree *tree = search->pool->tlks[0]->tree;
	BranchModel *bm = search->pool->tlks[0]->bm;
	int nNodes = Tree_node_count(tree);
    
    memcpy(search->best_indexes, bm->indicators, nNodes);
    
	Tree_heights_to_vector(tree, search->best_heights);
    memcpy(search->current_heights, search->best_heights, nNodes * sizeof(double));
    
	BranchModel_rates_to_vector( bm, search->best_rates );
    search->default_rate = dmean(search->best_rates, Parameters_count(bm->rates) );
	
	GAClockSearch *gasearch = (GAClockSearch*)search->pDerivedObj;
	
	
	// logs all the configurations with the likelihood value and the number of times the GA generated it
	// ex: 0:1:5:6:9: -2555.2 3
	bool logHash = false;
	FILE *hashLogFile = NULL;
	if ( logHash ) {
		StringBuffer *buff = new_StringBuffer(100);
		StringBuffer_append_strings(buff, 2, search->logFileName, ".log");
		char *hashLog = StringBuffer_tochar(buff);
		free_StringBuffer(buff);
		hashLogFile = fopen(hashLog, "w");
		free(hashLog);
		assert(hashLogFile);
	}
	
    //ga->mutation_rate = 1.75/(30*sqrt(i));
    fprintf(stdout, "Start GA with local clocks\n");
    
    
    GA *ga = new_GA( GA_CHC, GA_CHROMOSOME_UNSIGNED, gasearch->popSize, nNodes );
    
    ga->fitness     = localclock_fitness2;
    ga->mate        = localclock_mate_chc2;
    ga->mutate      = localclock_mutate_chc2;
    ga->diversity   = ga_diversity_fitness; //ga_diversity_atomic;
    ga->termination = ga_default_termination;
    
    ga->mutation_threshold =  0.0001;
    ga->mutation_rate = 0.15;
    
    ga_set_ngeneration(ga,gasearch->ngen);
    
    ga->data = search;
    
    ga->use_max_no_improvement = true;
    ga->max_no_improvement = gasearch->max_no_improvement;
    
    ga_set_n_threads(ga, search->nThreads);
    
    ga->feedback_computeall = true;
    ga->verbosity = gasearch->verbosity;
    ga->feedback_sampling = gasearch->feedback_sampling;
    
    ga->termination_callback = _localclock_termination_callback;
    
    gasearch->initGA(ga, search); // init the Ks but don't compute the fitness
    
    ga_init(ga);
    
    ga_evolve(ga);
    
    if ( logHash ) {
        fprintf(hashLogFile, "# %d\n", search->n_rate);
        print_ModelTally_Hashtable(hashLogFile, ga->lookup , search->n_rate);
        fprintf(hashLogFile, "\n");
    }
    
    free_GA(ga);
	
    
    
    // reset BranchModel
    LocalClock_set_number_of_clocks(bm, search->n_rate);
    bool *chr = bvector(nNodes);
	for ( int i = 0; i < nNodes; i++) {
        if(search->best_indexes[i] == 1){
            chr[i] = true;
        }
	}
    localclock_set_indicators(bm, chr);
    free(chr);
    
    BranchModel_vector_to_rates(bm, search->best_rates);
    // reset tree heights in order to get the root height
    Tree_vector_to_heights( search->best_heights, tree);
    Tree_constraint_heights(tree);
    SingleTreeLikelihood_update_all_nodes(search->pool->tlks[0]);
    
	
	if( search->logFile != NULL ){
		// we only write the header if it is on write mode.
		// If we are in append mode then the header will be printed by the caller
		if ( search->logFileMode[0] == 'w') {
			fprintf(search->logFile, "\nEnd;\n");
		}
		
		fclose(search->logFile);
	}
	
	if ( logHash ) {
		fclose(hashLogFile);
	}
    
	return search->best_lk;
}

// when the number of nodes is lower than the pop size it might cause the init of the GA go on forever
void _localclock_ga_init( GA *ga, ClockSearch *search ){
    
    GAClockSearch *gaSearch = (GAClockSearch*)search->pDerivedObj;
    
    int node_count = Tree_node_count(search->pool->tlks[0]->tree);
	
    Hashtable *hash    = new_Hashtable_string(ga->pop_size);
	StringBuffer *buff = new_StringBuffer( ga->chromosome_size *2 );
    
    unsigned int *chromosome = NULL;
    
    int i = 0;
	if ( search->n_rate == search->starting_n_rate ) {
        
        if ( gaSearch->pop_individuals == NULL ) {
            memcpy( ga->individuals[0]->chromosome, search->best_indexes, ga->chromosome_size * sizeof(unsigned) );
            if ( is_chromosome_ok( ga->individuals[i], node_count) ){
                qsort(ga->individuals[i]->chromosome, ga->individuals[i]->size, sizeof(unsigned int), qsort_asc_uivector );

                chromosome = ga->individuals[i]->chromosome;
                
                for ( int j = 0; j < ga->chromosome_size; j++) {
                    StringBuffer_append_format(buff, "%u:", chromosome[j]);
                }
                fprintf(stderr, "%d %s\n", 0, buff->c);
                Hashtable_add(hash, String_clone(buff->c), new_Double(1));
                ++i;
            }
        }
        else {
            for ( int i = 0; i < ga->pop_size; i++ ) {
                memcpy( ga->individuals[i]->chromosome, gaSearch->pop_individuals[i], ga->chromosome_size * sizeof(unsigned) );
            }
            free(gaSearch->pop_individuals);
            free(gaSearch->pop_heights);
            free(gaSearch->pop_rates);
            free_StringBuffer(buff);
            free_Hashtable(hash);
            
            return;
        }
    }
    else {
        for (int j = 0; j < search->pool->count; j++) {
			LocalClock_set_number_of_clocks(search->pool->tlks[j]->bm, search->n_rate-1);
		}
        memcpy( ga->individuals[0]->chromosome, search->best_indexes, (ga->individuals[0]->size-1) * sizeof(unsigned int) );
		int j = 0;
        
        chromosome = ga->individuals[0]->chromosome;
        
		while(  j != ga->individuals[0]->size-1 ){
			chromosome[ga->individuals[0]->size-1] = random_int(node_count-2);
			for ( j = 0; j < ga->individuals[0]->size-1; j++) {
				if ( chromosome[ga->individuals[0]->size-1] == chromosome[j] ) {
					break;
				}
			}
		}
		i++;
    }
    
	
	int counter = 0; // avoid infinite loops
	for ( ; i < ga->pop_size; i++ ) {
        chromosome = ga->individuals[i]->chromosome;
		counter = 0;
		while (1){
			StringBuffer_empty(buff);
			_random_individual( ga->individuals[i], node_count-1 );
			
			if ( !is_chromosome_ok( ga->individuals[i], node_count) ) continue;
			
			qsort(ga->individuals[i]->chromosome, ga->individuals[i]->size, sizeof(unsigned), qsort_asc_uivector );
			
			for ( int j = 0; j < ga->chromosome_size; j++) {
				StringBuffer_append_format(buff, "%u:", chromosome[j]);
			}
			
			if( !Hashtable_exists( hash, buff->c) || counter == 10 ){				
				Hashtable_add(hash, String_clone(buff->c), new_Double(1));
                fprintf(stderr, "%d %s\n", 0, buff->c);
				break;
			}
			counter++;
		}
	}
		
	free_StringBuffer(buff);
	free_Hashtable(hash);
}


void _random_individual( Individual *indiv, int n_nodes ){
	memset(indiv->chromosome, 0, indiv->size * sizeof(int) );
	
    unsigned int *chromosome = indiv->chromosome;
    
	int location = 0;
	for ( int i = 0; i < indiv->size; i++) {
		
		while( 1 ){
			location = random_int(n_nodes-1);
			int j = 0;
			for ( ; j < i; j++) {
				if ( location == chromosome[j] ) {
					break;
				}
			}
			if( i == j ){
				chromosome[i] = location;
				break;
			}
		}
	}
	
}


/* might revert to an existing clock after a second mutation since
 I do not keep the initial chromose. could create a temporary
 chromosome */
void localclock_mutate_chc( GA *ga, Individual *individual ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
	int nNodes = Tree_node_count(searcher->pool->tlks[0]->tree);
	
	unsigned *backup = uivector(individual->size);
    unsigned int *chromosome = individual->chromosome;
    
	memcpy(backup, chromosome, individual->size *sizeof(unsigned) );
	   
	int n_mutation = 0;
	while ( 1 ){
		n_mutation = 0;
		for (int ith = 0; ith < individual->size; ith++) {
			double rnum = random_double();
			if( rnum < ga->mutation_rate ){
				n_mutation++;
				int i = 0;
				int b = 0;
				while ( i != individual->size ) {
					b = random_int(nNodes-2);
					for ( i = 0; i < individual->size; i++) {
						if ( chromosome[i] == b ) break;
					}
					
				}
				chromosome[ith] = b;
				
			}
		}
		if ( is_chromosome_ok( individual, nNodes-1) ) break;
		memcpy( chromosome, backup, individual->size *sizeof(unsigned) );
	}
	
	qsort(chromosome, individual->size, sizeof(unsigned), qsort_asc_uivector );
	
	free(backup);
}


// Switch ~50% different local clock positions
//#define DEBUG_LOCALCLOCK_MATE
void  localclock_mate_chc( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
	int i = 0;
	int counter = 0;
	int nclocks = individual1->size;
	
	
	bool *mask1 = bvector(nclocks);
	bool *mask2 = bvector(nclocks);
	
	memset(mask1, false, nclocks*sizeof(bool) );
	memset(mask2, false, nclocks*sizeof(bool) );
    
    unsigned int *chromosome1 = individual1->chromosome;
    unsigned int *chromosome2 = individual2->chromosome;
    unsigned int *newchromosome = new1->chromosome;
	
	// check the number of clocks that can be swapped
	for ( ; i < nclocks; i++) {
		int j = 0;
		for ( ; j < nclocks; j++) {
			if ( chromosome1[i] == chromosome2[j] ) {
				break;
			}
		}
		if ( j == nclocks ){
			counter++;
		}
		else {
			mask1[i] = true;
			mask2[j] = true;
		}

	}
	
	memcpy(newchromosome, chromosome1, individual1->size * sizeof(unsigned int));
		
	int n_mate = 0;
	
	// same indivs or only 1 clock not in common (it is like swaping the chromosomes of the 2 individuals)
	if( counter == 0 || counter == 1 ){ 
		free(mask1);
		free(mask2);
		return;
	}
		
		
	for ( i = 0; i < counter; i++) {
		double rnum = random_double();
		if ( rnum < ga->mate_prob ) {
			// choose position of the first element of the array to swap
			int p1 = 0;
			int p1counter = 0;
			for ( ; p1 < nclocks; p1++) {
				if ( mask1[p1] == false ) {
					if (p1counter == i - n_mate ) break;
					else p1counter++;
				}
			}
			
			int p2 = random_int(nclocks-1);
			while ( mask2[p2] == true || chromosome1[p1] == chromosome2[p2] ) {
				p2 = random_int(nclocks-1);
			}

			mask1[p1] = true;
			mask2[p2] = true;
						
			newchromosome[p1] = chromosome2[p2];
			n_mate++;
		}
		
	}
	
	
	free(mask1);
	free(mask2);
	
	qsort(newchromosome, new1->size, sizeof(unsigned), qsort_asc_uivector );
}

static SingleTreeLikelihood * _localclock_reset_paramters_( ClockSearch *searcher, unsigned index ){
	Tree_vector_to_heights( searcher->best_heights, searcher->pool->tlks[index]->tree );
	
	for (int p = 0; p < Parameters_count(searcher->pool->tlks[index]->bm->rates); p++) {
		Parameters_set_value(searcher->pool->tlks[index]->bm->rates, p, searcher->default_rate);
	}
	
	Tree_constraint_heights( searcher->pool->tlks[index]->tree );
	SingleTreeLikelihood_update_all_nodes( searcher->pool->tlks[index] );
	
	return searcher->pool->tlks[index];
}

double localclock_fitness( GA *ga, Individual *individual ){
    
    unsigned int *chromosome = individual->chromosome;
    
	StringBuffer *buff = new_StringBuffer( individual->size*2 );
	int i = 0;
	for ( i = 0; i < individual->size; i++) {
		StringBuffer_append_format(buff, "%u:", chromosome[i]);
	}
    
	ModelTally *pFit = NULL;
    bool found = false;

	#pragma omp critical
	{
		
        found = Hashtable_exists(ga->lookup, buff->c);
        
        if ( found ){
            ModelTally *tally = Hashtable_get(ga->lookup, buff->c);
			individual->fitness = tally->lk; // at this point tally->lk can be NAN
			tally->count++;
			
			if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
				ga->feedback_count--;
				fprintf(stderr, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
			}
		}
		else {
			pFit = new_ModelTally(NAN);
			Hashtable_add( ga->lookup, String_clone(buff->c), pFit);
		}		
	}
	
	if ( !found ) {
        
		int tid = 0;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		
		ClockSearch *searcher = (ClockSearch*)ga->data;
		
		SingleTreeLikelihood *tlk = _localclock_reset_paramters_( searcher, tid );
		
		localclock_set_indicators2( tlk->bm, individual->chromosome );
				
		individual->fitness = optimize_singletreelikelihood(tlk);
		
		#pragma omp critical
		{
			if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
				ga->feedback_count--;
				fprintf(stderr, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
			}
			
			if( individual->fitness > ga->maxfitness ){
				Tree_heights_to_vector(tlk->tree, searcher->best_heights);
				BranchModel_rates_to_vector(tlk->bm, searcher->best_rates);
				memcpy(searcher->best_indexes, individual->chromosome, individual->size * sizeof(unsigned) );
				ga->maxfitness = individual->fitness;
			}
			pFit->lk = individual->fitness;
			if( searcher->logFile != NULL ){
				double aic = searcher->IC(searcher, pFit->lk, searcher->df);//AICc(pFit->lk, searcher->df, tlk->sp->nsites);
				fprintf(searcher->logFile, "tree TREE%d [&LnL=%f,IC=%f] = [&R] ", searcher->tree_count++, individual->fitness, aic);
				
                Node **nodes = Tree_nodes(tlk->tree);
                StringBuffer * buffer = new_StringBuffer(10);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buffer);
                        StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buffer->c);
                        
                        if ( tlk->bm->indicators[i] ) {
                            Node_set_annotation(nodes[i], "local", "1");
                        }
                    }
                }
                Tree_print_nexus_with_annotation(searcher->logFile, tlk->tree);
                free_StringBuffer(buffer);
                
                //print_tree_nexus_treelikelihood(searcher->logFile, tlk);
				fprintf(searcher->logFile, "\n");
			}
		}
	}
	
	free_StringBuffer( buff );	
	
	return individual->fitness;
}



// Check if a local cock is assigned at most once
static bool is_chromosome_ok( const Individual *individual, const int nLocalClock ){
    
    unsigned int *chromosome = individual->chromosome;
    
	unsigned char *map = (unsigned char *)calloc( nLocalClock, sizeof(unsigned char) );
	assert(map);
	int i = 0;
	for ( i = 0; i < individual->size; i++) { // ignore the root
		map[ chromosome[i] ]++;
	}
	
	for ( i = 0; i < nLocalClock; i++) {
		if ( map[i] > 1 ){
			free(map);
			return false;
		}
	}
	free(map);
	return true;
	
}

#pragma mark -

bool isOK(unsigned * chromosome, int length ){
    int zero_count = 0;
    for( int i = 0; i < length; i++ ){
        if( chromosome[i] == 0 ) zero_count++;
    }
    return zero_count!= length;
}

void _localclock_ga_init2( GA *ga, ClockSearch *search ){
    
    int node_count = Tree_node_count(search->pool->tlks[0]->tree);
	
    unsigned int *chromosome = NULL;
    
    
    for (int j = 0; j < ga->pop_size; j++) {
        chromosome = ga->individuals[j]->chromosome;
        
        memset(chromosome, 0, sizeof(unsigned)*node_count);
        int local_count = random_int2(1, 4);
        
        for ( int i = 0; i < local_count; i++ ) {
            int pos = random_int(node_count-1);
            while( chromosome[pos] == 1 ){
                pos = random_int(node_count-1);
            }
            chromosome[pos] = 1;
        }
    }
}

// the number of local clock is between 1 and the number of local clock in the current best model.
// It avoids getting too many rates
// Could cluster some rates in order to get fewer parameters
void localclock_mutate_chc4( GA *ga, Individual *individual ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
    unsigned int *chromosome = individual->chromosome;
    
    fprintf(searcher->logFile,"-\n");
    
    memset(chromosome, 0, sizeof(unsigned)*individual->size);
    int min = imax(searcher->n_rate-2, 1);
    int local_count = random_int2(min, searcher->n_rate);
    
    for ( int i = 0; i < local_count; i++ ) {
        int pos = random_int(individual->size-1);
        while( chromosome[pos] == 1 ){
            pos = random_int(individual->size-1);
        }
        chromosome[pos] = 1;
    }
}

void localclock_mutate_chc2( GA *ga, Individual *individual ){
    unsigned int *chromosome = individual->chromosome;
    
    int ones = 0;
    int *pos = ivector(individual->size);
    for ( int i = 0; i < individual->size; i++ ) {
        if(chromosome[i] == 1 ){
            pos[ones] = i;
            ones++;
        }
    }

    // choose number of indicators to change
    int n = imax(ga->mutation_rate * ones, dmin(ones,2));
   
    for ( int i = 0; i < n; i++ ) {
        // choose a position with a 1
        int rand = random_int(ones-1);
        int p1 = pos[rand];
        while ( p1 == -1 ) {
            rand = random_int(ones-1);
            p1 = pos[rand];
        }
        
        // choose a position with a 0
        int p2 = random_int(individual->size-1);
        while( p1 == p2 || chromosome[p2] == 1 ){
            p2 = random_int(individual->size-1);
        }
        // swap
        swap_uint(&chromosome[p1], &chromosome[p2]);
        pos[rand] = -1;
    }
    free(pos);
}

void localclock_mutate_chc3( GA *ga, Individual *individual ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
    unsigned int *chromosome = individual->chromosome;
    unsigned *backup = uivector(individual->size);
	memcpy(backup, individual->chromosome, individual->size *sizeof(unsigned) );
    int root_id = Node_id(Tree_root(searcher->pool->tlks[0]->tree));
    
	while ( 1 ){
		for (int i = 0; i < individual->size; i++) {
            // we don't need to mutate the root
            if( i == root_id ) continue;
            
			double rnum = random_double();
			if( rnum < ga->mutation_rate ){
				chromosome[i] = 1-chromosome[i];
			}
		}
        if ( root_id == 0 ) {
            chromosome[0] = chromosome[1];
        }
        else {
            chromosome[root_id] = chromosome[root_id-1];
        }
        
        if( isOK(chromosome, individual->size) ){
            break;
        }
        memcpy( chromosome, backup, individual->size *sizeof(unsigned) );
	}
    free(backup);
}

// echange indicators between adjacent lineages
void localclock_mutate_chc5( GA *ga, Individual *individual ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
    Node  **nodes = Tree_nodes(searcher->pool->tlks[0]->tree);
    unsigned int *chromosome = individual->chromosome;
    
    //fprintf(searcher->logFile,"-\n");
    
    int ones = 0;
    int *pos = ivector(individual->size);
    for ( int i = 0; i < individual->size; i++ ) {
        if(chromosome[i] == 1 ){
            pos[ones] = i;
            ones++;
        }
    }
    
    // choose number of indicators to change
    int n = imax(ga->mutation_rate * ones, dmin(ones,2));
    
    for ( int i = 0; i < n; i++ ) {
        // choose a position with a 1
        int rand = random_int(ones-1);
        int p1 = pos[rand];
        while ( p1 == -1 ) {
            rand = random_int(ones-1);
            p1 = pos[rand];
        }
        
        Node *node = nodes[p1];
        int r = random_int(1);
        // exchange with its parent
        if( r < 0.5 && !Node_isroot(Node_parent(node)) ){
            swap_uint(&chromosome[p1], &chromosome[Node_id(Node_parent(node))]);
        }
        // exchange with its sibling
        else {
            swap_uint(&chromosome[p1], &chromosome[Node_id(Node_sibling(node))]);
        }
        
        pos[rand] = -1;
    }
    free(pos);
}

void  localclock_mate_chc2( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
	
    //    IndividualData *new_data = (IndividualData*)new1->data;
    ClockSearch *searcher = (ClockSearch*)ga->data;
    int count_failure = 0;
    
    int root_id = Node_id(Tree_root(searcher->pool->tlks[0]->tree));
    
    Individual *best  = individual2;
    Individual *worst = individual1;
    if ( individual1->fitness > individual2->fitness ) {
        best  = individual1;
        worst = individual2;
    }
    
    unsigned int *newchromosome = new1->chromosome;
    unsigned int *bestchromosome = best->chromosome;
    unsigned int *worstchromosome = worst->chromosome;
    
    
    
	do {
		// we need to reset the K otherwise at each iteration it will look more and more like the best chromosome
        memcpy(newchromosome, bestchromosome, best->size * sizeof(unsigned int));
        
		for ( int i = 0; i < best->size; i++) {
            // we don't need to mate the root
            if( i == root_id ) continue;
            
			double rnum = random_double();
			if ( rnum < ga->mate_prob ) {
				newchromosome[i] = worstchromosome[i];
			}
		}
        
        if ( root_id == 0 ) {
            newchromosome[0] = newchromosome[1];
        }
        else {
            newchromosome[root_id] = newchromosome[root_id-1];
        }
        count_failure++;
        if(count_failure == 20) break;
	}
    while ( !isOK(newchromosome, new1->size) );
    
    // if individuals are not compatible for mating we just mutate the best one
    if( count_failure == 20 ){
        memcpy(newchromosome, bestchromosome, best->size * sizeof(unsigned int));
        localclock_mutate_chc2(ga, best);
    }
}

double localclock_fitness2( GA *ga, Individual *individual ){
    
    unsigned int *chromosome = individual->chromosome;
    bool *chr = bvector(individual->size);
    
	StringBuffer *buff = new_StringBuffer( individual->size*3 );
	int i = 0;
    int local_count = 0;
	for ( i = 0; i < individual->size; i++) {
		StringBuffer_append_format(buff, "%u:", chromosome[i]);
        if(chromosome[i] == 1){
            chr[i] = true;
            local_count++;
        }
	}
	
	ModelTally *pFit = NULL;
    bool found = false;
	
#pragma omp critical
	{
        found = Hashtable_exists(ga->lookup, buff->c);
		
		if ( found ){
            //fprintf(stderr, "fund\n");
            ModelTally *tally = Hashtable_get(ga->lookup, buff->c);
			individual->fitness = tally->lk; // at this point tally->lk can be NAN
			tally->count++;
			
			if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
				ga->feedback_count--;
				fprintf(stderr, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
			}
		}
		else {
			pFit = new_ModelTally(NAN);
			Hashtable_add( ga->lookup, String_clone(buff->c), pFit);
		}
        
	}
	
	if ( !found ) {
		
		int tid = 0;
#ifdef _OPENMP
		tid = omp_get_thread_num();
		//printf("tid %d %d\n",tid,omp_get_num_threads());
#endif
		
		ClockSearch *searcher = (ClockSearch*)ga->data;
        
        SingleTreeLikelihood *tlk = searcher->pool->tlks[tid];
        
        //IndividualData *data = (IndividualData*)individual->data;
        
        // heights
        //Tree_vector_to_heights( data->heights, tlk->tree );
        Tree_vector_to_heights( searcher->best_heights, tlk->tree );
        Tree_constraint_heights( tlk->tree );
        SingleTreeLikelihood_update_all_nodes( tlk );
        
        // rates
        LocalClock_set_number_of_clocks(tlk->bm, local_count);
        localclock_set_indicators(tlk->bm, chr);
        
        for (int p = 0; p < Parameters_count(tlk->bm->rates); p++) {
            Parameters_set_value(tlk->bm->rates, p, searcher->default_rate);
        }
		
        // optimize
		double lnl = optimize_singletreelikelihood(tlk);
        individual->fitness = -searcher->IC(searcher, lnl, SingleTreeLikelihood_df_count(tlk));
        
        //Tree_heights_to_vector(tlk->tree, data->heights);
        //BranchModel_rates_to_vector(tlk->bm, data->rates);
        
#pragma omp critical
		{
			if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
				ga->feedback_count--;
				fprintf(stderr, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
			}
			
            // We are comparing AICs not LnLs but we take the negative
			if( individual->fitness > ga->maxfitness ){
				Tree_heights_to_vector(tlk->tree, searcher->best_heights);
				BranchModel_rates_to_vector(tlk->bm, searcher->best_rates);
				memcpy(searcher->best_indexes, individual->chromosome, individual->size * sizeof(unsigned) );
				ga->maxfitness = individual->fitness;
                searcher->best_lk = lnl;
                searcher->n_rate = local_count;
			}
			pFit->lk = individual->fitness;
			if( searcher->logFile != NULL ){
				fprintf(searcher->logFile, "tree TREE%d [&LnL=%f,%s=%f,local=%d] = [&R] ", searcher->tree_count++, lnl, INFORMATION_CRITERION[searcher->ic], -individual->fitness, local_count);
                
                Node **nodes = Tree_nodes(tlk->tree);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buff);
                        StringBuffer_append_format(buff, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buff->c);
                        
                        if ( tlk->bm->indicators[i] ) {
                            Node_set_annotation(nodes[i], "local", "1");
                        }
                    }
                }
                
                Tree_print_nexus_with_annotation(searcher->logFile, tlk->tree);
				//print_tree_nexus_treelikelihood(searcher->logFile, tlk);
				fprintf(searcher->logFile, "\n");
			}
		}
	}
	
	free_StringBuffer( buff );
    free(chr);
	
	return individual->fitness;
}

void _localclock_termination_callback( GA *ga ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
	fprintf( stderr, "Number of local clocks: %d\n", searcher->n_rate );
    fprintf( stderr, "LnL: %f %s: %f\n", searcher->best_lk, INFORMATION_CRITERION[searcher->ic], -ga->maxfitness );
}

