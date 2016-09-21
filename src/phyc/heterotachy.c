/*
 *  heterotachy.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/17/11.
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

#include "heterotachy.h"

#include <math.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif

#include "tree.h"
#include "treeio.h"
#include "branchmodel.h"
#include "matrix.h"

#include "discreteclock.h"
#include "localclock.h"

#define DEBUG_HETEROTACHY 1

const char *HETEROTACHY_STRING[4] = {
	"Local clock exhaustive",
	"Local clock greedy",
	"Local clock GA",
	"Discrete clock GA"
};



#pragma mark -
#pragma mark ClockSearh


ClockSearch * new_ClockSearch( SingleTreeLikelihood *tlk, unsigned nThreads ){
    ClockSearch *search = (ClockSearch*)malloc( sizeof(ClockSearch) );
    assert(search);
    
    search->pDerivedObj = NULL;
    
	search->nThreads = 1;
#if defined (PTHREAD_ENABLED)
    search->nThreads = nThreads;
#elif defined (_OPENMP)
	search->nThreads = nThreads;
	omp_set_num_threads(nThreads);
#endif
	
   	search->pool = new_PooledTreeLikelihood(tlk, search->nThreads, true, false );
    
    search->n_rate = Parameters_count(tlk->bm->rates);
	search->starting_n_rate = Parameters_count(tlk->bm->rates);
	search->max_n_rate = (int)(Tree_node_count(tlk->tree)/2);
    
	search->alpha = 0.05;
	
	search->best_lk = tlk->lk;
	search->best_indexes = uivector(Tree_node_count(tlk->tree));
	search->best_heights = dvector(Tree_node_count(tlk->tree));
	search->best_rates   = dvector(Parameters_count(tlk->bm->rates)); // length(best_rates)-1 == length(best_indexes) for local clocks
	
	search->algorithm = 0;
	
	search->optimize =  NULL;
    search->free     = free_ClockSearch;
	search->init     = ClockSearch_init;
	
	search->starting_df = SingleTreeLikelihood_df_count(tlk);
	search->starting_lk = tlk->lk;
	
	// for resetting the parameters
	search->default_rate = 0.0001;
	
	search->verbosity = 1;
	
	search->mode = 0; // 0 all represented, 1 can be imcomplete assignment
    
    
    
    // THE VARIABLES BELOW SHOULD BE PRIVATE
	// to save parameters from a round to another
	search->current_indexes = uivector(Tree_node_count(tlk->tree));
	search->current_heights = dvector(Tree_node_count(tlk->tree));
	search->current_rates   = dvector(Parameters_count(tlk->bm->rates));
    
    search->df = SingleTreeLikelihood_df_count(tlk); //degrees of freedom for AIC, AICc...
    
	search->initialized = false;
    
    search->logFileName = NULL;
	search->logFileMode = "w";
    search->logFile = NULL;
	search->tree_count = 1; // when we restart a GA we could set this to last tree id+1
    
    search->ic = INFORMATION_CRITERION_AICc;
    search->IC = ClockSearch_IC;
    //search->ic_sample_size = search->pool->tlks[0]->sp->count;
    search->ic_sample_size = search->pool->tlks[0]->sp->nsites;
    
    return search;
}


void free_ClockSearch( ClockSearch *search ){
	search->pool->free(search->pool);
	if( search->logFileName != NULL ) free(search->logFileName);
    
    free(search->current_indexes);
    free(search->current_heights);
    free(search->current_rates);
    
    free(search->best_indexes);
    free(search->best_heights);
    free(search->best_rates);
    
	free(search);
}

void ClockSearch_init( ClockSearch *search ){
	
    if ( search->init == NULL ) {
        error("The clock search alogrithm function not initialzed\n");
    }
    
	if ( search->logFileName != NULL) {
		//fprintf(stderr, "heterotachy_init: file %s\n", search->logFileName);
		search->logFile = fopen(search->logFileName, search->logFileMode);
		assert(search->logFile);

		// we only write the header if it is on write mode.
		// If we are in append mode then the header will be printed by the caller
		if ( search->logFileMode[0] == 'w') {
            Tree_print_nexus_header_figtree(search->logFile, search->pool->tlks[0]->tree);
		}
	}
	
	search->best_lk = search->starting_lk;
	
	search->n_rate = search->starting_n_rate;
    
    search->df = search->starting_df;
    
	search->initialized = true;
}



void ClockSearch_set_max_n_rate( ClockSearch *search, unsigned max ){
    search->max_n_rate = max;
}

void ClockSearch_set_verbosity( ClockSearch *search, int verbosity ){
    search->verbosity = verbosity;
}

void ClockSearch_set_starting_lnl( ClockSearch *search, double lnl ){
    search->starting_lk = lnl;
}

void ClockSearch_set_starting_df( ClockSearch *search, int df ){
    search->starting_df = df;
}

void ClockSearch_set_alpha( ClockSearch *search, double alpha ){
    search->alpha = alpha;
}

void ClockSearch_set_logfile_name( ClockSearch *search, const char *name ){
    search->logFileName = String_clone(name);
}

void ClockSearch_set_logfile_mode( ClockSearch *search, char *mode ){
    search->logFileMode = mode;
}

void * ClockSearch_derived_object( ClockSearch *search ){
    return search->pDerivedObj;
}

void ClockSearch_set_mode( ClockSearch *search, int mode ){
    search->mode = mode;
}

void ClockSearch_set_ic( ClockSearch *search, information_criterion ic ){
    search->ic = ic;
}

double ClockSearch_IC( ClockSearch * search, double lnl, int k ){
    double ic = NAN;
    switch ( search->ic ) {
        case INFORMATION_CRITERION_AIC:{
            ic = AIC(lnl, k);
            break;
        }
        case INFORMATION_CRITERION_BIC:{
            ic = BIC(lnl, k, search->ic_sample_size);
            break;
        }
        case INFORMATION_CRITERION_AICc:{
            ic = AICc(lnl, k, search->ic_sample_size);
            break;
        }
        default:
            break;
    }
    return ic;
}

void ClockSearch_set_ic_sample_size( ClockSearch *search, int sample_size ){
    search->ic_sample_size = sample_size;
}

#pragma mark -
#pragma mark GAClockSearh



ClockSearch * new_GAClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize ){
	switch ( tlk->bm->name ) {
		case CLOCK_LOCAL:
			return new_GALocalClockSearch(tlk, nThreads, popSize);
			break;
		case CLOCK_DISCRETE:
			return new_GADiscreteClockSearch(tlk, nThreads, popSize);
			break;
		case CLOCK_RELAXED:
			error("GA for relaxed clocks is not yet implemented\n");
			break;
		case CLOCK_STRICT:
			error("that's dumb\n");
			break;
		default:
			assert(0);
	}
	return NULL;
}



void free_GAClockSearch( ClockSearch *search ){
	GAClockSearch *gasearch = (GAClockSearch*)search->pDerivedObj;
	
    // it should be freed before
    if (gasearch->pop_individuals != NULL ) {
        for ( int i = 0; i < gasearch->popSize; i++ ) {
            free(gasearch->pop_individuals[i]);
        }
        free(gasearch->pop_individuals);
        
        for ( int i = 0; i < gasearch->popSize; i++ ) {
            free(gasearch->pop_heights[i]);
        }
        free(gasearch->pop_heights);
        for ( int i = 0; i < gasearch->popSize; i++ ) {
            free(gasearch->pop_rates[i]);
        }
        free(gasearch->pop_rates);
    }
    
    free(gasearch); gasearch = NULL;
	free_ClockSearch(search);
}






#pragma mark -
#pragma mark GA stuff


Individual * clone_Individual_GA( Individual *indiv ){
    Individual *clone = new_Individual_unsigned(indiv->size);
    clone->data = (IndividualData *)malloc(sizeof(IndividualData));
    assert(clone->data);
    
    IndividualData *clone_data = (IndividualData*)clone->data;
    IndividualData *indiv_data = (IndividualData*)indiv->data;
    
    clone_data->heights = clone_dvector(indiv_data->heights, indiv->size);
    clone_data->rates   = clone_dvector(indiv_data->rates, indiv_data->n_rate);
    clone_data->n_rate  = indiv_data->n_rate;
    
    memcpy(clone->chromosome, indiv->chromosome, indiv->size*sizeof(unsigned));
    clone->fitness = indiv->fitness;
    clone->id = indiv->id;
    
    return clone;
}

void free_Individual_GA( Individual *indiv ){
    IndividualData *data = (IndividualData*)indiv->data;
	
	free(data->heights);
	free(data->rates);
    free(data);
    free_Individual(indiv);
}

ModelTally * new_ModelTally( const double value ){
	ModelTally *tally = (ModelTally*)malloc(sizeof(ModelTally));
    assert(tally);
	tally->lk = value;
	tally->count = 1;
	return tally;
}

void print_ModelTally_Hashtable( FILE *pf, Hashtable *hash, const int nClass ){
	Hashtable_init_iterator(hash);
	HashEntry *entry = NULL;
	while ( (entry = Hashtable_next(hash) ) != NULL ) {
		ModelTally *value = (ModelTally*)HashEntry_value(entry);
		char *key = (char*)HashEntry_key(entry);
		fprintf(pf, "%s %f %d\n", key, value->lk, value->count);
	}
}
