/*
 *  qsearch.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 20/08/2015.
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

#include "qsearch.h"

#include <assert.h>

#include "matrix.h"
#include "optimize.h"
#include "ga.h"
#include "random.h"
#include "mathconstant.h"

#include "classification.h"

typedef struct QSearchModelTally{
    double lk;
    unsigned count;
} QSearchModelTally;


QSearchModelTally * new_QSearchModelTally( const double value );

Individual * clone_Indiv_GA( Individual *indiv );

void free_Indiv_GA( Individual *indiv );


#pragma mark -
#pragma mark QSearh



void _free_QSearch( QSearch *search );

void _QSearch_init( QSearch *search );

static double _IC( QSearch *search, double lnl, int k );

QSearch * new_QSearch( SingleTreeLikelihood *tlk, unsigned nThreads ){
    QSearch *search = (QSearch*)malloc( sizeof(QSearch) );
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
    
    search->rateCount = Parameters_count(tlk->sm->m->rates);
    search->startingRateCount = search->rateCount;
    //TODO: change
    search->max_n_rate = tlk->sm->nstate*(tlk->sm->nstate-1)/2;
    
    search->alpha = 0.05;
    
    search->indexCount = tlk->sm->nstate*(tlk->sm->nstate-1);
    if(tlk->sm->m->reversible){
        search->indexCount /= 2;
    }
    search->best_lk = tlk->lk;
    search->best_indexes = uivector(search->indexCount);
    search->best_rates   = dvector(search->rateCount);
    
    search->optimize =  NULL;
    search->free     = _free_QSearch;
    search->init     = _QSearch_init;
    
    search->starting_df = SingleTreeLikelihood_df_count(tlk);
    search->starting_lk = tlk->lk;
    
    // for resetting the parameters
    search->default_rate = 1;
    
    search->verbosity = 1;
    
    // THE VARIABLES BELOW SHOULD BE PRIVATE
    // to save parameters from a round to another
    search->current_indexes = uivector(search->indexCount);
    search->current_rates   = dvector(search->rateCount);
    
    search->df = SingleTreeLikelihood_df_count(tlk); //degrees of freedom for AIC, AICc...
    
    search->initialized = false;
    
    search->logFileName = NULL;
    search->logFileMode = "w";
    search->logFile = NULL;
    
    search->ic = INFORMATION_CRITERION_AIC;
    search->IC = _IC;
    //search->ic_sample_size = search->pool->tlks[0]->sp->count;
    search->ic_sample_size = search->pool->tlks[0]->sp->nsites;
    
    search->mode = 0;
    
    return search;
}

void _free_QSearch( QSearch *search ){
    search->pool->free(search->pool);
    if( search->logFileName != NULL ) free(search->logFileName);
    
    free(search->current_indexes);
    free(search->current_rates);
    
    free(search->best_indexes);
    free(search->best_rates);
    
    free(search);
}

void _QSearch_init( QSearch *search ){
    
    if ( search->init == NULL ) {
        error("The clock search alogrithm function not initialzed\n");
    }
    
    if ( search->logFileName != NULL) {
        search->logFile = fopen(search->logFileName, search->logFileMode);
    }
    
    search->best_lk = search->starting_lk;
    
    search->rateCount = search->startingRateCount;
    
    search->df = search->starting_df;
    
    search->initialized = true;
}

double _IC( QSearch *search, double lnl, int k ){
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

#pragma mark -
#pragma mark GAQSearh

static double _ga_optimize( QSearch *search );
static void _ga_init( GA *ga, QSearch *search );


static void _free_GAQSearch( QSearch *search ){
    GAQSearch *gasearch = (GAQSearch*)search->pDerivedObj;
    
    // it should be freed before
    if (gasearch->pop_individuals != NULL ) {
        for ( int i = 0; i < gasearch->popSize; i++ ) {
            free(gasearch->pop_individuals[i]);
            free(gasearch->pop_rates[i]);
        }
        free(gasearch->pop_individuals);
        free(gasearch->pop_rates);
    }
    
    free(gasearch);
    _free_QSearch(search);
}


QSearch * new_GAQSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize ){
    QSearch *search = new_QSearch(tlk, nThreads);
    
    GAQSearch *gaSearch = (GAQSearch*)malloc( sizeof(GAQSearch) );
    assert(gaSearch);
    
    search->pDerivedObj = gaSearch;
    gaSearch->pBaseObj = search;
    
    search->optimize = _ga_optimize;
    search->free     = _free_GAQSearch;
    
    
    // GA stuff
    gaSearch->initGA = _ga_init;
    gaSearch->popSize = popSize;
    gaSearch->ngen = 500;
    gaSearch->max_no_improvement = 50;
    
    gaSearch->verbosity = 1;
    gaSearch->feedback_sampling = 1;
    
    gaSearch->pop_individuals = NULL;
    gaSearch->pop_rates       = NULL;
    
    return search;
}

void model_to_chromosome(unsigned *chromosome, const unsigned *model, int dim, bool reversible){
    int index = 0;
    if(reversible){
        for ( int i = 0; i < dim; i++ ) {
            for ( int j = i+1; j < dim; j++ ) {
                chromosome[index++] = model[i*dim+j];
            }
        }
    }
    else {
        for ( int i = 0; i < dim; i++ ) {
            for ( int j = 0; j < dim; j++ ) {
                if (i != j) {
                    chromosome[index++] = model[i*dim+j];
                }
            }
        }
    }
}

void chromosome_to_model(DiscreteParameter* model, const unsigned *chromosome, int dim, bool reversible){
    int index = 0;
    if(reversible){
        for ( int i = 0; i < dim; i++ ) {
            for ( int j = i+1; j < dim; j++ ) {
                //model[j*dim+i] = model[i*dim+j] = chromosome[index++];
				model->set_value(model, i*dim+j, chromosome[index++]);
            }
        }
    }
    else {
        for ( int i = 0; i < dim; i++ ) {
            for ( int j = 0; j < dim; j++ ) {
                if (i != j) {
                    //model[i*dim+j] = chromosome[index++];
					model->set_value(model, i*dim+j, chromosome[index++]);
                }
            }
        }
    }
}

static void _init_population_from_best(GA *ga, QSearch *search){
    Hashtable *hash = new_Hashtable_string(ga->pop_size);
    StringBuffer *buff = new_StringBuffer( ga->chromosome_size *3 );
    int *v =  ivector(search->rateCount);
    
    unsigned int *chromosome = ga->individuals[0]->chromosome;
    
    uivector_canonical(chromosome, ga->chromosome_size, v, search->rateCount);
    
    StringBuffer_empty(buff);
    for ( int j = 0; j < ga->chromosome_size; j++) {
        StringBuffer_append_format(buff, "%u:",chromosome[j]);
    }
    Hashtable_add(hash, String_clone(buff->c), new_Double(1));
    
    
    // We create an ordered population I0,I1,...,IN-1 where I0 is the best individual set up above and each individual starting
    // from I1 will be a copy of the previous individual with some randomization in the location of the classes
    for ( int i = 1; i < ga->pop_size; i++ ) {
        chromosome = ga->individuals[i]->chromosome;
        
        int counter = ga->pop_size;
        while ( counter-- ){
            memcpy( ga->individuals[i]->chromosome, ga->individuals[i-1]->chromosome, ga->chromosome_size * sizeof(unsigned) );
            // swap some classes. No need to worry about full allocation of the classes
            for ( int j = 0; j < ga->chromosome_size; j++) {
                
                double rnum = random_double();
                if ( rnum < 0.2 ) {
                    int rpos = random_int( ga->chromosome_size-1 );
                    swap_uint(&chromosome[j], &chromosome[rpos]);
                }
            }
            if(search->mode == 0){
                uivector_canonical(chromosome, ga->chromosome_size, v, search->rateCount);
            }
            
            StringBuffer_empty(buff);
            for ( int j = 0; j < ga->chromosome_size; j++) {
                StringBuffer_append_format(buff, "%u:",chromosome[j]);
            }
            
            if( !Hashtable_exists( hash, buff->c) ){
                break;
            }
        }
        
        if( !Hashtable_exists( hash, buff->c) ) {
            Hashtable_add(hash, String_clone(buff->c), new_Double(1));
        }
    }
    free(v);
    
    
    free_StringBuffer(buff);
    free_Hashtable(hash);
}

static void _init_population( GA *ga, QSearch *search ){
    Hashtable *hash = new_Hashtable_string(ga->pop_size);
    StringBuffer *buff = new_StringBuffer( ga->chromosome_size *3 );
    
    // The first chromosome contains the map of the branchmodel passed to the constructor
    unsigned int *chromosome = ga->individuals[0]->chromosome;
    
    memcpy( chromosome, search->best_indexes, ga->chromosome_size * sizeof(unsigned) );
    int *v =  ivector(search->rateCount);
    
    if ( search->rateCount != search->startingRateCount ) {
        char name[50] = "q.";
        int index = 0;
        
        // Add an additional parameter
        for ( int j = 0; j < search->pool->count; j++ ) {
			//TODO: check why +4 and not +2
            sprintf(name+4, "%d", (Parameters_count(search->pool->tlks[j]->sm->m->rates)+1) );
            double rate = Parameters_value(search->pool->tlks[j]->sm->m->rates, 0);
            Parameter *p = new_Parameter_with_postfix_and_ownership(name, "relativerate", rate, Parameters_constraint(search->pool->tlks[j]->sm->m->rates, index), false);
            Parameters_move(search->pool->tlks[j]->sm->m->rates, p);
        }
        
        int *map = ivector(search->rateCount);
        // find the first class represented twice (except if one of the nodes is the root) and assign the new class at this position
        // that would break down if K size == nClass
        int j = 0;
        for ( ; j < ga->chromosome_size; j++ ) {
            if ( map[ chromosome[j] ] == 1 ) {
                chromosome[j] = search->rateCount-1;
                break;
            }
            map[ chromosome[j] ]++;
        }
        uivector_canonical( chromosome, ga->chromosome_size, map, search->rateCount);
        free(map);
        
    }
    // The initial chromosome is not necessarily canonical
    else if(search->mode == 0){
        uivector_canonical(chromosome, ga->chromosome_size, v, search->rateCount);
    }
    
    if ( !check_complete_assignment( chromosome, search->rateCount, ga->chromosome_size) ){
        fprintf(stderr, "The initial individual is not ok\n");
    }
    
    
    StringBuffer_empty(buff);
    for ( int j = 0; j < ga->chromosome_size; j++) {
        StringBuffer_append_format(buff, "%u:",chromosome[j]);
    }
    Hashtable_add(hash, String_clone(buff->c), new_Double(1));
    
    if(0){
        // We create an ordered population I0,I1,...,IN-1 where I0 is the best individual set up above and each individual starting
        // from I1 will be a copy of the previous individual with some randomization in the location of the classes
        for ( int i = 1; i < ga->pop_size; i++ ) {
            chromosome = ga->individuals[i]->chromosome;
            
            int counter = ga->pop_size;
            while ( counter-- ){
                memcpy( ga->individuals[i]->chromosome, ga->individuals[i-1]->chromosome, ga->chromosome_size * sizeof(unsigned) );
                // swap some classes. No need to worry about full allocation of the classes
                for ( int j = 0; j < ga->chromosome_size; j++) {
                    
                    double rnum = random_double();
                    if ( rnum < 0.2 ) {
                        int rpos = random_int( ga->chromosome_size-1 );
                        swap_uint(&chromosome[j], &chromosome[rpos]);
                    }
                }
                if(search->mode == 0){
                    uivector_canonical(chromosome, ga->chromosome_size, v, search->rateCount);
                }
                
                StringBuffer_empty(buff);
                for ( int j = 0; j < ga->chromosome_size; j++) {
                    StringBuffer_append_format(buff, "%u:",chromosome[j]);
                }
                
                if( !Hashtable_exists( hash, buff->c) ){
                    break;
                }
            }
            
            if( !Hashtable_exists( hash, buff->c) ) {
                Hashtable_add(hash, String_clone(buff->c), new_Double(1));
            }
        }
    }
    else {
        for ( int i = 1; i < ga->pop_size; i++ ) {
            chromosome = ga->individuals[i]->chromosome;
            
            int counter = ga->pop_size;
            while ( counter-- ){
                memcpy( chromosome, ga->individuals[0]->chromosome, ga->chromosome_size * sizeof(unsigned) );
                // swap some classes. No need to worry about full allocation of the classes
                for ( int j = 0; j < ga->chromosome_size; j++) {
                    double rnum = random_double();
                    if ( rnum < 0.2 ) {
                        int rpos = random_int( ga->chromosome_size-1 );
                        swap_uint(&chromosome[j], &chromosome[rpos]);
                    }
                }
                if(search->mode == 0){
                    uivector_canonical(chromosome, ga->chromosome_size, v, search->rateCount);
                }
                
                StringBuffer_empty(buff);
                for ( int j = 0; j < ga->chromosome_size; j++) {
                    StringBuffer_append_format(buff, "%u:",chromosome[j]);
                }
                
                if( !Hashtable_exists( hash, buff->c) ){
                    break;
                }
            }
            
            if( !Hashtable_exists( hash, buff->c) ) {
                Hashtable_add(hash, String_clone(buff->c), new_Double(1));
            }
        }
    }
    free(v);
    
    
    free_StringBuffer(buff);
    free_Hashtable(hash);
}

void _ga_init( GA *ga, QSearch *search ){
    if( search->mode == 2 ){

    }
    else {
        _init_population(ga, search);
    }
}

static void _mutate_chc( GA *ga, Individual *individual );
static void _mutate2_chc( GA *ga, Individual *individual );
static double _fitness( GA *ga, Individual *individual );
#if defined (PTHREAD_ENABLED)
static void * _fitness_threads( void *threadpool );
#endif
static void  _mate( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );
static void  _mate2( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );


double _ga_optimize( QSearch *search ){
    
    if ( !search->initialized ) {
        search->init(search);
    }
    
    SingleTreeLikelihood *tlk = search->pool->tlks[0];
    SubstitutionModel *m = tlk->sm->m;
    
    model_to_chromosome(search->current_indexes, m->model->values, m->nstate, m->reversible);
    model_to_chromosome(search->best_indexes, m->model->values, m->nstate, m->reversible);
    //memcpy(search->current_indexes, m->model, search->indexCount * sizeof(unsigned));
    //memcpy(search->best_indexes, search->current_indexes, search->indexCount * sizeof(unsigned));
    
    Parameters_store_value(m->rates, search->current_rates);
    Parameters_store_value(m->rates, search->best_rates);
    
    GAQSearch *gasearch = (GAQSearch*)search->pDerivedObj;
    
    double previous_lk = search->starting_lk;
    int previous_df = search->starting_df;
    double previous_IC = search->IC(search, previous_lk, search->starting_df);
    double IC;
    
    int previousRateCount = SingleTreeLikelihood_df_count(search->pool->tlks[0])-previous_df;
    
    int searchMode = search->mode;
    
    for ( ; search->rateCount <= search->max_n_rate; search->rateCount++ ) {
        fprintf(stdout, "\nStart GA with %d rate classes\n", search->rateCount);
        
        GA *ga = new_GA( GA_CHC, GA_CHROMOSOME_UNSIGNED, gasearch->popSize, search->indexCount );
        
        ga->fitness    = _fitness;
        
        if(search->mode == 0){
            ga->mate   = _mate;
            ga->mutate = _mutate_chc;
        }
        else if( search->mode == 1 ) {
            ga->mate   = _mate2;
            ga->mutate = _mutate2_chc;
        }
        else if( search->mode == 2 ){
            ga->mate   = _mate2;
            ga->mutate = _mutate2_chc;
        }
        
        ga->termination = ga_default_termination;
        ga->diversity   = ga_diversity_fitness;
        //ga->termination_callback = _discreteclock_termination_callback;
        
        ga->free_individual = free_Indiv_GA;
        ga->clone_individual = clone_Indiv_GA;
        
        ga->mutation_threshold =  0.0001;
        ga->mutation_rate = 0.15;
        
        ga_set_ngeneration(ga, gasearch->ngen);
        
        ga->data = search;
        
        ga->use_max_no_improvement = true;
        ga->max_no_improvement = gasearch->max_no_improvement;
        
        ga_set_n_threads(ga, search->nThreads);
        
        ga->feedback_computeall = true;
        ga->verbosity = gasearch->verbosity;
        ga->feedback_sampling = gasearch->feedback_sampling;
        
        #if defined (PTHREAD_ENABLED)
            ga->thread_worker = _fitness_threads;
        #endif
        
        gasearch->initGA(ga, search); // init the Ks but don't compute the fitness
        
        search->df = SingleTreeLikelihood_df_count(search->pool->tlks[0]);
        
        ga_init(ga);
        
        if ( search->startingRateCount == search->rateCount) {
            fprintf(stdout, "Mating probability:   %f\n", ga->mate_prob);
            fprintf(stdout, "Mutation probability: %f\n", ga->mutation_rate);
            fprintf(stdout, "Mutation threshold:   %f\n", ga->mutation_threshold);
        }
        
        
        
        if(searchMode == 1){
            
            int maxGen = ga->max_generation;
            int chunk = maxGen/10;
            ga->max_generation = chunk;
            for (int i = 0; i < 10; i++ ) {
                ga_evolve(ga);
                
                if(search->logFile != NULL) fprintf(search->logFile, "\n");
                
                Hashtable_empty(ga->lookup);
                
                search->mode = 0;
                ga->fitness(ga, ga->individuals[0]);
                search->mode = 1;
                
                if(search->logFile != NULL) fprintf(search->logFile, "\n");
                ga->max_generation += chunk;
                
                // We need to update the other treelikelihood to get the new rates
                for ( int k = 1; k < search->nThreads; k++) {
                    for ( int j = 0; j < Parameters_count(search->pool->tlks[0]->sm->m->rates); j++ ) {
                        Parameters_set_value(search->pool->tlks[k]->sm->m->rates, j, Parameters_value(search->pool->tlks[0]->sm->m->rates,j));
                    }
                }
                
                if (maxGen <= ga->generation) {
                    break;
                }
            }
            
            ga->max_generation = maxGen;
            search->mode = 0;
            
            ga->mate        = _mate;
            ga->mutate      = _mutate_chc;
            ga->generation = 0;
            int *v =  ivector(search->rateCount);
            Hashtable_empty(ga->lookup);
            
            for ( int i = 0; i < ga->pop_size; i++ ) {
                unsigned * chromosome = ga->individuals[i]->chromosome;
                uivector_canonical(chromosome, ga->individuals[i]->size, v, search->rateCount);
            }
            free(v);
        }
        else if(searchMode == 2 && search->rateCount == search->startingRateCount){
            
            int rateCount = search->rateCount; // reset!!
            
            
            int n = 50;
            double *values = log_spaced_spaced_vector(0.00001, 500, n);
            
            for ( int i = 0; i < n; i++ ) {
                printf("%f\n",values[i]);
            }
            
            char name[50] = "q.";
            // Add an additional parameter
            for ( int j = 0; j < search->pool->count; j++ ) {
                int i = 0;
                for ( ; i < rateCount; i++ ) {
                    Parameters_set_value(search->pool->tlks[j]->sm->m->rates, i, values[i]);
                }
                
                for ( ; i < n; i++ ) {
                    sprintf(name+4, "%d", (Parameters_count(search->pool->tlks[j]->sm->m->rates)+1) );
                    Parameter *p = new_Parameter_with_postfix_and_ownership(name, "relativerate", values[i], Parameters_constraint(search->pool->tlks[j]->sm->m->rates, 0), false);
                    Parameters_move(search->pool->tlks[j]->sm->m->rates, p);
                }
            }
        
            
            for ( int j = 0; j < ga->pop_size; j++ ) {
                unsigned *chromosome = ga->individuals[j]->chromosome;
                for (int i = 0; i < ga->chromosome_size; i++) {
                    chromosome[i] = random_int( n-1 );
                }
            }
            free(values);
            
            search->rateCount = n;
            
            ga_evolve(ga);
            
            
            int *perm = ivector(search->indexCount);
            double *d = dvector(search->indexCount);
            unsigned * chromosome = ga->individuals[0]->chromosome;
            int count = 0;
            
            if( tlk->sm->m->reversible ){
                for ( int i = 0; i < tlk->sm->nstate; i++ ) {
                    for ( int j = i+1; j < tlk->sm->nstate; j++ ) {
                        d[count] = Parameters_value(search->pool->tlks[0]->sm->m->rates, chromosome[count]);
                        perm[count] = count;
                        count++;
                    }
                }
            }
            
            dvector_sort_track(d, count, perm);
            int *jenks_cats = classification_Jenks_breaks(d, count, rateCount);
            
            for ( int k = 0; k < rateCount; k++ ) {
                printf("%d \n", jenks_cats[k]);
            }
            
            count = 0;
            if( tlk->sm->m->reversible ){
                for ( int i = 0; i < tlk->sm->nstate; i++ ) {
                    for ( int j = i+1; j < tlk->sm->nstate; j++ ) {
                        for ( int k = 0; k < rateCount; k++ ) {
                            if( count <= jenks_cats[k] ){
                                printf("%d) %d %d %d %d value: %f\n",count,perm[count],chromosome[count], jenks_cats[k], k, d[count]);
                                chromosome[ perm[count] ] = k;
                                break;
                            }
                        }
                        count++;
                    }
                }
            }
            
            
            free(jenks_cats);
            free(perm);
            free(d);
            
            
            
            // Remove additional parameters
            for ( int j = 0; j < search->pool->count; j++ ) {
                for ( int i = n-1; i >= rateCount; i-- ) {
                    //Parameters_remove(search->pool->tlks[j]->sm->m->rates, i);
                    Parameters_pop(search->pool->tlks[j]->sm->m->rates);
                }
            }
            search->rateCount = rateCount;
            
            if(search->logFile != NULL) fprintf(search->logFile, "\n");
            
            Hashtable_empty(ga->lookup);
            
            
            // We need to update the other treelikelihood to get the new rates
            for ( int k = 1; k < search->nThreads; k++) {
                for ( int j = 0; j < Parameters_count(search->pool->tlks[0]->sm->m->rates); j++ ) {
                    Parameters_set_value(search->pool->tlks[k]->sm->m->rates, j, Parameters_value(search->pool->tlks[0]->sm->m->rates,j));
                }
            }
                
             
            search->mode = 0;
            
            ga->mate        = _mate;
            ga->mutate      = _mutate_chc;
            ga->generation = 0;
            
            Hashtable_empty(ga->lookup);
            
            _init_population_from_best(ga, search);
            
        }
        
        ga_evolve(ga);
        
        
        search->best_lk = ga->maxfitness;
        
        free_GA(ga);
        
        IC = search->IC(search, search->best_lk, search->df);
        fprintf(stdout, "\n");
        fprintf(stdout, "IC_%d = %f LnL = %f df = %d\n", previousRateCount, previous_IC, previous_lk, previous_df);
        fprintf(stdout, "IC_%d = %f LnL = %f df = %d\n", search->rateCount, IC, search->best_lk, search->df );
        fprintf(stdout, "==============================================\n\n");
        
        if ( IC > previous_IC ) {
            break;
        }
        
        previous_df = search->df;
        previous_lk = search->best_lk;
        previous_IC = IC;
        previousRateCount = search->rateCount;
        
        // copy best_indexes, best_rates
        memcpy(search->current_indexes, search->best_indexes, search->indexCount * sizeof(unsigned) );
        
        search->best_rates    = realloc( search->best_rates, (search->rateCount+1)*sizeof(double));
        search->current_rates = realloc( search->current_rates, (search->rateCount+1)*sizeof(double));
        
        search->best_rates[search->rateCount] = search->default_rate;
        memcpy(search->current_rates, search->best_rates, (search->rateCount+1) * sizeof(double) );
        
        if( search->logFile != NULL ){
            // seperate the models with different number of classes
            fprintf(search->logFile, "\n[#rates=%d LnL=%f]\n\n", search->rateCount, search->best_lk);
        }
    }
    
    // we save the previous model only if it did not reach the maximum number of clock
    if( search->rateCount != search->max_n_rate ){
        // save the best heights and indexes
        memcpy(search->best_indexes, search->current_indexes, search->indexCount * sizeof(unsigned) );
        
        search->rateCount--;
        search->best_rates = realloc( search->best_rates, search->rateCount * sizeof(double));
        memcpy(search->best_rates, search->current_rates, search->rateCount * sizeof(double) );
        search->best_lk = previous_lk;
        Parameters_pop(m->rates);//remove rate
    }
    
    for ( int i = 0; i < Parameters_count(m->rates); i++ ) {
        tlk->sm->m->set_rate(tlk->sm->m, search->best_rates[i], i);
    }
    
    chromosome_to_model(m->model, search->best_indexes, m->nstate, m->reversible);
    
    SingleTreeLikelihood_update_all_nodes( tlk );
    
    if( search->logFile != NULL ) {
        fclose(search->logFile);
    }
    
    return search->best_lk;
}

void _mutate_chc( GA *ga, Individual *individual ){
    unsigned int bit = 0;
    unsigned int *chromosome = individual->chromosome;
    QSearch *searcher = (QSearch*)ga->data;
    SingleTreeLikelihood *tlk = searcher->pool->tlks[0];
    int nRates = Parameters_count(tlk->sm->m->rates);
    
    unsigned *backup = uivector(individual->size);
    int *v =  ivector(nRates);
    memcpy(backup, chromosome, individual->size *sizeof(unsigned) );
    
    while (1) {
        for (int i = 0; i < individual->size; i++) {
            double rnum = random_double();
            if( rnum < ga->mutation_rate ){
                bit = chromosome[i];
                while ( bit == chromosome[i] ) {
                    bit = random_int( nRates-1 );
                }
                chromosome[i] = bit;
            }
        }
        if ( check_complete_assignment( chromosome, nRates, individual->size) ){
            break;
        }
        memcpy( chromosome, backup, individual->size *sizeof(unsigned) );
    }
    uivector_canonical(chromosome, individual->size, v, nRates);
    free(backup);
    free(v);
}

void  _mate( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
    int count_failure = 0;
    Individual *best  = individual2;
    Individual *worst = individual1;
    QSearch *searcher = (QSearch*)ga->data;
    SingleTreeLikelihood *tlk = searcher->pool->tlks[0];
    int nRates = Parameters_count(tlk->sm->m->rates);
    
    if ( individual1->fitness > individual2->fitness ) {
        best  = individual1;
        worst = individual2;
    }
    
    unsigned int *newchromosome = new1->chromosome;
    unsigned int *bestchromosome = best->chromosome;
    unsigned int *worstchromosome = worst->chromosome;
    
    do {
        // we need to reset the K otherwise at each iteration it will look more and more like the best chromosome
        memcpy(newchromosome, bestchromosome, best->size * sizeof(unsigned));
        
        for ( int i = 0; i < best->size; i++) {
            double rnum = random_double();
            if ( rnum < ga->mate_prob ) {
                newchromosome[i] = worstchromosome[i];
            }
        }
        count_failure++;
        if(count_failure == 20) break;
    }
    while ( !check_complete_assignment(newchromosome, nRates, new1->size) );
    
    // if individuals are not compatible for mating we just mutate the best one
    if( count_failure == 20 ){
        memcpy(newchromosome, bestchromosome, best->size * sizeof(unsigned));
        _mutate_chc(ga, best);
    }
    else {
        int *v =  ivector(nRates);
        uivector_canonical(new1->chromosome, new1->size, v, nRates);
        free(v);
    }
}

double _fitness( GA *ga, Individual *individual ){
    
    unsigned int *chromosome = individual->chromosome;
    
    StringBuffer *buff = new_StringBuffer( individual->size*3 );
    int i = 0;
    for ( i = 0; i < individual->size; i++) {
        StringBuffer_append_format(buff, "%u:", chromosome[i]);
    }
    
    QSearchModelTally *pFit = NULL;
    bool found = false;
    
#pragma omp critical
    {
        found = Hashtable_exists(ga->lookup, buff->c);
        
        if ( found ){
            QSearchModelTally *tally = Hashtable_get(ga->lookup, buff->c);
            individual->fitness = tally->lk; // at this point tally->lk can be NAN
            tally->count++;
            
            if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
                ga->feedback_count--;
                fprintf(stdout, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
            }
        }
        else {
            pFit = new_QSearchModelTally(NAN);
            Hashtable_add( ga->lookup, String_clone(buff->c), pFit);
        }
        
    }
    
    if ( !found ) {
        
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        
        QSearch *searcher = (QSearch*)ga->data;
        
        SingleTreeLikelihood *tlk = searcher->pool->tlks[tid];
        
        //IndivData *data = (IndivData*)individual->data;
        
//        for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++ ) {
//            tlk->sm->m->set_rate(tlk->sm->m, data->rates[i], i);
//        }
        chromosome_to_model(tlk->sm->m->model, individual->chromosome, tlk->sm->m->nstate, tlk->sm->m->reversible);
        
        SingleTreeLikelihood_update_all_nodes( tlk );
        
        if(searcher->mode == 0 ) {
            Parameters_set_all_value(tlk->sm->m->rates, 1.0);
            individual->fitness = optimize_singletreelikelihood(tlk);
        }
        else {
            tlk->sm->m->need_update = true;
            individual->fitness = tlk->calculate(tlk);
        }
        
//        Parameters_store_value(tlk->sm->m->rates, data->rates);
        
        if( searcher->logFile != NULL ){
            double ic = searcher->IC(searcher, individual->fitness, searcher->df);
            
            StringBuffer_empty(buff);
            StringBuffer_append_format(buff, "%f", individual->fitness);
            StringBuffer_append_format(buff, ",%f", ic);
            for ( int i = 0; i < tlk->sm->m->nstate * tlk->sm->m->nstate; i++ ) {
                StringBuffer_append_format(buff, ",%d", tlk->sm->m->model[i]);
            }
            
            for ( int i = 0; i < searcher->rateCount; i++ ) {
                StringBuffer_append_format(buff, ",%f", Parameters_value(tlk->sm->m->rates, i));
            }
        }
        
#pragma omp critical
        {
            if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
                ga->feedback_count--;
                fprintf(stdout, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
                fflush(stdout);
            }
            
            if( individual->fitness > ga->maxfitness ){
                memcpy(searcher->best_indexes, individual->chromosome, individual->size * sizeof(unsigned) );
                
                if(searcher->mode < 2) Parameters_store_value(tlk->sm->m->rates, searcher->best_rates);
                ga->maxfitness = individual->fitness;
            }
            pFit->lk = individual->fitness;
            if( searcher->logFile != NULL ){
                fprintf(searcher->logFile, "%s\n", buff->c);
                //fflush(searcher->logFile);
            }
        }
    }
    
    free_StringBuffer( buff );
    
    return individual->fitness;
}

#if defined (PTHREAD_ENABLED)
void * _fitness_threads( void *threadpool ){
    threadpool_ga_t *pool = (threadpool_ga_t *)threadpool;
    GA *ga = pool->ga;
    
    int index_treelikelihood = 0;
    pthread_mutex_lock(&(pool->lock));
    index_treelikelihood = pool->index_treelikelihood++;
    pthread_mutex_unlock(&(pool->lock));
    
    StringBuffer *buff = new_StringBuffer( ga->chromosome_size*3 );
    
    QSearch *searcher = (QSearch*)ga->data;
    
    SingleTreeLikelihood *tlk = searcher->pool->tlks[index_treelikelihood];
    
    Individual *individual = NULL;
    
    
    while( 1 ) {
        
        pthread_mutex_lock(&(pool->lock));
        
        if( pool->count == ga->pop_size ){
            pthread_mutex_unlock(&(pool->lock));
            break;
        }
        //printf("Thread %d individual %d\n",index_treelikelihood, pool->count );
        
        individual = ga->individuals[pool->count]; // count contains the indexes of the next individual to evaluate
        
        unsigned int *chromosome = individual->chromosome;
        
        StringBuffer_empty(buff);
        int i = 0;
        for ( i = 0; i < ga->chromosome_size; i++) {
            StringBuffer_append_format(buff, "%u:", chromosome[i]);
        }
        
        QSearchModelTally *pFit = NULL;
        bool found = Hashtable_exists(ga->lookup, buff->c);
        
        if ( found ){
            QSearchModelTally *tally = Hashtable_get(ga->lookup, buff->c);
            individual->fitness = tally->lk; // at this point tally->lk can be NAN
            tally->count++;
            
            if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
                ga->feedback_count--;
                fprintf(stdout, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
                fflush(stdout);
            }
        }
        else {
            pFit = new_QSearchModelTally(NAN);
            Hashtable_add( ga->lookup, String_clone(buff->c), pFit);
        }
        
        
        pool->count++;
        pthread_mutex_unlock(&(pool->lock));
        
        if ( !found ) {
            
//            IndivData *data = (IndivData*)individual->data;
//            
//            for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++ ) {
//                tlk->sm->m->set_rate(tlk->sm->m, data->rates[i], i);
//            }
            chromosome_to_model(tlk->sm->m->model, individual->chromosome, tlk->sm->m->nstate, tlk->sm->m->reversible);
            
            SingleTreeLikelihood_update_all_nodes( tlk );
            
            
            if(searcher->mode == 0 ) {
                Parameters_set_all_value(tlk->sm->m->rates, 1.0);
                individual->fitness = optimize_singletreelikelihood(tlk);
            }
            else {
                tlk->sm->m->need_update = true;
                individual->fitness = tlk->calculate(tlk);
            }
//            Parameters_store_value(tlk->sm->m->rates, data->rates);
            
            if( searcher->logFile != NULL ){
                double ic = searcher->IC(searcher, pFit->lk, searcher->df);
                
                StringBuffer_empty(buff);
                StringBuffer_append_format(buff, "%f", individual->fitness);
                StringBuffer_append_format(buff, ",%f", ic);
                for ( int i = 0; i < tlk->sm->m->nstate * tlk->sm->m->nstate; i++ ) {
                    StringBuffer_append_format(buff, ",%d", tlk->sm->m->model[i]);
                }
                
                for ( int i = 0; i < searcher->rateCount; i++ ) {
                    StringBuffer_append_format(buff, ",%f", Parameters_value(tlk->sm->m->rates, i));
                }
            }
            
            pthread_mutex_lock(&(pool->lock));
            {
                if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
                    ga->feedback_count--;
                    fprintf(stdout, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
                    fflush(stdout);
                }
                
                if( individual->fitness > ga->maxfitness ){
                    memcpy(searcher->best_indexes, individual->chromosome, individual->size * sizeof(unsigned) );
                    Parameters_store_value(tlk->sm->m->rates, searcher->best_rates);
                    ga->maxfitness = individual->fitness;
                }
                pFit->lk = individual->fitness;
                if( searcher->logFile != NULL ){
                    fprintf(searcher->logFile, "%s\n", buff->c);
                    //fflush(searcher->logFile);
                }
            }
            pthread_mutex_unlock(&(pool->lock));
        }
    }
    free_StringBuffer( buff );
    
    pthread_exit(NULL);
    return NULL;
}
#endif

QSearchModelTally * new_QSearchModelTally( const double value ){
    QSearchModelTally *tally = (QSearchModelTally*)malloc(sizeof(QSearchModelTally));
    assert(tally);
    tally->lk = value;
    tally->count = 1;
    return tally;
}

Individual * clone_Indiv_GA( Individual *indiv ){
    Individual *clone = new_Individual_unsigned(indiv->size);
    
    memcpy(clone->chromosome, indiv->chromosome, indiv->size*sizeof(unsigned));
    clone->fitness = indiv->fitness;
    clone->id = indiv->id;
    
    return clone;
}

void free_Indiv_GA( Individual *indiv ){
    
    free_Individual(indiv);
}


void _mutate2_chc( GA *ga, Individual *individual ){
    unsigned int bit = 0;
    unsigned int *chromosome = individual->chromosome;
    QSearch *searcher = (QSearch*)ga->data;
    SingleTreeLikelihood *tlk = searcher->pool->tlks[0];
    int nRates = Parameters_count(tlk->sm->m->rates);
    
    for (int i = 0; i < individual->size; i++) {
        double rnum = random_double();
        if( rnum < ga->mutation_rate ){
            bit = chromosome[i];
            while ( bit == chromosome[i] ) {
                bit = random_int( nRates-1 );
            }
            chromosome[i] = bit;
        }
    }
    
    
}

void  _mate2( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
    Individual *best  = individual2;
    Individual *worst = individual1;
    
    if ( individual1->fitness > individual2->fitness ) {
        best  = individual1;
        worst = individual2;
    }
    
    unsigned int *newchromosome = new1->chromosome;
    unsigned int *worstchromosome = worst->chromosome;
    
    memcpy(newchromosome,  best->chromosome, best->size * sizeof(unsigned));
    
    for ( int i = 0; i < best->size; i++) {
        double rnum = random_double();
        if ( rnum < ga->mate_prob ) {
            newchromosome[i] = worstchromosome[i];
        }
    }
    
}


#pragma mark -
#pragma mark LASSO


double l1norm( void *data ){
    SingleTreeLikelihood *tlk = ((BrentData*)data)->tlk;
    
    double penalty = 0;
    for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++) {
        if( Parameters_estimate( tlk->sm->m->rates, i ) ){
			penalty += Parameters_value(tlk->sm->m->rates,i);
		}
    }
    return penalty;
}

double penalyzed_loglikelihood( void *data ){
    SingleTreeLikelihood *tlk = ((BrentData*)data)->tlk;
    
    double lambda = *((BrentData*)data)->backup;
    double penalty = l1norm(data);
    //printf("%f %f %f ", tlk->calculate(tlk), lambda, penalty);
    return tlk->calculate(tlk) - lambda * penalty;
}

double _brent_optimize_relative_rate( Parameters *params, double *grad, void *data ){
    BrentData *mydata = (BrentData*)data;
    SingleTreeLikelihood *stlk = mydata->tlk;
    
    stlk->sm->m->set_rate( stlk->sm->m, Parameters_value(params, 0), mydata->index_param );
    SingleTreeLikelihood_update_all_nodes(stlk);
    double lnl = mydata->f(mydata);
    //printf("  %f %f\n",Parameters_value(params, 0),lnl);
    return -lnl;
}

double _cg_optimize_rel_rates( Parameters *params, double *grad, void *data ){
    MultivariateData *mydata = (MultivariateData*)data;
    SingleTreeLikelihood *stlk = mydata->tlk;
    
    for (int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
        stlk->sm->m->set_rate( stlk->sm->m, Parameters_value(params, i), i );
    }
    SingleTreeLikelihood_update_all_nodes(stlk);
    //double lk = stlk->calculate(stlk);
    double lk = mydata->f(mydata);
    fprintf(stderr, "PLnL: %f LnL: %f lnorm: %f\n",lk, stlk->calculate(stlk), -lk+stlk->calculate(stlk));
    
    if ( grad != NULL ) {
        for (int i = 0; i < Parameters_count(stlk->sm->m->rates); i++) {
            double h = SQRT_EPS*(fabs(Parameters_value(params,i)) + 1.0);
            
            double oldx = Parameters_value(params,i);
            
            stlk->sm->m->set_rate( stlk->sm->m, oldx + h, i );
            SingleTreeLikelihood_update_all_nodes(stlk);
            double fxplus = -mydata->f(mydata);
            
            stlk->sm->m->set_rate( stlk->sm->m, oldx - h, i );
            SingleTreeLikelihood_update_all_nodes(stlk);
            double fxminus = -mydata->f(mydata);
            
            stlk->sm->m->set_rate( stlk->sm->m, oldx, i );
            SingleTreeLikelihood_update_all_nodes(stlk);
            
            // Centered first derivative
            grad[i] = (fxplus-fxminus)/(2.0*h);
            
            fprintf(stderr, "cg_optimize_rates %s value=%f fxplus=%f fxminus=%f grad=%f h=%e\n", Parameters_name(stlk->sm->m->rates, i), oldx, fxplus, fxminus, grad[i], h);
        }
        fprintf(stderr, "\n");
    }
    //fprintf(stderr, "LnL %f\n",lk);
    
    return -lk;
}

int df( SingleTreeLikelihood *tlk ){
    int count = 0;
    for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++) {
        if( Parameters_estimate( tlk->sm->m->rates, i ) &&
			Parameters_value(tlk->sm->m->rates,i) <= Constraint_lower(Parameters_at(tlk->sm->m->rates, i)->cnstr)*1.2 ){
            count++;
        }
    }
    return count;
}

void qsearch_lasso(SingleTreeLikelihood *tlk ){
    double lambda = 1;
    
    Optimizer *opt_rel_rate = new_Optimizer( tlk->opt.relative_rates.method );
    opt_set_max_iteration(opt_rel_rate, tlk->opt.relative_rates.max_iteration );
    
    MultivariateData *data = new_MultivariateData( tlk, NULL);
    data->f = penalyzed_loglikelihood;
    data->model = &lambda;
    opt_set_data(opt_rel_rate, data);
    opt_set_objective_function(opt_rel_rate, _cg_optimize_rel_rates);
    opt_set_tolfx(opt_rel_rate, 0.001);
    
    
    Parameters *oneparameter = new_Parameters(1);
    FILE *pFile = fopen("/Users/mathieu/Desktop/pml.log", "w");
    printf("Number of parameters %d\n", Parameters_count(tlk->sm->m->rates));
    
    for ( lambda = 0.0001; lambda < 1000; lambda *= 10 ) {
        double lnl = penalyzed_loglikelihood(data);
        double fret = lnl;
        
        double status = opt_optimize( opt_rel_rate, tlk->sm->m->rates, &fret);
        if( status == OPT_ERROR ) error("OPT.relative_rates No SUCCESS!!!!!!!!!!!!\n");
        fret = -fret;
        if( status != OPT_SUCCESS ){
            tlk->sm->m->need_update = true;
            SingleTreeLikelihood_update_all_nodes(tlk);
        }
        printf("-------------------------------\n");
        printf("Lambda %f PML %f LnL %f AIC %f df %d\n", lambda, -lnl, tlk->calculate(tlk), AIC(tlk->calculate(tlk), df(tlk)), df(tlk));
        
        fprintf(pFile,"Lambda %f PML %f LnL %f AIC %f df %d\n", lambda, -lnl, tlk->calculate(tlk), AIC(tlk->calculate(tlk), df(tlk)), df(tlk));
        
        for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++) {
            if( Parameters_estimate( tlk->sm->m->rates, i ) ){
				fprintf(pFile," %e", Parameters_value(tlk->sm->m->rates,i));
			}
			
        }
        fprintf(pFile, "\n");
        fflush(pFile);
    }
    fclose(pFile);
    free_Optimizer(opt_rel_rate);
    free(data);
    free_Parameters(oneparameter);
}

void qsearch_lasso2(SingleTreeLikelihood *tlk ){
    double lambda = 1;
    BrentData *data = (BrentData*)malloc(sizeof(BrentData));
    assert(data);
    data->tlk = tlk;
    data->index_param = 0;
    data->f = penalyzed_loglikelihood;
    data->backup = &lambda;
    
    Optimizer *opt_rates = new_Optimizer( OPT_BRENT );
    opt_set_max_iteration(opt_rates, 1000);
    opt_set_data(opt_rates, data );
    opt_set_objective_function(opt_rates, _brent_optimize_relative_rate );
    
    Parameters *oneparameter = new_Parameters(1);
    FILE *pFile = fopen("/Users/mathieu/Desktop/pml.log", "w");
    printf("Number of parameters %d\n", Parameters_count(tlk->sm->m->rates));
    double status;
    for ( lambda = 0.001; lambda < 1000; lambda *= 10 ) {
        double lnl = penalyzed_loglikelihood(data);
        double fret = lnl;
        
        for ( int rounds = 0; rounds < 1000; rounds++ ) {
            
            for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++) {
                if( Parameters_estimate( tlk->sm->m->rates, i ) ){
					data->index_param = i;
					Parameters_set(oneparameter, 0, Parameters_at(tlk->sm->m->rates, i));
					
					status = opt_maximize( opt_rates, oneparameter, &fret);
					if( status == OPT_ERROR ) error("OPT.RATES No SUCCESS!!!!!!!!!!!!\n");
					fret = -fret;
					//printf("%d LnL: %f\n", i, fret);
				}
            }
            
            if ( fret - lnl < 0.01 ) {
                lnl = fret;
                break;
            }
            lnl = fret;
        }
        printf("-------------------------------\n");
        printf("Lambda %f PML %f LnL %f AIC %f df %d\n", lambda, lnl, tlk->calculate(tlk), AIC(tlk->calculate(tlk), df(tlk)), df(tlk));
        
        fprintf(pFile,"Lambda %f PML %f LnL %f AIC %f df %d\n", lambda, lnl, tlk->calculate(tlk), AIC(tlk->calculate(tlk), df(tlk)), df(tlk));
        
        for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++) {
            if( Parameters_estimate( tlk->sm->m->rates, i ) ){
				fprintf(pFile," %e", Parameters_value(tlk->sm->m->rates,i));
			}
			
        }
        fprintf(pFile, "\n");
    }
    fclose(pFile);
    free_Optimizer(opt_rates);
    free(data);
    free_Parameters(oneparameter);
}

