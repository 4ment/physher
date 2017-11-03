/*
 *  discreteclock.c
 *  PhyC
 *
 *  Created by Mathieu Fourment on 2/18/11.
 *  Copyright (C) 2010 Mathieu Fourment. All rights reserved.
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


#include "discreteclock.h"

#include <math.h>
#include <assert.h>
#include <unistd.h> // for sleep

#ifdef _OPENMP
#include <omp.h>
#endif

//#define PTHREAD_ENABLED 1
#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif

#include "treelikelihood.h"
#include "pooledtreelikelihood.h"
#include "utils.h"
#include "ga.h"
#include "optimize.h"
#include "matrix.h"
#include "treeio.h"
#include "modelselection.h"
#include "random.h"
#include "descriptivestats.h"
#include "classification.h"

#include "heterotachy.h"

#pragma mark Function declaration


#if defined (PTHREAD_ENABLED)
static void * _discreteclock_fitness_threads( void *threadpool );
#endif
    
static void discreteclock_mutate_chc( GA *ga, Individual *individual );
void  _DiscreteClock_mate( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );
void  _DiscreteClock_mate_clade( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );
static double discreteclock_fitness( GA *ga, Individual *individual );

static double _discreteclock_ga_optimize( ClockSearch *search  );

static void _discreteclock_ga_init( GA *ga, ClockSearch *search );


static double _discreteclock_fitness_order( GA *ga, Individual *individual );
static void _discreteclock_mutate_order( GA *ga, Individual *individual );
static void _discreteclock_mate_order( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );

static void _discreteclock_termination_callback( GA *ga );


#pragma mark -
// MARK: ClockSearch

ClockSearch * new_GADiscreteClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize ){
    ClockSearch *search = new_ClockSearch(tlk, nThreads);
	assert(search);
	
	GAClockSearch *gaSearch = (GAClockSearch*)malloc( sizeof(GAClockSearch) );
	assert(gaSearch);
	
	search->pDerivedObj = gaSearch;
    gaSearch->pBaseObj = search;

	search->optimize = _discreteclock_ga_optimize;
	search->free     = free_GAClockSearch;
	
    search->algorithm = HETEROTACHY_DISCRETE_GA;
	
	// GA stuff
	gaSearch->initGA = _discreteclock_ga_init;
	gaSearch->popSize = popSize;
	gaSearch->ngen = 500;
	gaSearch->max_no_improvement = 50;
	
	gaSearch->tree_count = 0;
	gaSearch->verbosity = 1;
	gaSearch->feedback_sampling = 1;
    
    gaSearch->pop_individuals = NULL;
    gaSearch->pop_heights     = NULL;
    gaSearch->pop_rates       = NULL;
	
	return search;
}

static void _sort_increasing_rate( double *rates, int *map, int size ){
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( rates[i] > rates[i+1] ) {
				done = false;
                dswap(&rates[i], &rates[i+1]);
                swap_int(&map[i], &map[i+1]);
			}
		}
		size--;
	}
}

static double _discreteclock_ga_optimize( ClockSearch *search  ){

	if ( !search->initialized ) {
		search->init(search);
	}
	
	Tree *tree = search->pool->tlks[0]->tree;
	BranchModel *bm = search->pool->tlks[0]->bm;
	int nNodes = Tree_node_count(tree);
    
    BranchModel_rates_to_vector(bm, search->current_rates);
    
    if ( search->mode == 1 ) {
        
        int *map = ivector(Parameters_count(bm->rates));
        for ( int i = 0; i < Parameters_count(bm->rates); i++ ) {
            map[i] = i;
        }
        _sort_increasing_rate(search->current_rates, map, Parameters_count(bm->rates));
        for ( int i = 0; i < nNodes; i++ ) {
            bm->map->values[i] = map[ bm->map->values[i] ];
        }
        
        Parameters_sort_from_ivector(bm->rates, map);
        
        free(map);
    }
    
    memcpy(search->current_indexes, bm->map->values, nNodes * sizeof(unsigned));
    Tree_heights_to_vector(search->pool->tlks[0]->tree, search->current_heights);

    memcpy(search->best_indexes, bm->map->values, nNodes * sizeof(unsigned)); // should set indexes here instead of the init_population
    Tree_heights_to_vector(search->pool->tlks[0]->tree, search->best_heights);
    BranchModel_rates_to_vector(bm, search->best_rates);
    
    
	search->default_rate = dmean(search->best_rates, Parameters_count(bm->rates) );
	
	GAClockSearch *gasearch = (GAClockSearch*)search->pDerivedObj;

	
	// logs all the configurations with the likelihood value and the number of times the GA generated it
	// ex: 0:1:2:0:3: -2555.2 3
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
		
	double previous_lk = search->starting_lk;
    int previous_df = search->starting_df;
	double previous_IC = search->IC(search, previous_lk, search->starting_df);
	double IC = previous_IC;
    
    int previous_n_rate = SingleTreeLikelihood_df_count(search->pool->tlks[0])-previous_df;
	
	double precision = search->pool->tlks[0]->opt.precision;
	for ( int j = 0; j < search->pool->count; j++ ) {
		//search->pool->tlks[j]->opt.precision = 0.01;
        if ( search->mode == 1 ) {
            search->pool->tlks[j]->opt.rates.optimize = false;
        }
        
        //search->pool->tlks[j]->opt.verbosity =  2;
	}
	
//    for( int i = 0; i < BranchModel_n_rate(bm); i++ ){
//        if( Parameters_value(bm->rates, i) <= BRANCHMODEL_RATE_MIN*1.01 ){
//            Parameters_set_fixed(bm->rates, true, i);
//        }
//    }

	for ( ; search->n_rate <= search->max_n_rate; search->n_rate++ ) {
		//ga->mutation_rate = 1.75/(30*sqrt(i));
		fprintf(stderr, "Start GA with %d discrete rate classes [%d]\n", search->n_rate, search->mode);
        
		GA *ga = new_GA( GA_CHC, GA_CHROMOSOME_UNSIGNED, gasearch->popSize, nNodes );

		if ( search->mode == 0 ) {
			ga->fitness     = discreteclock_fitness;
			ga->mate        = _DiscreteClock_mate;
			ga->mutate      = discreteclock_mutate_chc;
			
		}
		else {
			ga->fitness     = _discreteclock_fitness_order;
			ga->mate        = _discreteclock_mate_order;
			ga->mutate      = _discreteclock_mutate_order;
		}

		ga->termination = ga_default_termination;
		ga->diversity   = ga_diversity_fitness;
		ga->termination_callback = _discreteclock_termination_callback;
        
        ga->free_individual = free_Individual_GA;
        ga->clone_individual = clone_Individual_GA;
		
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
        ga->thread_worker = _discreteclock_fitness_threads;
#endif
		
		gasearch->initGA(ga, search); // init the Ks but don't compute the fitness
        
        search->df = SingleTreeLikelihood_df_count(search->pool->tlks[0]);
		
		ga_init(ga);
		
		if ( search->starting_n_rate == search->n_rate) {
			fprintf(stderr, "Mating probability:   %f\n", ga->mate_prob);
			fprintf(stderr, "Mutation probability: %f\n", ga->mutation_rate);
			fprintf(stderr, "Mutation threshold:   %f\n", ga->mutation_threshold);
		}
		
		ga_evolve(ga);
		
		search->best_lk = ga->maxfitness;
		
        // when we restart the GA from a log file with n rates and there was no n-1 rates in the file
        // when we restart from a (single) tree file
        // when we start a GA from a local clock
        // when we start a run from scratch without running a strict clock (that should not be possible)
		if ( search->df == search->starting_df && search->n_rate == search->starting_n_rate ){
            
        }
        else {
            
            IC = search->IC(search, search->best_lk, search->df);
			fprintf(stderr, "\n");
			fprintf(stderr, "IC_%d = %f LnL = %f df = %d\n", previous_n_rate, previous_IC, previous_lk, previous_df);
			fprintf(stderr, "IC_%d = %f LnL = %f df = %d", search->n_rate, IC, search->best_lk, search->df );
            
            if( search->pool->tlks[0]->approx != TREELIKELIHOOD_APPROXIMATION_NONE ){
                
                // reset BranchModel
                DiscreteClock_set_classes(bm, search->best_indexes);
                BranchModel_vector_to_rates(bm, search->best_rates);
                // reset tree heights
                Tree_vector_to_heights(search->best_heights, tree);
                Tree_constraint_heights(tree);
                SingleTreeLikelihood_update_all_nodes(search->pool->tlks[0]);
                
                fprintf(stderr, " (%f)\n", search->pool->tlks[0]->calculate(search->pool->tlks[0]) );
                SingleTreeLikelihood_update_all_nodes(search->pool->tlks[0]);
            }
            else {
                fprintf(stderr, "\n");
                
            }
			fprintf(stderr, "==============================================\n\n");
			
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

		previous_df = search->df;
		previous_lk = search->best_lk;
		previous_IC = IC;
        previous_n_rate = search->n_rate;

        
		// copy best_heights, best_indexes, best_rates
		memcpy(search->current_indexes, search->best_indexes, nNodes * sizeof(unsigned) );
		memcpy(search->current_heights, search->best_heights, nNodes * sizeof(double) );
		
		search->best_rates    = realloc( search->best_rates, (search->n_rate+1)*sizeof(double));
		search->current_rates = realloc( search->current_rates, (search->n_rate+1)*sizeof(double));
		
		search->best_rates[search->n_rate] = search->default_rate;
		memcpy(search->current_rates, search->best_rates, (search->n_rate+1) * sizeof(double) );
				
		if ( logHash ) {
			fprintf(hashLogFile, "# %d\n", search->n_rate);
			print_ModelTally_Hashtable(hashLogFile, ga->lookup , search->n_rate);
			fprintf(hashLogFile, "\n");
		}
		
		if( search->logFile != NULL ){
			fprintf(search->logFile, "\n");// seperate the trees with different number of classes
			fprintf(search->logFile, "[nclasses=%d, LnL=%f]\n\n", search->n_rate, search->best_lk);
		}
		
		free_GA(ga);
		
		for ( int j = 0; j < search->pool->count; j++ ) {
			search->pool->tlks[j]->opt.precision = precision;
		}
	}    
    
    // we save the previous model only if it did not reach the maximum number of clock
    // in thi
    if( search->n_rate != search->max_n_rate ){
        // save the best heights and indexes
        memcpy(search->best_indexes, search->current_indexes, nNodes * sizeof(unsigned) );
        memcpy(search->best_heights, search->current_heights, nNodes * sizeof(double) );
        
        search->n_rate--;
        search->best_rates = realloc( search->best_rates, search->n_rate * sizeof(double));
        memcpy(search->best_rates, search->current_rates, search->n_rate * sizeof(double) );
        search->best_lk = previous_lk;
    }
    
    // reset BranchModel
    DiscreteClock_set_number_of_rate_classes(bm, search->n_rate);
    DiscreteClock_set_classes(bm, search->best_indexes);
    BranchModel_vector_to_rates(bm, search->best_rates);
    // reset tree heights
    Tree_vector_to_heights(search->best_heights, tree);
    Tree_constraint_heights(tree);
    SingleTreeLikelihood_update_all_nodes(search->pool->tlks[0]);
	
	if( search->logFile != NULL ) {
        fprintf(search->logFile, "\nEnd;\n");	
		fclose(search->logFile);
	}

	if ( logHash ) {
		fclose(hashLogFile);
	}
		   
	return search->best_lk;
}

static void _init_population( GA *ga, ClockSearch *search );


static void _discreteclock_ga_init( GA *ga, ClockSearch *search ){
	    
    GAClockSearch *gaSearch = (GAClockSearch*)search->pDerivedObj;

    IndividualData *data = NULL;
    
    for ( int i = 0; i < gaSearch->popSize; i++ ) {
        ga->individuals[i]->data = (IndividualData *)malloc(sizeof(IndividualData));
        assert(ga->individuals[i]->data);
        data = (IndividualData*)ga->individuals[i]->data;
        data->heights = dvector(ga->chromosome_size);
        data->rates   = dvector(search->n_rate);
        data->n_rate  = search->n_rate;
    }
    
//    if ( search->n_rate == search->starting_n_rate ) {
    
        if ( gaSearch->pop_individuals == NULL ) {

            _init_population(ga, search);
            
            for ( int i = 0; i < gaSearch->popSize; i++ ) {
                data = (IndividualData*)ga->individuals[i]->data;
                memcpy(data->heights, search->best_heights, ga->chromosome_size * sizeof(double) );
            }
            
            if ( search->mode == 0 ) {
                data = (IndividualData*)ga->individuals[0]->data;
                memcpy(data->rates, search->best_rates, search->n_rate * sizeof(double) );
            
                // should find a better default rate, especially when we have several rates from a local clock
                // for now it is the mean rate sum r_i/n
                for ( int i = 1; i < gaSearch->popSize; i++ ) {
                    data = (IndividualData*)ga->individuals[i]->data;
                    for (int p = 0; p < search->n_rate; p++) {
                        data->rates[p] = search->default_rate;
                    }
                }
            }
            else {
                for ( int i = 0; i < gaSearch->popSize; i++ ) {
                    data = (IndividualData*)ga->individuals[i]->data;
                    memcpy(data->rates, search->best_rates, search->n_rate * sizeof(double) );
                }
            }
            
        }
        else {
            for ( int i = 0; i < ga->pop_size; i++ ) {
                
                int index = Node_id( Tree_root(search->pool->tlks[0]->tree));
                if ( !check_complete_assignment_ignore_one( gaSearch->pop_individuals[i], search->n_rate, ga->individuals[i]->size, index) ){// don't include the root
                    fprintf(stderr, "The initial individual is not ok\n");
                }
                else {
                    memcpy( ga->individuals[i]->chromosome, gaSearch->pop_individuals[i], ga->chromosome_size * sizeof(unsigned) );
                }
                
                data = (IndividualData*)ga->individuals[i]->data;
                memcpy(data->heights, gaSearch->pop_heights[i], ga->chromosome_size * sizeof(double) );
                memcpy(data->rates, gaSearch->pop_rates[i], search->n_rate * sizeof(double) );
            }
            free(gaSearch->pop_individuals);
            free(gaSearch->pop_heights);
            free(gaSearch->pop_rates);
        }
        
//    }
//    else {
//        
//        for ( int j = 0; j < search->pool->count; j++ ) {
//            DiscreteClock_set_number_of_rate_classes( search->pool->tlks[j]->bm, search->n_rate );
//        }
//        
//        int *map = ivector(search->n_rate);
//        
//        // use the previous population and add a class rate to each individual
//        for ( int i = 0; i < gaSearch->popSize; i++ ) {
//            data = (IndividualData*)ga->individuals[i]->data;
//            data->rates  = realloc(data->rates, search->n_rate * sizeof(double));
//            data->n_rate = search->n_rate;
//            
//            memset(map, 0, search->n_rate*sizeof(int));
//            
//            // find the first class represented twice and assign the new class at this position
//			// that would break down if K size == nClass
//            int j = 0;
//			for ( ; j < ga->chromosome_size; j++ ) {
//				if ( map[ ga->individuals[i]->chromosome[j] ] == 1 ) {
//                    data->rates[ search->n_rate-1 ]   = data->rates[ ga->individuals[i]->chromosome[j] ];
//					ga->individuals[i]->chromosome[j] = search->n_rate-1;
//					break;
//				}
//				map[ ga->individuals[i]->chromosome[j] ]++;
//			}
//            
//            uivector_canonical(ga->individuals[i]->chromosome, ga->chromosome_size-1, map, search->n_rate);
//            
//            dvector_sort_from_ivector(data->rates, map, search->n_rate);            
//        }
//        free(map);
//    }
}


void _init_population( GA *ga, ClockSearch *search ){
    Hashtable *hash = new_Hashtable_string(ga->pop_size);
	StringBuffer *buff = new_StringBuffer( ga->chromosome_size *3 );
    bool use_jenks = false;
    
    int root_id = Node_id(Tree_root(search->pool->tlks[0]->tree));
    
    // The first chromosome contains the map of the branchmodel passed to the constructor
    unsigned int *chromosome = ga->individuals[0]->chromosome;
    
    if ( search->n_rate != search->starting_n_rate ) {
        if( use_jenks ){
            //TODO: should be indexed by id not postorder_idx
            fprintf(stderr, "Replace postorder_idx with id %s %d\n",__FILE__,__LINE__);
            exit(1);
            BranchModel *bm = search->pool->tlks[0]->bm;
            Tree_heights_to_vector(bm->tree, search->best_heights);
            BranchModel_rates_to_vector(bm, search->best_rates);
            // do not count the root
            int n = ga->chromosome_size-1;
            int *perm = ivector(n);
            double *d = dvector(n);
            Node **nodes = Tree_get_nodes(bm->tree, POSTORDER);
            for ( int i = 0; i < n; i++ ) {
                d[i] = Node_distance(nodes[i]) - (Node_time_elapsed(nodes[i]) * bm->get(bm,nodes[i]));
                perm[i] = i;
            }
            dvector_sort_track(d, n, perm);
            int *jenks_cats = classification_Jenks_breaks(d, n, search->n_rate);
            
            for ( int i = 0; i < n; i++ ) {                
                for ( int j = 0; j < search->n_rate; j++ ) {
                    if( perm[i] <= jenks_cats[j] ){
                        chromosome[i] = j;
                        break;
                    }
                }
            }
            free(d);
            free(perm);
            free(jenks_cats);
        }

        // Allocate parameters
        for ( int j = 0; j < search->pool->count; j++ ) {
            DiscreteClock_set_number_of_rate_classes( search->pool->tlks[j]->bm, search->n_rate );
        }
    }
	
	if ( search->mode > 0) {
		memcpy( chromosome, search->current_indexes, ga->chromosome_size * sizeof(unsigned) );
	}
	else {
        memcpy( chromosome, search->best_indexes, ga->chromosome_size * sizeof(unsigned) );
        
        
        if( search->n_rate != search->starting_n_rate ){
            
            
            if( !use_jenks ){
                int *map = ivector(search->n_rate);
                // find the first class represented twice (except if one of the nodes is the root) and assign the new class at this position
                // that would break down if K size == nClass
                int j = 0;
                for ( ; j < ga->chromosome_size; j++ ) {
                    if ( map[ chromosome[j] ] == 1 ) {
                        chromosome[j] = search->n_rate-1;
                        break;
                    }
                    if( j != root_id ) map[ chromosome[j] ]++;
                }
                
                if( root_id == 0 ) chromosome[0] = chromosome[1];
                else chromosome[root_id] = chromosome[root_id-1];
                
                uivector_canonical( chromosome, ga->chromosome_size, map, search->n_rate);
                free(map);
            }
        }
        
        if ( !check_complete_assignment( chromosome, search->n_rate, ga->chromosome_size) ){
            fprintf(stderr, "The initial individual is not ok\n");
        }
	}
	
	StringBuffer_empty(buff);
	for ( int j = 0; j < ga->chromosome_size; j++) {
        StringBuffer_append_format(buff, "%u:",chromosome[j]);
	}
	Hashtable_add(hash, String_clone(buff->c), new_Double(1));
	
	// We create an ordered population I0,I1,...,IN-1 where I0 is the best individual set up above and each individual starting
	// from I1 will be a copy of the previous individual with some randomization in the location of the classes

    int *v =  ivector(search->n_rate);
    for ( int i = 1; i < ga->pop_size; i++ ) {
        chromosome = ga->individuals[i]->chromosome;
        
        int counter = ga->pop_size;
        while ( counter-- ){
            memcpy( ga->individuals[i]->chromosome, ga->individuals[i-1]->chromosome, ga->chromosome_size * sizeof(unsigned) );
            // swap some classes. No need to worry about full allocation of the classes
            for ( int j = 0; j < ga->chromosome_size; j++) {
                if( j == root_id ) continue;
                
                double rnum = random_double();
                if ( rnum < 0.2 ) {
                    int rpos = random_int( ga->chromosome_size-1 );
                    swap_uint(&chromosome[j], &chromosome[rpos]);
                }
            }
            
            if( root_id == 0 ) chromosome[0] = chromosome[1];
            else chromosome[root_id] = chromosome[root_id-1];
            
            if ( search->mode == 0 ) {					
                uivector_canonical(chromosome, ga->chromosome_size, v, search->n_rate);
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

// We create a population for random individuals starting from 0 to N-1
// We avoid identical individuals unless it takes too long.
void init_random_population( GA *ga, ClockSearch *search ){
    
    Hashtable *hash = new_Hashtable_string(ga->pop_size);
	StringBuffer *buff = new_StringBuffer( ga->chromosome_size *2 );
    int *v =  ivector(search->n_rate);
    
    unsigned int *chromosome = NULL;
    
    for ( int i = 0; i < ga->pop_size; i++ ) {
        chromosome = ga->individuals[i]->chromosome;
        
        int counter = ga->pop_size;
        while ( counter-- ){
            ga_random_individual_unsigned( ga->individuals[i], search->n_rate ); // no need to check it since this function returns an indiv with at least 1 assignment of each class
            
            if ( search->mode == 0 ) {
                uivector_canonical(chromosome, ga->chromosome_size-1, v, search->n_rate);
            }
            
            StringBuffer_empty(buff);
            for ( int j = 0; j < ga->chromosome_size-1; j++) {
                StringBuffer_append_format(buff, "%u:", chromosome[j]);
            }
            
            if( !Hashtable_exists( hash, buff->c) ){
                break;
            }
        }
        Hashtable_add(hash, String_clone(buff->c), new_Double(1));
    }
    
	free_StringBuffer(buff);
	free_Hashtable(hash);
    free(v);
}

// Avoid the situtation where the root's children have their own rates.
// At least one of the two children must share its rate with another node
bool check_roots_kids_assignment( BranchModel *bm ){
    Node *left  = Tree_root(bm->tree)->left;
    Node *right = Tree_root(bm->tree)->right;
    int class1 = bm->map->values[ Node_id(left)  ];
    int class2 = bm->map->values[ Node_id(right) ];
    
    for ( int i = 0; i < Tree_node_count(bm->tree); i++ ) {
        if( Node_id(Tree_root(bm->tree)) == i ) continue;
        
        //if ( (i != left->postorder_idx && bm->map[i] == class1) ^ (i != right->postorder_idx && bm->map[i] == class2) ) {
        if ( (i != Node_id(left) && i != Node_id(right)) && (bm->map->values[i] == class1 || bm->map->values[i] == class2) ) {
            return true;
        }
    }
    return false;
}

// Avoid the situtation where the root's children have the same rate
bool check_roots_kids_assignment2( Tree *tree, unsigned *k ){
    Node *left  = Tree_root(tree)->left;
    Node *right = Tree_root(tree)->right;
    return ( k[ Node_id(left) ] != k[ Node_id(right) ] );
}


static void _sort_decreasing_lnls( int *map, double *lnls, int size ){
	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( lnls[i] < lnls[i+1] ) {
				done = false;
                dswap(&lnls[i], &lnls[i+1]);
                swap_int(&map[i], &map[i+1]);
			}
		}
		size--;
	}
}


// BranchModel should have 2 rates already
void DiscreteClock_greedy( SingleTreeLikelihood *tlk ){
    
    tlk->opt.verbosity = 0;
    Node **nodes = Tree_nodes(tlk->tree);
    double *lnls = dvector(Tree_node_count(tlk->tree)-1);
    int *map = ivector(Tree_node_count(tlk->tree)-1);
    SingleTreeLikelihood_update_all_nodes(tlk);
    
    double rate = Parameters_value(tlk->bm->rates, 0);
    
    int root_id = Node_id(Tree_root(tlk->tree));
    
    lnls[ root_id ] = -INFINITY;
    map[ root_id ] = root_id;
    
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        if( root_id == i  ) continue;
        
        memset(tlk->bm->map->values, 0, sizeof(unsigned)*Tree_node_count(tlk->tree));
		tlk->bm->map->set_value(tlk->bm->map, i, 1);
        Parameters_set_value(tlk->bm->rates, 0, rate);
        Parameters_set_value(tlk->bm->rates, 1, rate);
        
        SingleTreeLikelihood_update_all_nodes(tlk);
        tlk->bm->need_update = true;
        
        lnls[i] = optimize_singletreelikelihood(tlk);
        map[i] = i;
        
        printf("LnL %f %s %f %f\n", lnls[i], Node_name(nodes[i]), Parameters_value(tlk->bm->rates, 0), Parameters_value(tlk->bm->rates, 1));
    }
    
    _sort_decreasing_lnls(map, lnls, Tree_node_count(tlk->tree)-1);
    
    printf("\n-----------------------------------------\n");
    
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        if( root_id == map[i]  ) continue;
        printf("LnL %f %s (%d)\n", lnls[i], Node_name(nodes[map[i]]), Node_tip_count(nodes[map[i]]));
    }
    
    printf("\n-----------------------------------------\n");
    
    
    
    int i = 1;
    for ( ; i < Tree_node_count(tlk->tree)-1; i++ ) {
        
        DiscreteClock_set_number_of_rate_classes(tlk->bm, i+2);
        
        memset(tlk->bm->map->values, 0, sizeof(unsigned)*Tree_node_count(tlk->tree)-1);
        for ( int j = 0 ; j <= i; j++ ) {
			tlk->bm->map->set_value(tlk->bm->map, map[j], j+1);
        }
        
        SingleTreeLikelihood_update_all_nodes(tlk);
        
        double lnl = optimize_singletreelikelihood(tlk);
        printf("LnL %f [%d]\n", lnl, i);
    }
    
    printf("finished\n");
    free(lnls);
    free(map);
}

static int _discreteclock_greedy_search2( SingleTreeLikelihood *tlk, const int current, int *current_indexes, double *best_lk ) {
    
    int n_rate = Parameters_count(tlk->bm->rates);
    int tot = Tree_node_count(tlk->tree)-n_rate+1;
    
    int best_pos = 0;
    tlk->opt.verbosity = 0;
    
	for ( int i = 0; i < tot; i++ ) {
        
        memset(tlk->bm->map->values, 0, sizeof(unsigned)*Tree_node_count(tlk->tree)-1);
        
        for ( int j = 0; j < current; j++ ) {
			tlk->bm->map->set_value(tlk->bm->map, current_indexes[j], j+1);
        }
        
        int idx = 0;
        for ( ; idx < Tree_node_count(tlk->tree)-1; idx++ ) {
            int k = 0;
            for ( ; k <= current; k++ ) {
                if( current_indexes[k] == idx ){
                    break;
                }
            }
            //printf("idx %d current %d k %d\n", idx, current, k);
            if( k == current+1 && idx > current_indexes[current] ){
                break;
            }
		}
		tlk->bm->map->set_value(tlk->bm->map, idx, current+1);
		current_indexes[current] = idx;
        //print_ivector(current_indexes, current+1);
        
        SingleTreeLikelihood_update_all_nodes(tlk);
        
        double lk = optimize_singletreelikelihood(tlk);
        
        //printf("LnL %f %s %f %f\n", lk, Node_name(nodes[i]), Parameters_value(tlk->bm->rates, 0), Parameters_value(tlk->bm->rates, 1));
        
        
        if ( lk > *best_lk ){
            best_pos = idx;
            *best_lk = lk;
            //Tree_heights_to_vector(tlk->tree, search->best_heights);
            //BranchModel_rates_to_vector(tlk->bm, search->best_rates);
            //memcpy(search->best_indexes, indexes, (search->n_rate-1) * sizeof(unsigned) );
        }
	}
	return best_pos;
}

void DiscreteClock_greedy2( SingleTreeLikelihood *tlk ){
    int *indexes = ivector(Tree_node_count(tlk->tree));
    double lnl = -INFINITY;
    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
    
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        //    int i = 1;
        //    indexes[0] =3;
        if( i != 0 ) DiscreteClock_set_number_of_rate_classes(tlk->bm, i+2);
        indexes[i] = -1;
        indexes[i] = _discreteclock_greedy_search2(tlk, i, indexes, &lnl);
        
        printf("\n-----------------------------------------\n");
        
        printf("%d LnL %f %s\n", (i+2),lnl, Node_name(nodes[indexes[i]]));
    }
    
    free(indexes);
}


static int _discreteclock_greedy_search3( SingleTreeLikelihood *tlk, const int current, int *current_indexes, double *best_lk ) {
    
    int n_rate = Parameters_count(tlk->bm->rates);
    int tot = Tree_node_count(tlk->tree)-n_rate+1;
  
    int best_pos = 0;
    tlk->opt.verbosity = 0;
    
	for ( int i = 0; i < tot; i++ ) {
        
        memset(tlk->bm->map->values, 0, sizeof(unsigned)*Tree_node_count(tlk->tree)-1);
        
        for ( int j = 0; j < current; j++ ) {
			tlk->bm->map->set_value(tlk->bm->map, current_indexes[j], j+1);
        }
        
        int idx = 0;
        for ( ; idx < Tree_node_count(tlk->tree)-1; idx++ ) {
            int k = 0;
            for ( ; k <= current; k++ ) {
                if( current_indexes[k] == idx ){
                    break;
                }
            }
            //printf("idx %d current %d k %d\n", idx, current, k);
            if( k == current+1 && idx > current_indexes[current] ){
                break;
            }
        }
        
        for ( int j = 1; j <= current+1; j++ ) {
			tlk->bm->map->set_value(tlk->bm->map, idx, j);
            current_indexes[current] = idx;
            //print_ivector(current_indexes, current+1);
            
            SingleTreeLikelihood_update_all_nodes(tlk);
            
            double lk = optimize_singletreelikelihood(tlk);
            
            //printf("LnL %f %s %f %f\n", lk, Node_name(nodes[i]), Parameters_value(tlk->bm->rates, 0), Parameters_value(tlk->bm->rates, 1));
            
            
            if ( lk > *best_lk ){
                best_pos = idx;
                *best_lk = lk;
                //Tree_heights_to_vector(tlk->tree, search->best_heights);
                //BranchModel_rates_to_vector(tlk->bm, search->best_rates);
                //memcpy(search->best_indexes, indexes, (search->n_rate-1) * sizeof(unsigned) );
            }
        }
	}
	return best_pos;
}

void DiscreteClock_greedy3( SingleTreeLikelihood *tlk ){
    int *indexes = ivector(Tree_node_count(tlk->tree));
    double lnl = -INFINITY;
    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
    
    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
        //    int i = 1;
        //    indexes[0] =3;
        if( i != 0 ) DiscreteClock_set_number_of_rate_classes(tlk->bm, i+2);
        indexes[i] = -1;
        indexes[i] = _discreteclock_greedy_search3(tlk, i, indexes, &lnl);
        
        printf("\n-----------------------------------------\n");
        
        printf("%d LnL %f %s\n", (i+2),lnl, Node_name(nodes[indexes[i]]));
    }
    
    free(indexes);
}

#pragma mark -
#pragma mark Implementation of GA interface

void  _DiscreteClock_mate( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
	
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
    IndividualData *best_data  = (IndividualData*)best->data;
    
	int classCount = best_data->n_rate;
    
    //    memcpy(new_data->heights, best_data->heights, new1->size * sizeof(double));
    //    memcpy(new_data->rates, best_data->rates, classCount * sizeof(double));
    
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
    while ( !check_complete_assignment(newchromosome, classCount, new1->size) );
    
    // if individuals are not compatible for mating we just mutate the best one
    if( count_failure == 20 ){
        memcpy(newchromosome, bestchromosome, best->size * sizeof(unsigned int));
        discreteclock_mutate_chc(ga, best);
    }
    else {
        int *v =  ivector(classCount);
        uivector_canonical(new1->chromosome, new1->size, v, classCount);
        //dvector_sort_from_ivector(new_data->rates, v, classCount);
        free(v);
    }
}

static void _change_clade( Node *node, unsigned *parent, unsigned *progeny){
    if(node == NULL)return;
    _change_clade(Node_left(node),  parent, progeny);
    _change_clade(Node_right(node), parent, progeny);
    progeny[Node_id(node)] = parent[Node_id(node)];
}

void  _DiscreteClock_mate_clade( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){

    if( random_double() > 0.1 ){
        _DiscreteClock_mate(ga, individual1, individual2, new1);
        return;
    }
    
	ClockSearch *searcher = (ClockSearch*)ga->data;
    int count_failure = 0;
    int nNodes = individual1->size;
    int root_id = Node_id(Tree_root(searcher->pool->tlks[0]->tree));
    IndividualData *best_data  = (IndividualData*)individual1->data;
    int classCount = best_data->n_rate;
    
    unsigned *best  = individual2->chromosome;
    unsigned *worst = individual1->chromosome;
    unsigned *chromosome = new1->chromosome;
    
    if ( individual1->fitness > individual2->fitness ) {
        best  = individual1->chromosome;
        worst = individual2->chromosome;
    }
    
    while ( count_failure < 5 ) {
        int index = root_id;
        memcpy(chromosome, best, nNodes * sizeof(unsigned));
        
        while ( index == root_id) {
            index = random_int(nNodes-1);
        }
        _change_clade( Tree_node(searcher->pool->tlks[0]->tree, index), worst, chromosome);
        
        if ( root_id == 0 ) {
            chromosome[0] = chromosome[1];
        }
        else {
            chromosome[root_id] = chromosome[root_id-1];
        }
        
        if( check_complete_assignment(chromosome, classCount, nNodes ) ) break;
        count_failure++;
    }
    
    if( count_failure == 5 ){
        _DiscreteClock_mate(ga, individual1, individual2, new1);
    }
    
    int *v =  ivector(classCount);
    uivector_canonical(chromosome, nNodes, v, classCount);
    free(v);
}


void discreteclock_mutate_chc( GA *ga, Individual *individual ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
	//GAClockSearch *gaSearch = (GAClockSearch*)searcher->pDerivedObj;
//    IndividualData *data = (IndividualData*)individual->data;

    //memcpy(gaSearch->pop_individuals[individual->id], individual->chromosome, individual->size *sizeof(unsigned) );
    
    unsigned int *chromosome = individual->chromosome;
    unsigned *backup = uivector(individual->size);
	memcpy(backup, individual->chromosome, individual->size *sizeof(unsigned) );
    int nRates = Parameters_count(searcher->pool->tlks[0]->bm->rates);
    int *v =  ivector(nRates);
    int root_id = Node_id(Tree_root(searcher->pool->tlks[0]->tree));
    
    if( memcmp( chromosome, ga->best_chromosome, sizeof(unsigned)*individual->size) == 0 ){
        ClockSearch *searcher = (ClockSearch*)ga->data;
        int left = Node_id(Node_left(Tree_root(searcher->pool->tlks[0]->tree)));
        int right = Node_id(Node_right(Tree_root(searcher->pool->tlks[0]->tree)));
        
        if( chromosome[left] != chromosome[right] ){
            chromosome[left] = chromosome[right];
            if ( root_id == 0 ) {
                chromosome[0] = chromosome[1];
            }
            else {
                chromosome[root_id] = chromosome[root_id-1];
            }
            
            if ( check_complete_assignment( chromosome, nRates, individual->size) ){
                uivector_canonical(chromosome, individual->size, v, nRates);
                free(v);
                free(backup);
                return;
            }
            memcpy( chromosome, backup, individual->size *sizeof(unsigned) );
            
            chromosome[right] = chromosome[left];
            if ( root_id == 0 ) {
                chromosome[0] = chromosome[1];
            }
            else {
                chromosome[root_id] = chromosome[root_id-1];
            }
            if ( check_complete_assignment( chromosome, nRates, individual->size) ){
                uivector_canonical(chromosome, individual->size, v, nRates);
                free(v);
                free(backup);
                return;
            }
        }
    }

    memcpy( chromosome, backup, individual->size *sizeof(unsigned) );
    
	unsigned int bit = 0;
	while ( 1 ){
		for (int i = 0; i < individual->size; i++) {
            // we don't need to mutate the root
            if( i == root_id ) continue;
            
			double rnum = random_double();
			if( rnum < ga->mutation_rate ){
				bit = chromosome[i];
				while ( bit == chromosome[i] ) {
					bit = random_int( nRates-1 );
				}
				chromosome[i] = bit;
			}
		}
        if ( root_id == 0 ) {
            chromosome[0] = chromosome[1];
        }
        else {
            chromosome[root_id] = chromosome[root_id-1];
        }
        
		if ( check_complete_assignment( individual->chromosome, nRates, individual->size) ){
            break;
        }
        memcpy( chromosome, backup, individual->size *sizeof(unsigned) );
	}
    
    
	uivector_canonical(chromosome, individual->size, v, nRates);
	//dvector_sort_from_ivector(data->rates, v, nRates);
	free(v);
    free(backup);
}


double discreteclock_fitness( GA *ga, Individual *individual ){
    
    unsigned int *chromosome = individual->chromosome;
    
	StringBuffer *buff = new_StringBuffer( individual->size*3 );
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
        
        IndividualData *data = (IndividualData*)individual->data;
        
        Tree_vector_to_heights( data->heights, tlk->tree );
        
        BranchModel_vector_to_rates(tlk->bm, data->rates);
        
        Tree_constraint_heights( tlk->tree );
        SingleTreeLikelihood_update_all_nodes( tlk );
        
        DiscreteClock_set_classes( tlk->bm, individual->chromosome );
		
		individual->fitness = optimize_singletreelikelihood(tlk);
        
        Tree_heights_to_vector(tlk->tree, data->heights);
        BranchModel_rates_to_vector(tlk->bm, data->rates);
        
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
				double ic = searcher->IC(searcher, pFit->lk, searcher->df);
				fprintf(searcher->logFile, "tree TREE%d [&LnL=%f,%s=%f] = [&R] ", searcher->tree_count++, individual->fitness, INFORMATION_CRITERION[searcher->ic], ic);
                
                Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buff);
                        StringBuffer_append_format(buff, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buff->c);
                        
                        StringBuffer_empty(buff);
                        StringBuffer_append_format(buff, "%d", tlk->bm->map->values[ Node_id(nodes[i]) ]);
                        Node_set_annotation(nodes[i], "class", buff->c);
                    }
                }
                Tree_print_nexus_with_annotation(searcher->logFile, tlk->tree);
				//print_tree_nexus_treelikelihood(searcher->logFile, tlk);
				fprintf(searcher->logFile, "\n");
			}
		}
	}
	
	free_StringBuffer( buff );
	
	return individual->fitness;
}

#if defined (PTHREAD_ENABLED)
void * _discreteclock_fitness_threads( void *threadpool ){
    threadpool_ga_t *pool = (threadpool_ga_t *)threadpool;
    GA *ga = pool->ga;
    
    int index_treelikelihood = 0;
    pthread_mutex_lock(&(pool->lock));
    index_treelikelihood = pool->index_treelikelihood++;
    pthread_mutex_unlock(&(pool->lock));
    
    StringBuffer *buff = new_StringBuffer( ga->chromosome_size*3 );
    
    ClockSearch *searcher = (ClockSearch*)ga->data;
    
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
        
        ModelTally *pFit = NULL;
        bool found = Hashtable_exists(ga->lookup, buff->c);
		
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
        
        
        pool->count++;
        pthread_mutex_unlock(&(pool->lock));
        
        if ( !found ) {
            
            IndividualData *data = (IndividualData*)individual->data;
            
            Tree_vector_to_heights( data->heights, tlk->tree );
            
            BranchModel_vector_to_rates(tlk->bm, data->rates);
            
            Tree_constraint_heights( tlk->tree );
            SingleTreeLikelihood_update_all_nodes( tlk );
            
            DiscreteClock_set_classes( tlk->bm, chromosome );
            
            individual->fitness = optimize_singletreelikelihood(tlk);
            
            Tree_heights_to_vector(tlk->tree, data->heights);
            BranchModel_rates_to_vector(tlk->bm, data->rates);
            
            pthread_mutex_lock(&(pool->lock));
            {
                if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
                    ga->feedback_count--;
                    fprintf(stderr, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
                }
                
                if( individual->fitness > ga->maxfitness ){
                    Tree_heights_to_vector(tlk->tree, searcher->best_heights);
                    BranchModel_rates_to_vector(tlk->bm, searcher->best_rates);
                    memcpy(searcher->best_indexes, chromosome, individual->size * sizeof(unsigned) );
                    ga->maxfitness = individual->fitness;
                }
                pFit->lk = individual->fitness;
                if( searcher->logFile != NULL ){
                    double ic = searcher->IC(searcher, pFit->lk, searcher->df);
                    fprintf(searcher->logFile, "tree TREE%d [&LnL=%f,%s=%f] = [&R] ", searcher->tree_count++, individual->fitness, INFORMATION_CRITERION[searcher->ic], ic);
                    
                    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
                    for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                        Node_empty_annotation(nodes[i]);
                        if( !Node_isroot(nodes[i]) ){
                            StringBuffer_empty(buff);
                            StringBuffer_append_format(buff, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                            Node_set_annotation(nodes[i], "rate", buff->c);
                            
                            StringBuffer_empty(buff);
                            StringBuffer_append_format(buff, "%d", tlk->bm->map->values[ Node_id(nodes[i]) ]);
                            Node_set_annotation(nodes[i], "class", buff->c);
                        }
                    }
                    Tree_print_nexus_with_annotation(searcher->logFile, tlk->tree);
                    fprintf(searcher->logFile, "\n");
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



#pragma mark -
// MARK: without optimization of rates

// Rates and their order should never change, only indexes.
// All the individuals have the same rates

// Same as discreteclock_mutate_chc but does not transform the string into canonical form
void _discreteclock_mutate_order( GA *ga, Individual *individual ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
    
    int root_id = Node_id(Tree_root(searcher->pool->tlks[0]->tree));
    unsigned int *chromosome = individual->chromosome;
	
    if ( /* DISABLES CODE */ (1) == 0 ) {
        unsigned *backup = uivector(individual->size);
        memcpy(backup, individual->chromosome, individual->size *sizeof(unsigned) );
        
        int nRates = Parameters_count(searcher->pool->tlks[0]->bm->rates);
        
        unsigned int bit = 0;
        while ( 1 ){
            for (int i = 0; i < individual->size; i++) {
                // we don't need to mutate the root
                if( i == root_id ) continue;
                
                double rnum = random_double();
                if( rnum < ga->mutation_rate ){
                    bit = chromosome[i];
                    while ( bit == chromosome[i] ) {
                        bit = random_int( nRates-1 );
                    }
                    chromosome[i] = bit;
                }
            }
            
            if ( root_id == 0 ) {
                chromosome[0] = chromosome[1];
            }
            else {
                chromosome[root_id] = chromosome[root_id-1];
            }
            
            if ( check_complete_assignment( chromosome, nRates, individual->size) ) break;
            memcpy( chromosome, backup, individual->size *sizeof(unsigned) );
        }
        free(backup);
    }
    else {
        IndividualData *data = (IndividualData*)individual->data;
        
        Tree *tree = searcher->pool->tlks[0]->tree;
        Node **nodes = Tree_nodes(tree);
        
        int *indexes = ivector(Tree_node_count(tree));
        double earliest = 0;
        int count = 0;
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            if ( Node_isleaf(nodes[i]) ) {
                earliest = dmax(earliest, Node_height(nodes[i]));
            }
        }
        
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            if( i == root_id ) continue;
            
            if ( Node_height(nodes[i]) >= earliest || (Node_height(nodes[i]) <= earliest && Node_height(Node_parent(nodes[i])) > earliest) ) {
                indexes[count++] = Node_id(nodes[i]);
                //printf("%s\n", nodes[i]->name);
                if( chromosome[Node_id(nodes[i])] < data->n_rate-2 ){
                    chromosome[Node_id(nodes[i])] += 2;
                }
            }
        }
        
        free(indexes);
        
    }
}

void  _discreteclock_mate_order( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
	
    ClockSearch *searcher = (ClockSearch*)ga->data;
    
    IndividualData *new_data = (IndividualData*)new1->data;
    
    Individual *best  = individual2;
    Individual *worst = individual1;
    if ( individual1->fitness > individual2->fitness ) {
        best  = individual1;
        worst = individual2;
    }
    IndividualData *best_data  = (IndividualData*)best->data;
    
    
    unsigned int *newchromosome   = new1->chromosome;
    unsigned int *bestchromosome  = best->chromosome;
    unsigned int *worstchromosome = worst->chromosome;
    
	int classCount = best_data->n_rate;
    
    memcpy(new_data->heights, best_data->heights, new1->size * sizeof(double));
    
    int root_id = Node_id(Tree_root(searcher->pool->tlks[0]->tree));
    
	do {
		// we need to reset the K otherwise at each iteration it will look more and more like the best chromosome
        memcpy(newchromosome, bestchromosome, best->size * sizeof(unsigned int));
        
		for ( int i = 0; i < best->size; i++) {
            // we don't need to mutate the root
            if( i == root_id ) continue;
            
			double rnum = random_double();
			if ( rnum < ga->mate_prob ) {
				newchromosome[i] = worstchromosome[i];
			}
		}
	}
    while ( !check_complete_assignment(new1->chromosome, classCount, new1->size) );
}

double _discreteclock_fitness_order( GA *ga, Individual *individual ){
    
    unsigned int *chromosome = individual->chromosome;
    
	StringBuffer *buff = new_StringBuffer( individual->size*3 );
	int i = 0;
	for ( i = 0; i < individual->size-1; i++) {
		StringBuffer_append_format(buff, "%u:", chromosome[i]);
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
        
        IndividualData *data = (IndividualData*)individual->data;
        
        Tree_vector_to_heights( data->heights, tlk->tree );
        
        BranchModel_vector_to_rates(tlk->bm, data->rates);
        
        Tree_constraint_heights( tlk->tree );
        SingleTreeLikelihood_update_all_nodes( tlk );
        
        DiscreteClock_set_classes( tlk->bm, individual->chromosome );
		
		individual->fitness = optimize_singletreelikelihood(tlk);
        
        Tree_heights_to_vector(tlk->tree, data->heights);
        //BranchModel_rates_to_vector(tlk->bm, data->rates);
        
        #pragma omp critical
		{
            //printf("%d %f  %s\n", individual->id, individual->fitness, buff->c);
			if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ){
				ga->feedback_count--;
				fprintf(stderr, "%d/%d\r", (ga->pop_size-ga->feedback_count), ga->pop_size);
			}
			
			if( individual->fitness > ga->maxfitness ){
				Tree_heights_to_vector(tlk->tree, searcher->best_heights);
				//BranchModel_rates_to_vector(tlk->bm, searcher->best_rates);
				memcpy(searcher->best_indexes, individual->chromosome, individual->size * sizeof(unsigned) );
				ga->maxfitness = individual->fitness;
			}
			pFit->lk = individual->fitness;
			if( searcher->logFile != NULL ){
				double ic = searcher->IC(searcher, pFit->lk, searcher->df);
				fprintf(searcher->logFile, "tree TREE%d [&LnL=%f,IC=%f] = [&R] ", searcher->tree_count++, individual->fitness, ic);
                
                Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
                for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
                    Node_empty_annotation(nodes[i]);
                    if( !Node_isroot(nodes[i]) ){
                        StringBuffer_empty(buff);
                        StringBuffer_append_format(buff, "%e", tlk->bm->get(tlk->bm,nodes[i]));
                        Node_set_annotation(nodes[i], "rate", buff->c);
                        
                        StringBuffer_empty(buff);
                        StringBuffer_append_format(buff, "%d", tlk->bm->map->values[ Node_id(nodes[i]) ]);
                        Node_set_annotation(nodes[i], "class", buff->c);
                    }
                }
                Tree_print_nexus_with_annotation(searcher->logFile, tlk->tree);
                
				//print_tree_nexus_treelikelihood(searcher->logFile, tlk);
				fprintf(searcher->logFile, "\n");
			}
		}
	}
	
	free_StringBuffer( buff );
	
	return individual->fitness;
}



void _discreteclock_termination_callback( GA *ga ){
	ClockSearch *searcher = (ClockSearch*)ga->data;
	fprintf( stderr, "Root height: %f\n", searcher->best_heights[ Node_id(Tree_root(searcher->pool->tlks[0]->tree)) ] );
}


