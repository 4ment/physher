/*
 *  ga.c
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

#include "ga.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <signal.h>

#include "utils.h"
#include "matrix.h"
#include "random.h"
#include "mstring.h"
#include "utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif

#define GA_DEFAULT_MUTATION_PROB 0.15
#define GA_DEFAULT_MATING_PROB 0.5
#define GA_DEFAULT_DIVERSITY_THRESHOLD 0.0001

//#define DEBUG_GA_MUTATE 1
//#define DEBUG_CHC_SELECTION

static void ga_core( GA *ga );
static void ga_chc( GA *ga );

static int sort_compare_individuals( const void *a, const void *b );
static bool ga_are_chromosomes_different_unsigned( const unsigned *k1, const unsigned *k2, const unsigned size );


static void _copy_chromosome_unsigned( void *dst, const void *src, int size );

static void _copy_chromosome_double( void *dst, const void *src, int size );

static void _copy_chromosome_bool( void *dst, const void *src, int size );

void ga_void( GA *ga ){}

#pragma mark -
#pragma mark GA
/**********************************************************************************************************************************
 ********************************************************* GA *********************************************************************
 **********************************************************************************************************************************/

GA * new_GA( ga_selection_algorithm name, ga_chromosome_type chromosome_type, const unsigned int pop_size, const unsigned int chromosome_size ){
	GA *ga = (GA*)malloc( sizeof(GA) );
	assert(ga);
	
	if ( name == GA_SA ) {
		
	}

	ga->mutation_rate = GA_DEFAULT_MUTATION_PROB;
    ga->min_mutation_rate = 0.01;
	ga->max_generation = 100;
	ga->selection = name;
	ga->pop_size = pop_size;
	ga->tol = 0.001;
	ga->use_tol = false;
	ga->use_annealing = false;
	ga->generation = 0;
	ga->cooling_factor = 0.95;
	ga->annealing_phase = 0.1 * ga->max_generation;
    
    ga->nelites = pop_size/2;
    
	ga->chromosome_size = chromosome_size;
    ga->chromosome_type = chromosome_type;
	
	ga->max_no_improvement     = 20;
	ga->n_without_improvement  = 0;
	ga->use_max_no_improvement = false;
	
	ga->maxfitness = -INFINITY;
	ga->best_chromosome = NULL;
	
	ga->individuals = NULL;
	
	ga->mutate  = NULL;
	ga->mate    = NULL;
	ga->fitness = NULL;
	
	ga->termination = ga_default_termination;
	ga->termination_callback = NULL;
	
	ga->log = ga_void;
	
	// dumping is not used for now
	ga->dump     = false;
	ga->output   = NULL;
	ga->filename = NULL;
	ga->dump_frequency = 0;
	
	ga->lookup = new_Hashtable_string(100);
	
	ga->data = NULL;
	
	if ( name == GA_CHC) {
		ga->mate_prob = GA_DEFAULT_MATING_PROB;
	}
	else {
		ga->mate_prob = 0.01;
	}
	
	// CHC specific
	ga->mutation_threshold = GA_DEFAULT_DIVERSITY_THRESHOLD;
	ga->mate_prob = GA_DEFAULT_MATING_PROB;
	
	ga->nThreads = 1;
	
	ga->initialized = false;
	
	ga->start_time   = 0;
	ga->current_time = 0;
	
	ga->eval = 0;
	
	ga->feedback_computeall = false;
	ga->feedback_count = 0;
	
	ga->verbosity = 1;
	ga->feedback_sampling = 1;
	
	ga->interrupted   = false;
	ga->interruptible = true;
    
    ga->free_individual  = free_Individual;

    ga->individuals = (Individual**)calloc( pop_size, sizeof(Individual*) );
    assert(ga->individuals);
    
    switch (chromosome_type) {
        case GA_CHROMOSOME_UNSIGNED:{
            ga->best_chromosome = uivector(chromosome_size);
            ga->copy_chromosome = _copy_chromosome_unsigned;

            ga->new_individual = new_Individual_unsigned;
            ga->clone_individual = clone_Individual_unsigned;
            
            ga->diversity = ga_diversity_atomic_unsigned;
            ga->mate    = ga_default_mate_unsigned;
            break;
        }
        case GA_CHROMOSOME_DOUBLE:{
            ga->best_chromosome = dvector(chromosome_size);
            ga->copy_chromosome = _copy_chromosome_double;
            
            ga->new_individual = new_Individual_double;
            ga->clone_individual = clone_Individual_double;
            
            ga->diversity = ga_diversity_atomic_double;
            ga->mate    = ga_default_mate_double;
            break;
        }
        case GA_CHROMOSOME_BOOL:{
            ga->best_chromosome = bvector(chromosome_size);
            ga->copy_chromosome = _copy_chromosome_bool;
            
            ga->new_individual = new_Individual_bool;
            ga->clone_individual = clone_Individual_bool;
            
            ga->diversity = ga_diversity_atomic_bool;
            ga->mate    = ga_default_mate_bool;
            break;
        }
        default:{
            error("new_GA not implemented yet\n");
            break;
        }
    }
    
    
    for (int i = 0; i < pop_size; i++) {
        ga->individuals[i]     = ga->new_individual(ga->chromosome_size);
        ga->individuals[i]->id = i;
    }
	
	return ga;
}

void free_GA( GA *ga ){
	for (int i = 0; i < ga->pop_size; i++) {
		ga->free_individual(ga->individuals[i]);
	}
	free(ga->individuals);
	ga->individuals = NULL;
	if ( ga->lookup != NULL ) free_Hashtable( ga->lookup);
	if ( ga->filename != NULL ) free(ga->filename);
	if ( ga->output != NULL ) fclose(ga->output);
	if( ga->best_chromosome != NULL )free(ga->best_chromosome);
	free(ga);
	signal(SIGINT, SIG_DFL);// restore the default handler
}

void ga_evolve( GA *ga ){
	if ( !ga->initialized ) {
		ga_init(ga);
	}
	switch ( ga->selection ) {
		case GA_CHC:
			ga_chc( ga );
			break;
		case GA_ROULETTE:
			ga_core(ga);
			break;
		default:
			fprintf(stderr, "No GA algorithm specified\n");
			exit(1);
			break;
	}
}

bool ga_interrupted = false;

void signal_callback_handler( int signum ) {
	printf("Caught signal %d\n",signum);
	if ( ga_interrupted == true ) {
		exit(SIGINT);
	}
	ga_interrupted = true;
}


void ga_init( GA *ga ){
	if( ga->dump && ga->filename != NULL ){
		ga->output = fopen( ga->filename, "w");
		assert( ga->output );
	}

	time(&ga->start_time);
	
	if ( ga->interruptible ) {
		signal(SIGINT, signal_callback_handler);
		ga_interrupted = ga->interrupted = false;
	}
    ga->initialized = true;
    ga_set_n_threads(ga, ga->nThreads);
}

#pragma mark -
#pragma mark Individual

/**********************************************************************************************************************************
 **************************************************** INDIVIDUAL ******************************************************************
 **********************************************************************************************************************************/

Individual * new_Individual_unsigned( const unsigned int chromosome_size ){
    Individual *indiv = (Individual*)malloc( sizeof(Individual) );
    assert(indiv);
	indiv->chromosome = uivector(chromosome_size);
	indiv->size = chromosome_size;
	indiv->fitness = 0;
	indiv->data = NULL;
	return indiv;
}

void free_Individual( Individual *indiv ){
	if( indiv->chromosome != NULL ) free(indiv->chromosome);
	free(indiv);
	indiv = NULL;
}

Individual * clone_Individual_unsigned( Individual *individual ){
	Individual *indiv = new_Individual_unsigned( individual->size );
    unsigned int *chromosome = individual->chromosome;
    unsigned int *newchromosome = indiv->chromosome;
    
	memcpy(newchromosome, chromosome, individual->size * sizeof(unsigned int));
	indiv->fitness = individual->fitness;
	indiv->size = individual->size;
	indiv->id = individual->id;
	indiv->data = NULL;
	return indiv;
}

void _copy_chromosome_unsigned( void *dst, const void *src, int size ){
    memcpy(dst, src, size*sizeof(unsigned) );
}

Individual * new_Individual_bool( const unsigned int chromosome_size ){
	Individual *indiv = (Individual*)malloc( sizeof(Individual) );
	indiv->chromosome = bvector(chromosome_size);
	assert(indiv->chromosome);
	indiv->size = chromosome_size;
	indiv->fitness = 0;
	indiv->data = NULL;
	return indiv;
}

Individual * clone_Individual_bool( Individual *individual ){
	Individual *indiv = new_Individual_double( individual->size );
    bool *chromosome = individual->chromosome;
    bool *newchromosome = indiv->chromosome;
    
	memcpy(newchromosome, chromosome, individual->size * sizeof(bool));
	indiv->fitness = individual->fitness;
	indiv->size = individual->size;
	indiv->id = individual->id;
	indiv->data = NULL;
	return indiv;
}

void _copy_chromosome_bool( void *dst, const void *src, int size ){
    memcpy(dst, src, size*sizeof(bool) );
}

Individual * new_Individual_double( const unsigned int chromosome_size ){
	Individual *indiv = (Individual*)malloc( sizeof(Individual) );
    assert(indiv);
	indiv->chromosome = dvector(chromosome_size);
	indiv->size = chromosome_size;
	indiv->fitness = 0;
	indiv->data = NULL;
	return indiv;
}

Individual * clone_Individual_double( Individual *individual ){
	Individual *indiv = new_Individual_double( individual->size );
    double *chromosome = individual->chromosome;
    double *newchromosome = indiv->chromosome;
    
	memcpy(newchromosome, chromosome, individual->size * sizeof(double));
	indiv->fitness = individual->fitness;
	indiv->size = individual->size;
	indiv->id = individual->id;
	indiv->data = NULL;
	return indiv;
}

void _copy_chromosome_double( void *dst, const void *src, int size ){
    memcpy(dst, src, size*sizeof(double) );
}

void Individual_swap( Individual **indiv1, Individual **indiv2 ){
    Individual *indiv = *indiv1;
    *indiv1 = *indiv2;
    *indiv2 = indiv;
}

#pragma mark -
#pragma mark Roulette wheel selection

/**********************************************************************************************************************************
 ************************************** Stochastic - roulette wheel selection *****************************************************
 **********************************************************************************************************************************/

static void ga_roulette_wheel_selection( GA *ga, Individual **sel_parent ){
	double total_fitness = 0;
	int i = 0;
	int isNegative = false;
	double minf = INFINITY;
    int popSize = ga->pop_size;
    assert(popSize > 0);
	
	for ( i = 0; i < popSize; i++) {
		total_fitness += ga->individuals[i]->fitness;
		if( ga->individuals[i]->fitness < minf ) minf = ga->individuals[i]->fitness;
	}
	
	if( minf < 0){
		total_fitness = 0;
		for ( i = 0; i < popSize; i++) {
			total_fitness += ga->individuals[i]->fitness - minf;
		}
		isNegative = true;
		//fprintf(stderr, "minf %f\n",minf);
	}
	
	double rnum;
	
    assert(ga->nelites > 0);
    int count = 0;
    for( ; count < ga->nelites; count++ ){
        ga->copy_chromosome(sel_parent[count]->chromosome, ga->individuals[0]->chromosome, ga->chromosome_size);
		
		sel_parent[count]->fitness = ga->individuals[0]->fitness;
		sel_parent[count]->id = count;
    }
    
	for( ; count < popSize; count++ ){
		double fitnessSoFar = 0.;
		rnum = random_double2(total_fitness);
		
		for ( i = 0; i < popSize; i++){
			fitnessSoFar += ga->individuals[i]->fitness;
			//fprintf(stderr, "total_fitness %f rnum %f  fitnessSoFar %f %f\n",total_fitness,rnum,fitnessSoFar, (fitnessSoFar-minf ));
			if( isNegative ) fitnessSoFar -= minf;
			
			if( fitnessSoFar >= rnum ) break;
		}
        ga->copy_chromosome(sel_parent[count]->chromosome, ga->individuals[i]->chromosome, ga->chromosome_size);
		
		sel_parent[count]->fitness = ga->individuals[i]->fitness;
		sel_parent[count]->id = count;
		
	}
}

#pragma mark -
#pragma mark GA core
/**********************************************************************************************************************************
 ********************************************************** GA CORE ***************************************************************
 **********************************************************************************************************************************/
void ga_core( GA *ga ){
    int popSize = ga->pop_size;
    assert(popSize > 0);
    if( ga->verbosity > 0 ){
		fprintf(stdout, "\nEvaluate initial (random) individuals [GA core]\n" );
	}
    
	// Calculate initial fitnesses
	compute_all_fitnesses(ga);
	
	for ( int i = 0; i < popSize; i++ ) {
		if( ga->individuals[i]->fitness > ga->maxfitness ){
			ga->maxfitness = ga->individuals[i]->fitness;
		}
	}
    
    if ( ga->lookup != NULL ) {
        for ( int i = 0; i < popSize; i++) {
            if ( isnan(ga->individuals[i]->fitness) ) {
                ga->fitness( ga, ga->individuals[i] );
            }
        }
    }
	
    qsort( ga->individuals, popSize, sizeof(Individual*), sort_compare_individuals );
	
	ga->previous_maxfitness = ga->maxfitness = ga->individuals[0]->fitness;
	
	if( ga->verbosity > 0 ){
		fprintf(stdout, "\nInitial best fit model LnL = %f\n\n", ga->maxfitness );
	}
	
    ga->copy_chromosome(ga->best_chromosome, ga->individuals[0]->chromosome, ga->chromosome_size);
    
	// Create temporary population
    Individual ** sel_pop = (Individual **)malloc( popSize*sizeof(Individual*) );
	assert(sel_pop);
	for ( int i = 0; i < popSize; i++) {
		sel_pop[i] = ga->clone_individual(ga->individuals[i]);
	}
	
	//printf("GA start\n");
	
	while( !ga->termination(ga) ){
		
		if( ga->use_annealing ){
			//if( ga->generation % ga->annealing_phase == 0 ){
            if( ga->n_without_improvement > ga->max_no_improvement/2 ){
				ga->mutation_rate = fmax(ga->min_mutation_rate, ga->mutation_rate * ga->cooling_factor );
			}
		}
		
		// Selection step
		switch ( ga->selection ) {
			case GA_ROULETTE:
				ga_roulette_wheel_selection( ga, sel_pop );
				break;
			default:
				assert(0);
		}
		
		//printf("After selection\n");
		//print_individuals( sel_pop, ga->pop_size);
		
        if ( ga->mate != NULL ) {
            if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
                fprintf(stdout, "Mating step\n");
            }
            for ( int i = 0; i < popSize/2; i+=2 ) {
                ga->mate( ga, sel_pop[i], sel_pop[i+1], ga->individuals[i] );
            }
        }
        else {
            for ( int i = 0; i < popSize; i++ ) {
                ga->copy_chromosome( ga->individuals[i]->chromosome, sel_pop[i]->chromosome, ga->chromosome_size);
            }
        }
        
		if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
            fprintf(stdout, "Mutation step\n");
        }
		// Mutation step
		for ( int i = 0; i < popSize; i++) {
            ga->mutate( ga, ga->individuals[i] );
		}
		
		//printf("After mutate %f\n",ga->mutation_rate);
		//print_individuals( sel_pop, ga->pop_size);
				
		if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
            fprintf(stdout, "Evaluate mutated individuals\n" );
        }
        
        compute_all_fitnesses(ga);
        
        if ( ga->lookup != NULL ) {
            for ( int i = 0; i < popSize; i++) {
                if ( isnan(ga->individuals[i]->fitness) ) {
                    ga->fitness( ga, ga->individuals[i] );
                }
            }
        }
        
        qsort( ga->individuals, ga->pop_size, sizeof(Individual*), sort_compare_individuals );
        
        //ga->maxfitness = ga->individuals[0]->fitness;
        
        //ga->copy_chromosome(ga->best_chromosome, ga->individuals[0]->chromosome, ga->chromosome_size);
        
		ga->generation++;
	}
    
    for ( int i = 0; i < popSize; i++ ) {
        ga->free_individual(sel_pop[i]);
	}
	free(sel_pop);
}


#pragma mark -
#pragma mark CHC
/**********************************************************************************************************************************
 ********************************************************** GA CHC  ***************************************************************
 **********************************************************************************************************************************/

// Use this function only for unsigned chromosomes
static bool chc_selection( GA *ga, Individual **sel_parent, int **family_count ){

	if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
		fprintf(stdout, "Mating step... ");
	}

	
	int i = 0;
	double minf = INFINITY;
	
	double *fitnesses = dvector( ga->pop_size);
	//fprintf(stdout, "fitnesses\n");
	for ( i = 0; i < ga->pop_size; i++ ) {
		if( ga->individuals[i]->fitness < minf ) minf = ga->individuals[i]->fitness;
		fitnesses[i] = ga->individuals[i]->fitness;
	}
	
	// rescale fitnesses if there are negative values
	if( minf < 0){
		for ( i = 0; i < ga->pop_size; i++) {
			fitnesses[i] -= minf;
		}
	}
	//fprintf(stderr, "probs\n");
	double *probs = dvector( ga->pop_size );
	double fitness_sum = probs[0] = 1.;
	for ( i = 1; i < ga->pop_size; i++ ) {
		probs[i] = exp( (fitnesses[i] - fitnesses[0])/2.0 );
		fitness_sum += probs[i];
	}
	// rescale to 1
	//fprintf(stderr, "rescale to 1\n");
	for ( i = 0; i < ga->pop_size; i++ ) {
		probs[i] /= fitness_sum;
	}
	
	bool incestuous = false;
	
	//fprintf(stderr, "memset\n");
	for ( i = 0; i < ga->pop_size; i++) {
		memset(family_count[i], 0, ga->pop_size*sizeof(int) );
	}
	
	for( int count = 1; count < ga->pop_size; count++ ){
		
		//fprintf(stderr, "\n--------------\ncount %d\n", count);
		
		double fitnessSoFar = 0.0;
		double rnum = random_double();
		//fprintf(stderr, "rnum %f\n", rnum);
		
		// choose first parent
		int idx_p1 = 0;
		for ( ; idx_p1< ga->pop_size; idx_p1++ ){
			fitnessSoFar += probs[idx_p1];
			if( fitnessSoFar >= rnum ) break;
		}
		//fprintf(stderr, "idx_p1 %d [%f]\n", idx_p1, fitnessSoFar);
		// choose second parent
		int max_iter = ga->pop_size*10;
		int idx_p2 = 0;
		int try = 0;
		bool diff = false;
		for ( ; try < max_iter; try++) {
			idx_p2 = idx_p1;
			// the best solution is way better than the others, would fail in the next loop
			// it happens in the early generations of the GA
			if ( probs[idx_p1] == 1.0 ){
				//fprintf(stderr, "probs[idx_p1] == 1\n");
				while ( idx_p2 == idx_p1 ) {
					idx_p2 = random_int(ga->pop_size-1);
				}
				//fprintf(stderr, "end\n");
			}
			else {
				//fprintf(stderr, "probs[idx_p1] != 1\n");
				// should be able to get out of this loop unless probs[idx_p1] is very high
				while ( idx_p1 == idx_p2 ) {
					rnum = random_double2( 1.0-probs[idx_p1] );
					//fprintf(stderr, "rnum %f\n", rnum);
					//fprintf(stderr, "%d %f %f\n", idx_p1,probs[idx_p1],rnum);
					/*fitnessSoFar = ( idx_p1 == 0 ? 0 : probs[0] );
					
					for ( idx_p2 = 0; idx_p2< ga->pop_size; idx_p2++ ){
						if ( idx_p1 != idx_p2 ) fitnessSoFar += probs[idx_p2];
						if( fitnessSoFar >= rnum ) break;
					}*/
					
					//NEW
					fitnessSoFar = ( idx_p1 == 0 ? 0 : probs[ga->pop_size-1] );
					
					for ( idx_p2 = ga->pop_size-1; idx_p2 >= 0; idx_p2-- ){
						if ( idx_p1 != idx_p2 ) fitnessSoFar += probs[idx_p2];
						if( fitnessSoFar >= rnum ) break;
					}
					//fprintf(stderr, "idx_p2 %d [%f]\n", idx_p2, fitnessSoFar);
					
					// HACK HACK HACK
					// I get a -1 for idx_p2 sometimes from the previous loop
					if ( idx_p2 == -1 ) {
						idx_p2 = idx_p1;
						while ( idx_p2 == idx_p1 ) {
							idx_p2 = random_int(ga->pop_size-1);
						}
					}
				}
				//fprintf(stderr, "end\n");
			}
			
			// check if there are too many identical individuals
			while ( family_count[idx_p1][idx_p2] > 5 ) {
				idx_p2 = random_int(ga->pop_size-1);
			}
			//fprintf(stderr, "final idx_p2 %d\n", idx_p2);
			
			diff  = ga_are_chromosomes_different_unsigned(ga->individuals[idx_p1]->chromosome, ga->individuals[idx_p2]->chromosome, ga->chromosome_size);
			//fprintf(stderr, "diff %d\n", diff);
			
			if( diff ){
				//fprintf(stderr, "try %d id1=%d id2=%d  index1=%d index2=%d, fitness1=%f fitness2=%f\n",try, ga->individuals[idx_p1]->id, ga->individuals[idx_p2]->id, idx_p1,idx_p2,ga->individuals[idx_p1]->fitness, ga->individuals[idx_p2]->fitness );
				family_count[idx_p1][idx_p2]++;
				family_count[idx_p2][idx_p1]++;
				//fprintf(stderr, "break\n");
				break;
			}
		}
		
		// failed to find an appropriate partner for individual 1
		if( try == max_iter ){
			incestuous = true;
			break;
#ifdef DEBUG_CHC_SELECTION
			fprintf(stderr, "Incest! %d %d\n", idx_p1, idx_p2 );
			print_individual_unsigned(ga->individuals[idx_p1]);
			print_individual_unsigned(ga->individuals[idx_p2]);
#endif
		}
		else{
			//fprintf(stderr, "\n------------------------------\n");exit(0);
			//fprintf(stderr, "id1=%d id2=%d  index1=%d idndex2=%d, fitness1=%f fitness2=%f\n",ga->individuals[idx_p1]->id, ga->individuals[idx_p2]->id, idx_p1,idx_p2,ga->individuals[idx_p1]->fitness, ga->individuals[idx_p2]->fitness );
			//fprintf(stderr, "mate\n");
			ga->mate( ga, ga->individuals[idx_p1], ga->individuals[idx_p2], sel_parent[count] );
			//fprintf(stderr, "end mate\n");
			//ga->fitness( ga, sel_parent[count] ); this is calculated later in a threaded for loop

#ifdef DEBUG_CHC_SELECTION

			int n_mate = 0;//compareChromosomes( ga->individuals[idx_p1]->chromosome, sel_parent[count]->chromosome, ga->chromosome_size );
			for (int ii = 0; ii < ga->chromosome_size; ii++) {
				if ( ga->individuals[idx_p1]->chromosome[ii] != sel_parent[count]->chromosome[ii] ) {
					n_mate++;
				}
			}
			fprintf(stderr, "RECOMBINATION (%d/%d)\n", n_mate, ga->individuals[idx_p1]->size);
			fprintf(stderr, "Mating %d[p=%f] %d[p=%f] - %d mismatches\n",idx_p1, probs[idx_p1], idx_p2, probs[idx_p2], diff);
			print_individual( ga->individuals[idx_p1] );
			print_individual( ga->individuals[idx_p2] );
			if ( n_mate != 0 ) {
				fprintf(stderr, "----------------------------------------------------\n");
				print_individual( sel_parent[count] );
			}
			fprintf(stderr, "====================================================\n");

#endif	
		}
	}
	
	if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
		fprintf(stdout, "done\n");
	}
	
	free(probs);
	free(fitnesses);
	return incestuous;
}



void ga_chc( GA *ga ){
    int popSize  = ga->pop_size;
    assert(popSize > 0);
	int i = 0;	
		
	if( ga->verbosity > 0 ){
		fprintf(stdout, "\nEvaluate initial (random) individuals\n" );
	}
	
	// Calculate initial fitnesses
	compute_all_fitnesses(ga);
    
    if ( ga->lookup != NULL ) {
        for ( i = 0; i < popSize; i++) {
            if ( isnan(ga->individuals[i]->fitness) ) {
                ga->fitness( ga, ga->individuals[i] );
            }
        }
    }
	
	qsort( ga->individuals, popSize, sizeof(Individual*), sort_compare_individuals );
	
	ga->previous_maxfitness = ga->maxfitness = ga->individuals[0]->fitness;
	
	if( ga->verbosity > 0 ){
		fprintf(stdout, "\nInitial best fit model LnL = %f\n\n", ga->maxfitness );
	}
	
    ga->copy_chromosome(ga->best_chromosome, ga->individuals[0]->chromosome, ga->chromosome_size);
	
	
	double diversity = ga->diversity( ga );
	bool is_incestuous = (diversity < ga->mutation_threshold);
	
	double *probs = dvector( popSize);	// probability vector
	
	int **family_count = imatrix(popSize, popSize); // control the number of identical individuals
    
    int previous_unique_model = Hashtable_length(ga->lookup);
    int previous_unique_model_counter = 0;
    
    // Create temporary population
    Individual ** sel_pop = (Individual **)malloc( popSize*sizeof(Individual*) );
	assert(sel_pop);
	for ( i = 0; i < popSize; i++ ) {
		sel_pop[i] = ga->clone_individual(ga->individuals[i]);
	}
	
	while( !ga->termination(ga) ){
		
		if( ga->use_annealing ){
			if( ga->generation % ga->annealing_phase == 0 ){
				ga->mutation_rate = fmax(ga->min_mutation_rate, ga->mutation_rate * ga->cooling_factor );
			}
		}
		
		// Muate, except the best individual
		if( is_incestuous ){
			if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
				fprintf(stdout, "Trigger hypermutation: threshold = %f\n", ga->mutation_threshold);
                double accum = 0;
                for ( int i = 0; i < popSize; i++ ) {
                    accum += ga->individuals[i]->fitness;
                }
                fprintf(stdout, "Minimum: %f Median: %f Mean: %f ... ", ga->individuals[popSize-1]->fitness, ga->individuals[popSize/2]->fitness, (accum/popSize) );
			}
			
			#pragma omp parallel for
			for ( i = 1; i < popSize-1; i++) {
#ifdef DEBUG_GA_MUTATE
				
				#pragma omp critical
				{
					unsigned *old = clone_uivector(ga->individuals[i]->chromosome, ga->individuals[i]->size );
					ga->mutate( ga, ga->individuals[i] );
					int n_mutation = 0;
					for (int ii = 0; ii < ga->chromosome_size; ii++) {
						if ( ga->individuals[i]->chromosome[ii] != old[ii] ) {
							n_mutation++;
						}
					}
					fprintf(stderr, "MUTATION (%d/%d): rate %f\n", n_mutation, ga->individuals[i]->size, ga->mutation_rate );
					print_individual(ga->individuals[i]);
					fprintf(stderr, "%3d) fitness %f K =", ga->individuals[i]->id, ga->individuals[i]->fitness);
					for (int ii = 0; ii < ga->individuals[i]->size; ii++) {
						fprintf(stderr, " %d", old[ii]);
					}
					fprintf(stderr, "\n============================================================\n");
					free(old);
					old = NULL;
				}
#else
				ga->mutate( ga, ga->individuals[i] );
				
#endif
			}
			
			// replaced the worst individual with the best one and mutate it
			memcpy(ga->individuals[popSize-1]->chromosome, ga->individuals[0]->chromosome, ga->chromosome_size * sizeof(unsigned));
			ga->mutate( ga, ga->individuals[popSize-1] );
			
            if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
				fprintf(stdout, "done\n");
				fprintf(stdout, "Evaluate mutated individuals\n" );
			}
            
			compute_fitnesses(ga, 1);
			
			if ( ga->lookup != NULL ) {
				for ( i = 1; i < popSize; i++) {
					if ( isnan(ga->individuals[i]->fitness) ) {
						ga->fitness( ga, ga->individuals[i] );
					}
				}
			}
			
			qsort( ga->individuals, popSize, sizeof(Individual*), sort_compare_individuals );
			
			ga->maxfitness = ga->individuals[0]->fitness;
            
            if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
				fprintf(stdout, "Diversity = %f\n", ga->diversity( ga ) );
			}

#ifdef DEBUG_GA_MUTATE
            ga_print_population_unsigned(ga);
			fprintf(stderr, "++++++++++++++++++++++++++++++++++++++++++\n");
#endif
		}
		
		// Selection
		is_incestuous = chc_selection( ga, sel_pop, family_count );
		
		if( !is_incestuous ){
			// assign new population to current population
			// keep the fittest one
			for ( i = 1; i < popSize; i++) {
                Individual_swap(&ga->individuals[i], &sel_pop[i]);
			}
			
			if( ga->verbosity > 0 && ga->generation % ga->feedback_sampling == 0 ){
				fprintf(stdout, "Evaluate next generation of individuals\n" );
			}
			compute_fitnesses(ga, 1);
			
			// In case there are more than 2 individuals in the population that have the same chromosome
			// We would have to evaluate several times the same model
			// use a second pass to fetch from the hastable the fitness
			if ( ga->lookup != NULL ) {
				for ( i = 1; i < popSize; i++) {
					if ( isnan(ga->individuals[i]->fitness) ) {
						ga->fitness( ga, ga->individuals[i] );
					}
				}
			}
						
			qsort( ga->individuals, popSize, sizeof(Individual*), sort_compare_individuals );
			
			diversity = ga->diversity( ga );
			
			if(diversity < ga->mutation_threshold ){
				is_incestuous = true;
				//ga->generation--;
				
				// Even if the diversity is low we could have found a good candidate during the selection process
				//continue;
			}
            
            // check if the current population has already been evaluated
            if( Hashtable_length(ga->lookup) == previous_unique_model ){
                previous_unique_model_counter++;
                if( previous_unique_model_counter == 2 ){
                    is_incestuous = true;
                    previous_unique_model_counter = 0;
                    ga->generation++;
                    continue;
                }
            }
		
		}
        else {
            previous_unique_model_counter = 0;
        }
//		else {
//			ga->generation--;
//			continue;
//		}
		
		
		previous_unique_model = Hashtable_length(ga->lookup);
		
        ga->copy_chromosome(ga->best_chromosome, ga->individuals[0]->chromosome, ga->chromosome_size);
		
		ga->generation++;
		
	}
	
	for ( i = 0; i < popSize; i++ ) {
        ga->free_individual(sel_pop[i]);
	}
	free(sel_pop);
	free(probs);
	
	free_imatrix(family_count, popSize);
}

#pragma mark -
#pragma mark Fitness functions

void compute_all_fitnesses( GA *ga ){
	compute_fitnesses( ga, 0 );
}

void compute_fitnesses( GA *ga, int start_index ){
	if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ) {
		ga->feedback_count = ga->pop_size - start_index;
		fprintf(stdout, "%d/%d\r", ( ga->pop_size-ga->feedback_count), ga->pop_size);
	}
    if(ga->nThreads > 1){
#if defined (PTHREAD_ENABLED)
    // No need to create more threads than jobs
    //int n_threads = imin(ga->nThreads, (ga->pop_size - start_index));
    int n_threads = ga->nThreads;
    
    threadpool_ga_t *threadpool = malloc(sizeof(threadpool_ga_t));
    assert(threadpool);
    threadpool->count = start_index;
    threadpool->ga = ga;
    threadpool->index_treelikelihood = 0;
    threadpool->threads = malloc(n_threads*sizeof(pthread_t));
    assert(threadpool->threads);
    
    pthread_mutex_init(&(threadpool->lock), NULL);
    
    for ( int i = 0; i < n_threads; i++ ) {
        pthread_create( &(threadpool->threads[i]), NULL, ga->thread_worker, threadpool );
    }
    for ( int i = 0; i < n_threads; i++ ) {
        pthread_join(threadpool->threads[i], NULL);
    }
    free(threadpool->threads);
    pthread_mutex_destroy(&(threadpool->lock));
    free(threadpool);
#else
    #pragma omp parallel for
	for ( int i = start_index; i < ga->pop_size; i++) {
		ga->fitness( ga, ga->individuals[i] );
	}
#endif
    }
    else {
        for ( int i = start_index; i < ga->pop_size; i++) {
            ga->fitness( ga, ga->individuals[i] );
        }
    }
    
	if ( ga->verbosity > 0 && ga->feedback_computeall && ga->generation % ga->feedback_sampling == 0 ) {
		fprintf(stdout, "\n\n");
	}
	ga->eval += (ga->pop_size - start_index);//exit(1);
}

#pragma mark -
#pragma mark get/setters

void ga_set_n_threads( GA *ga, const unsigned int nThreads ){
#if defined (PTHREAD_ENABLED)
    ga->nThreads = nThreads;
#elif defined (_OPENMP)
		ga->nThreads = nThreads;
		omp_set_num_threads(nThreads);
#endif
}

unsigned int ga_get_n_threads( const GA *ga ){
	return ga->nThreads;
}

void ga_set_filename( GA *ga, const char *filename ){
	if( filename != NULL ){
		ga->filename = String_clone(filename);
		ga->dump = true;
		ga->dump_frequency = 1;
	}
}

void ga_set_mutation_rate( GA *ga, const double rate ){
	ga->mutation_rate = rate;
}

void ga_set_ngeneration( GA *ga, int ngen){
	ga->max_generation = ngen;
}

void ga_set_nelites( GA *ga, int nelites){
    if( nelites >= ga->pop_size){
        fprintf(stderr, "The number of offspring of the elite (%d) should be lower than the population size (%d)\n", nelites, ga->pop_size);
        exit(1);
    }
	ga->nelites = nelites;
}

#pragma mark -
//MARK: Termination functions


bool ga_default_termination( GA *ga ){
	
	// Little improvement but never worse since we always keep the best individual
	if( ga->generation != 0 && fabs( ga->individuals[0]->fitness - ga->previous_maxfitness) <= 0.001 ) {
		ga->n_without_improvement++;
	}
	else{
		ga->n_without_improvement = 0;
	}
	
	if ( ga->verbosity > 0 ) {
		if ( ga->generation % ga->feedback_sampling == 0 ){
			time(&ga->current_time);
			double t = difftime(ga->current_time, ga->start_time);
			fprintf(stdout, "-----------------------------------------------------------------\n");
			fprintf(stdout, "Generation = %d\n", ga->generation);
			fprintf(stdout, "Number of generation without improvement = %d\n", ga->n_without_improvement);
			fprintf(stdout, "Max fitness = %f\n", ga->maxfitness);
			fprintf(stdout, "Total runtime ");
			print_pretty_time(stdout,t);
			if( ga->lookup != NULL )fprintf(stdout, "Number of unique models evaluated %u (total %u)\n", Hashtable_length(ga->lookup), ga->eval);
			fprintf(stdout, "Diversity = %f\n", ga->diversity(ga) );
			if ( ga->termination_callback != NULL ) {
				ga->termination_callback(ga);
			}
			
			//print_population(ga);
			fprintf(stdout, "\n");
		}
	}
	if( ga->use_tol && fabs( ga->individuals[0]->fitness - ga->previous_maxfitness) < ga->tol ) return true;
	
    //FIXME: IT used to work like this
	//ga->previous_maxfitness = ga->maxfitness = ga->individuals[0]->fitness;
    ga->previous_maxfitness = ga->maxfitness;
	
	ga->interrupted = ga_interrupted;
	if ( ga->interruptible && ga->interrupted ) {
		fprintf(stdout, "Ctrl-c\n\n");
		fprintf(stdout, "0) To exit the program\n");
		fprintf(stdout, "1) Continue the GA\n");
		fprintf(stdout, "2) To stop the GA\n");
		fprintf(stdout, "3) Change maximum number of generations\n");
		fprintf(stdout, "4) Change maximum number of generations without improvement\n");
		//fprintf(stdout, "5) Dump individuals for restart (TODO)\n");
		int answer = -1;
		
		while ( answer < 0 || answer > 4 ) {
			fprintf(stdout, "Choice: ");
			int r = scanf("%d",  &answer);
            if ( r == 0 ) {
                answer = -1;
            }
		}
		
		switch (answer) {
			case 0:{
				exit(SIGINT);
				break;
            }
			case 1:{
				break;
            }
			case 2:{
				return true;
				break;
            }
			case 3:{
				answer = -2;
                int r = 0;
				while ( answer <= 0 ) {
					fprintf(stdout, "Number of generations [%d]:\n", ga->max_generation);
					r = scanf("%d",  &answer);
                    if ( r == 0 ) {
                        answer = -2;
                    }
				}
				ga->max_generation = answer;
				break;
            }
			case 4:{
				answer = -2;
                int r = 0;
				while ( answer <= 0 ) {
					fprintf(stdout, "Number of generations [%d]:\n", ga->max_no_improvement);
					r = scanf("%d",  &answer);
                    if ( r == 0 ) {
                        answer = -2;
                    }
				}
				ga->max_no_improvement = answer;
				break;
            }
			default:{
				assert(0);
            }
		}
	}
	ga->interrupted = ga_interrupted = false;
	
	if( ga->max_generation <= ga->generation) return true;
	if( ga->use_max_no_improvement && ga->max_no_improvement == ga->n_without_improvement ) return true;
	return false;
}

#pragma mark -
#pragma mark Diversity functions

// Genotypic diversity
double ga_diversity_atomic_unsigned( const GA *ga ){
	int j = 0;
	int k = 0;
	double count = 0.;
	int tot = 0;
    unsigned int *chromosome1 = NULL;
    unsigned int *chromosome2 = NULL;
	for (int i = 0; i < ga->pop_size; i++) {
        chromosome1 = ga->individuals[i]->chromosome;
		for ( j = i + 1; j < ga->pop_size; j++) {
            chromosome2 = ga->individuals[j]->chromosome;
			for ( k = 0; k < ga->chromosome_size; k++) {
				if ( chromosome1[k] != chromosome2[k]) count++;
				tot++;
			}
		}
	}
	return count/tot;
}

// Phenotypic diversity

// Requires the individual to be sorted (individual[0] highest fitness)
double ga_diversity_fitness( const GA *ga ){
	return fabs( (ga->individuals[ga->pop_size-1]->fitness - ga->individuals[0]->fitness)/ga->individuals[ga->pop_size-1]->fitness);
}

// Coefficient of variation
double ga_diversity_CV( const GA *ga ){
	double mean = 0;
	for (int i = 0; i < ga->pop_size; i++) {
		mean += ga->individuals[i]->fitness;
	}
	mean /= ga->pop_size;
	double sigma = 0;
	for (int i = 0; i < ga->pop_size; i++) {
		sigma += pow( ga->individuals[i]->fitness - mean ,2);
	}
	sigma = sqrt( sigma/ga->pop_size );
	return sigma/fabs(mean);
}

#pragma mark -
#pragma mark Functions for unsigned chromosomes


// check if a class is assigned at least once but does not check if maximum class in chromosome is really the number of class intended by the user
// make the individual canonical
// Should only be used if classses are assigned at least once
void  ga_default_mate_unsigned( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
	
    unsigned *chromosome1   = individual1->chromosome;
    unsigned *chromosome2   = individual2->chromosome;
    unsigned *newchromosome = new1->chromosome;
    
	// always assume that the parent is ok
	int classCount = uimax_vector(chromosome1, individual1->size)+1;
	
	bool ok = false;
	while ( !ok ) {
		// we need to reset the K otherwise at each iteration it will look more and more like indiv2
		memcpy(newchromosome, chromosome1, individual1->size * sizeof(unsigned int));
		// ignore the root
		for ( int i = 0; i < individual1->size-1; i++) {
			double rnum = random_double();
			if ( rnum < ga->mate_prob ) {
				newchromosome[i] = chromosome2[i];
			}
		}
		ok = check_complete_assignment(newchromosome, classCount, new1->size-1);// don't include the root
	}
	
	MakeStringCanonical2(newchromosome, classCount, new1->size-1); // don't include the root
}

// does not check if a class is assigned at least once
void  ga_mate_imcomplete_unsigned( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
    
    unsigned *chromosome1   = individual1->chromosome;
    unsigned *chromosome2   = individual2->chromosome;
    unsigned *newchromosome = new1->chromosome;
    
	memcpy(newchromosome, chromosome1, individual1->size * sizeof(unsigned int));
	// ignore the root
	for ( int i = 0; i < individual1->size-1; i++) {
		double rnum = random_double();
		if ( rnum < ga->mate_prob ) {
			newchromosome[i] = chromosome2[i];
		}
	}
	
	
}

// make the individual canonical
void  ga_mate_break_point_unsigned( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
    
    unsigned *chromosome1   = individual1->chromosome;
    unsigned *chromosome2   = individual2->chromosome;
    unsigned *newchromosome = new1->chromosome;
    
	// always assume that the parent is ok
	int classCount = uimax_vector(chromosome1, individual1->size)+1;
	bool ok = false;
	while ( !ok ) {
		int break_point = random_int(individual1->size-1);
		if ( break_point != 0 && break_point != individual2->size-1 ) {
			memcpy(newchromosome, chromosome1, break_point * sizeof(unsigned int));
			memcpy(newchromosome+break_point, chromosome2, individual2->size-break_point * sizeof(unsigned int));
		}
		ok = check_complete_assignment(newchromosome, classCount, new1->size-1);// don't include the root
	}
	MakeStringCanonical2(newchromosome, classCount, new1->size-1);	// don't include the root
}

void ga_random_individual_unsigned( Individual *indiv, const int n_classes ){
	
	int i = 0;
    unsigned *chromosome = indiv->chromosome;
	// ensure that each rate class is assigned at least once
	for ( ; i < n_classes; i++ ) {
		chromosome[i] = i;
	}
	
	for ( ; i < indiv->size; i++ ) {
		chromosome[i] = random_int(n_classes-1);
	}
	
	int rpos;
	for ( i = 0; i < n_classes; i++ ) {
		rpos = random_int( indiv->size-2 ); // we don't want to swap the clock of the root
		swap_uint( &chromosome[rpos], &chromosome[i]);
	}
}

int ga_compare_chromosomes_unsigned( const unsigned *k1, const unsigned *k2, const unsigned size ){
	int count = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if ( k1[i] == k2[j] ) {
				count++;
				break;
			}
		}
	}
	return size-count;
}

bool ga_are_chromosomes_different_unsigned( const unsigned *k1, const unsigned *k2, const unsigned size ){
	for ( int i = 0; i < size; i++ ) {
		if ( k1[i] != k2[i] ) {
			return true;
		}		
	}
	return false;
}

static int cmpFunc(const void *p1, const void *p2){
	int **a = (int **)p1;
	int **b = (int **)p2;
	return ( a[0][1] - b[0][1]);
}

// return FASLE if vector = {0,1,3,4} (missing class 2)
bool check_complete_assignment( unsigned *randomModel, int classCount, int dim ){
	
	int *compressedString = ivector( classCount);
	int h = 0;
	for ( h = 0; h < dim; h++){
		compressedString[ randomModel[h] ] = 1;
	}
	for ( h = 0; h < classCount; h++){
		if ( compressedString[h] == 0 ) {
			break;
		}
	}
	free(compressedString);
	compressedString = NULL;
	return ( h == classCount );
}

// same above but we ignore 1 index. This used to screen the root
bool check_complete_assignment_ignore_one( unsigned *randomModel, int classCount, int dim, int ignore ){
	
	int *compressedString = ivector( classCount);
	int h = 0;
	for ( h = 0; h < dim; h++){
        if ( ignore == h) continue;
        
		compressedString[ randomModel[h] ] = 1;
	}
	for ( h = 0; h < classCount; h++){
		if ( compressedString[h] == 0 ) {
			break;
		}
	}
	free(compressedString);
	return ( h == classCount );
}

// Ported to C from HYPHY language
// return false if a class is not represented in the vector
bool MakeStringCanonical ( unsigned *randomModel, const int nClass, const int stateVectorDimension ){
	int classCount = nClass;
	if ( nClass < 0 ) {
		classCount = uimax_vector(randomModel, stateVectorDimension)+1;
	}
	
	int *compressedString = ivector( classCount);
	unsigned v,l;
	int h = 0;
	
	for ( h = 0; h < stateVectorDimension; h++){
		compressedString[ randomModel[h] ] = 1;
	}
	
	compressedString[0] = 0;
	for ( h = 1; h < classCount; h++){
		compressedString[h] = compressedString[h] + compressedString[h-1];
	}
	
	for ( h = 0; h < stateVectorDimension; h++){
		v = compressedString[ randomModel[h] ];
		randomModel[h] = v;
	}
	
	l = v = compressedString[classCount-1]+1;

	if ( v > 1 ){
		int **sortedOrder = imatrix( v,2) ;
		for ( h = 0; h < v; h++){
			sortedOrder[h][0] = -1;
		}
		int cc = 0;
		for ( h = 0; h < stateVectorDimension; h++){
			int hshift = randomModel[h];
			if (sortedOrder[hshift][0] < 0){
				sortedOrder[hshift][0] = cc;
				sortedOrder[hshift][1] = hshift;
				cc++;
			}
		}
		
		qsort(sortedOrder, v, sizeof(sortedOrder[0]), cmpFunc);
		
		for ( h = 0; h < stateVectorDimension; h++){
			v = randomModel[h];
			randomModel[h] = sortedOrder[v][0];
		}
		free_imatrix(sortedOrder, v);
		
	}
	
	free(compressedString);
	compressedString = NULL;
	return (l == stateVectorDimension);
}

void MakeStringCanonical2 ( unsigned *randomModel, const int classCount, const int stateVectorDimension ){
	
	int h = 0;
	int **sortedOrder = imatrix( stateVectorDimension, 2) ;
	for ( h = 0; h < stateVectorDimension; h++){
		sortedOrder[h][0] = -1;
	}
	int cc = 0;
	for ( h = 0; h < stateVectorDimension; h++){
		int hshift = randomModel[h];
		if (sortedOrder[hshift][0] < 0){
			sortedOrder[hshift][0] = cc;
			sortedOrder[hshift][1] = hshift;
			cc++;
		}
	}
	
	qsort(sortedOrder, classCount, sizeof(sortedOrder[0]), cmpFunc);
	
	unsigned v;
	for ( h = 0; h < stateVectorDimension; h++){
		v = randomModel[h];
		randomModel[h] = sortedOrder[v][0];
	}
	free_imatrix(sortedOrder, stateVectorDimension);
}

// take a chromosome with nclasses-1 classes, pick one class and divide it into 2 classes
void ga_nested_individual_unsigned( Individual *indiv, const int nclasses ){
    
    unsigned int *chromosome = indiv->chromosome;
    
	unsigned char *map = (unsigned char *)calloc( nclasses, sizeof(unsigned char) );
	assert(map);
	int i = 0;
	for ( i = 0; i < indiv->size-1; i++) { // ignore the root
		map[ chromosome[i] ]++;
	}
	assert(map[nclasses-1] == 0);
	
	int pick = nclasses-1;
	
	while ( map[pick] <= 1 ) {
		pick = random_int(nclasses-2);
	}
	
	int count =	0;
	for ( i = 0; i < indiv->size-1; i++) {
		if ( chromosome[i] == pick && count < map[pick]/2  ) {
			chromosome[i] = nclasses-1;
			count++;
		}
	}
	
	free(map);
}

// take a chromosome with nclasses-1 classes, pick one class and divide it into 2 classes
void ga_nested_individual2_unsigned( Individual *indiv, const int nclasses ){
	unsigned char *map = (unsigned char *)calloc( nclasses, sizeof(unsigned char) );
	assert(map);
	int i = 0;
	
	int class = 0;
	int max = 0;
	
    unsigned int *chromosome = indiv->chromosome;
    
	for ( i = 0; i < indiv->size-1; i++) { // ignore the root
		map[ chromosome[i] ]++;
		if ( map[ chromosome[i] ] > max ) {
			class = chromosome[i];
		}
	}
	assert(map[nclasses-1] == 0);
	
	int count =	0;
	for ( i = 0; i < indiv->size-1; i++) {
		if ( chromosome[i] == class && count < map[class]/2  ) {
			chromosome[i] = nclasses-1;
			count++;
		}
	}
	
	free(map);
}

void ga_print_individual_unsigned( Individual *indiv ){
    unsigned int *chromosome = indiv->chromosome;
    
	fprintf(stdout, "%3d) fitness %f K =", indiv->id, indiv->fitness);
	for (int i = 0; i < indiv->size; i++) {
		fprintf(stdout, " %u", chromosome[i]);
	}
	//fprintf(stderr, "  [%d]\n", bool2dec(indiv->chromosome, indiv->size) );
	fprintf(stdout, "\n");
}

void ga_print_population_unsigned( GA *ga ){
	for (int i = 0; i < ga->pop_size; i++) {
		ga_print_individual_unsigned( ga->individuals[i] );
	}
	fprintf(stdout, "\n");
}

void ga_print_individuals_unsigned( Individual **individuals, const int n ){
	for (int i = 0; i < n; i++) {
		ga_print_individual_unsigned( individuals[i] );
	}
}

#pragma mark -
#pragma mark Functions for double chromosomes

// Genotypic diversity
double ga_diversity_atomic_double( const GA *ga ){
	int j = 0;
	int k = 0;
	double count = 0.;
	int tot = 0;
    double *chromosome1 = NULL;
    double *chromosome2 = NULL;
	for (int i = 0; i < ga->pop_size; i++) {
        chromosome1 = ga->individuals[i]->chromosome;
		for ( j = i + 1; j < ga->pop_size; j++) {
            chromosome2 = ga->individuals[j]->chromosome;
			for ( k = 0; k < ga->chromosome_size; k++) {
				if ( chromosome1[k] != chromosome2[k]) count++;
				tot++;
			}
		}
	}
	return count/tot;
}

// check if a class is assigned at least once but does not check if maximum class in chromosome is really the number of class intended by the user
// make the individual canonical
// Should only be used if classses are assigned at least once
void  ga_default_mate_double( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
	
    double *chromosome1   = individual1->chromosome;
    double *chromosome2   = individual2->chromosome;
    double *newchromosome = new1->chromosome;
    // we need to reset the K otherwise at each iteration it will look more and more like indiv2
    memcpy(newchromosome, chromosome1, individual1->size * sizeof(double));
    // ignore the root
    for ( int i = 0; i < individual1->size; i++) {
        double rnum = random_double();
        if ( rnum < ga->mate_prob ) {
            newchromosome[i] = chromosome2[i];
        }
    }
}

// Genotypic diversity
double ga_diversity_atomic_bool( const GA *ga ){
	int j = 0;
	int k = 0;
	double count = 0.;
	int tot = 0;
    bool *chromosome1 = NULL;
    bool *chromosome2 = NULL;
	for (int i = 0; i < ga->pop_size; i++) {
        chromosome1 = ga->individuals[i]->chromosome;
		for ( j = i + 1; j < ga->pop_size; j++) {
            chromosome2 = ga->individuals[j]->chromosome;
			for ( k = 0; k < ga->chromosome_size; k++) {
				if ( chromosome1[k] != chromosome2[k]) count++;
				tot++;
			}
		}
	}
	return count/tot;
}

void  ga_default_mate_bool( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 ){
	
    double *chromosome1   = individual1->chromosome;
    double *chromosome2   = individual2->chromosome;
    double *newchromosome = new1->chromosome;
    // we need to reset the K otherwise at each iteration it will look more and more like indiv2
    memcpy(newchromosome, chromosome1, individual1->size * sizeof(bool));
    // ignore the root
    for ( int i = 0; i < individual1->size; i++) {
        double rnum = random_double();
        if ( rnum < ga->mate_prob ) {
            newchromosome[i] = chromosome2[i];
        }
    }
}

#pragma mark -
#pragma mark Misc

int sort_compare_individuals( const void *a, const void *b ){
	const Individual *const *aa = a;
	const Individual *const *bb = b;
	if( (*aa)->fitness < (*bb)->fitness ) return 1;
	else if( (*aa)->fitness == (*bb)->fitness ) return 0;
	return -1;
}

int bool2dec( bool *bin, int len){
	int b, k;
	int sum = 0; 
	
	for(k = 0; k < len; k++){
		b = 1;
		b = b<<(len-k);
		// sum it up
		sum += (bin[k] * b);
		//printf("%d*%d + ",n,b);  // uncomment to show the way this works
	}
	return(sum);
}

int bin2dec(char *bin){
	int  b , k, n;
	int  sum = 0; 
	
	size_t len = strlen(bin) - 1;
	for(k = 0; k <= len; k++) {
		b = 1;
		n = (bin[k] - '0'); // char to numeric value
		if ((n > 1) || (n < 0)) 
		{
			puts("\n\n ERROR! BINARY has only 1 and 0!\n");
			return (0);
		}
		b = b<<(len-k);
		// sum it up
		sum = sum + n * b;
		//printf("%d*%d + ",n,b);  // uncomment to show the way this works
	}
	return(sum);
}

