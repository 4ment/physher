/*
 *  ga.h
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


#ifndef _GENETIC_ALGORITHM_H_
#define _GENETIC_ALGORITHM_H_

#include <time.h>


#include "utils.h"
#include "hashtable.h"

//#define PTHREAD_ENABLED 1
#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif


typedef struct _GA GA;

typedef struct Individual{
	unsigned int id;
	//unsigned int *chromosome;
    void *chromosome;
	unsigned int size;
	double fitness;
	
	void *data;
}Individual;


typedef void (*ga_mutate)( GA *ga, Individual *individual );

typedef double (*ga_fitness)( GA *ga, Individual *individual );

typedef void (*ga_mate)( GA *ga, Individual *individual1, Individual *individual2, Individual *offspring );

typedef bool (*ga_termination)( GA *ga );

typedef double (*ga_diversity)( const GA *ga );

typedef void (*ga_log)( GA *ga );

typedef void (*ga_termination_callback)( GA *ga );


typedef Individual * (*ga_new_individual)( const unsigned int chromosome_size );

typedef Individual * (*ga_clone_individual)( Individual *indiv );

typedef void (*ga_free_individual)( Individual *indiv );

typedef void (*ga_copy_chromosome)( void *dst, const void *src, int size );


typedef enum ga_selection_algorithm{
	GA_ROULETTE,
	GA_CHC,
	
	GA_SA
}ga_selection_algorithm;

typedef enum ga_chromosome_type{
	GA_CHROMOSOME_UNSIGNED,
	GA_CHROMOSOME_DOUBLE,
	GA_CHROMOSOME_FLOAT,
	GA_CHROMOSOME_BOOL,
    GA_CHROMOSOME_COMPLEX
}ga_chromosome_type;

struct _GA{
	ga_selection_algorithm selection;
	ga_chromosome_type chromosome_type;
    
	unsigned int chromosome_size;
	unsigned int pop_size;
	Individual **individuals;
	
	unsigned int max_generation;
	unsigned int generation;
	
	double mutation_rate;
	
	double cooling_factor;
	unsigned int annealing_phase;
	bool use_annealing;
	double min_mutation_rate;
    
    int nelites; // The number of offsprings of the best individual that is passed on
	
	double tol;
	bool use_tol;
	
	double maxfitness;
	double previous_maxfitness;
	//unsigned int *best_chromosome;
    void *best_chromosome;
	
	bool use_max_no_improvement;
	int max_no_improvement;
	int n_without_improvement;
	
	ga_mutate  mutate;
	ga_fitness fitness;
	ga_mate    mate;
	
	ga_termination termination;
	ga_termination_callback termination_callback;

	ga_diversity diversity;
	
	ga_log log;
    
	
	FILE *output;
	char *filename;
	bool dump;
	int dump_frequency;
	
	Hashtable *lookup;
	
	void *data;
	
	bool initialized;
	
	// CHC specific
	double mutation_threshold;
	double mate_prob;
	
	unsigned int nThreads; // for OPENMP
	
	time_t start_time;
	time_t current_time;
	
	bool feedback_computeall;
	int feedback_count;
	
	unsigned eval;
	
	int verbosity;
	int feedback_sampling;
	bool interruptible;
	bool interrupted;
    
    
    
    ga_new_individual new_individual;

    ga_clone_individual clone_individual;

    ga_free_individual free_individual;
    
    ga_copy_chromosome copy_chromosome;

#if defined (PTHREAD_ENABLED)
    void * (*thread_worker)(void*);
#endif
};


#if defined (PTHREAD_ENABLED)
typedef struct threadpool_ga_t{
    pthread_t *threads;
    pthread_mutex_t lock;
    GA *ga;      // read-write
    int index_treelikelihood; // read
    int count; // read-write
}threadpool_ga_t;
#endif

#pragma mark -
// MARK: Individual

Individual * new_Individual_unsigned( const unsigned int chromosome_size );

Individual * clone_Individual_unsigned( Individual *individual );

Individual * new_Individual_double( const unsigned int chromosome_size );

Individual * clone_Individual_double( Individual *individual );

Individual * new_Individual_bool( const unsigned int chromosome_size );

Individual * clone_Individual_bool( Individual *individual );


void free_Individual( Individual *individual );

void ga_print_individual_unsigned( Individual *indiv );

void ga_print_individuals_unsigned( Individual **individuals, const int n );


#pragma mark -
// MARK: GA

GA * new_GA( ga_selection_algorithm name, ga_chromosome_type chromosome_type, const unsigned int pop_size, const unsigned int chromosome_size );

void free_GA( GA *ga );

void ga_init( GA *ga );

void ga_evolve( GA *ga );

#pragma mark -
// MARK: getter/setter
void ga_set_filename( GA *ga, const char *filename );

void ga_set_n_threads( GA *ga, const unsigned int nThreads );

unsigned int ga_get_n_threads( const GA *ga );

void ga_set_mutation_rate( GA *ga, const double rate );

void ga_set_ngeneration( GA *ga, int ngen);

void ga_set_nelites( GA *ga, int nelites);


#pragma mark -
// MARK: Diversity calculation

double ga_diversity_fitness( const GA *ga );
double ga_diversity_CV( const GA *ga );

#pragma mark -
#pragma mark Termination function

bool ga_default_termination( GA *ga );

#pragma mark -
#pragma mark Fitness functions

void compute_all_fitnesses( GA *ga );
void compute_fitnesses( GA *ga, int start_index );


#pragma mark -
#pragma mark Functions for unsigned chromosomes

double ga_diversity_atomic_unsigned( const GA *ga );

void  ga_default_mate_unsigned( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );

void  ga_mate_imcomplete_unsigned( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );

void  ga_mate_break_point_unsigned( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );

void ga_print_population_unsigned( GA *ga );

void ga_random_individual_unsigned( Individual *indiv, const int n_classes );

void ga_nested_individual_unsigned( Individual *indiv, const int nclasses );

void ga_nested_individual2_unsigned( Individual *indiv, const int nclasses );


int ga_compare_chromosomes_unsigned( const unsigned *k1, const unsigned *k2, const unsigned size );

bool check_complete_assignment( unsigned *randomModel, int classCount, int dim );

bool check_complete_assignment_ignore_one( unsigned *randomModel, int classCount, int dim, int ignore );

bool MakeStringCanonical ( unsigned *randomModel, const int classCount, const int stateVectorDimension );
void MakeStringCanonical2 ( unsigned *randomModel, const int classCount, const int stateVectorDimension );

#pragma mark -
#pragma mark Functions for double chromosomes

double ga_diversity_atomic_double( const GA *ga );

double ga_diversity_atomic_bool( const GA *ga );

void  ga_default_mate_double( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );

void  ga_default_mate_bool( GA *ga, Individual *individual1, Individual *individual2, Individual *new1 );

#pragma mark -
// MARK: Miscellaneous

int bin2dec(char *bin);
int bool2dec( bool *bin, int len);





#endif
