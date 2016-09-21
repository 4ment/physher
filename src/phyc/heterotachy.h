/*
 *  heterotachy.h
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


#ifndef _HETEROTACHY_H_
#define _HETEROTACHY_H_

#include "treelikelihood.h"
#include "pooledtreelikelihood.h"
#include "ga.h"
#include "modelselection.h"




typedef enum heterotachy_algorithm {
	HETEROTACHY_LOCALCLOCK_EXHAUSTIVE,
	HETEROTACHY_LOCALCLOCK_GREEDY,
	HETEROTACHY_LOCALCLOCK_GA,
	HETEROTACHY_DISCRETE_GA,
	HETEROTACHY_RELAXED_GA
}heterotachy_algorithm;




//MARK: struct definitions

struct _ClockSearch;

typedef struct _ClockSearch ClockSearch;

typedef double (*clocksearch_optimize_f)( ClockSearch* );
typedef void   (*clocksearch_free_f)( ClockSearch* );
typedef void   (*clocksearch_init_f)(ClockSearch*);

struct _ClockSearch {
	void *pDerivedObj; // pointer to subclass (GAClockSearch)

	unsigned max_n_rate;
    
	double alpha;	
	
	char *logFileName;
	char *logFileMode;
	
	int    starting_df;
	double starting_lk;	
	
	unsigned verbosity; 
	
    
     /////////////////////////////////////////////
    // THE VARIABLES BELOW SHOULD BE PROTECTED //
   /////////////////////////////////////////////
    
    heterotachy_algorithm algorithm;
	
	clocksearch_optimize_f optimize;
	clocksearch_free_f     free;
	clocksearch_init_f     init;

   	PooledTreeLikelihood *pool;
    
    unsigned int nThreads;

    double    best_lk;
	unsigned *best_indexes;
	double   *best_heights;
	double   *best_rates; // length(best_rates)-1 == length(best_indexes) for local clocks
    
    int       df; //degrees of freedom for AIC, AICc...
    
    unsigned  n_rate;
    
	int starting_n_rate;
    
	// to save parameters from a round to another
	unsigned *current_indexes;
	double   *current_heights;
	double   *current_rates;
    
    double  default_rate; // for resetting the parameters, not great way to do that
    
	bool initialized;
    
    FILE *logFile;
	unsigned tree_count;
    
	int mode; // 0 all represented, 1 can be imcomplete assignment
    
    double (*IC)(ClockSearch*, double, int);
    information_criterion ic;
    int ic_sample_size;
    
};


typedef void (*pInitGA)(GA *, ClockSearch*);

typedef struct GAClockSearch {
	ClockSearch *pBaseObj;
	unsigned popSize;
	unsigned ngen;
	//unsigned nGeneration;
	pInitGA initGA;
	int tree_count; // number of trees, used for logging the trees
	int max_no_improvement;
	int verbosity;
	int feedback_sampling;
    
    // used ONLY to seed the GA population, removed after the first GA is created
    unsigned **pop_individuals;
    double   **pop_heights;
    double   **pop_rates;
} GAClockSearch;

typedef struct IndividualData {
    double *heights;
    double *rates;
    int     n_rate;
} IndividualData;

#pragma mark -
#pragma mark ClockSearch

double ClockSearch_IC( ClockSearch * search, double lnl, int k );

ClockSearch * new_ClockSearch( SingleTreeLikelihood *tlk, unsigned nThreads );


void free_ClockSearch( ClockSearch *search );

void ClockSearch_init( ClockSearch *search );


void ClockSearch_set_max_n_rate( ClockSearch *search, unsigned max );

void ClockSearch_set_verbosity( ClockSearch *search, int verbosity );

void ClockSearch_set_starting_lnl( ClockSearch *search, double lnl );

void ClockSearch_set_starting_df( ClockSearch *search, int df );

void ClockSearch_set_alpha( ClockSearch *search, double alpha );

void ClockSearch_set_logfile_name( ClockSearch *search, const char *name );

void ClockSearch_set_logfile_mode( ClockSearch *search, char *mode );

void * ClockSearch_derived_object( ClockSearch *search );

void ClockSearch_set_mode( ClockSearch *search, int mode );

void ClockSearch_set_ic( ClockSearch *search, information_criterion ic );

void ClockSearch_set_ic_sample_size( ClockSearch *search, int sample_size );

#pragma mark -
#pragma mark GAClockSearch

//ClockSearch * new_GAClockSearch( TreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize );
ClockSearch * new_GAClockSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize );

void free_GAClockSearch( ClockSearch *search );


#pragma mark -

void free_Individual_GA( Individual *indiv );

Individual * clone_Individual_GA( Individual *indiv );


#pragma mark -

typedef struct ModelTally{
	double lk;
	unsigned count;
} ModelTally;


ModelTally * new_ModelTally( const double value );

void print_ModelTally_Hashtable( FILE *pf, Hashtable *hash, const int nClass );

#endif
