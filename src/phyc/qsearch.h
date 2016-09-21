/*
 *  qsearch.h
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


#ifndef __PhyC__qsearch__
#define __PhyC__qsearch__

#include <stdio.h>


#include "treelikelihood.h"
#include "pooledtreelikelihood.h"
#include "modelselection.h"
#include "ga.h"

struct _QSearch;

typedef struct _QSearch QSearch;

struct _QSearch {
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
    
    
    double (*optimize)( QSearch* );
    void (*free)( QSearch* );
    void (*init)( QSearch* );
    
   	PooledTreeLikelihood *pool;
    
    unsigned int nThreads;
    
    double    best_lk;
    unsigned *best_indexes;
    double   *best_rates;
    
    int       df; //degrees of freedom for AIC, AICc...
    
    unsigned  rateCount;
    int startingRateCount;
    
    int indexCount; // length of the chromosome
    
    // to save parameters from a round to another
    unsigned *current_indexes;
    double   *current_rates;
    
    double  default_rate; // for resetting the parameters, not great way to do that
    
    bool initialized;
    
    FILE *logFile;
    
    double (*IC)(QSearch*, double, int);
    information_criterion ic;
    int ic_sample_size;
    
    int mode;
    
};


QSearch * new_GAQSearch( SingleTreeLikelihood *tlk, const unsigned nThreads, const unsigned popSize );


typedef struct GAQSearch {
    QSearch *pBaseObj;
    unsigned popSize;
    unsigned ngen;
    void (*initGA)(GA *, QSearch*);
    int max_no_improvement;
    int verbosity;
    int feedback_sampling;
    
    // used ONLY to seed the GA population, removed after the first GA is created
    unsigned **pop_individuals;
    double   **pop_rates;
} GAQSearch;

#pragma mark -
#pragma mark LASSO

void qsearch_lasso(SingleTreeLikelihood *tlk );

#endif /* defined(__PhyC__qsearch__) */
