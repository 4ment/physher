/*
 *  phyboot.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 29/10/12.
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

#ifndef PhyC_phyboot_h
#define PhyC_phyboot_h

#include "treelikelihood.h"
#include "tree.h"
#include "distancematrix.h"
#include "phyclustering.h"
#include "phyresampling.h"


#pragma mark *** Maximum likelihood ***


void SingleTreeLikelihood_bootstrap( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads );

void SingleTreeLikelihood_jackknife( const SingleTreeLikelihood *tlk, const char *output, bool save_patterns, int nthreads );


void SingleTreeLikelihood_bootstrap_openmp( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads );

void SingleTreeLikelihood_jackknife_openmp( const SingleTreeLikelihood *tlk, const char *output, bool save_patterns, int nthreads );

void SingleTreeLikelihood_resampling_openmp( const SingleTreeLikelihood *tlk, resampling_scheme scheme, int count, const char *output, bool save_patterns, int nthreads );


#if defined (PTHREAD_ENABLED)
void SingleTreeLikelihood_bootstrap_threads( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads );

void SingleTreeLikelihood_jackknife_threads( const SingleTreeLikelihood *tlk, const char *output, bool save_patterns, int nthreads );

void SingleTreeLikelihood_resampling_threads( const SingleTreeLikelihood *tlk, resampling_scheme scheme, int count, const char *output, bool save_patterns, int nthreads );
#endif

#pragma mark -

void SingleTreeLikelihood_bootstrap_strict_local_openmp( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads );

void SingleTreeLikelihood_bootstrap_greedy( const SingleTreeLikelihood *tlk, int count, const char *output, const char *prefix, const char *postfix, double lk_strict, double rate_strict, double significance_level, unsigned nthreads );


double ** Phyboot_read_rate_and_heights( const char *filename, int *b, branchmodel *clock );


void Phyboot_annotate_fixed_topology( Tree *tree, const char *file, double ci, const char *jackfile );

void Phyboot_confidence_intervals( const char *boot_file, double ci, const char *jack_file );

void Phyboot_annotate( Tree *tree, const char *filename);


#pragma mark -
#pragma mark *** Distance ***

void Distance_bootstrap( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads );

void Distance_jackknife( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, const char *output, bool save_patterns, int nthreads );


void Distance_resampling_openmp( const Sequences *sequences, resampling_scheme scheme, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads );

void Distance_bootstrap_openmp( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads );

void Distance_jackknife_openmp( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, const char *output, bool save_patterns, int nthreads );


#if defined (PTHREAD_ENABLED)
void Distance_bootstrap_threads( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads );

void Distance_jackknife_threads( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, const char *output, bool save_patterns, int nthreads );

void Distance_resampling_threads( const Sequences *sequences, resampling_scheme scheme, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads );
#endif

#endif
