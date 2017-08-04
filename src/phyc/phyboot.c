/*
 *  phyboot.c
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

#include "phyboot.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <strings.h>

#ifdef _OPENMP
#include <omp.h>
#endif


#if defined (PTHREAD_ENABLED)
#include <pthread.h>
#endif

#include "treelikelihood.h"
#include "optimize.h"
#include "treeio.h"
#include "tree.h"
#include "heterotachy.h"
#include "localclock.h"
#include "descriptivestats.h"
#include "matrix.h"
#include "sequence.h"
#include "sequenceio.h"
#include "statistics.h"
#include "gaussian.h"
#include "boot.h"
#include "nj.h"
#include "upgma.h"

#include "splitsystem.h"

#include "topologyopt.h"


#pragma mark -
#pragma mark *** Maximum likelihood ***

#ifdef PTHREAD_ENABLED
typedef struct threadpool_resampling_t{
    pthread_t *threads;
    pthread_mutex_t lock;
    SingleTreeLikelihood *tlk;// read-only
    resampling_scheme scheme;
    FILE *pfile_trees;
    FILE *pfile_params;
    bool save_patterns;
    char *output;// to save data
    int count; // read-write
    int total; // read-only
}threadpool_resampling_t;

static void * _resampling_thread_worker( void *threadpool  ){
    threadpool_resampling_t *pool = (threadpool_resampling_t*)threadpool;
    SingleTreeLikelihood *tlk = pool->tlk;
    
    StringBuffer *buffer = new_StringBuffer(10);
    
    while ( 1 ) {
        int index;
        pthread_mutex_lock(&(pool->lock));
        if( pool->count == pool->total ){
            pthread_mutex_unlock(&(pool->lock));
            break;
        }
        index = pool->count++;
        pthread_mutex_unlock(&(pool->lock));
        
        SitePattern *sp = NULL;
        if( pool->scheme == RESAMPLING_BOOTSTRAP ){
            sp = SitePattern_bootstrap(pool->tlk->sp);
        }
        else if( pool->scheme == RESAMPLING_JACKKNIFE ){
            sp = SitePattern_jackknife(pool->tlk->sp, index);
        }
        else {
            sp = SitePattern_jackknife_n(pool->tlk->sp, pool->tlk->sp->count/2);
            
        }
        
        SiteModel *sm = clone_SiteModel(pool->tlk->sm);
        Tree *tree = clone_Tree(pool->tlk->tree);
        BranchModel *bm = ( pool->tlk->bm == NULL ? NULL : clone_BranchModel(pool->tlk->bm, tree) );
        
        SingleTreeLikelihood *tlk2 = new_SingleTreeLikelihood(tree, sm, sp, bm);
        OptConfig_copy(&pool->tlk->opt, &tlk2->opt);
        
        if( tlk2->opt.topology_optimize ){
            TopologyOptimizer *topology = new_TopologyOptimizer( tlk2, TREE_SEARCH_PARSIMONY_SPR );
            topology->optimize(topology);
            free_TopologyOptimizer(topology);
            SingleTreeLikelihood_update_all_nodes(tlk2);
        }
        
        double lk = optimize_singletreelikelihood(tlk2);
        
        StringBuffer *buffer_local = new_StringBuffer(100);
        
        if( bm != NULL ){
            Node **nodes = Tree_nodes(tlk2->tree);
            for ( int i = 0; i < Tree_node_count(tlk2->tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    StringBuffer_empty(buffer_local);
                    StringBuffer_append_format(buffer_local, "%e", tlk2->bm->get(tlk2->bm,nodes[i]));
                    Node_set_annotation(nodes[i], "rate", buffer_local->c);
                }
            }
        }
        
        StringBuffer_empty(buffer_local);
        if( tlk->opt.freqs.optimize ){
            
            tlk2->sm->m->update_frequencies(tlk2->sm->m);
            StringBuffer_append_format(buffer_local, "%e,%e,%e,%e,", tlk2->sm->m->_freqs[0], tlk2->sm->m->_freqs[1], tlk2->sm->m->_freqs[2], tlk2->sm->m->_freqs[3]);
            //            for ( int j = 0; j < Parameters_count(tlk2->sm->m->freqs); j++ ) {
            //                StringBuffer_append_format(buffer_local, "%e,", Parameters_value(tlk2->sm->m->freqs, j));
            //            }
        }
        if( tlk->opt.relative_rates.optimize ){
            for ( int j = 0; j < Parameters_count(tlk2->sm->m->rates); j++ ) {
                StringBuffer_append_format(buffer_local, "%e,", Parameters_value(tlk2->sm->m->rates, j));
            }
        }
        for (int i = 0; i < Parameters_count(tlk->sm->rates); i++) {
            if(!Parameters_fixed(tlk->sm->rates, i)){
                StringBuffer_append_format(buffer_local, "%e,", Parameters_value(tlk->sm->rates, i));
            }
        }
        if( buffer_local->length != 0 ){
            StringBuffer_chop(buffer_local);
        }
        
        int n = (pool->scheme == RESAMPLING_JACKKNIFE ? pool->tlk->sp->weights[index] : 1);
        
        pthread_mutex_lock(&(pool->lock));
        {
            for ( int j = 0; j < n; j++) {
                fprintf(pool->pfile_trees, "tree TREE%i [&lk=%f] = [&R] ", index, lk);
                Tree_print_nexus_with_annotation2(pool->pfile_trees, tlk2->tree, bm!=NULL);
                fprintf(pool->pfile_trees, "\n");
            }
            
            if( pool->pfile_params != NULL ){
                for ( int j = 0; j < n; j++) {
                    fprintf(pool->pfile_params, "%s\n", buffer_local->c);
                }
                fflush(pool->pfile_params);
            }
            
            fprintf(stderr, "Replicate %d/%d\r", pool->count, pool->total);
            
        }
        pthread_mutex_unlock(&(pool->lock));
        if ( pool->save_patterns ) {
            StringBuffer_set_string(buffer_local, pool->output);
            StringBuffer_append_format(buffer_local, ".%d", index);
            StringBuffer_append_string(buffer_local, ".fasta");
            SitePattern_save(sp, buffer_local->c);
        }
        free_SingleTreeLikelihood(tlk2);
        free_StringBuffer(buffer_local);
    }
    free_StringBuffer(buffer);
    pthread_exit(NULL);
    return NULL;
}

void SingleTreeLikelihood_bootstrap_threads( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads ){
    SingleTreeLikelihood_resampling_threads(tlk, RESAMPLING_BOOTSTRAP, count, output, save_patterns, nthreads);
}

void SingleTreeLikelihood_jackknife_threads( const SingleTreeLikelihood *tlk, const char *output, bool save_patterns, int nthreads ){
    SingleTreeLikelihood_resampling_threads(tlk, RESAMPLING_JACKKNIFE, tlk->sp->count, output, save_patterns, nthreads);
}

void SingleTreeLikelihood_jackknife_n_threads( const SingleTreeLikelihood *tlk, int replicates, const char *output, bool save_patterns, int nthreads ){
    SingleTreeLikelihood_resampling_threads(tlk, RESAMPLING_JACKKNIFE, replicates, output, save_patterns, nthreads);
}

void SingleTreeLikelihood_resampling_threads( const SingleTreeLikelihood *tlk, resampling_scheme scheme, int count, const char *output, bool save_patterns, int nthreads ){
    
    nthreads = imin(nthreads, count);
    
    StringBuffer *buffer = new_StringBuffer(10);
    
    StringBuffer_append_string(buffer, output);
    StringBuffer_append_string(buffer, ".trees");
    FILE *pfile_trees = fopen(buffer->c,"w");
    FILE *pfile_params = NULL;
    
    if( tlk->opt.freqs.optimize || tlk->opt.relative_rates.optimize || tlk->opt.gamma.optimize || tlk->opt.pinv.optimize ){
        StringBuffer *buffer2 = new_StringBuffer(10);
        StringBuffer_empty(buffer);
        StringBuffer_append_string(buffer, output);
        StringBuffer_append_string(buffer, ".params");
        pfile_params = fopen(buffer->c,"w");
        
        StringBuffer_empty(buffer);
        
        if( tlk->opt.freqs.optimize ){
            tlk->sm->m->update_frequencies(tlk->sm->m);
            StringBuffer_append_string(buffer, "A,C,G,T,");
            StringBuffer_append_format(buffer2, "%e,%e,%e,%e,", tlk->sm->m->_freqs[0], tlk->sm->m->_freqs[1], tlk->sm->m->_freqs[2], tlk->sm->m->_freqs[3]);
            //            for ( int j = 0; j < Parameters_count(tlk->sm->m->freqs); j++ ) {
            //                StringBuffer_append_format(buffer, "%s,", Parameters_name(tlk->sm->m->freqs, j));
            //                StringBuffer_append_format(buffer2, "%e,", Parameters_value(tlk->sm->m->freqs, j));
            //            }
        }
        if ( tlk->opt.relative_rates.optimize) {
            for ( int j = 0; j < Parameters_count(tlk->sm->m->rates); j++ ) {
                StringBuffer_append_format(buffer, "%s,", Parameters_name(tlk->sm->m->rates, j));
                StringBuffer_append_format(buffer2, "%e,", Parameters_value(tlk->sm->m->rates, j));
            }
        }
        for (int i = 0; i < Parameters_count(tlk->sm->rates); i++) {
            if(!Parameters_fixed(tlk->sm->rates, i)){
                StringBuffer_append_format(buffer, "%s,", Parameters_name(tlk->sm->rates, i));
                StringBuffer_append_format(buffer2, "%e,", Parameters_value(tlk->sm->rates, i));
            }
        }
        if( buffer->length != 0 ){
            StringBuffer_chop(buffer);
            StringBuffer_chop(buffer2);
            fprintf(pfile_params, "#%s\n", buffer->c);
            fprintf(pfile_params, "#%s\n", buffer2->c);
            fflush(pfile_params);
        }
        free_StringBuffer(buffer2);
    }
    
    threadpool_resampling_t *threadpool = malloc(sizeof(threadpool_resampling_t));
    assert(threadpool);
    
    threadpool->tlk = (SingleTreeLikelihood*)tlk;
    
    threadpool->save_patterns = save_patterns;
    threadpool->output = (char *)output;
    
    threadpool->count = 0;
    threadpool->total = count;
    threadpool->threads = malloc(nthreads*sizeof(pthread_t));
    assert(threadpool->threads);
    threadpool->scheme = scheme;
    threadpool->pfile_trees = pfile_trees;
    threadpool->pfile_params = pfile_params;
    
    pthread_mutex_init(&(threadpool->lock), NULL);
    
    Tree_print_nexus_header_figtree_Taxa(pfile_trees, tlk->tree);
    Tree_print_nexus_header_figtree_BeginTrees(pfile_trees, tlk->tree);
    
    
    for ( int i = 0; i < nthreads; i++ ) {
        pthread_create( &(threadpool->threads[i]), NULL, _resampling_thread_worker, threadpool );
    }
    for ( int i = 0; i < nthreads; i++ ) {
        pthread_join(threadpool->threads[i], NULL);
    }
    
    fprintf(stderr, "\n");
    fprintf(pfile_trees, "End;");
    fclose(pfile_trees);
    free_StringBuffer(buffer);
    if( pfile_params != NULL ){
        fclose(pfile_params);
    }
    
    free(threadpool->threads);
    pthread_mutex_destroy(&(threadpool->lock));
    free(threadpool);
    
}
#endif

// can be used with any clock.
void SingleTreeLikelihood_bootstrap_openmp( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads ){
    SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_BOOTSTRAP, count, output, save_patterns, nthreads);
}

void SingleTreeLikelihood_jackknife_openmp( const SingleTreeLikelihood *tlk, const char *output, bool save_patterns, int nthreads ){
    SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_JACKKNIFE, tlk->sp->count, output, save_patterns, nthreads);
}

void SingleTreeLikelihood_jackknife_n_openmp( const SingleTreeLikelihood *tlk, int replicates,const char *output, bool save_patterns, int nthreads ){
    SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_JACKKNIFE_PROPORTION, replicates, output, save_patterns, nthreads);
}

void SingleTreeLikelihood_resampling_openmp( const SingleTreeLikelihood *tlk, resampling_scheme scheme, int count, const char *output, bool save_patterns, int nthreads ){
    
#if defined (_OPENMP)
    nthreads = imin(nthreads, count);
#else
    nthreads = 1;
    
#endif
    
    StringBuffer *buffer = new_StringBuffer(10);
    
    StringBuffer_append_string(buffer, output);
    StringBuffer_append_string(buffer, ".trees");
	FILE *pfile_trees = fopen(buffer->c,"w");
    FILE *pfile_params = NULL;
    
    bool optimize_sitemmodel = false;
    for (int i = 0; i < Parameters_count(tlk->sm->rates); i++) {
        if(Parameters_fixed(tlk->sm->rates, i) == false){
            optimize_sitemmodel = true;
        }
    }
    
    
    if( tlk->opt.freqs.optimize || tlk->opt.relative_rates.optimize || optimize_sitemmodel ){
        StringBuffer *buffer2 = new_StringBuffer(10);
        StringBuffer_empty(buffer);
        StringBuffer_append_string(buffer, output);
        StringBuffer_append_string(buffer, ".params");
        pfile_params = fopen(buffer->c,"w");
        
        StringBuffer_empty(buffer);
        
        if( tlk->opt.freqs.optimize ){
            tlk->sm->m->update_frequencies(tlk->sm->m);
            StringBuffer_append_string(buffer, "A,C,G,T,");
            StringBuffer_append_format(buffer2, "%e,%e,%e,%e,", tlk->sm->m->_freqs[0], tlk->sm->m->_freqs[1], tlk->sm->m->_freqs[2], tlk->sm->m->_freqs[3]);
//            for ( int j = 0; j < Parameters_count(tlk->sm->m->freqs); j++ ) {
//                StringBuffer_append_format(buffer, "%s,", Parameters_name(tlk->sm->m->freqs, j));
//                StringBuffer_append_format(buffer2, "%e,", Parameters_value(tlk->sm->m->freqs, j));
//            }
        }
        if ( tlk->opt.relative_rates.optimize) {
            for ( int j = 0; j < Parameters_count(tlk->sm->m->rates); j++ ) {
                StringBuffer_append_format(buffer, "%s,", Parameters_name(tlk->sm->m->rates, j));
                StringBuffer_append_format(buffer2, "%e,", Parameters_value(tlk->sm->m->rates, j));
            }
        }
        for (int i = 0; i < Parameters_count(tlk->sm->rates); i++) {
            if(!Parameters_fixed(tlk->sm->rates, i)){
                StringBuffer_append_format(buffer, "%s,", Parameters_name(tlk->sm->rates, i));
                StringBuffer_append_format(buffer2, "%e,", Parameters_value(tlk->sm->rates, i));
            }
        }
        if( buffer->length != 0 ){
            StringBuffer_chop(buffer);
            StringBuffer_chop(buffer2);
            fprintf(pfile_params, "#%s\n", buffer->c);
            fprintf(pfile_params, "#%s\n", buffer2->c);
            fflush(pfile_params);
        }
        free_StringBuffer(buffer2);
    }
    
    
	Tree_print_nexus_header_figtree_Taxa(pfile_trees, tlk->tree);
	Tree_print_nexus_header_figtree_BeginTrees(pfile_trees, tlk->tree);
    
	int rep = 1;
	
#pragma omp parallel for num_threads(nthreads)
	for ( int i = 0; i < count; i++ ) {
        SitePattern *sp = NULL;
        if( scheme == RESAMPLING_BOOTSTRAP ){
            sp = SitePattern_bootstrap(tlk->sp);
        }
        else if( scheme == RESAMPLING_JACKKNIFE ){
            sp = SitePattern_jackknife(tlk->sp,i);
            
        }
        else {
            sp = SitePattern_jackknife_n(tlk->sp, tlk->sp->count/2);
            
        }
        
        SiteModel *sm = clone_SiteModel(tlk->sm);
        Tree *tree = clone_Tree(tlk->tree);
        BranchModel *bm = ( tlk->bm == NULL ? NULL : clone_BranchModel(tlk->bm, tree) );
        
        SingleTreeLikelihood *tlk2 = new_SingleTreeLikelihood(tree, sm, sp, bm);
		OptConfig_copy(&tlk->opt, &tlk2->opt);
        
        if( tlk2->opt.topology_optimize ){
            TopologyOptimizer *topology = new_TopologyOptimizer( tlk2, TREE_SEARCH_PARSIMONY_SPR );            
            topology->optimize(topology);
            free_TopologyOptimizer(topology);
            SingleTreeLikelihood_update_all_nodes(tlk2);
        }
    
        double lk = optimize_singletreelikelihood(tlk2);
        
        StringBuffer *buffer_local = new_StringBuffer(100);
		
        if( bm != NULL ){
            Node **nodes = Tree_nodes(tlk2->tree);
            for ( int i = 0; i < Tree_node_count(tlk2->tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    StringBuffer_empty(buffer_local);
                    StringBuffer_append_format(buffer_local, "%e", tlk2->bm->get(tlk2->bm,nodes[i]));
                    Node_set_annotation(nodes[i], "rate", buffer_local->c);
                }
            }
        }
        
        StringBuffer_empty(buffer_local);
        
        if( tlk->opt.freqs.optimize ){
            
            tlk2->sm->m->update_frequencies(tlk2->sm->m);
            StringBuffer_append_format(buffer_local, "%e,%e,%e,%e,", tlk2->sm->m->_freqs[0], tlk2->sm->m->_freqs[1], tlk2->sm->m->_freqs[2], tlk2->sm->m->_freqs[3]);
//            for ( int j = 0; j < Parameters_count(tlk2->sm->m->freqs); j++ ) {
//                StringBuffer_append_format(buffer_local, "%e,", Parameters_value(tlk2->sm->m->freqs, j));
//            }
        }
        if( tlk->opt.relative_rates.optimize ){
            for ( int j = 0; j < Parameters_count(tlk2->sm->m->rates); j++ ) {
                StringBuffer_append_format(buffer_local, "%e,", Parameters_value(tlk2->sm->m->rates, j));
            }
        }
        for (int i = 0; i < Parameters_count(tlk->sm->rates); i++) {
            if(!Parameters_fixed(tlk->sm->rates, i)){
                StringBuffer_append_format(buffer_local, "%e,", Parameters_name(tlk->sm->rates, i));
            }
        }
        if( buffer_local->length != 0 ){
            StringBuffer_chop(buffer_local);
        }
        
        int n = (scheme == RESAMPLING_JACKKNIFE ? tlk->sp->weights[i] : 1);
#pragma omp critical
		{
            for ( int j = 0; j < n; j++) {
                fprintf(pfile_trees, "tree TREE%i [&lk=%f] = [&R] ", i, lk);
                Tree_print_nexus_with_annotation2(pfile_trees, tlk2->tree, bm!=NULL);
                fprintf(pfile_trees, "\n");
            }
            
            if( pfile_params != NULL ){
                for ( int j = 0; j < n; j++) {
                    fprintf(pfile_params, "%s\n", buffer_local->c);
                }
                fflush(pfile_params);
            }
            
			fprintf(stderr, "Replicate %d/%d\r", rep++, count);
		}
        
        if ( save_patterns ) {
            StringBuffer_set_string(buffer_local, output);
            StringBuffer_append_format(buffer_local, ".%d",i);
            StringBuffer_append_string(buffer_local, ".fasta");
            SitePattern_save(sp, buffer_local->c);
        }
        free_SingleTreeLikelihood(tlk2);
        free_StringBuffer(buffer_local);
	}
    fprintf(stderr, "\n");
	fprintf(pfile_trees, "End;");
	fclose(pfile_trees);
	
	free_StringBuffer(buffer);
    if( pfile_params != NULL ){
        fclose(pfile_params);
    }
}

void SingleTreeLikelihood_bootstrap( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads ){
#if defined (PTHREAD_ENABLED)
    if( nthreads == 1 ){
        SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_BOOTSTRAP, count, output, save_patterns, nthreads);
    }
    else {
        printf("Using PThread\n");
        SingleTreeLikelihood_resampling_threads(tlk, RESAMPLING_BOOTSTRAP, count, output, save_patterns, nthreads);
    }
#else
#ifdef _OPENMP
    printf("Using OpenMP\n");
#endif
    SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_BOOTSTRAP, count, output, save_patterns, nthreads);
#endif
}

void SingleTreeLikelihood_jackknife( const SingleTreeLikelihood *tlk, const char *output, bool save_patterns, int nthreads ){
#if defined (PTHREAD_ENABLED)
    if( nthreads == 1 ){
        SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_JACKKNIFE, tlk->sp->count, output, save_patterns, nthreads);
    }
    else {
        printf("Using PThread\n");
        SingleTreeLikelihood_resampling_threads(tlk, RESAMPLING_JACKKNIFE, tlk->sp->count, output, save_patterns, nthreads);
    }
#else
#ifdef _OPENMP
    printf("Using OpenMP\n");
#endif
    SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_JACKKNIFE, tlk->sp->count, output, save_patterns, nthreads);
#endif
}

void SingleTreeLikelihood_jackknife_n( const SingleTreeLikelihood *tlk, int replicates, const char *output, bool save_patterns, int nthreads ){
#if defined (PTHREAD_ENABLED)
    if( nthreads == 1 ){
        SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_JACKKNIFE_PROPORTION, replicates, output, save_patterns, nthreads);
    }
    else {
        printf("Using PThread\n");
        SingleTreeLikelihood_resampling_threads(tlk, RESAMPLING_JACKKNIFE_PROPORTION, replicates, output, save_patterns, nthreads);
    }
#else
#ifdef _OPENMP
    printf("Using OpenMP\n");
#endif
    SingleTreeLikelihood_resampling_openmp(tlk, RESAMPLING_JACKKNIFE_PROPORTION, replicates, output, save_patterns, nthreads);
#endif
}


#pragma mark -

// used for boostrapping local clock
// for each bootstrap a strict clock calculated first and then the greedy algorithm is applied. It starts from 1 local clock
// tThe input branchmode can be either strict or local (NOTHING ELSE)
void SingleTreeLikelihood_bootstrap_strict_local_openmp( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads ){
    
#if defined (_OPENMP)
    nthreads = imin(nthreads, count);
#else
    nthreads = 1;
    
#endif
    
    StringBuffer *buffer = new_StringBuffer(30);
    StringBuffer_append_string(buffer, output);
    StringBuffer_append_string(buffer, ".strict.boot.trees");
    
	FILE *pfile = fopen(buffer->c,"w");
    
    StringBuffer_empty(buffer);
    StringBuffer_append_string(buffer, output);
    StringBuffer_append_string(buffer, ".greedy.boot.trees");
    
    FILE *pfilelocal = fopen(buffer->c,"w");
	
	SiteModel **models = (SiteModel**)malloc(nthreads * sizeof(SiteModel*));
    assert(models);
	Tree **trees = (Tree**)malloc(nthreads * sizeof(Tree*));
    assert(trees);
    BranchModel **branchmodels = (BranchModel**)malloc(nthreads * sizeof(BranchModel*));
    assert(branchmodels);
	
    BranchModel *old_bm = tlk->bm;
    
    double rate = BranchModel_mean_rate_scaled(tlk->bm);
    
	for ( int i = 0; i < nthreads; i++ ) {
		models[i] = clone_SiteModel(tlk->sm);
		trees[i]  = clone_Tree(tlk->tree);
        
        if(old_bm->name != CLOCK_STRICT){
            branchmodels[i] = new_StrictClock(trees[i]);
            Parameters_set_value(branchmodels[i]->rates, 0, rate);
        }
        else {
            branchmodels[i] = clone_BranchModel(old_bm, trees[i]);
        }
	}
    
    double *heights = dvector(Tree_node_count(tlk->tree));
    Tree_heights_to_vector(tlk->tree, heights);

    
	Tree_print_nexus_header_figtree_Taxa(pfile, tlk->tree);
	Tree_print_nexus_header_figtree_BeginTrees(pfile, tlk->tree);
    
    Tree_print_nexus_header_figtree_Taxa(pfilelocal, tlk->tree);
	Tree_print_nexus_header_figtree_BeginTrees(pfilelocal, tlk->tree);
    
    
	int rep = 1;
    
    fprintf(stdout, "Replicate 0/%d\r", count);
    fflush(stdout);

#pragma omp parallel for num_threads(nthreads)
	for ( int i = 0; i < count; i++ ) {
		SitePattern *sp = SitePattern_bootstrap(tlk->sp);
		
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
        //printf("tid %d %d\n",tid,omp_get_num_threads());
#endif
        
		BranchModel *bm = branchmodels[tid];
		SingleTreeLikelihood *tlk2 = new_SingleTreeLikelihood(trees[tid], models[tid], sp, bm);
		OptConfig_copy(&tlk->opt, &tlk2->opt);
		
        Tree_vector_to_heights(heights, trees[tid]);
        Tree_constraint_heights(trees[tid]);
    
        bm->set(bm, 0, rate);
        
        
		double clock_strict_lnl = optimize_singletreelikelihood(tlk2);
        
        double clock_strict_rate = Parameters_value(tlk->bm->rates, 0);
        
        #pragma omp critical
		{
			fprintf(pfile, "tree TREE%d [&lk=%f] = [&R] ", i, clock_strict_lnl);
            Node **nodes = Tree_nodes(tlk2->tree);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", tlk2->bm->get(tlk2->bm,nodes[0]));
            for ( int i = 0; i < Tree_node_count(tlk2->tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    Node_set_annotation(nodes[i], "rate", buffer->c);
                }
            }
            Tree_print_nexus_with_annotation(pfile, tlk2->tree);
			fprintf(pfile, "\n");
            fflush(pfile);
			
			if ( save_patterns ) {
				StringBuffer_set_string(buffer, output);
				StringBuffer_append_format(buffer, ".%d",i);
				StringBuffer_append_string(buffer, ".fasta");
				SitePattern_save(sp, buffer->c);
			}
			
		}
        
        BranchModel *bm_local = new_LocalClock( tlk2->tree, 1 );
        Parameters_set_value(bm_local->rates, 0, clock_strict_rate);
        Parameters_set_value(bm_local->rates, 1, clock_strict_rate);
        
        //SingleTreeLikelihood_set_BranchModel(tlk2, bm_local, false);
        tlk2->bm = bm_local; // dirty
        SingleTreeLikelihood_update_all_nodes(tlk2);
        
        ClockSearch *greedy_local = new_LocalClockSearch( tlk2, 1 );
        ClockSearch_set_starting_lnl(greedy_local, clock_strict_lnl);
        
//        StringBuffer_set_string(buffer, output);
//        StringBuffer_append_format(buffer, ".%d",i);
//        StringBuffer_append_string(buffer, ".greedy.trees");
//        ClockSearch_set_logfile_name(greedy_local, buffer->c);

        ClockSearch_set_verbosity(greedy_local, 0);
        //ClockSearch_set_max_n_rate(greedy_local, 3);
        double clock_local_lnl = greedy_local->optimize(greedy_local);

        #pragma omp critical
		{
			fprintf(pfilelocal, "tree TREE%d [&lk=%f] = [&R] ", i, clock_local_lnl);
            Node **nodes = Tree_get_nodes(tlk2->tree, POSTORDER);
            
            for ( int i = 0; i < Tree_node_count(tlk2->tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    StringBuffer_empty(buffer);
                    StringBuffer_append_format(buffer, "%e", tlk2->bm->get(tlk2->bm,nodes[i]));
                    Node_set_annotation(nodes[i], "rate", buffer->c);
                    
                    if ( tlk2->bm->indicators[i] ) {
                        Node_set_annotation(nodes[i], "local", "1");
                    }
                }
            }
            Tree_print_nexus_with_annotation(pfilelocal, tlk2->tree);
			fprintf(pfilelocal, "\n");
            fflush(pfilelocal);
			
			fprintf(stdout, "Replicate %d/%d\r", rep++, count);
            fflush(stdout);
		}
        greedy_local->free(greedy_local);
		bm_local->free(bm_local, false);
        
        free_SingleTreeLikelihood_share2(tlk2, true, true, false, true);
	}
    fprintf(stdout, "\n");
	fprintf(pfile, "End;");
	fprintf(pfilelocal, "End;");
	fclose(pfile);
	fclose(pfilelocal);
	
	for ( int i = 0; i < nthreads; i++ ) {
		free_SiteModel( models[i] );
		free_Tree(trees[i]);
        branchmodels[i]->free(branchmodels[i], false);
	}
	free(models);
	free(trees);
	free(branchmodels);
	free_StringBuffer(buffer);
    free(heights);
    
}


// optimize topology and other parameters for each bootstraps and calculate strict clock
// output is a stem
void Bootstrap_strict_openmp( const SingleTreeLikelihood *tlk, int count, const char *output, bool save_patterns, int nthreads ){

#ifndef _OPENMP
    nthreads = 1;
#endif
    
    // save the best tree topology and the distances
    StringBuffer *buffer = new_StringBuffer(100);
    Tree_to_newick(buffer, Tree_root(tlk->tree));
    char *newick = StringBuffer_tochar(buffer);
    double *distances = dvector(Tree_node_count(tlk->tree));
    Tree_branch_length_to_vector(tlk->tree, distances);
    
    StringBuffer_empty(buffer);
    StringBuffer_append_string(buffer, ".free.tree");
	FILE *dfile = fopen(buffer->c,"w");
    
    StringBuffer_set_string(buffer, ".strict.tree");
	FILE *tfile = fopen(buffer->c,"w");
	
	SiteModel **models = (SiteModel**)malloc(nthreads * sizeof(SiteModel*));
    assert(models);
	for ( int i = 0; i < nthreads; i++ ) {
		models[i] = clone_SiteModel(tlk->sm);
	}
    double rate = tlk->bm->get(tlk->bm, 0);
	
	Tree_print_nexus_header_figtree_Taxa(dfile, tlk->tree);
	Tree_print_nexus_header_figtree_BeginTrees(dfile, tlk->tree);
	
	Tree_print_nexus_header_figtree_Taxa(tfile, tlk->tree);
	Tree_print_nexus_header_figtree_BeginTrees(tfile, tlk->tree);
    
	StringBuffer_empty(buffer);
    
    OptConfig opt;
    OptConfig_copy(&tlk->opt, &opt);
	
    int rep = 1;
	
#pragma omp parallel for num_threads(nthreads)
	for ( int i = 0; i < count; i++ ) {
		SitePattern *sp = SitePattern_bootstrap(tlk->sp);
        int tid = 0;
		#ifdef _OPENMP
            tid = omp_get_thread_num();
            //printf("tid %d %d\n",tid,omp_get_num_threads());
        #endif
        
		
        Tree *tree = new_Tree(newick, true);
        Tree_vector_to_branch_length(tree, distances);
		SingleTreeLikelihood *tlk2 = new_SingleTreeLikelihood(tree, models[tid], sp, NULL);
        OptConfig_copy(&tlk->opt, &opt);
		
		double lnl1 = optimize_singletreelikelihood(tlk2);
        
        #pragma omp critical
		{
			fprintf(dfile, "tree TREE0 [&LnL=%f] = [&U] ", lnl1);
			//print_tree_nexus_treelikelihood(dfile, tlk2);
            
            Node **nodes = Tree_get_nodes(tlk2->tree, POSTORDER);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", tlk2->bm->get(tlk2->bm,0));
            for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                Node_empty_annotation(nodes[i]);
            }
            Tree_print_nexus_with_annotation(dfile, tlk2->tree);
			fprintf(dfile, "\n");
			
			//fprintf(stderr, "Replicate %d/%d LnL = %f\n", rep++, count, lnl );
			
			if ( save_patterns ) {
				StringBuffer_set_string(buffer, output);
				StringBuffer_append_format(buffer, ".%d",i);
				StringBuffer_append_string(buffer, ".fasta");
				SitePattern_save(sp, buffer->c);
			}
			
		}
        
        // need to root the tree
        
        BranchModel *bm = new_BranchModel(tree, CLOCK_STRICT);
        bm->set(bm, 0, rate);
        
        SingleTreeLikelihood_set_BranchModel(tlk2, bm, false);
        
        double lnl2 = optimize_singletreelikelihood(tlk2);
		
#pragma omp critical
		{
			fprintf(tfile, "tree TREE0 [&LnL=%f] = [&R] ", lnl2);
            Node **nodes = Tree_get_nodes(tlk2->tree, POSTORDER);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", tlk2->bm->get(tlk2->bm,0));
            for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    Node_set_annotation(nodes[i], "rate", buffer->c);
                }
            }
            Tree_print_nexus_with_annotation(tfile, tlk2->tree);
            
			//print_tree_nexus_treelikelihood(tfile, tlk2);
			fprintf(tfile, "\n");
			
			fprintf(stderr, "Replicate %d/%d LnL = %f  Strict LnL: %f\n", rep++, count, lnl1, lnl2 );
			
		}
		
        free_SingleTreeLikelihood_share2(tlk2, false, true, false, false);
	}
	fprintf(dfile, "End;");
	fclose(dfile);
	fclose(tfile);
	
	for ( int i = 0; i < nthreads; i++ ) {
		free_SiteModel( models[i] );
	}
	free(models);
	free_StringBuffer(buffer);
    free(newick);
    free(distances);
}



// The parallelization is done within the greedy search.
// If there is less branches than  threads it would be better to parallelize the bootstrap (#branches < threads < bootstrap)
// we should reset the original Treelikelihood to its initial parameter values

void SingleTreeLikelihood_bootstrap_greedy( const SingleTreeLikelihood *tlk, int count, const char *output, const char *prefix, const char *postfix, double lk_strict, double rate_strict, double significance_level, unsigned nthreads ){
	FILE *pfile = fopen(output,"w");
	Tree_print_nexus_header_figtree_Taxa(pfile, tlk->tree);
	Tree_print_nexus_header_figtree_BeginTrees(pfile, tlk->tree);
	
	BranchModel *bm = new_LocalClock( tlk->tree, 1 );
	
	StringBuffer *buff = new_StringBuffer(50);
	
	fprintf(stderr, "\nBootstrapping: %d replicates\n", count);
	
	for ( int i = 0; i < count; i++ ) {
		StringBuffer_empty(buff);
		StringBuffer_append_format(buff, "%s.%d.%s", prefix, i, postfix);
		Sequences *seqs = readSequences(buff->c);
		SitePattern *sp = new_SitePattern(seqs);
		free_Sequences(seqs);
		
		
		// rate free
		SingleTreeLikelihood *tlk2 = new_SingleTreeLikelihood(tlk->tree, tlk->sm, sp, NULL);
#if defined (SSE3_ENABLED) || (AVX_ENABLED)
		tlk2->use_SIMD = tlk->use_SIMD;
#endif
		OptConfig_copy(&tlk->opt, &tlk2->opt);
		tlk2->opt.bl.optimize              = true;
		tlk2->opt.freqs.optimize           = true;
		tlk2->opt.relative_rates.optimize  = true;
		tlk2->opt.heights.optimize         = false;
		tlk2->opt.rates.optimize           = false;
		
		
		double lk_ratefree = optimize_singletreelikelihood(tlk2);
		fprintf(stderr, "Rate free %f\n", lk_ratefree);
		
		// strict
		tlk2->opt.bl.optimize             = false;
		tlk2->opt.freqs.optimize          = false;
		tlk2->opt.relative_rates.optimize = false;
		tlk2->opt.heights.optimize        = true;
		tlk2->opt.rates.optimize          = true;
		
        for (int i = 0; i < Parameters_count(tlk->sm->rates); i++) {
            Parameters_set_fixed(tlk->sm->rates, true, i);
        }
        
		BranchModel *bm_strict = new_StrictClock( tlk->tree );
		SingleTreeLikelihood_set_BranchModel(tlk2, bm_strict, false);
		bm_strict->set( bm_strict, 0, rate_strict );		
		
		lk_strict = optimize_singletreelikelihood(tlk2);
		
		fprintf(stderr, "Strict clock %f\n", lk_strict);
		
		// local
		LocalClock_set_number_of_clocks(bm, 1);
		Parameters_set_all_value(bm->rates, rate_strict);
		
		SingleTreeLikelihood_set_BranchModel(tlk2, bm, false);
		OptConfig_copy(&tlk->opt, &tlk2->opt);	
		
		//fprintf(stderr, "Not optimized LnL = %f (strict rate %f)\n", tlk2->calculate(tlk2), rate_strict );
		
		ClockSearch *greedy_local = new_LocalClockSearch( tlk2, nthreads );
		greedy_local->starting_lk = lk_strict;
		greedy_local->algorithm = HETEROTACHY_LOCALCLOCK_GREEDY;
		greedy_local->alpha = significance_level;
		
		double lk = greedy_local->optimize(greedy_local);
		
		fprintf(stderr, "Replicate %d/%d LnL = %f\n", i, count, lk );
		
		// reset rates and heights with the best model
		LocalClock_set_number_of_clocks(bm, greedy_local->n_rate);
		localclock_set_indicators2(bm, greedy_local->best_indexes);
		BranchModel_vector_to_rates(bm, greedy_local->best_rates);

		Tree_vector_to_heights(greedy_local->best_heights, tlk2->tree);
		Tree_constraint_heights(tlk2->tree);
		
		fprintf(pfile, "tree TREE0 [&lk=%f] = [&R] ", lk);
        Node **nodes = Tree_get_nodes(tlk2->tree, POSTORDER);
        StringBuffer *buffer = new_StringBuffer(10);
        for ( int i = 0; i < Tree_node_count(tlk2->tree); i++ ) {
            Node_empty_annotation(nodes[i]);
            if( !Node_isroot(nodes[i]) ){
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "%e", tlk2->bm->get(tlk2->bm, nodes[i]));
                Node_set_annotation(nodes[i], "rate", buffer->c);
                
                if( tlk2->bm->indicators[i] ){
                    Node_set_annotation(nodes[i], "local", "1");
                }
            }
        }
        Tree_print_nexus_with_annotation(pfile, tlk2->tree);
		//print_tree_nexus_treelikelihood(pfile, tlk2);
		fprintf(pfile, "\n");
        free_StringBuffer(buffer);

		greedy_local->free(greedy_local);
		
        free_SingleTreeLikelihood_share2(tlk2, true, true, false, true);
	}
	fprintf(pfile, "End;");
	fclose(pfile);
	
	bm->free(bm, false);
	free_StringBuffer(buff);
	fprintf(stderr, "\n\n");
}

double ** Phyboot_read_rate_and_heights( const char *filename, int *b, branchmodel *clock ){
	
	FileReader *reader = new_FileReader(filename,100);
	
	char *ptr = NULL;
	int count = 0;
	
    double **matrix = NULL;
    int matrix_capacity = 100;
    
    StringBuffer *buffer = new_StringBuffer(10);
    StringBuffer *buffer2 = new_StringBuffer(10);
    
    const char *NEXUS_TREE_TRANSLATE_BLOCK = "TRANSLATE";
	char **names = NULL;
	unsigned names_count = 0;
    
    int dim = 0;
	
	while ( reader->read_line(reader) ) {
		ptr = reader->line;
		
        if ( String_i_start_with(reader->line, NEXUS_TREE_TRANSLATE_BLOCK, true) ) {
			StringBuffer_empty(buffer);
			names = parseNexusTranslateBlock(reader, buffer, &names_count);
		}
        
		if ( String_start_with(reader->line, "tree", true) ) {
            assert(names);
            
			while ( *ptr != '(' ) ptr++;
			char *pch = strrchr(ptr, ';');
			if ( pch != NULL ) {
				ptr[pch-ptr] = '\0';
			}
			else error("missing ; in tree file or spaces at the end\n");
			
            StringBuffer_empty(buffer);
			pch = ptr;
			while ( *pch != '\0' ) {
				// taxa id
				if ( (*pch >= 48 && *pch <= 57) && ( *(pch-1) == '(' || *(pch-1) == ',') ) {
					StringBuffer_empty(buffer2);
					while ( *pch != ':' ) {
						StringBuffer_append_char(buffer2, *pch);
						pch++;
					}
					int index = atoi(buffer2->c);
					//fprintf(stderr, "index %d/%d %s %s\n", index, names_count, buffer2->c,names[index-1]);
					StringBuffer_append_string(buffer, names[index-1]);
					StringBuffer_append_char(buffer, ':');
				}
				else {
					StringBuffer_append_char(buffer, *pch);
				}
				pch++;
			}
			
			Tree *t = new_Tree(buffer->c, false);
            
            if( count == 0 ){
                dim = Tree_node_count(t); // space for the node heights
                
                if( String_contains_str(buffer->c, "local=1") ) *clock = CLOCK_LOCAL;
                else if( String_contains_str(buffer->c, "class=") ) *clock = CLOCK_DISCRETE;
                else *clock = CLOCK_STRICT;
                
                if( String_contains_str(buffer->c, "local=1") ){
                    dim += Tree_node_count(t)*2; // extra space for the rates and local cock assignment
                }
                else if (String_contains_str(buffer->c, "class=") ){
                    dim += Tree_node_count(t); // extra space for the rates
                }
                else {
                    dim++; // for strict clock we only need one rate
                }
                matrix = dmatrix(dim, matrix_capacity);
            }
            
            if( matrix_capacity == count ){
                matrix_capacity *= 2;
                for ( int i = 0; i < dim; i++ ) {
                    matrix[i] = realloc(matrix[i], matrix_capacity*sizeof(double));
                }
            }
            
            Node **nodes = Tree_get_nodes(t, POSTORDER);
            
            if( *clock == CLOCK_STRICT ){
                matrix[Tree_node_count(t)][count] = Node_get_double_from_info(nodes[0], "rate=");
            }
            else {
                for ( int i = 0; i < Tree_node_count(t); i++ ) {
                    matrix[Tree_node_count(t)+i][count] = Node_get_double_from_info(nodes[i], "rate=");
                    
                    if( String_contains_str(nodes[i]->info, "local=1") ){
                        matrix[Tree_node_count(t)*2+i][count] = 1;
                    }
                }
            }
            
			for ( int i = 0; i < Tree_node_count(t); i++ ) {
				matrix[i][count] = Node_height(nodes[i]);
			}
            
			free_Tree(t);
            
			count++;
		}
	}
	free_FileReader(reader);
	free_StringBuffer(buffer);
	free_StringBuffer(buffer2);
	free_cmatrix(names, names_count);
    
    *b = count;
    
    if(count == 0 ) return NULL;
    
    if( matrix_capacity != count ){
        for ( int i = 0; i <= dim; i++ ) {
            matrix[i] = realloc(matrix[i], count*sizeof(double));
        }
    }
    
	return matrix;
}



// Annotate tree with confidence intervals calculated from different bootstrap scheme (normal, percentile, BCa)
// tree: MLE tree
void Phyboot_annotate_fixed_topology( Tree *tree, const char *boot_file, double ci, const char *jack_file ){
    
    branchmodel clock_type;
    
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    
    int nboot = 0;
    
    printf("\nCalculating confidence intervals using percentile and normal bootstrap methods\n");
    
    fprintf(stderr, "\nReading trees... ");
    
    double **matrix = Phyboot_read_rate_and_heights(boot_file, &nboot, &clock_type);
    fprintf(stderr, "done [%d]\n\n", clock_type);
    
    double qlb = (1.0-ci)/2.0;
    double qub = (1.0-ci)/2.0+ci;
    
    StringBuffer *buffer  = new_StringBuffer(10);
    StringBuffer *buffer2 = new_StringBuffer(10);
    
    int nNodes = Tree_node_count(tree);
    
    double rate_median, rate_mean, rate_p_low, rate_p_high, rate_n_low, rate_n_high, rate_bca_low, rate_bca_high, rate_mle;
    double height_median, height_mean, height_p_low, height_p_high, height_bca_low, height_bca_high, height_mle;
    
    if ( clock_type == CLOCK_STRICT ){
        qsort(matrix[nNodes], nboot, sizeof(double), qsort_asc_dvector);
        
        rate_median = dmedian_ordered(matrix[nNodes], nboot);
        rate_mean   = dmean(matrix[nNodes], nboot);
        
        if( Node_annotation(nodes[0], "rate") != NULL ){
            rate_mle = atof(Node_annotation(nodes[0], "rate"));
        }
        else {
            rate_mle = Node_get_double_from_info(nodes[0], "rate=");
        }
        
        // Percentile confidence intervals
        rate_p_low  = dpercentile_ordered(matrix[nNodes], nboot, qlb);
        rate_p_high = dpercentile_ordered(matrix[nNodes], nboot, qub);
        
        // Normal confidence intervals
        rate_n_low  = boot_ci_norm2(matrix[nNodes], nboot, rate_mle, qlb);
        rate_n_high = boot_ci_norm2(matrix[nNodes], nboot, rate_mle, qub);
    }
    
    for ( int i = 0; i < nNodes; i++ ) {
        
        //height
        qsort(matrix[i], nboot, sizeof(double), qsort_asc_dvector);
        
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "%e", Node_height(nodes[i]));
        Node_set_annotation(nodes[i], "height", buffer->c);
        
        if( !Node_isleaf(nodes[i]) ){
            height_median = dmedian_ordered(matrix[i], nboot);
            height_mean   = dmean(matrix[i], nboot);
            height_p_low  = dpercentile_ordered(matrix[i], nboot, qlb);
            height_p_high = dpercentile_ordered(matrix[i], nboot, qub);
            height_mle = Node_height(nodes[i]);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", height_mean);
            Node_set_annotation(nodes[i], "height_mean", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", height_median);
            Node_set_annotation(nodes[i], "height_median", buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", height_p_low, height_p_high);
            StringBuffer_append_format(buffer2, "height_perc_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", boot_ci_norm2(matrix[i], nboot, height_mle, qlb), boot_ci_norm2(matrix[i], nboot, height_mle, qub) );
            StringBuffer_append_format(buffer2, "height_norm_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            if( Node_isroot(nodes[i]) ){
                
                fprintf(stdout, ". Root age\n");
                fprintf(stdout, "    Mean:    %e\n",  height_mean);
                fprintf(stdout, "    Median:  %e\n", height_median );
                fprintf(stdout, "    Range:  %e %e\n", matrix[i][0],matrix[i][nboot-1] );
                fprintf(stdout, "    MLE:  %e\n", height_mle );
                fprintf(stdout, "    Percentile %.0f%% CI: [%e %e]\n", (ci*100), height_p_low, height_p_high);
                fprintf(stdout, "    Normal     %.0f%% CI: [%e %e]\n", (ci*100), boot_ci_norm2(matrix[i], nboot, height_mle, qlb), boot_ci_norm2(matrix[i], nboot, height_mle,qub));
            }
        }
        
        // rates
        if( !Node_isroot(nodes[i]) ){
            
            if ( clock_type != CLOCK_STRICT ){
                qsort(matrix[nNodes+i], nboot, sizeof(double), qsort_asc_dvector);
                
                rate_median = dmedian_ordered(matrix[nNodes+i], nboot);
                rate_mean   = dmean(matrix[nNodes+i], nboot);
                
                rate_mle = Node_get_double_from_info(nodes[i], "rate=");
                
                rate_p_low  = dpercentile_ordered(matrix[nNodes+i], nboot, qlb);
                rate_p_high = dpercentile_ordered(matrix[nNodes+i], nboot, qub);
                
                rate_n_low  = boot_ci_norm2(matrix[nNodes+i], nboot, rate_mle, qlb);
                rate_n_high = boot_ci_norm2(matrix[nNodes+i], nboot, rate_mle, qub);
            }
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", rate_mean);
            Node_set_annotation(nodes[i], "rate_mean", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", rate_median);
            Node_set_annotation(nodes[i], "rate_median", buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", rate_p_low, rate_p_high);
            StringBuffer_append_format(buffer2, "rate_perc_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", rate_n_low, rate_n_high );
            StringBuffer_append_format(buffer2, "rate_norm_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            // count the number of times a local clock was assigned to this node
            if( clock_type == CLOCK_LOCAL ){
                int nlocal = 0;
                for ( int j = 0; j < nboot; j++ ) {
                    nlocal += matrix[nNodes*2+i][j];
                }
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "%d", nlocal);
                Node_set_annotation(nodes[i], "local_count", buffer->c);
            }
            
        }
    }
    
    if ( clock_type == CLOCK_STRICT ){
        fprintf(stdout, ". Substitution rate\n");
        fprintf(stdout, "    Mean:    %e\n",  rate_mean);
        fprintf(stdout, "    Median:  %e\n", rate_median );
        fprintf(stdout, "    MLE:  %e\n", rate_mle );
        fprintf(stdout, "    Percentile %.0f%% CI: [%e %e]\n", (ci*100), rate_p_low, rate_p_high);
        fprintf(stdout, "    Normal     %.0f%% CI: [%e %e]\n", (ci*100), rate_n_low, rate_n_high);
    }
    
    if( jack_file != NULL && file_exists(jack_file) ){
        printf("\nCalculating confidence intervals using bias-corrected and accelerated bootstrap method (BCa)\n");
        fprintf(stdout, "\nReading trees from %s... ", jack_file);
        int count = 0;// number of sitepatterns or trees in the file
        branchmodel clock_type2;
        double ** jack1 = Phyboot_read_rate_and_heights(jack_file, &count, &clock_type2);
        fprintf(stdout, "done\n\n");
        
        if( clock_type == CLOCK_STRICT ){
            rate_mle = Node_get_double_from_info(nodes[0], "rate=");
            
            rate_bca_low  = bootci_BCa(jack1[nNodes], count, matrix[nNodes], nboot, rate_mle, qlb);
            rate_bca_high = bootci_BCa(jack1[nNodes], count, matrix[nNodes], nboot, rate_mle, qub);
        }
        
        for ( int i = 0; i < nNodes; i++ ) {
            
            //height
            if( !Node_isleaf(nodes[i]) ){
                height_mle = Node_height(nodes[i]);
                
                height_bca_low  = bootci_BCa(jack1[i], count, matrix[i], nboot, height_mle, qlb);
                height_bca_high = bootci_BCa(jack1[i], count, matrix[i], nboot, height_mle, qub);
                
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", height_bca_low, height_bca_high);
                StringBuffer_append_format(buffer2, "height_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
                
                if( Node_isroot(nodes[i]) ){
                    fprintf(stdout, ". Root age\n");
                    fprintf(stdout, "    BCa %.0f%% CI: [%e %e]\n", (ci*100), height_bca_low, height_bca_high);
                }
            }
            
            // rate
            if( !Node_isroot(nodes[i]) ){
                if ( clock_type != CLOCK_STRICT ){
                    rate_mle = Node_get_double_from_info(nodes[i], "rate=");
                    
                    rate_bca_low  = bootci_BCa(jack1[nNodes+i], count, matrix[nNodes+i], nboot, rate_mle, qlb);
                    rate_bca_high = bootci_BCa(jack1[nNodes+i], count, matrix[nNodes+i], nboot, rate_mle, qub);
                }
                
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", rate_bca_low, rate_bca_high);
                StringBuffer_append_format(buffer2, "rate_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            }
            
        }
        
        if( clock_type == CLOCK_STRICT ){
            
            fprintf(stdout, ". Substitution rate\n");
            fprintf(stdout, "    BCa %.0f%% CI: [%e %e]\n", (ci*100), rate_bca_low, rate_bca_high);
        }
        
        free_dmatrix(jack1, Tree_node_count(tree)+1);
    }
    
    if( clock_type == CLOCK_STRICT) free_dmatrix(matrix, Tree_node_count(tree)+1);
    else free_dmatrix(matrix, Tree_node_count(tree)*3);
    
    
    free_StringBuffer(buffer);
    free_StringBuffer(buffer2);
    
    fprintf(stdout, "\n");
    
}

// Annotate tree with confidence intervals calculated from different bootstrap scheme (normal, percentile, BCa)
// tree: MLE tree
void Phyboot_annotate_fixed_topology2( Tree *tree, const char *boot_file, double ci, const char *jack_file ){
    
    branchmodel clock_type;
    
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    
    int nboot = 0;
    
    printf("\nCalculating confidence intervals using percentile and normal bootstrap methods\n");
    
    fprintf(stderr, "\nReading trees... ");
    
    double **matrix = Phyboot_read_rate_and_heights(boot_file, &nboot, &clock_type);
    fprintf(stderr, "done [%d]\n\n", clock_type);
    
    double qlb = (1.0-ci)/2.0;
    double qub = (1.0-ci)/2.0+ci;
    
    StringBuffer *buffer  = new_StringBuffer(10);
    StringBuffer *buffer2 = new_StringBuffer(10);
    
    int nNodes = Tree_node_count(tree);
    
    double rate_median, rate_mean, rate_p_low, rate_p_high, rate_n_low, rate_n_high, rate_bca_low, rate_bca_high, rate_mle;
    double height_median, height_mean, height_p_low, height_p_high, height_bca_low, height_bca_high, height_mle;
    
    if ( clock_type == CLOCK_STRICT ){
        qsort(matrix[nNodes], nboot, sizeof(double), qsort_asc_dvector);
        
        rate_median = dmedian_ordered(matrix[nNodes], nboot);
        rate_mean   = dmean(matrix[nNodes], nboot);
        
        if( Node_annotation(nodes[0], "rate") != NULL ){
            rate_mle = atof(Node_annotation(nodes[0], "rate"));
        }
        else {
            rate_mle = Node_get_double_from_info(nodes[0], "rate=");
        }
        
        // Percentile confidence intervals
        rate_p_low  = dpercentile_ordered(matrix[nNodes], nboot, qlb);
        rate_p_high = dpercentile_ordered(matrix[nNodes], nboot, qub);
        
        // Normal confidence intervals
        rate_n_low  = boot_ci_norm2(matrix[nNodes], nboot, rate_mle, qlb);
        rate_n_high = boot_ci_norm2(matrix[nNodes], nboot, rate_mle, qub);
    }
    
    for ( int i = 0; i < nNodes; i++ ) {
        
        //height
        qsort(matrix[i], nboot, sizeof(double), qsort_asc_dvector);
        
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "%e", Node_height(nodes[i]));
        Node_set_annotation(nodes[i], "height", buffer->c);
        
        if( !Node_isleaf(nodes[i]) ){
            height_median = dmedian_ordered(matrix[i], nboot);
            height_mean   = dmean(matrix[i], nboot);
            height_p_low  = dpercentile_ordered(matrix[i], nboot, qlb);
            height_p_high = dpercentile_ordered(matrix[i], nboot, qub);
            height_mle = Node_height(nodes[i]);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", height_mean);
            Node_set_annotation(nodes[i], "height_mean", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", height_median);
            Node_set_annotation(nodes[i], "height_median", buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", height_p_low, height_p_high);
            StringBuffer_append_format(buffer2, "height_perc_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", boot_ci_norm2(matrix[i], nboot, height_mle, qlb), boot_ci_norm2(matrix[i], nboot, height_mle, qub) );
            StringBuffer_append_format(buffer2, "height_norm_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            if( Node_isroot(nodes[i]) ){
                
                fprintf(stdout, ". Root age\n");
                fprintf(stdout, "    Mean:    %e\n",  height_mean);
                fprintf(stdout, "    Median:  %e\n", height_median );
                fprintf(stdout, "    Range:  %e %e\n", matrix[i][0],matrix[i][nboot-1] );
                fprintf(stdout, "    MLE:  %e\n", height_mle );
                fprintf(stdout, "    Percentile %.0f%% CI: [%e %e]\n", (ci*100), height_p_low, height_p_high);
                fprintf(stdout, "    Normal     %.0f%% CI: [%e %e]\n", (ci*100), boot_ci_norm2(matrix[i], nboot, height_mle, qlb), boot_ci_norm2(matrix[i], nboot, height_mle,qub));
            }
        }
        
        // rates
        if( !Node_isroot(nodes[i]) ){
            
            if ( clock_type != CLOCK_STRICT ){
                qsort(matrix[nNodes+i], nboot, sizeof(double), qsort_asc_dvector);
                
                rate_median = dmedian_ordered(matrix[nNodes+i], nboot);
                rate_mean   = dmean(matrix[nNodes+i], nboot);
                
                rate_mle = Node_get_double_from_info(nodes[i], "rate=");
                
                rate_p_low  = dpercentile_ordered(matrix[nNodes+i], nboot, qlb);
                rate_p_high = dpercentile_ordered(matrix[nNodes+i], nboot, qub);
                
                rate_n_low  = boot_ci_norm2(matrix[nNodes+i], nboot, rate_mle, qlb);
                rate_n_high = boot_ci_norm2(matrix[nNodes+i], nboot, rate_mle, qub);
            }
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", rate_mean);
            Node_set_annotation(nodes[i], "rate_mean", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", rate_median);
            Node_set_annotation(nodes[i], "rate_median", buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", rate_p_low, rate_p_high);
            StringBuffer_append_format(buffer2, "rate_perc_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", rate_n_low, rate_n_high );
            StringBuffer_append_format(buffer2, "rate_norm_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            // count the number of times a local clock was assigned to this node
            if( clock_type == CLOCK_LOCAL ){
                int nlocal = 0;
                for ( int j = 0; j < nboot; j++ ) {
                    nlocal += matrix[nNodes*2+i][j];
                }
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "%d", nlocal);
                Node_set_annotation(nodes[i], "local_count", buffer->c);
            }
            
        }
    }
    
    if ( clock_type == CLOCK_STRICT ){
        fprintf(stdout, ". Substitution rate\n");
        fprintf(stdout, "    Mean:    %e\n",  rate_mean);
        fprintf(stdout, "    Median:  %e\n", rate_median );
        fprintf(stdout, "    MLE:  %e\n", rate_mle );
        fprintf(stdout, "    Percentile %.0f%% CI: [%e %e]\n", (ci*100), rate_p_low, rate_p_high);
        fprintf(stdout, "    Normal     %.0f%% CI: [%e %e]\n", (ci*100), rate_n_low, rate_n_high);
    }
    
    if( jack_file != NULL && file_exists(jack_file) ){
        printf("\nCalculating confidence intervals using bias-corrected and accelerated bootstrap method (BCa)\n");
        fprintf(stdout, "\nReading trees from %s... ", jack_file);
        int count = 0;// number of sitepatterns or trees in the file
        branchmodel clock_type2;
        double ** jack1 = Phyboot_read_rate_and_heights(jack_file, &count, &clock_type2);
        fprintf(stdout, "done\n\n");
        
        if( clock_type == CLOCK_STRICT ){
            rate_mle = Node_get_double_from_info(nodes[0], "rate=");
            
            rate_bca_low  = bootci_BCa(jack1[nNodes], count, matrix[nNodes], nboot, rate_mle, qlb);
            rate_bca_high = bootci_BCa(jack1[nNodes], count, matrix[nNodes], nboot, rate_mle, qub);
        }
        
        for ( int i = 0; i < nNodes; i++ ) {
            
            //height
            if( !Node_isleaf(nodes[i]) ){
                height_mle = Node_height(nodes[i]);
                
                height_bca_low  = bootci_BCa(jack1[i], count, matrix[i], nboot, height_mle, qlb);
                height_bca_high = bootci_BCa(jack1[i], count, matrix[i], nboot, height_mle, qub);
                
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", height_bca_low, height_bca_high);
                StringBuffer_append_format(buffer2, "height_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
                
                if( Node_isroot(nodes[i]) ){
                    fprintf(stdout, ". Root age\n");
                    fprintf(stdout, "    BCa %.0f%% CI: [%e %e]\n", (ci*100), height_bca_low, height_bca_high);
                }
            }
            
            // rate
            if( !Node_isroot(nodes[i]) ){
                if ( clock_type != CLOCK_STRICT ){
                    rate_mle = Node_get_double_from_info(nodes[i], "rate=");
                    
                    rate_bca_low  = bootci_BCa(jack1[nNodes+i], count, matrix[nNodes+i], nboot, rate_mle, qlb);
                    rate_bca_high = bootci_BCa(jack1[nNodes+i], count, matrix[nNodes+i], nboot, rate_mle, qub);
                }
                
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", rate_bca_low, rate_bca_high);
                StringBuffer_append_format(buffer2, "rate_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            }
            
        }
        
        if( clock_type == CLOCK_STRICT ){
            
            fprintf(stdout, ". Substitution rate\n");
            fprintf(stdout, "    BCa %.0f%% CI: [%e %e]\n", (ci*100), rate_bca_low, rate_bca_high);
        }
        
        free_dmatrix(jack1, Tree_node_count(tree)+1);
    }
    
    if( clock_type == CLOCK_STRICT) free_dmatrix(matrix, Tree_node_count(tree)+1);
    else free_dmatrix(matrix, Tree_node_count(tree)*3);
    
    
    free_StringBuffer(buffer);
    free_StringBuffer(buffer2);
    
    fprintf(stdout, "\n");
    
}

void Phyboot_confidence_intervals( const char *boot_file, double ci, const char *jack_file ){
    
    int nboot = 0;
    int nboot_capcity = 50;
    int count = 0; // number of parameters
    int i = 0;
    
    double **matrix = NULL;
    double *mles    = NULL;
    char **params_str = NULL;
    
    printf("\nCalculating confidence intervals using percentile and normal bootstrap methods\n");
    
    fprintf(stderr, "\nReading parameters... ");

    FileReader *reader = new_FileReader(boot_file,100);
    
    // First line is the parameter description
    reader->read_line(reader);
    params_str = String_split_char(&reader->line[1], ',', &count);
    assert(count > 0);
    matrix = (double**)malloc(count*sizeof(double*));
    assert(matrix);
    for ( i = 0; i < count; i++ ) {
        matrix[i] = dvector(nboot_capcity);
    }
    
    // Second line contains MLEs
    reader->read_line(reader);
    int mle_count;
    mles = String_split_char_double(&reader->line[1], ',', &mle_count);
    assert(count == mle_count);
    
    while ( reader->read_line(reader) ) {
        if( String_is_empty(reader->line) || reader->line[0] == '#') continue;
        int count2;
        double *temp = String_split_char_double(reader->line, ',', &count2);
        if(nboot == nboot_capcity){
            nboot_capcity += 50;
            for ( i = 0; i < count; i++ ) {
                matrix[i] = realloc(matrix[i], nboot_capcity*sizeof(double));
            }
        }
        assert(count == count2);
        for ( i = 0; i < count; i++ ) {
            matrix[i][nboot] = temp[i];
        }
        free(temp);
        nboot++;
    }
    assert(nboot > 0);
    free_FileReader(reader);
    
    fprintf(stderr, "done\n\n");
    
    double qlb = (1.0-ci)/2.0;
    double qub = (1.0-ci)/2.0+ci;
    
    printf(", MLE, Percentile %.0f%% CI, Normal %.0f%% CI",(ci*100),(ci*100));
    if( jack_file ){
        printf(" BCa %.0f%% CI",(ci*100));
    }
    printf("\n");
    
    for ( i = 0; i < count; i++ ) {
        qsort(matrix[i], nboot, sizeof(double), qsort_asc_dvector);
        
        fprintf(stdout, "%s, %e", params_str[i], mles[i]);
        
        // Percentile confidence intervals
        double p_low  = dpercentile_ordered(matrix[i], nboot, qlb);
        double p_high = dpercentile_ordered(matrix[i], nboot, qub);
        
        // Normal confidence intervals
        double n_low  = boot_ci_norm2(matrix[i], nboot, mles[i], qlb);
        double n_high = boot_ci_norm2(matrix[i], nboot, mles[i], qub);
        
        fprintf(stdout, ", [%e %e], [%e %e]", p_low, p_high, n_low, n_high);
        
        if( jack_file ){
            fprintf(stdout, ", [%e %e]", p_low, p_high);
        }
        printf("\n");
        
    }
    free(mles);
    for ( i = 0; i < count; i++ ) {
        free(params_str[i]);
    }
    free(params_str);
    free_dmatrix(matrix, count);
    
    fprintf(stdout, "\n");
    
}

void Phyboot_annotate( Tree *tree, const char *filename){
    
    Hashtable *hash = new_Hashtable_string(Tree_tip_count(tree));
    int count = 0;
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    
    int nNodes = Tree_node_count(tree);
    int splitCount = nNodes - Tree_tip_count(tree) - 1;
    
    int *counts = ivector(splitCount);
    
    for ( int i = 0; i < nNodes; i++ ) {
        if ( Node_isleaf(nodes[i]) ) {
            Hashtable_add(hash, String_clone(Node_name(nodes[i])), new_Int(count++));
        }
    }
    bool **splitTree = getSplitSystem( hash, tree );
    
    FileReader *reader = new_FileReader(filename,100);
    
    bool nexus = false;
    while ( reader->read_line(reader) ) {
        if( String_is_empty(reader->line)) continue;
        if( strcasecmp(reader->line, "#nexus") == 0 ){
            nexus = true;
        }
        break;
        
    }
    free_FileReader(reader);
    
    reader = new_FileReader(filename,100);
    
    if( nexus ){
        TreeFileIterator *iter = new_TreeFileIterator(filename);
        char *line = NULL;
        while ( (line = TreeFileIterator_next_tree(iter)) != NULL  ) {
            Tree *t = new_Tree(line, false);
            bool **split = getSplitSystem( hash, t );
            
            int j = 0;
            for ( int i = 0; i < nNodes-1; i++ ) {
                if ( !Node_isleaf(nodes[i]) ) {
                    if ( hasSplit(split, splitTree[j], splitCount, nNodes) ) counts[j]++;
                    j++;
                }
            }
            free(line);
            free_Tree(t);
            free_bmatrix(split, splitCount);
        }
        free_TreeFileIterator(iter);
    }
    else{
        while ( reader->read_line(reader) ) {
            //printf("%s\n", reader->line);
            if( String_is_empty(reader->line)) continue;
            Tree *t = new_Tree(reader->line, false);
            bool **split = getSplitSystem( hash, t );
            
            int j = 0;
            for ( int i = 0; i < nNodes-1; i++ ) {
                if ( !Node_isleaf(nodes[i]) ) {
                    if ( hasSplit(split, splitTree[j], splitCount, nNodes) ) counts[j]++;
                    j++;
                }
            }
            free_Tree(t);
            free_bmatrix(split, splitCount);
        }
    }
    
    int j = 0;
    StringBuffer *buffer = new_StringBuffer(10);
    for ( int i = 0; i < nNodes-1; i++ ) {
        if ( !Node_isleaf(nodes[i]) && !(Node_isroot(Node_parent(nodes[i])) && Node_isleaf(Node_sibling(nodes[i]))) ) {
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%d", counts[j]);
            Node_set_annotation(nodes[i], "boot", buffer->c);
            j++;
        }
    }
    free_StringBuffer(buffer);
    free_bmatrix(splitTree, splitCount);
    free(counts);
    free_FileReader(reader);
    free_Hashtable(hash);
}

#pragma mark -
#pragma mark *** Distance ***

#if defined (PTHREAD_ENABLED)

typedef struct threadpool_resampling_distance_t{
    pthread_t *threads;
    pthread_mutex_t lock;
    const Sequences *sequences;
    distancematrix_model model;
    clustering_algorithm algorithm;
    resampling_scheme scheme;
    FILE *pfile;
    bool save_patterns;
    char *output;// to save data
    int count; // read-write
    int total; // read-only
}threadpool_resampling_distance_t;

static void * _resampling_thread_worker_distance( void *threadpool  ){
    threadpool_resampling_distance_t *pool = (threadpool_resampling_distance_t*)threadpool;
    
    StringBuffer *buffer = new_StringBuffer(10);
    
    while ( 1 ) {
        
        int count;
        pthread_mutex_lock(&(pool->lock));
        if( pool->count == pool->total ){
            pthread_mutex_unlock(&(pool->lock));
            break;
        }
        count = pool->count++;
        pthread_mutex_unlock(&(pool->lock));
        
        Sequences *resampled  = NULL;
        if ( pool->scheme == RESAMPLING_BOOTSTRAP ) {
            resampled = Sequences_bootstrap(pool->sequences);
        }
        else if( pool->scheme == RESAMPLING_JACKKNIFE){
            resampled = Sequences_jackknife(pool->sequences, count);
        }
        else if( pool->scheme == RESAMPLING_JACKKNIFE_PROPORTION){
            resampled = Sequences_jackknife_n(pool->sequences, pool->sequences->length/2);
        }
        else {
            error("Resampling scheme not recognized\n");
        }
        
        double **matrix = Sequences_distance(resampled, pool->model);
        Tree *tree = NULL;
        if ( pool->algorithm == CLUSTERING_NJ ) {
            tree = new_NJ(resampled, matrix);
        }
        else if( pool->algorithm == CLUSTERING_UPGMA){
            tree = new_UPGMA(resampled, matrix);
        }
        else {
            error("clustering algorithm not yet implemented\n");
        }
        
        free_dmatrix(matrix, resampled->size);
        
        pthread_mutex_lock(&(pool->lock));
        {
            Tree_print_newick(pool->pfile, tree, false);
            //fprintf(pool->pfile, "\n");
#ifdef FEEDBACK
            printf("#%d\n", (count+1));
#else
            fprintf(stderr, "Replicate %d/%d\r", (count+1), pool->total);
#endif
        }
        pthread_mutex_unlock(&(pool->lock));
        
        free_Tree(tree);
        
        if ( pool->save_patterns ) {
            StringBuffer_set_string(buffer, pool->output);
            StringBuffer_append_format(buffer, ".%d", count);
            StringBuffer_append_string(buffer, ".fasta");
            Sequences_save_fasta(resampled, buffer->c);
        }
        free_Sequences(resampled);
    }
    free_StringBuffer(buffer);
    pthread_exit(NULL);
    return NULL;
}

void Distance_resampling_threads( const Sequences *sequences, resampling_scheme scheme, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads ){
#ifdef FEEDBACK
    printf("#=%d",count);
#endif
    
    nthreads = imin(nthreads, count);
    
    FILE *pfile = fopen(output,"w");
    nthreads = imin(nthreads, count);
    threadpool_resampling_distance_t *threadpool = malloc(sizeof(threadpool_resampling_distance_t));
    assert(threadpool);
    threadpool->pfile = pfile;
    
    threadpool->save_patterns = save_patterns;
    threadpool->output = (char *)output;
    
    threadpool->count = 0;
    threadpool->total = count;
    threadpool->sequences = sequences;
    threadpool->model = model;
    threadpool->algorithm = algorithm;
    threadpool->scheme = scheme;
    threadpool->threads = malloc(nthreads*sizeof(pthread_t));
    assert(threadpool->threads);
    
    pthread_mutex_init(&(threadpool->lock), NULL);
    
    for ( int i = 0; i < nthreads; i++ ) {
        pthread_create( &(threadpool->threads[i]), NULL, _resampling_thread_worker_distance, threadpool );
    }
    for ( int i = 0; i < nthreads; i++ ) {
        pthread_join(threadpool->threads[i], NULL);
    }
    
    free(threadpool->threads);
    pthread_mutex_destroy(&(threadpool->lock));
    free(threadpool);
    fclose(pfile);
    
#ifdef FEEDBACK
    
#else
    fprintf(stderr, "\n");
#endif
}

void Distance_bootstrap_threads( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads ){
    Distance_resampling_threads(sequences, RESAMPLING_BOOTSTRAP, model, algorithm, count, output, save_patterns, nthreads );
}

void Distance_jackknife_threads( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, const char *output, bool save_patterns, int nthreads ){
    Distance_resampling_threads(sequences, RESAMPLING_JACKKNIFE, model, algorithm, sequences->length, output, save_patterns, nthreads );
}

void Distance_jackknife_n_threads( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int replicates, const char *output, bool save_patterns, int nthreads ){
    Distance_resampling_threads(sequences, RESAMPLING_JACKKNIFE_PROPORTION, model, algorithm, replicates, output, save_patterns, nthreads );
}

#endif

void Distance_resampling_openmp( const Sequences *sequences, resampling_scheme scheme, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads ){
#ifdef FEEDBACK
    printf("#=%d",count);
#endif
    
    nthreads = imin(nthreads, count);
    
    FILE *pfile = fopen(output,"w");
    
    int rep = 1;
    
    #pragma omp parallel for num_threads(nthreads)
    for ( int i = 0; i < count; i++ ) {
        Sequences *resampled  = NULL;
        if ( scheme == RESAMPLING_BOOTSTRAP ) {
            resampled = Sequences_bootstrap(sequences);
        }
        else if( scheme == RESAMPLING_JACKKNIFE){
            resampled = Sequences_jackknife(sequences, i);
        }
        else if( scheme == RESAMPLING_JACKKNIFE_PROPORTION){
            resampled = Sequences_jackknife(sequences, sequences->length/2);
        }
        else {
            error("Resampling scheme not recognized\n");
        }
        
        double **matrix = Sequences_distance(resampled, model);
        Tree *tree = NULL;
        if ( algorithm == CLUSTERING_NJ ) {
            tree = new_NJ(sequences, matrix);
        }
        else if( algorithm == CLUSTERING_UPGMA){
            tree = new_UPGMA(resampled, matrix);
        }
        else {
            error("clustering algorithm not yet implemented\n");
        }
        free_dmatrix(matrix, sequences->size);
        
        #pragma omp critical
        {
            
            Tree_print_newick(pfile, tree, false);
            //fprintf(pfile, "\n");
#ifdef FEEDBACK
            printf("#%d\n", rep);
#else
            fprintf(stderr, "Replicate %d/%d\r", rep++, count);
#endif
        }
        if ( save_patterns ) {
            StringBuffer *buffer = new_StringBuffer(10);
            StringBuffer_set_string(buffer, output);
            StringBuffer_append_format(buffer, ".%d", i);
            StringBuffer_append_string(buffer, ".fasta");
            Sequences_save_fasta(resampled, buffer->c);
            free_StringBuffer(buffer);
        }
        
        free_Sequences(resampled);
        free_Tree(tree);
#ifdef FEEDBACK
        printf("#%d\n", (i+1));
#endif
    }
#ifdef FEEDBACK

#else
    fprintf(stderr, "\n");
#endif
    fclose(pfile);
}


void Distance_bootstrap_openmp( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads ){
    Distance_resampling_openmp(sequences, RESAMPLING_BOOTSTRAP, model, algorithm, count, output, save_patterns, nthreads );
}

void Distance_jackknife_openmp( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, const char *output, bool save_patterns, int nthreads ){
    Distance_resampling_openmp(sequences, RESAMPLING_JACKKNIFE, model, algorithm, sequences->length, output, save_patterns, nthreads );
}


void Distance_jackknife_n_openmp( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int replicates, const char *output, bool save_patterns, int nthreads ){
    Distance_resampling_openmp(sequences, RESAMPLING_JACKKNIFE_PROPORTION, model, algorithm, replicates, output, save_patterns, nthreads );
}

void Distance_bootstrap( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int count, const char *output, bool save_patterns, int nthreads ){
#if defined (PTHREAD_ENABLED)
    if( nthreads == 1 ){
        Distance_bootstrap_openmp(sequences, model, algorithm, count, output, save_patterns, nthreads);
    }
    else {
        printf("Using PThread\n");
        Distance_bootstrap_threads(sequences, model, algorithm, count, output, save_patterns, nthreads);
    }
#else
#ifdef _OPENMP
    printf("Using OpenMP\n");
#endif
    Distance_bootstrap_openmp(sequences, model, algorithm, count, output, save_patterns, nthreads);
#endif
}

void Distance_jackknife( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, const char *output, bool save_patterns, int nthreads ){
#if defined (PTHREAD_ENABLED)
    if( nthreads == 1 ){
        Distance_jackknife_openmp(sequences, model, algorithm, output, save_patterns, nthreads);
    }
    else {
        printf("Using PThread\n");
        Distance_jackknife_threads(sequences, model, algorithm, output, save_patterns, nthreads);
    }
#else
#ifdef _OPENMP
    printf("Using OpenMP\n");
#endif
    Distance_jackknife_openmp(sequences, model, algorithm, output, save_patterns, nthreads);
#endif
}

void Distance_jackknife_n( const Sequences *sequences, distancematrix_model model, clustering_algorithm algorithm, int replicates, const char *output, bool save_patterns, int nthreads ){
#if defined (PTHREAD_ENABLED)
    if( nthreads == 1 ){
        Distance_jackknife_n_openmp(sequences, model, algorithm, replicates, output, save_patterns, nthreads);
    }
    else {
        printf("Using PThread\n");
        Distance_jackknife_n_threads(sequences, model, algorithm, replicates, output, save_patterns, nthreads);
    }
#else
#ifdef _OPENMP
    printf("Using OpenMP\n");
#endif
    Distance_jackknife_n_openmp(sequences, model, algorithm, replicates, output, save_patterns, nthreads);
#endif
}


