/*
 *  physher.c
 *  physher
 *
 *  Created by Mathieu Fourment on 11/10/10.
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

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h> // for sleep
#include <time.h>
#include <ctype.h>
#include <strings.h>

#include "phyc/substmodel.h"
#include "phyc/jc69.h"
#include "phyc/K80.h"
#include "phyc/f81.h"
#include "phyc/hky.h"
#include "phyc/gtr.h"
#include "phyc/gy94.h"
#include "phyc/mg94.h"
#include "phyc/nucsubst.h"
#include "phyc/gensubst.h"
#include "phyc/unrest.h"
#include "phyc/nonstat.h"
#include "phyc/dayhoff.h"
#include "phyc/lg.h"
#include "phyc/wag.h"

#include "phyc/matrix.h"
#include "phyc/sequence.h"
#include "phyc/sequenceio.h"
#include "phyc/sitepattern.h"
#include "phyc/tree.h"
#include "phyc/node.h"
#include "phyc/treeio.h"
#include "phyc/branchmodel.h"
#include "phyc/treelikelihood.h"
#include "phyc/optimize.h"
#include "phyc/mstring.h"
#include "phyc/localclock.h"
#include "phyc/discreteclock.h"
#include "phyc/heterotachy.h"
#include "phyc/optimizer.h"
#include "phyc/mathconstant.h"
#include "phyc/geneticcode.h"
#include "phyc/roottotip.h"
#include "phyc/args.h"
#include "phyc/random.h"
#include "phyc/descriptivestats.h"

#include "phyc/phyboot.h"
#include "phyc/phyci.h"
#include "phyc/boot.h"

#include "phyc/datatype.h"

#include "phyc/nj.h"
#include "phyc/upgma.h"
#include "phyc/distancematrix.h"
#include "phyc/phyclustering.h"

#include "phyc/topologyopt.h"

#include "phyc/treelikelihoodX.h"
#include "phyc/treelikelihood4.h"

#include "phyc/parsimony.h"

#include "phyc/solve.h"


#include "phyc/asr.h"
#include "phyc/qsearch.h"

#ifndef DISABLED_CONFIG_HEADER
#include "phyc/PhyCConfig.h"
#endif




// should only be used when the tree is dated and dates are forward in time
double get_earliest_date( Tree *tree ){
    double earliest_date = 0;
    int i = 0;
    if ( Tree_dated(tree) ) {
        Node **nodes = Tree_get_nodes(tree, POSTORDER);
        for ( ; i < Tree_node_count(tree); i++ ) {
            if ( Node_isleaf(nodes[i])) {
                earliest_date = dmax(earliest_date, nodes[i]->time);
            }
        }
    }
    return earliest_date;
}


// Return the rate
double init_times( SingleTreeLikelihood *tlk, bool forward, double rate_guess_user, bool unRooted ){
    Tree *tree = tlk->tree;
    
    bool isDated = parse_dates( tree );
    Tree_set_dated(tree, isDated);
    
    if ( !isDated) {
        fprintf(stdout, "Isochronous sequences\n");
    }
    
    double rate_guess = 0.03;
    
    
    if ( rate_guess_user > 0 ) {
        rate_guess = rate_guess_user;
        if(rate_guess > 0 ) fprintf(stdout, "User rate: %f\n", rate_guess);
    }
    
    
    // if there is no starting tree then it is an NJ tree
    
    if ( unRooted ) {
        double **regres = max_lm_tree(tree, forward, false);
        Node **nodes = Tree_nodes(tree);
        int count = 0;
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            if( Node_isroot(nodes[i]) ) continue;
            printf("%s slope: %f r: %f residuals: %e intercept: %f\n", Node_name(nodes[(int)regres[count][0]]), regres[count][1], regres[count][10], regres[count][12], regres[count][7]);
            count++;
        }
        rate_guess = regres[0][1]; // use the slope of the top regression
        Tree_reroot(tree, nodes[(int)regres[0][0]]);
        free(regres);
    }
    
    if ( !Tree_dated(tree) ) rate_guess = -1;
    
    // use linear regression or the calibrations to guess the mean rate
    if ( !unRooted && rate_guess <= 0 ) {
        if ( Tree_dated(tree) ) {
            
            //lm_tree_cluster(tree, forward, 2);
            
            double *lm = lm_tree(tree, forward, false);
            rate_guess = lm[1];
            
            if ( rate_guess <= 0 ) {
                fprintf(stderr, "Rate approximation using linear regression failed: %f\n", rate_guess);
                fprintf(stderr, "The sequences do no appear to be clock like!\n");
                fprintf(stderr, "or something is wrong with the dates provided\n");
                fprintf(stderr, "Provide a rate that you think is appropriate\n\n");
                exit(0);
            }
            else {
                fprintf(stdout, "Rate approximation using linear regression: %f\n", rate_guess);
                fprintf(stdout, "Correlation coefficient: %f RSS: %e Intercept: %f\n", lm[10], lm[12], lm[7]);
            }
            free(lm);
        }
        else {
            if ( !parse_calibrations(tree) ){
                error("Homochronous tree has no calibration points\n");
            }
            rate_guess = mean_rate_calibrations(tree);
            fprintf(stdout, "Rate approximation using calibrations: %f\n", rate_guess);
        }
    }
    
    
    if ( Tree_dated(tree) ) {
        Tree_init_heights_heterochronous(tree, rate_guess, forward);
    }
    else {
        Tree_init_heights_homochronous(tree, rate_guess);
        //summarize_calibrations(tree);
    }
    return rate_guess;
}

void calculate_posteriors_sites( SingleTreeLikelihood *tlk, const char *filename ){
    int i,s,c,best_index;
    double conditionalMean,best_prob;
	FILE *file = fopen(filename, "w");
    
	if( tlk->sm->shape != NULL ) fprintf(file, "Alpha %e\n", Parameter_value(tlk->sm->shape));
    if( tlk->sm->pinv  != NULL ) fprintf(file, "Pinv %e\n", Parameter_value(tlk->sm->pinv));
    
    fprintf(file, "Proportions\n");
	for ( i = 0; i < tlk->sm->cat_count; i++ ) {
		fprintf(file, "%f ", tlk->sm->get_proportion(tlk->sm, i ));
	}
	fprintf(file, "\nRates\n");
	for ( i = 0; i < tlk->sm->cat_count; i++ ) {
		fprintf(file, "%e ",tlk->sm->get_rate(tlk->sm, i));
	}
    fprintf(file, "\n");
	
	
	double **posteriors = SingleTreeLikelihood_posterior_sites(tlk);
	for ( s = 0; s < tlk->sp->nsites; s++ ) {
		i = tlk->sp->indexes[s];
		fprintf(file, "%d", s);
		conditionalMean = 0.0;
        
        best_index = 0;
        best_prob = posteriors[i][0];
		for ( c = 0; c < tlk->sm->cat_count; c++ ) {
            if( posteriors[i][c] > best_prob ){
                best_index = c;
                best_prob = posteriors[i][c];
            }
			fprintf(file, ",%e", posteriors[i][c]);
			conditionalMean += posteriors[i][c] * tlk->sm->get_rate(tlk->sm, c);
		}
		fprintf(file, ",%e,%d\n", conditionalMean, best_index);
	}
    
    fprintf(file, "\n");
    
    free_dmatrix(posteriors, tlk->sp->count);
    fclose(file);
}


void init_likelihood_approximation( SingleTreeLikelihood *tlk, treelikelihood_approximation approximation ){
    Tree *tree = tlk->tree;
    
    int dim = Tree_node_count(tree);
    if( approximation == TREELIKELIHOOD_APPROXIMATION_HESSIAN_DIAGONAL ){
        fprintf(stdout, "\nCalculating diagonal Hessian\n");
        tlk->hessian = dvector(dim);
        tlk->hessian_length = dim;
        SingleTreeLikelihood_Hessian_diag(tlk, tlk->hessian, &tlk->hessian_length);
//        int verbo = tlk->opt.verbosity;
//        tlk->opt.verbosity = 0;
//        while (!success) {
//            tlk->opt.bl.tolx = 0.000001;// needed to get negative second derivative
//            lnl = optimize_singletreelikelihood(tlk);
//            success = SingleTreeLikelihood_Hessian_diag(tlk, tlk->hessian, &tlk->hessian_length);
//        }
//        tlk->opt.verbosity = verbo;
        
        Node **nodes = Tree_get_nodes(tree, POSTORDER);
        int count = 0;
        for( int i = 0; i < dim; i++ ){
            //printf("%s %f %f\n", Node_name(nodes[i]), Node_distance(nodes[i]), -1.0/tlk->hessian[i]);
            if( Node_distance(nodes[i]) > 1e-7 ){
                printf("%s %e %e %s\n", Node_name(nodes[i]), Node_distance(nodes[i]), tlk->hessian[count], (tlk->hessian[count] >0?"*":""));
                count++;
            }
            else{
                printf("%s %e ?\n", Node_name(nodes[i]), Node_distance(nodes[i]) );
            }
        }
    }
    else if( approximation == TREELIKELIHOOD_APPROXIMATION_HESSIAN ){
        fprintf(stdout, "\nCalculating Hessian\n");
        tlk->hessian = dvector(dim*dim);
        tlk->hessian_length = dim*dim;
        bool success = SingleTreeLikelihood_Hessian(tlk, tlk->hessian, NULL);
//        while (!success) {
//            tlk->opt.bl.tolx = 0.000001;// needed to get negative second derivative
//            optimize_singletreelikelihood(tlk);
//            success = SingleTreeLikelihood_Hessian(tlk, tlk->hessian, NULL);
//        }
        
        if(true){
            Node **nodes = Tree_get_nodes(tree, POSTORDER);
            int count = 0;
            for( int i = 0; i < dim; i++ ){
                for( int j = i; j < dim; j++ ){
                    printf("%s %s %e %e %e\n", Node_name(nodes[i]), Node_name(nodes[j]), Node_distance(nodes[i]), Node_distance(nodes[j]),tlk->hessian[i*dim+j] );
                    count++;
                }
            }
        }
    }
    else if( approximation == TREELIKELIHOOD_APPROXIMATION_DERIVATIVE_DIAGONAL){
        fprintf(stdout, "\nCalculating derivatives\n");
        tlk->hessian = dvector(dim);
        tlk->hessian_length = dim;
        SingleTreeLikelihood_update_all_nodes(tlk);
        tlk->lnl_bl = tlk->calculate(tlk);
        tlk->sm->m->need_update = true;
        calculate_all_d2lnl_d2t(tlk, tlk->hessian);
        
        int a = 0;
        double *d2lnl2 = dvector(Tree_node_count(tree));
        SingleTreeLikelihood_Hessian_diag(tlk, d2lnl2, &a);
        
        Node **nodes = Tree_get_nodes(tree, POSTORDER);
        for( int i = 0; i < dim; i++ ){
            //tlk->hessian[i] = -1/tlk->hessian[i];
            //printf("%s %e %e\n", Node_name(nodes[i]), Node_distance(nodes[i]), tlk->hessian[i]);
            printf("%s %e %e %e %s\n", Node_name(nodes[i]), Node_distance(nodes[i]), tlk->hessian[i], d2lnl2[i], (tlk->hessian[i] >= 0?"*":""));
            //printf("%s %e %e %s\n", Node_name(nodes[i]), Node_distance(nodes[i]), tlk->hessian[i], (tlk->hessian[i] >0?"*":""));
        }
        free(d2lnl2);
    }
    
    fprintf(stdout, " \n\n");
}

// same topology if clock
void run_bootstrap( SingleTreeLikelihood *tlk, const char *output_stem, int nthreads_bootstrap, bool bootstrap, int samples, bool bca ){
    tlk->opt.topology_threads = 1;
    
    StringBuffer *buffer  = new_StringBuffer(10);
    Node **nodes = NULL;
    Tree *tree = tlk->tree;
    double ci_alpha = 0.05; // should be an option
    bool save_patterns = true; // should be an option
    char *jackfile_trees = NULL;
    char *jackfile_params = NULL;
    
    tlk->opt.verbosity = 0;
    
    if( tlk->bm != NULL ){
        nodes = Tree_get_nodes(tree, POSTORDER);
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            Node_empty_annotation(nodes[i]);
            
            if( !Node_isroot(nodes[i]) ){
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm, nodes[i]));
                Node_set_annotation(nodes[i], "rate",  buffer->c);
            }
        }
    }
    
    StringBuffer_set_string(buffer, output_stem);
    if( tlk->bm == NULL ){
        StringBuffer_append_string(buffer, ".freerate");
    }
    else {
        switch ( tlk->bm->name ) {
            case CLOCK_STRICT:{
                StringBuffer_append_string(buffer, ".strict");
                break;
            }
            case CLOCK_DISCRETE:{
                StringBuffer_append_string(buffer, ".discrete");
                break;
            }
            case CLOCK_LOCAL:{
                StringBuffer_append_string(buffer, ".local");
                break;
            }
    //        case CLOCK_RELAXED:{
    //            StringBuffer_append_string(buffer, ".strict.boot.trees");
    //            break;
    //        }
            default:
                break;
        }
    }

    
    if( bootstrap ){
        
        StringBuffer_append_string(buffer, ".boot");
    
        fprintf(stdout, "\nRunning %d bootstrap replicates using %d threads\n", samples, nthreads_bootstrap);
        SingleTreeLikelihood_bootstrap(tlk, samples, buffer->c, save_patterns, nthreads_bootstrap);
        
        
        if( bca ){
            StringBuffer *buffer2 = new_StringBuffer(10);
            StringBuffer_set_string(buffer2, buffer->c);
            StringBuffer_append_string(buffer2, ".jackknife");
            
            SingleTreeLikelihood_jackknife(tlk, buffer2->c, false, nthreads_bootstrap);
            
            StringBuffer_append_string(buffer2, ".trees");
            jackfile_trees = StringBuffer_tochar(buffer2);
            
            StringBuffer_set_string(buffer2, buffer->c);
            StringBuffer_append_string(buffer2, ".params");
            jackfile_params = StringBuffer_tochar(buffer2);
            
            free_StringBuffer(buffer2);
        }
    }
    else {
        StringBuffer_append_string(buffer, ".jackknife");
        
        fprintf(stdout, "\nRunning jackknife using %d threads\n", nthreads_bootstrap);
        SingleTreeLikelihood_jackknife(tlk, buffer->c, save_patterns, nthreads_bootstrap);
    }
    
    char *stem = StringBuffer_tochar(buffer);
    
    StringBuffer_append_string(buffer, ".trees");
    
    if( tlk->opt.topology_optimize ){
        Phyboot_annotate(tlk->tree, buffer->c);
    }
    else {
        Phyboot_annotate_fixed_topology( tree, buffer->c, 1-ci_alpha, jackfile_trees );
    }
    
    StringBuffer_chop(buffer);
    
    FILE *pfile_tree = fopen(buffer->c,"w");
    assert(pfile_tree);
    Tree_print_nexus_header_figtree_Taxa(pfile_tree, tree);
    Tree_print_nexus_header_figtree_BeginTrees(pfile_tree, tree);
    fprintf(pfile_tree, "tree TREE0 = [&R] ");
    Tree_print_nexus_with_annotation2( pfile_tree, tree, tlk->bm!=NULL );
    fprintf(pfile_tree, "\nEnd;");
    fclose(pfile_tree);
    
    
    StringBuffer_set_string(buffer, stem);
    StringBuffer_append_string(buffer, ".params");
    
    if( file_exists(buffer->c) ){
        Phyboot_confidence_intervals(buffer->c, 1-ci_alpha, jackfile_params);
    }
    
    
    free_StringBuffer(buffer);
    if(jackfile_params != NULL )free(jackfile_params);
    if(jackfile_trees != NULL )free(jackfile_trees);
    free(stem);
}

void run_custom_likelihood( SingleTreeLikelihood *tlk, const char *output_stem, int nthreads, int bootstrap, int jackknife, bool bca, bool forward ){
    Tree *tree = tlk->tree;
    StringBuffer *buffer = NULL;
    int i = 0;
    double lk;
    bool local_clock =false ,discrete_clock = false;
    Node **nodes = NULL;
    FILE *customFile = NULL;
    BranchModel *bm_custom = NULL;
    
    buffer = new_StringBuffer(10);
    StringBuffer_set_string(buffer, output_stem);
	StringBuffer_append_string(buffer, ".custom.tree");
    
    customFile = fopen(buffer->c,"w");
    assert(customFile);
    
    
    nodes = Tree_get_nodes(tree, POSTORDER);

    for( i = 0; i < Tree_node_count(tree); i++ ){
        if( String_contains_str(nodes[i]->info, "local") ){
            local_clock = true;
            break;
        }
        else if( String_contains_str(nodes[i]->info, "class") ){
            discrete_clock = true;
            break;
        }
    }
    
    if( local_clock ){
        printf("Local clock\n\n");
        bm_custom = new_LocalClock_from_tree(tree);
    }
    else if( discrete_clock ){
        printf("Discrete clock\n\n");
        bm_custom = new_DiscreteClock_from_tree(tree);
    }
    else {
        error("Could not read the clock type in the input tree\n");
    }
    
    // take the rate from the previous branchmodel
    for ( i = 0; i < Parameters_count(tlk->bm->rates); i++ ) {
        Parameters_set_all_value(bm_custom->rates, Parameters_value(tlk->bm->rates, 0));
    }
    
    SingleTreeLikelihood_set_BranchModel(tlk, bm_custom, false);
    
    lk = optimize_singletreelikelihood(tlk);
    
    fprintf(stdout, "\nLnL = %f\n", lk );
    fprintf(stdout, "Root height = %f", Node_height(Tree_root(tree)));
    if ( forward && Tree_dated(tree) ) {
        fprintf(stdout, " (%f)", (get_earliest_date(tree) - Node_height(Tree_root(tree))) );
    }
    fprintf(stdout, "\n");
    
    
    
    Tree_print_nexus_header_figtree_Taxa(customFile, tlk->tree);
    Tree_print_nexus_header_figtree_BeginTrees(customFile, tlk->tree);
    
    fprintf(customFile, "tree TREE0 [&lk=%f] = [&R] ", lk);
    
    for ( i = 0; i < Tree_node_count(tree); i++ ) {
        Node_empty_annotation(nodes[i]);
        if( !Node_isroot(nodes[i]) ){
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", bm_custom->get(bm_custom, nodes[i]));
            Node_set_annotation(nodes[i], "rate", buffer->c);
        }
    }
    Tree_print_nexus_with_annotation(customFile, tlk->tree);
    fprintf(customFile, "\nEnd;");
    fclose(customFile);
    
    if ( bootstrap > 0 || jackknife ) {
        run_bootstrap(tlk, output_stem, nthreads, bootstrap>0, (bootstrap>0? bootstrap: jackknife), bca);
    }
    
    free_StringBuffer(buffer);
}

// The branchmodel should be a strict clock and tlk should be optimized
void run_local_greedy_likelihood( SingleTreeLikelihood *tlk, const char* output_stem, const char* ic, int nThreads, bool forward, int sample_size, int verbosity ){
    BranchModel *bm_local = NULL;
    Tree *tree = tlk->tree;
    double rate = Parameters_value(tlk->bm->rates, 0);
    StringBuffer *buffer = new_StringBuffer(10);
    double clock_strict_lnl = tlk->lk;
    
    double clock_local_greedy_lnl;
	
    bm_local = new_LocalClock( tree, 1 );
    Parameters_set_value(bm_local->rates, 0, rate);
    Parameters_set_value(bm_local->rates, 1, rate);
    
    SingleTreeLikelihood_set_BranchModel(tlk, bm_local, false);
    
    ClockSearch *greedy_local = new_LocalClockSearch( tlk, nThreads );
    fprintf(stdout, "Threads: %d\n", greedy_local->nThreads );
    fflush(stdout);
    ClockSearch_set_starting_lnl(greedy_local, clock_strict_lnl);
    StringBuffer_append_strings(buffer, 2, output_stem, ".greedy.trees");
    ClockSearch_set_logfile_name(greedy_local, buffer->c);
    
    greedy_local->verbosity = verbosity;
    
    
    if( ic != NULL ){
        if( strcasecmp(ic, "AIC") == 0 ){
            ClockSearch_set_ic(greedy_local, INFORMATION_CRITERION_AIC);
        }
        else if( strcasecmp(ic, "AICc") == 0 ){
            ClockSearch_set_ic(greedy_local, INFORMATION_CRITERION_AICc);
        }
        else if( strcasecmp(ic, "BIC") == 0 ){
            ClockSearch_set_ic(greedy_local, INFORMATION_CRITERION_BIC);
        }
        else {
            fprintf(stderr, "IC not supported %s\n", ic);
            exit(1);
        }
    }
    fprintf(stdout, "Information criterion: %s\n", INFORMATION_CRITERION[greedy_local->ic] );
    
    
    ClockSearch_set_starting_df(greedy_local, SingleTreeLikelihood_df_count(tlk)-1); // -1 if the previous clock was strict
    
    if ( sample_size > 0 ) {
        ClockSearch_set_ic_sample_size(greedy_local, sample_size );
    }
    
    
    clock_local_greedy_lnl = greedy_local->optimize(greedy_local);
    
    int clock_local_greedy_n_clock = BranchModel_n_rate(bm_local)-1;
    
    if ( clock_local_greedy_n_clock >= 1 ) {
        
        double   *clock_local_greedy_rates = dvector( clock_local_greedy_n_clock+1 );
        BranchModel_rates_to_vector(bm_local, clock_local_greedy_rates);
        
        double max = dmax_vector(clock_local_greedy_rates, clock_local_greedy_n_clock+1);
        double min = dmin_vector(clock_local_greedy_rates, clock_local_greedy_n_clock+1);
        double meanRate_greedy = dmean(clock_local_greedy_rates, clock_local_greedy_n_clock+1);
        double meanRateScaled_greedy = BranchModel_mean_rate_scaled(bm_local);
        
        
        fprintf(stdout, "\nNumber of local clocks: %d\n", clock_local_greedy_n_clock );
        fprintf(stdout, "LnL = %f\n", clock_local_greedy_lnl);
        fprintf(stdout, "Rate\n");
        fprintf(stdout, "Median = %f Min = %f Max = %f\n", dmedian(clock_local_greedy_rates, clock_local_greedy_n_clock+1), min, max );
        fprintf(stdout, "Mean   = %f\n", meanRate_greedy);
        fprintf(stdout, "Mean   = %f (weighted average of branch-specific rates, i.e. sum r_i*t_i /sum t_i)\n", meanRateScaled_greedy );
        fprintf(stdout, "Root height = %f", Node_height(Tree_root(tree)));
        if ( forward && Tree_dated(tree) ) {
            fprintf(stdout, " (%f)", (get_earliest_date(tree) - Node_height(Tree_root(tree))) );
        }
        fprintf(stdout, "\n");
        
        double aicc = AICc(clock_local_greedy_lnl, SingleTreeLikelihood_df_count(tlk), tlk->sp->nsites);
        double aic  = AIC(clock_local_greedy_lnl, SingleTreeLikelihood_df_count(tlk));
        double bic  = BIC(clock_local_greedy_lnl, SingleTreeLikelihood_df_count(tlk), tlk->sp->nsites);
        
        fprintf(stdout, "AIC:  %f df: %d\n", aic, SingleTreeLikelihood_df_count(tlk));
        fprintf(stdout, "AICc: %f df: %d\n", aicc, SingleTreeLikelihood_df_count(tlk));
        fprintf(stdout, "BIC:  %f df: %d\n\n", bic, SingleTreeLikelihood_df_count(tlk));

        
        free(clock_local_greedy_rates);
    }
    else {
        fprintf(stdout, "\nThe strict clock model is significantly bettern than the local clock models\n" );
        fprintf(stdout, "LnL = %f\n", clock_local_greedy_lnl);
    }
    
    greedy_local->free(greedy_local);
    
    StringBuffer_empty(buffer);
    StringBuffer_append_strings(buffer, 2, output_stem, ".greedy.tree");
    FILE *greedyBestFile = fopen(buffer->c,"w");
    assert(greedyBestFile);
    Tree_print_nexus_header_figtree_Taxa(greedyBestFile, tree);
    Tree_print_nexus_header_figtree_BeginTrees(greedyBestFile, tree);
    
    fprintf(greedyBestFile, "tree TREE0 [&LnL=%f] = [&R] ", clock_local_greedy_lnl);
    
    Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
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
    Tree_print_nexus_with_annotation(greedyBestFile, tlk->tree);
    fprintf(greedyBestFile, "\nEnd;");
    fclose(greedyBestFile);
}

void run_searchq( SingleTreeLikelihood *tlk ,  const char* output_stem, int popSize, int ngen, int max_no_improvement, int nThreads ){
    QSearch *search_q = NULL;
    GAQSearch *ga_q = NULL;
    double lnl;
    StringBuffer *buffer = new_StringBuffer(10);
    
    // wont work with non reversible
    if( Parameters_count(tlk->sm->m->rates) == 1 ){
        Parameters_add(tlk->sm->m->rates, clone_Parameter(Parameters_at(tlk->sm->m->rates, 0), true));
        for (int i = 0; i <  tlk->sm->nstate*(tlk->sm->nstate-1)/2; i+=2 ) {
            tlk->sm->m->model[i] = 1;
        }
    }
//    else {
//        int index = 0;
//        // Find a parameter that is not fixed
//        for ( ; index < Parameters_count(tlk->sm->m->rates); index++ ) {
//            if( Parameters_fixed(tlk->sm->m->rates, index)) break;
//        }
//        for ( int i = 0; i < Parameters_count(tlk->sm->m->rates); i++ ) {
//            if( Parameters_value(tlk->sm->m->rates, i) < 0.001){
//                Parameters_set_value(tlk->sm->m->rates, i, 0);
//                Parameters_set_fixed(tlk->sm->m->rates, true, i);
//                break;
//            }
//        }
//        
//    }
    tlk->opt.relative_rates.method = OPT_BRENT;
    search_q = new_GAQSearch(tlk, nThreads, popSize);
    
    search_q->mode = 0;
    
    StringBuffer_append_strings(buffer, 2, output_stem, ".qsearch.log");
    search_q->logFileName = String_clone( buffer->c);
    
    
    ga_q = (GAQSearch*)search_q->pDerivedObj;
    ga_q->ngen = ngen;
    ga_q->max_no_improvement = max_no_improvement;
    
    
    lnl = search_q->optimize(search_q);
    
    fprintf(stdout, "LnL = %f\n", lnl);
    
    search_q->free(search_q);
    free_StringBuffer(buffer);
}

void run_discrete_ga_likelihood( SingleTreeLikelihood *tlk , const char* output_stem, int popSize, int ngen, int max_no_improvement,
                                int nThreads, const char* ic, bool forward, int clock_categories, int approximation ){
    
    double previous_lnl  = tlk->lk;
    int    previous_df   = SingleTreeLikelihood_df_count(tlk);
    
    int nNodes;
    Tree *tree = tlk->tree;
    BranchModel *bm_discrete = NULL;
    ClockSearch *search_discrete = NULL;
    StringBuffer *buffer = NULL;
    
    buffer = new_StringBuffer(10);
    
    nNodes = Tree_node_count(tree);
    
    tlk->approx = approximation;
    
    if( tlk->bm->name == CLOCK_DISCRETE){
        bm_discrete = tlk->bm;
    }
    else if( tlk->bm->name == CLOCK_LOCAL){
        bm_discrete = new_DiscreteClock_from_LocalClock(tlk->bm);
        SingleTreeLikelihood_set_BranchModel(tlk, bm_discrete, false);
    }
    // it is a strict clock
    else {
        if ( clock_categories > 2 ) {
            bm_discrete = new_DiscreteClock(tree, clock_categories);
        }
        else {
            bm_discrete = new_DiscreteClock(tree, 2);
            
            Node **nodes = Tree_get_nodes(tree, POSTORDER);
            bool first = ( Node_distance(nodes[0]) >= Node_time_elapsed(nodes[0])*Parameters_value(tlk->bm->rates, 0) );
            for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
                int index = Node_id(nodes[i]);
                bm_discrete->map[index] = ( Node_distance(nodes[index]) >= Node_time_elapsed(nodes[index])*Parameters_value(tlk->bm->rates, 0) );
                bm_discrete->map[index] ^= first;
            }
            
            if( Node_id(Tree_root(tree)) == 0 ) bm_discrete->map[0] = bm_discrete->map[1];
            else bm_discrete->map[Node_id(Tree_root(tree))] = bm_discrete->map[Node_id(Tree_root(tree))-1];
            
            Parameters_set_value(bm_discrete->rates, 0, Parameters_value(tlk->bm->rates, 0));
            Parameters_set_value(bm_discrete->rates, 1, Parameters_value(tlk->bm->rates, 0));
            //FIXME: should check this
            //Parameters_set_lower(bm_discrete->rates, 0, Parameters_value(tlk->bm->rates, 0)*0.01);
            //Parameters_set_lower(bm_discrete->rates, 1, Parameters_value(tlk->bm->rates, 0)*0.01);
        }
        
        SingleTreeLikelihood_set_BranchModel(tlk, bm_discrete, false);
    }
    
    search_discrete = new_GADiscreteClockSearch(tlk, nThreads, popSize);
    
    ClockSearch_set_starting_lnl(search_discrete, previous_lnl);
    ClockSearch_set_starting_df(search_discrete, previous_df);
    
    
    //search_discrete->mode = 1;
    //ClockSearch_set_alpha(search_discrete, significance_level);
    
    StringBuffer_empty(buffer);
    StringBuffer_append_strings(buffer, 2, output_stem, ".ga.discrete.trees");
    ClockSearch_set_logfile_name(search_discrete, buffer->c);
    
    if( ic != NULL ){
        if( strcasecmp(ic, "AIC") == 0 ){
            ClockSearch_set_ic(search_discrete, INFORMATION_CRITERION_AIC);
        }
        else if( strcasecmp(ic, "AICc") == 0 ){
            ClockSearch_set_ic(search_discrete, INFORMATION_CRITERION_AICc);
        }
        else if( strcasecmp(ic, "BIC") == 0 ){
            ClockSearch_set_ic(search_discrete, INFORMATION_CRITERION_BIC);
        }
        else {
            fprintf(stderr, "IC not supported %s\n", ic);
            exit(1);
        }
    }
    
    GAClockSearch *ga_discrete = (GAClockSearch*)search_discrete->pDerivedObj;
    
    ga_discrete->ngen = ngen;
    ga_discrete->max_no_improvement = max_no_improvement;
    
    
    
    double lk = search_discrete->optimize(search_discrete);
    
    if(approximation != TREELIKELIHOOD_APPROXIMATION_NONE ){
        
        fprintf(stdout, "\nApproximated LnL = %f\n", lk );
        SingleTreeLikelihood_update_all_nodes(tlk);
        fprintf(stdout, "LnL recalculated = %f\n", tlk->calculate(tlk) );
        tlk->approx = TREELIKELIHOOD_APPROXIMATION_NONE;
        tlk->opt.verbosity = 1;
        lk = optimize_singletreelikelihood(tlk);
        tlk->approx = approximation;
    }
    
    int clock_discrete_n_rate = BranchModel_n_rate(bm_discrete);
    
    // check if we get a better a model
    
    unsigned *clock_discrete_indexes = clone_uivector( bm_discrete->map, nNodes );
    double   *clock_discrete_heights = dvector( nNodes );
    Tree_heights_to_vector(tree, clock_discrete_heights);
    double   *clock_discrete_rates   = dvector( clock_discrete_n_rate );
    BranchModel_rates_to_vector(bm_discrete, clock_discrete_rates);
    
    double max = dmax_vector(clock_discrete_rates, clock_discrete_n_rate);
    double min = dmin_vector(clock_discrete_rates, clock_discrete_n_rate);
    double meanRate_ga_discrete = dmean(clock_discrete_rates, clock_discrete_n_rate);
    double meanRateScaled_ga_discrete = BranchModel_mean_rate_scaled(bm_discrete);
    
    fprintf(stdout, "\nNumber of discrete rate classes clocks: %d\n", clock_discrete_n_rate );
    fprintf(stdout, "LnL = %f\n", lk);
    fprintf(stdout, "Mean Rate   = %f Min = %f Max = %f\n", meanRate_ga_discrete, min, max );
    fprintf(stdout, "Mean Rate   = %f (weighted average of branch-specific rates, i.e. sum r_i*t_i /sum t_i)\n", meanRateScaled_ga_discrete );
    fprintf(stdout, "Coefficient of correlation = %f\n", BranchModel_correlation(bm_discrete));
    fprintf(stdout, "Root height = %f", clock_discrete_heights[nNodes-1]);
    if ( forward && Tree_dated(tree) ) {
        fprintf(stdout, " (%f)", (get_earliest_date(tree) - clock_discrete_heights[nNodes-1]) );
    }
    fprintf(stdout, "\n");
    
    fprintf(stdout, "Tips     mean rate   = %f\n", BranchModel_mean_rate_tips_scaled(bm_discrete) );
    fprintf(stdout, "Internal mean rate   = %f\n\n", BranchModel_mean_rate_internal_scaled(bm_discrete) );
    
    fprintf(stdout, "Correlation with branch length %f\n\n", BranchModel_correlation_distance(bm_discrete) );
    
    SingleTreeLikelihood_update_all_nodes(tlk);
    fprintf(stdout, "recalculated %f\n", tlk->calculate(tlk));
    
    StringBuffer_empty(buffer);
    StringBuffer_append_strings(buffer, 2, output_stem, ".ga.discrete.tree");
    
    FILE *gaDiscreteBestFile = fopen( buffer->c,"w");
    assert(gaDiscreteBestFile);
    Tree_print_nexus_header_figtree_Taxa(gaDiscreteBestFile, tree);
    Tree_print_nexus_header_figtree_BeginTrees(gaDiscreteBestFile, tree);
    
    fprintf(gaDiscreteBestFile, "tree TREE0 [&LnL=%f] = [&R] ", lk);
    
    Node **nodes = Tree_nodes(tlk->tree);
    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        Node_empty_annotation(nodes[i]);
        if( !Node_isroot(nodes[i]) ){
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", tlk->bm->get(tlk->bm,nodes[i]));
            Node_set_annotation(nodes[i], "rate", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%d", tlk->bm->map[ Node_id(nodes[i]) ]);
            Node_set_annotation(nodes[i], "class", buffer->c);
        }
    }
    Tree_print_nexus_with_annotation(gaDiscreteBestFile, tlk->tree);
    
    fprintf(gaDiscreteBestFile, "\nEnd;");
    fclose(gaDiscreteBestFile);
    
    
    free(clock_discrete_indexes);
    free(clock_discrete_heights);
    free(clock_discrete_rates);
    
    
    search_discrete->free(search_discrete);
    free_StringBuffer(buffer);
	
}

// The branchmodel is already a LocalClock with already several local clocks
void run_local_ga_likelihood( SingleTreeLikelihood *tlk, const char* output_stem, int popSize, int ngen, int max_no_improvement,
                             int nThreads, const char* ic, bool forward ){
    Tree *tree = tlk->tree;
    time_t start_time, end_time, diff_time;
    double max,min, meanRate, meanRateScaled;
    double *rates = NULL;
    double lnl;
    FILE *outputFile = NULL;

    int nClock = 0;
    
    StringBuffer *buffer = NULL;
    
    ClockSearch *search_local = NULL;
    GAClockSearch *ga_local = NULL;
	   
    buffer = new_StringBuffer(10);
    
    StringBuffer_set_string(buffer, output_stem);
	StringBuffer_append_string(buffer, ".ga.local.tree");
    
    outputFile = fopen(buffer->c,"w");
    assert(outputFile);
    
    time(&start_time);
    
    BranchModel *bm_local = new_LocalClock(tlk->tree, 2);
    SingleTreeLikelihood_set_BranchModel(tlk, bm_local, false);
    
    search_local = new_GALocalClockSearch2(tlk, nThreads, popSize);
    ga_local = (GAClockSearch*)search_local->pDerivedObj;
    
    if( ic != NULL ){
        if( strcasecmp(ic, "AIC") == 0 ){
            ClockSearch_set_ic(search_local, INFORMATION_CRITERION_AIC);
        }
        else if( strcasecmp(ic, "AICc") == 0 ){
            ClockSearch_set_ic(search_local, INFORMATION_CRITERION_AICc);
        }
        else if( strcasecmp(ic, "BIC") == 0 ){
            ClockSearch_set_ic(search_local, INFORMATION_CRITERION_BIC);
        }
        else {
            fprintf(stderr, "IC not supported %s\n", ic);
            exit(1);
        }
    }
    
    StringBuffer_set_string(buffer, output_stem);
	StringBuffer_append_string(buffer, ".ga.local.trees");
    ClockSearch_set_logfile_name(search_local, buffer->c);
    
    ga_local->max_no_improvement = max_no_improvement;
    ga_local->ngen = ngen;
    
     lnl = search_local->optimize(search_local);
    
    nClock = BranchModel_n_rate(tlk->bm)-1;
    
    rates  = dvector( nClock+1 );
    BranchModel_rates_to_vector(tlk->bm, rates);
    
    max = dmax_vector(rates, nClock+1);
    min = dmin_vector(rates, nClock+1);
    meanRate = dmean(rates, nClock+1);
    meanRateScaled = BranchModel_mean_rate_scaled(tlk->bm);
    
    
    fprintf(stdout, "\nNumber of local clocks: %d\n", nClock );
    fprintf(stdout, "LnL = %f\n", lnl);
    fprintf(stdout, "Mean Rate   = %f Min = %f Max = %f\n", meanRate, min, max );
    fprintf(stdout, "Mean Rate   = %f (weighted average of branch-specific rates, i.e. sum r_i*t_i /sum t_i)\n", meanRateScaled );
    fprintf(stdout, "Root height = %f", Node_height(Tree_root(tree)));
    if ( forward && Tree_dated(tree) ) {
        fprintf(stdout, " (%f)", (get_earliest_date(tree) - Node_height(Tree_root(tree))) );
    }
    fprintf(stdout, "\n");
    
    free(rates);
    
    search_local->free(search_local);
    
    
    time(&end_time);
    diff_time = difftime(end_time, start_time);
    fprintf(stdout, "\nTotal runtime ");
    print_pretty_time(stdout, diff_time);
    
    SingleTreeLikelihood_update_all_nodes(tlk);
    fprintf(stdout, "recalculated %f\n", tlk->calculate(tlk));
    
    Tree_print_nexus_header_figtree_Taxa(outputFile, tree);
    Tree_print_nexus_header_figtree_BeginTrees(outputFile, tree);
    
    fprintf(outputFile, "tree TREE0 [&LnL=%f] = [&R] ", lnl);
    
    Node **nodes = Tree_nodes(tlk->tree);
    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
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
    Tree_print_nexus_with_annotation(outputFile, tlk->tree);
    fprintf(outputFile, "\nEnd;");
    fclose(outputFile);
    free_StringBuffer(buffer);
}

void append_simultron( SingleTreeLikelihood* tlk, const char* model_string, const char* outout_stem, StringBuffer *info ){
    SubstitutionModel *mod = tlk->sm->m;
    StringBuffer *buffer = new_StringBuffer(10);
    
    StringBuffer_set_string(buffer, outout_stem);
	StringBuffer_append_string(buffer, ".strict.tree");
    
    if( mod->dtype == DATA_TYPE_NUCLEOTIDE ){
        StringBuffer_append_string(info, "\n\n");
        StringBuffer_append_string(info, "simultron -m ");
        StringBuffer_append_string(info, model_string);
        
        if( mod->_freqs != NULL ){
            StringBuffer_append_format(info, " -f %f,%f,%f,%f", mod->_freqs[0], mod->_freqs[1], mod->_freqs[2], mod->_freqs[3] );
        }
        else {
            StringBuffer_append_string(info, " -f 0.25,0.25,0.25,0.25" );
        }
        if( Parameters_count(mod->rates) != 0 ){
            StringBuffer_append_string(info, " -r " );
            
            for ( int i = 0; i < Parameters_count(mod->rates); i++) {
                StringBuffer_append_format(info, "%f,", Parameters_value(mod->rates, i));
            }
            StringBuffer_chop(info);
        }
        if(tlk->sm->shape != NULL ){
            int cat = tlk->sm->cat_count;
            if ( tlk->sm->pinv != NULL ) {
                cat--;
            }
            StringBuffer_append_format(info, " -a %f -c %d",  Parameter_value(tlk->sm->shape), cat);
        }
        StringBuffer_append_format(info, " -l %d", tlk->sp->nsites);
        StringBuffer_append_format(info, " -s %f", Parameters_value( tlk->bm->rates, 0) );
        StringBuffer_append_format(info, " -i %s", buffer->c);
        StringBuffer_append_format(info, " -o %s", "sim.nex");
    }
    free_StringBuffer(buffer);
}

//gcc -Wall -O2 -std=c99 -o mldater test.c -L/Users/mathieu/CProjects/Math/build/Release/ -lMath -I/Users/mathieu/CProjects/Math/ -fopenmp

// valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-fds=yes --log-file=valgrind.txt --dsymutil=yes -v ./mldater -overwrite ~/Desktop/MLDate/GADate/2localclocks.config

/*************************************************************************************************
 *************************************** Parse config file ***************************************
 *************************************************************************************************/

Hashtable * parse_config_file( const char *configfile ){
	Hashtable *hash = new_Hashtable_string( 50 );
	FileReader *reader = new_FileReader(configfile, 1000);
	while ( reader->read_line(reader) ) {
				
		StringBuffer_trim(reader->buffer);
		if ( reader->buffer->length == 0 || String_start_with(reader->line, "#", true)) continue;
		char *p = reader->buffer->c;
		while ( *p != '\0') {
			if ( *p == '#') {
				*p = '\0';
				break;
			}
			p++;
		}
		reader->buffer->length = strlen(reader->buffer->c);
		StringBuffer_trim(reader->buffer);
		
		
		int l = 0;
		char **temp = String_split_char( reader->buffer->c, '=', &l );
		
		if ( l != 2 ) {
			fprintf(stderr, "Could not parse: %s . There should be only 1 '='(%d)\n", reader->line, l);
			continue;
		}
		char *key   = String_rtrim( temp[0] );
		char *value = String_trim( temp[1] );
		free(temp);
		temp = NULL;
		
        if ( Hashtable_exists(hash, key) ) {
			fprintf(stderr, "Option already set: %s\n", reader->buffer->c);
            exit(1);
		}
        
		Hashtable_add(hash, key, value);
		
	}
	free_FileReader(reader);
	return hash;
}


// Check if the names in the tree match the names in the sequence file
void check_aln_tree( Tree *tree, SitePattern *sp ){
	Node **nodes = Tree_get_nodes(tree, POSTORDER);

	if( Tree_tip_count(tree) != sp->size ) {
		fprintf(stdout, "Number of tips %d\n", Tree_tip_count(tree));
		fprintf(stdout, "Number of sequences %d\n\n", sp->size);
		fprintf(stdout, "The tree and the alignment don't have the same number of sequences\n\n");
		exit(1);
	}

	int flag = 0;
	for ( int i = 0; i < Tree_node_count(tree); i++ ) {
		if( !Node_isleaf(nodes[i]) ) continue;
		int j = 0;
		for ( ; j < sp->size; j++ ) {
			if( strcmp( Node_name(nodes[i]), sp->names[j] ) == 0 ){
				break;
			}
		}
		if( j == sp->size ){
			fprintf(stderr, "%s tree file not found in sequence file\n", Node_name(nodes[i]));
			
//			for ( j=0 ; j < sp->size; j++ ) {
//				fprintf(stderr, "  %s\n", sp->names[j]);
//			}
//			exit(0);
			
			flag++;
		}
	}
	if( flag > 0 ){
		fprintf(stderr, "Another program could have modified the name of the sequence\n");
		exit(1);
	}
}

void check_alignment( const Sequences *seqs){
    for ( int i = 0; i < seqs->size; i++ ) {
        char *name1 = seqs->seqs[i]->name;
        for ( int k = 0; k <strlen(name1); k++) {
            if( name1[k] == ',' || name1[k] == '(' || name1[k] == ')' || name1[k] == ':' || name1[k] == '[' || name1[k] == ']' || name1[k] == ';'){
                fprintf(stderr, "Illegal sequence name: %s\n", seqs->seqs[i]->name);
                fprintf(stderr, "Sequence name cannot contain the following characters ,():[];\n\n");
                exit(1);
                
            }
        }
        for ( int j = i+1; j < seqs->size; j++ ){
            if( strcmp(name1, seqs->seqs[j]->name) == 0 ){
                fprintf(stderr, "Duplicate sequence name: %s\n", seqs->seqs[i]->name);
                exit(1);
            }
        }
    }
    
    int length = seqs->seqs[0]->length;
    for ( int i = 1; i < seqs->size; i++ ) {
        if( seqs->seqs[i]->length != length ){
            fprintf(stderr, "Sequences of different length\n");
            exit(1);
        }
    }
}

void run_Distance( const char* seq_file, const char* output_stem, const char* model_string, int nthreads,
                  const char* method, int nexus_index, long seed, int bootstrap, int jackknife ){
    Sequences * seqs = NULL;
    if( nexus_index >= 0 ){
        seqs = readNexus(seq_file, nexus_index);
    }
    else {
        seqs = readSequences(seq_file);
    }
    
    check_alignment(seqs);
    
    fprintf(stdout, "Number of sequences: %d\n",seqs->size );
    fprintf(stdout, "Alignment length: %d\n",seqs->length );
    
    distancematrix_model model = DISTANCE_MATRIX_UNCORRECTED;
    if ( model_string != NULL ) {
        if( strcasecmp(model_string, "K2P") == 0 ){
            seqs->datatype = new_NucleotideDataType();
            model = DISTANCE_MATRIX_K2P;
            fprintf(stdout, "Model: K2P\n");
        }
        else if( strcasecmp(model_string, "JC69") == 0 ){
            seqs->datatype = new_NucleotideDataType();
            model = DISTANCE_MATRIX_JC69;
            fprintf(stdout, "Model: JC69\n");
        }
        else if( strcasecmp(model_string, "K83") == 0 ){
            seqs->datatype = new_AminoAcidDataType();
            model = DISTANCE_MATRIX_KIMURA;
            fprintf(stdout, "Model: K83 (Amino acid)\n");
        }
        else {
            seqs->datatype = new_NucleotideDataType();
            fprintf(stdout, "Model: raw\n");
        }
    }
    else {
        seqs->datatype = new_NucleotideDataType();
    }
    
    clustering_algorithm algorithm = CLUSTERING_NJ;
    
    StringBuffer *buffer = new_StringBuffer(10);
    
    StringBuffer_set_string(buffer, output_stem);
    
    char rooted = 'R';
    
    if( method != NULL  ){
        if( strcasecmp(method, "NJ") == 0 ){
            algorithm = CLUSTERING_NJ;
            StringBuffer_append_string(buffer, ".nj.tree");
        }
        else if( strcasecmp(method, "UPGMA") == 0 ){
            algorithm = CLUSTERING_UPGMA;
            StringBuffer_append_string(buffer, ".upgma.tree");
        }
        else {
            error("Cannot recognize algorithm\n");
        }
        fprintf(stdout, "Algorithm: %s\n", method);
    }
    time_t t0;
    time_t t1;
    time_t t2;
    double diff_time;
    
    
    time(&t0);
    
    Tree *tree = NULL;
    
    double double_precision = false;
    
    if(double_precision){
        fprintf(stdout, "Building double precision distance matrix... ");
        fflush(stdout);
        double **matrix = Sequences_distance(seqs, model);
        time(&t1);
        
        diff_time = difftime(t1, t0);
        print_pretty_time(stdout, diff_time);
        
        fprintf(stdout, "\n");
        fprintf(stdout, "Building tree... ");
        fflush(stdout);
        
        switch (algorithm) {
            case CLUSTERING_NJ:
                tree = new_NJ(seqs, matrix);
                rooted = 'U';
                break;
            case CLUSTERING_UPGMA:
                tree = new_UPGMA(seqs, matrix);
                break;
            default:
                assert(0);
        }
        
        free_dmatrix(matrix, seqs->size);
    }
    else {
        fprintf(stdout, "Building single precision distance matrix... ");
        fflush(stdout);
        float **matrix = Sequences_distance_float(seqs, model);
        time(&t1);
        
        diff_time = difftime(t1, t0);
        print_pretty_time(stdout, diff_time);
        
        fprintf(stdout, "\n");
        fprintf(stdout, "Building tree... ");
        fflush(stdout);
        
        switch (algorithm) {
            case CLUSTERING_NJ:
                tree = new_NJ_float(seqs, matrix);
                rooted = 'U';
                break;
            case CLUSTERING_UPGMA:
                tree = new_UPGMA_float(seqs, matrix);
                break;
            default:
                break;
        }
        free_fmatrix(matrix, seqs->size);
    }
    time(&t2);
    diff_time = difftime(t2, t1);
    print_pretty_time(stdout, diff_time);
    
    fprintf(stdout, "\n");
    fflush(stdout);
    
    FILE *pfile = fopen(buffer->c,"w");
    
    if( bootstrap > 0 || jackknife > 0 ){
        
        init_genrand(seed);
        
        fprintf(stdout, "Random seed: %lu\n\n", seed);
        
        bool save_patterns = false; // should be an option
        StringBuffer_set_string(buffer, output_stem);
        
        switch (algorithm) {
            case CLUSTERING_NJ:
                StringBuffer_append_string(buffer, ".nj.trees");
                break;
            case CLUSTERING_UPGMA:
                StringBuffer_append_string(buffer, ".upgma.trees");
                break;
            default:
                assert(0);
        }
        
        if( bootstrap > 0 ){
            Distance_bootstrap( seqs, model, algorithm, bootstrap, buffer->c, save_patterns, nthreads );
        }
        else {
            Distance_jackknife( seqs, model, algorithm, buffer->c, save_patterns, nthreads );
        }
        
        Phyboot_annotate(tree, buffer->c);
    }
    
    Tree_print_nexus_header_figtree_Taxa(pfile, tree);
    Tree_print_nexus_header_figtree_BeginTrees(pfile, tree);
    fprintf(pfile, "tree TREE0 = [&%c] ", rooted);
    Tree_print_nexus_with_annotation2(pfile, tree,false);
    fprintf(pfile, "\nEnd;");
    
    time(&t2);
    diff_time = difftime(t2, t0);
    fprintf(stdout, "\nTotal runtime ");
    print_pretty_time(stdout, diff_time);
    fprintf(stdout, "\n");
    
    fclose(pfile);
    
    free_StringBuffer(buffer);
    free_Tree(tree);
    free_Sequences(seqs);
}

bool array_of_string_contains(const char *str, const char *array[], int count){
    for ( int i = 0; i < count; i++) {
        if( strcasecmp(str, array[i]) == 0 ){
            return true;
        }
    }
    return false;
}

char* get_program_name(char* argv[]){
    char *name = argv[0]+strlen(argv[0]);
    while( name != argv[0] ){
        if ( *name == '/' || *name == '\\' ) {
            name++;
            break;
        }
        name--;
    }
    return name;
}

int main(int argc, char* argv[]){
    
    bool overwrite = false;
    char *seq_file  = NULL;
    char *tree_file = NULL;
    char *output_stem_user = NULL;
    
    char* clock = NULL;
    char* clock_algorithm = NULL;
    int nexus_index = -1;
    int genetic_code = 0;
    bool ambiguity = false;
    long seed = time(NULL);
    
    char* markov_states = NULL;
    char* model_string_user = NULL;
    char* rates_user = NULL;
    bool rates_fixed = false;
    char* frequencies_string_user = NULL;
    bool frequencies_unknown = false;
    bool frequencies_fixed = false;
    bool normalize_q = true;
    
    int rate_category_count = 1;
    double alpha = 0.5;
    bool alpha_fixed = false;
    bool use_pinv = false;
    double pinv = 0;
    bool pinv_fixed = false;
    double tree_scaler = -1;
    bool tree_unrooted = false;
    bool forward  = false;
    bool use_parsimony = true;
    char* topology_optimization_algorithm = NULL;
    
    char* qsearch = NULL;
    
    int ga_population_size = 30;
    int ga_generations = 500;
    int ga_max_no_improvement = 50;
    
    int clock_categories = 2;
    
    double clock_rate_guess = -1;
    char* ic = String_clone("AIC");
    int ic_sample_size = -1;
    
    bool use_double_distance = false;
    
    
    int approximation = TREELIKELIHOOD_APPROXIMATION_NONE;
    bool use_sse   = true;
    bool use_upper = true;
    int tlk_threads = 1; //for treelikleihood
    
    int nthreads =  1; // ga...
    
    int verbosity = 1;
    
    int bootstrap = 0;
    int jackknife = 0;
    bool bca = false;
    
    bool posterior_sites = false;
    bool asr = false;
    bool bl_fixed = false;
    
    char* distance = NULL;
    
    bool hessian = false;
    bool batch = false;
    
    bool help = false;
    bool help2 = false;
    
    char* fix = NULL;
    char* sitemodel_string = NULL;
    
    struct argsparser_option options[] = {
        {ARGS_OPTION_STRING,  'i', "sequences",    "input.sequences", &seq_file, "Input alignment file"},
        {ARGS_OPTION_STRING,  't', "tree",         "input.tree", &tree_file, "Input tree file"},
        {ARGS_OPTION_STRING,  'o', "stem",         "output.stem", &output_stem_user, "Output stem file"},
        
        {ARGS_OPTION_INTEGER, 0,   "gc",           "sequences.geneticcode", &genetic_code, "Genetic Code"},
        
        {ARGS_OPTION_STRING,  'm', "model",        "substmodel.type", &model_string_user, "Susbtitution model"},
        {ARGS_OPTION_STRING,  0,   "states",       "substmodel.states", &markov_states, "State space of Markov process. For nucleotide --states A,C,G,T"},
        {ARGS_OPTION_STRING,  'f', "frequencies",  "substmodel.freqs", &frequencies_string_user, "Frequencies as an array or 'e' for equal frequencies. For nucleotide -f 0.2,0.3,0.4,0.1"},
        {ARGS_OPTION_STRING,  'r', "rates",        "substmodel.rates", &rates_user, "Relative rates of the susbtitution matrix"},
        {ARGS_OPTION_BOOLEAN, 0,   "q-normalize",  "substmodel.matrix.normalize", &normalize_q, "Normalize rate matrix"},
        {ARGS_OPTION_FLAG,    0,   "f-unknown",    "substmodel.root.freqs.unknown", &frequencies_unknown, "Set frequencies to 1.0"},
        {ARGS_OPTION_STRING,  0,   "q-search",     "substmodel.qsearch", &qsearch, "Find best rate matrix (ga)"},
        
        {ARGS_OPTION_INTEGER, 'c', "cat",          "sitemodel.heterogeneity.gamma.cat", &rate_category_count, "Number of rate categories for gamma distribution"},
        {ARGS_OPTION_STRING,  'H', "het",          "sitemodel.heterogeneity.type", &sitemodel_string, "discrete or gammaquad"},
        {ARGS_OPTION_DOUBLE,  'a', "alpha",        "sitemodel.heterogeneity.gamma.alpha", &alpha, "Value of the alpha parameter of the gamma distribution"},
        {ARGS_OPTION_FLAG,    'I', "invariant",    "sitemodel.heterogeneity.pinv", &use_pinv, "Switch on a proportion of invariant sites"},
        {ARGS_OPTION_DOUBLE,  0,   "I-value",      "sitemodel.heterogeneity.pinv.value", &pinv, "Value of the proportion sites"},
        {ARGS_OPTION_FLAG,    0,   "ps",           "sitemodel.heterogeneity.posterior", &posterior_sites, "Caclulate posterior estimates of rates at each site"},
        
        {ARGS_OPTION_FLAG,    0,   "ambiguity",    "sequences.ambiguity", &ambiguity, "Use ambiguity for likelihood calculation"},
        
        {ARGS_OPTION_STRING,  'F', "fix",          "treelikelihood.fix", &fix, "Fix d: branch length, i: invariant, a: alpha, f: frequencies, r: rates"},
        {ARGS_OPTION_BOOLEAN, 0,   "sse",          "treelikelihood.sse", &use_sse, "Use SSE [default true]"},
        {ARGS_OPTION_BOOLEAN, 0,   "upper",        "treelikelihood.upper", &use_upper, "Use upper likelihood"},
        {ARGS_OPTION_INTEGER, 0,   "lk-threads",   "treelikelihood.threads", &tlk_threads, "Number of threads for likelihood calculation"},
        {ARGS_OPTION_INTEGER, 'A', "approx",       "treelikelihood.approximation", &approximation, "Input alignment file"},
        
        {ARGS_OPTION_FLAG,    'U', "unrooted",     "tree.unrooted", &tree_unrooted, "Input tree is rooted"},
        {ARGS_OPTION_DOUBLE,  's', "scaler",       "tree.scaler", &tree_scaler, "Scale input tree"},
        {ARGS_OPTION_STRING,  'O', "treeopt",      "tree.topolology.optimize", &topology_optimization_algorithm, "Optimize topology nni or spr (experimental)"},
        {ARGS_OPTION_FLAG,    0,   "parsimony",    "tree.topolology.optimize.parsimony", &use_parsimony, "Quick optimizaton of tree topology using parsimony before ML optimization"},
        
        {ARGS_OPTION_STRING,  'C', "clock",        "clock", &clock, "Clock type: strict, local, discrete"},
        {ARGS_OPTION_FLAG,    0,   "forward",      "clock.forward", &forward, "Time is forward"},
        {ARGS_OPTION_DOUBLE,  0,   "clock-rate",   "clock.rate", &clock_rate_guess, "A rate guess"},
        {ARGS_OPTION_STRING,  'S', "clock-search", "clock.algorithm", &clock_algorithm, "Algorithm for local and discrete clock: ga or greedy or exhaustive"},
        {ARGS_OPTION_INTEGER, 0,   "clock-cat",    "clock.discrete.cat", &clock_categories, "Number of discrete rate categories along phylogeny"},
        
        // GA
        {ARGS_OPTION_INTEGER, 0,   "ga-pop",       "ga.popsize", &ga_population_size, "Genetic algorithm population size"},
        {ARGS_OPTION_INTEGER, 0,   "ga-gen",       "ga.ngen", &ga_generations, "Genetic algorithm number of generations"},
        {ARGS_OPTION_INTEGER, 0,   "ga-no-improv", "ga.maxnoimprovement", &ga_max_no_improvement, "Genetic algorithm number of generation without improvment before stopping"},
        
        {ARGS_OPTION_STRING,  0,   "ic",           "ic", &ic, "Information criterion (AIC, AICc, BIC)"},
        {ARGS_OPTION_INTEGER, 0,   "ic-ss",        "ic.samplesize", &ic_sample_size, "Sample size for information criterion"},

        {ARGS_OPTION_INTEGER, 'b', "bootstrap",    "resampling.bootstrap", &bootstrap, "Number of bootstrap replicates"},
        {ARGS_OPTION_FLAG,    0,   "bca",          "resampling.bootstrap.bca", &bca, "Use BCA bootstrap"},
        {ARGS_OPTION_INTEGER, 'j', "jackknife",    "resampling.jackknife", &jackknife, "Jackknife"},
        
        // for GA, greedy
        {ARGS_OPTION_INTEGER, 'T', "nthreads",     "nthreads", &nthreads, "Number of threads for GA and bootstrap algorithms"},

        {ARGS_OPTION_INTEGER, 'V', "verbose",      "verbosity", &verbosity, "Verbosity"},
        {ARGS_OPTION_LONG,    'R', "seed",         "random.seed", &seed, "Random seed"},
        {ARGS_OPTION_FLAG,    0,   "asr",          "asr", &asr, "Ancestral sequence reconstruction"},
        {ARGS_OPTION_FLAG,    0,   "double-matrix","distance.matrix.double.precision", &use_double_distance, "Use double precision for distance matrix"},
        {ARGS_OPTION_FLAG,    0,   "batch",        "batch", &batch, "Use batch mode"},
        {ARGS_OPTION_FLAG,    0,   "hessian",      "derivative.hessian", &hessian, "Calculate Hessian"},
        
        {ARGS_OPTION_STRING,  'D', "distance",      "distance", &distance, "NJ or UPGMA"},
        
        
        {ARGS_OPTION_INTEGER, 0,   "nexus-index",   "nexus.index", &nexus_index, "Index of tree in a multi-tree nexus file"},
        
        {ARGS_OPTION_FLAG,  '!',   "overwrite",     "", &overwrite, "Overwrite output files"},
        
        {ARGS_OPTION_FLAG,  'h',   "help",          "", &help, "Print help"},
        {ARGS_OPTION_FLAG,  0,     "hh",            "", &help2, "Print more help"},
    };
    
    args_parser* argsparser = argsparser_init(options, sizeof(options)/sizeof(struct argsparser_option));
    
    
    if (argc > 1 && argc <= 3) {
        if(!file_exists(argv[argc-1])){
            fprintf(stderr, "Cannot read input file: %s\n", argv[argc-1]);
            exit(1);
        }
        argsparser_parse_file(argsparser, argv[argc-1]);
    }
    else{
        argsparser_parse(argsparser, argv, argc);
    }
    
    char *program_name = get_program_name(argv);
    
    if (help || help2 || argc == 1) {
        
        fprintf(stdout, "\n%s:\n\n", program_name);
        printf("Fourment M and Holmes EC. Novel non-parametric models to estimate evolutionary rates and divergence times from heterochronous sequence data.\n");
        printf("BMC Evolutionary Biology 14:163, 2014\n\n");
        
        printf("Command options\n\n");
        
        argsparser_help(argsparser, (help2?1:0));
        
#ifndef DISABLED_CONFIG_HEADER
        fprintf(stdout, "\nLibrary used:\n");
        fprintf(stdout, "PhyC v%d.%d\n", PHYC_VERSION_MAJOR, PHYC_VERSION_MINOR );
        if(PHYC_SSE_ENABLED){
            fprintf(stdout, "  SSE     support: %s\n", PHYC_SSE_LEVEL);
        }
        else{
            fprintf(stdout, "  SSE     support: disabled\n" );
        }
        //fprintf(stdout, "  AVX     support: %s\n", (PHYC_AVX_ENABLED ? "enabled" : "disabled") );
        fprintf(stdout, "  OpenMP  support: %s\n", (PHYC_OPENMP_ENABLED ? "enabled" : "disabled") );
        fprintf(stdout, "  PThread support: %s\n", (PHYC_PTHREAD_ENABLED ? "enabled" : "disabled") );
        fprintf(stdout, "\n\n");
#endif
        printf("\n");
        printf("Example\n\n");
        printf("%s -i alignment.fa -o stem -m GTR -c 4\n\n", program_name);
        exit(0);
    }
	
    if (fix != NULL) {
        alpha_fixed = String_contains(fix, 'a');
        pinv_fixed = String_contains(fix, 'i');
        bl_fixed = String_contains(fix, 'd');
        frequencies_fixed = String_contains(fix, 'f');
        rates_fixed = String_contains(fix, 'r');
    }
    
    char* output_stem = NULL;
	if (output_stem_user != NULL) {
		output_stem = output_stem_user;
	}
	else {
		output_stem = String_clone(seq_file);
	}
    
    if (distance != NULL) {
        run_Distance( seq_file, output_stem, model_string_user, nthreads, distance, nexus_index, seed, bootstrap, jackknife );
        free(output_stem);
        free(distance);
        exit(0);
    }
    
    
    StringBuffer *buffer = new_StringBuffer(100);
	StringBuffer_append_strings(buffer, 2, output_stem, ".freerate.tree");
	char *freerate_filename = String_clone(buffer->c);
	
	StringBuffer_empty(buffer);
	StringBuffer_append_strings(buffer, 2, output_stem, ".strict.tree");
	char *strict_filename = String_clone(buffer->c);
    
	StringBuffer_empty(buffer);
	
	// can have crashed too
	bool freerate_done         = file_exists(freerate_filename);
	bool strict_done           = file_exists(strict_filename);
	
	bool run_strict = false;
	bool run_greedy_local = false;
	bool run_ga_local = false;
	bool run_ga_discrete = false;
	bool run_custom = false;
	
    if( clock != NULL ){
        run_strict = true;
        if( strcmp(clock, "strict") == 0 ){

        }
        else if( strcmp(clock, "local") == 0 ){
            run_greedy_local = true;
            
            if( clock_algorithm != NULL ){
                if( strcasecmp(clock_algorithm, "ga") == 0 ){
                    run_ga_local = true;
                    run_greedy_local = false;
                }
                else if( strcasecmp(clock_algorithm, "greedy") == 0 ){
                    run_greedy_local = true;
                }
                else if( strcasecmp(clock_algorithm, "exhaustive") == 0 ){
                    run_greedy_local = false;
                }
                else {
                    fprintf(stdout, "Algorithm type not recognized: %s\n\n", clock_algorithm);
                    exit(1);
                }
            }
            
        }
        else if( strcmp(clock, "discrete") == 0 ){
            run_ga_discrete = true;
        }
        else if( strcmp(clock, "custom") == 0 ){
            run_custom = true;
        }
        else {
            fprintf(stdout, "Clock type not recognized: %s\n\n", clock);
            exit(1);
        }
    }
	
	if ( !overwrite && (freerate_done || strict_done) ) {
		fprintf(stderr, "Some files that %s will create are already present in the directory\n", program_name);
		fprintf(stderr, "Use --overwrite to... overwrite\n");
		exit(1);
	}
	
	fprintf(stdout, "Methods selected:\n");
	fprintf(stdout, "  * No clock\n");
	if ( run_strict ) fprintf(stdout, "  * Strict clock\n");
	if ( run_greedy_local ) fprintf(stdout, "  * Greedy local clock\n");
	if ( run_ga_local )     fprintf(stdout, "  * GA local clock\n");
	if ( run_ga_discrete )  fprintf(stdout, "  * GA discrete clock\n");
	if ( run_custom )  fprintf(stdout, "  * GA custom clock\n");
	fprintf(stdout, "\n");
	

	init_genrand(seed);
	
	fprintf(stdout, "Random seed: %lu\n\n", seed);
	

	StringBuffer *info = new_StringBuffer(100);
	StringBuffer_append_string(info, "[!Physher\n\n");
	StringBuffer_append_strings(info, 3, "Alignment: ", seq_file, "\n");
    if( tree_file != NULL ){
        StringBuffer_append_strings(info, 3, "Tree:      ", tree_file, "\n");
    }
    else {
        StringBuffer_append_string(info, "Tree:      NJ\n");
    }
	
	fprintf(stdout, "Sequence file: %s\n", seq_file);
    if( tree_file != NULL ){
        fprintf(stdout, "Tree file:     %s\n\n", tree_file);
    }
    else {
        fprintf(stdout, "Tree file: NJ\n\n");
    }
	Sequences* seqs = NULL;
	if( nexus_index >= 0 ){
		seqs = readNexus(seq_file, nexus_index);
	}
	else {
		seqs = readSequences(seq_file);
	}
	
    check_alignment(seqs);
    
    assert( genetic_code >= 0 && genetic_code < 15 );
	
	fprintf(stdout, "Number of sequences: %d\n",seqs->size );
    
	fprintf(stdout, "Alignment length: %d\n",seqs->length );
    
    const char* nucleotide_models[] = {"JC69", "K80", "F81", "HKY", "GTR", "UREV", "NONSTAT"};
    const char* amino_acid_models[] = {"DAYHOFF", "LG", "WAG"};
    const char* codon_models[]      = {"GY94", "MG94"};
	
    /*************************************************************************************************
     ***************************************** Create DataType ***************************************
     *************************************************************************************************/
    
    // Determine the data type from model or number of states
    unsigned matrixDimension = 0;
    
    DataType* dataType = NULL;
    char* model_string = NULL;
    
    if( markov_states != NULL ){
        char ** states = String_to_string_array( markov_states, ',', &matrixDimension );
        dataType = new_GenericDataType(matrixDimension, (const char**)states);
        for ( int i = 0; i < matrixDimension; i++ ) {
            free(states[i]);
        }
        free(states);
        
        if( model_string_user != NULL ){
            const char* temp = model_string_user;
            if( isInt(temp) ){
                int classCount = atoi(temp);
                StringBuffer_empty(buffer);
                assert(classCount < matrixDimension*matrixDimension);
                
                
                int count = 0;
                for ( int i = 0; i < matrixDimension; i++) {
                    for ( int j = 0; j < matrixDimension; j++) {
                        if(i == j ){
                            if( i == 0 ) StringBuffer_append_string(buffer, ",0");
                            else StringBuffer_append_string(buffer, ",0");
                        }
                        else if( count < classCount ) {
                            StringBuffer_append_format(buffer, ",%d", count);
                            count++;
                        }
                        else {
                            StringBuffer_append_format(buffer, ",%d", random_int(classCount-1));
                        }
                    }
                }
                
                model_string = StringBuffer_tochar(buffer);
                
            }
            else {
                model_string = String_clone(model_string_user);
            }
        }
        else{
            model_string = String_clone("ER");
        }
    }
	else if( model_string_user != NULL ){
        model_string = String_clone(model_string_user);
        
        
        if(strlen(model_string) == 5 && model_string[0] == '0'){
            dataType = new_NucleotideDataType();
        }
        else if( array_of_string_contains(model_string, nucleotide_models, sizeof(nucleotide_models) / sizeof(nucleotide_models[0])) ){
            dataType = new_NucleotideDataType();
        }
        else if( array_of_string_contains(model_string, amino_acid_models, sizeof(amino_acid_models) / sizeof(amino_acid_models[0])) ){
            dataType = new_AminoAcidDataType();
            fprintf(stdout, "Genetic code: %s (%d)\n", GENETIC_CODE_NAMES[dataType->genetic_code], dataType->genetic_code );
        }
        else if( array_of_string_contains(model_string, codon_models, sizeof(codon_models) / sizeof(codon_models[0])) ){
            dataType = new_CodonDataType(genetic_code);
        }
        else {
            fprintf(stderr,"Does not recognize the model type (%s)\n", model_string);
            exit(1);
        }
        
        matrixDimension = dataType->state_count(dataType);
    }
    else {
        error("No substitution model\n");
    }
    
    seqs->datatype = clone_DataType(dataType);
	
	fprintf(stdout, "Sequence type: %s (State count: %d)\n", dataType->desc, dataType->state_count(dataType) );
    StringBuffer_append_format(info, "Sequence type: %s (State count: %d)\n", dataType->desc, dataType->state_count(dataType));
    
    if( markov_states != NULL ){
        StringBuffer_append_format(info, "States: %s\n\n", markov_states);
    }

	
	StringBuffer_append_format(info, "Number of sequences: %d\n", seqs->size);
	StringBuffer_append_format(info, "Alignment length: %d\n", seqs->length);
    
    /*************************************************************************************************
     **************************************** Create SitePattern *************************************
     *************************************************************************************************/
    
	SitePattern* sp = new_SitePattern( seqs );
    
    if( sp->nsites == 1 ){
        int* test = calloc(dataType->state_count(dataType), sizeof(int));
        for (int i = 0; i < sp->size; i++) {
            test[sp->patterns[0][i]]++;
        }
        bool flag = false;
        for (int i = 0; i < dataType->state_count(dataType); i++) {
            if( test[i] == 0 ){
                fprintf(stderr, "State %s not assigned\n", dataType->state_string(dataType,i));
                flag = true;
            }
        }
        if(flag){
            exit(1);
        }
        free(test);
    }
	
	double unc_lk = unconstrained_lk(sp, seqs->length);
	
    int polymorphisms = SitePattern_polymorphic_count(sp);
    
	fprintf(stdout, "Number of polymorphic sites: %d/%d (%f)\n", polymorphisms, sp->nsites,((double)polymorphisms/sp->nsites) );
    
	StringBuffer_append_format(info, "Number of polymorphic site: %d/%d (%f)\n\n", polymorphisms, sp->nsites,((double)polymorphisms/sp->nsites) );
    
	fprintf(stdout, "Number of patterns: %d\n", sp->count );
	fprintf(stdout, "Unconstrained likelihood: %f\n\n", unc_lk );
	
	StringBuffer_append_format(info, "Number of patterns: %d\n", sp->count);
	StringBuffer_append_format(info, "Unconstrained likelihood: %f\n", unc_lk);
	
    
    StringBuffer_append_format(info, "Random seed: %lu\n\n", seed);
    
    
	/*************************************************************************************************
	 ****************************************** Create tree ******************************************
	 *************************************************************************************************/
	Tree* tree = NULL;

    if( tree_file != NULL ){
        char *treestring =  readTree( tree_file );
        
        tree = new_Tree( treestring, true );
        
        free(treestring);
	}
    else {
        distancematrix_model model = DISTANCE_MATRIX_UNCORRECTED;
        if( dataType->type == DATA_TYPE_AMINO_ACID ){
            model = DISTANCE_MATRIX_KIMURA;
        }
        
        if(use_double_distance){
            double** matrix = Sequences_distance(seqs,model);
            tree = new_NJ(seqs, matrix);
            free_dmatrix(matrix, seqs->size);
        }
        else {
            float** matrix = Sequences_distance_float(seqs,model);
            tree = new_NJ_float(seqs, matrix);
            free_fmatrix(matrix, seqs->size);
        }
    }
    
	if ( tree_scaler > 0 ) {
		Tree_scale_distance(tree, tree_scaler);
	}
	
	check_aln_tree(tree, sp);
    
	/*************************************************************************************************
	 *********************************** Create substitution model ***********************************
	 *************************************************************************************************/
	
	SubstitutionModel* mod = NULL;
    
    double *frequencies_user = NULL;
    if ( frequencies_string_user != NULL ) {
        if(strlen(frequencies_string_user) == 1 && tolower(frequencies_string_user[0]) == 'e'){
            frequencies_user = dvector(matrixDimension);
            for (int i = 0; i < matrixDimension; i ++) {
                frequencies_user[i] = 1.0/matrixDimension;
            }
        }
        else{
            unsigned nfreqs = 0;
            frequencies_user = String_to_double_array( frequencies_string_user, ',', &nfreqs);
            assert(nfreqs==matrixDimension);
        }
    }
    
    if ( dataType->type == DATA_TYPE_NUCLEOTIDE ) {
        bool unequal_frequencies = false;
        if(frequencies_user != NULL){
            for (int i = 0; i < 4; i ++) {
                if(fabs(frequencies_user[i] - 0.25) > 0.00001){
                    unequal_frequencies = true;
                    break;
                }
            }
        }
        
        if( strcasecmp("JC69", model_string) == 0 ){
            mod = new_JC69();
        }
        else if( strcasecmp("K80", model_string) == 0 ){
            mod = new_K80();
        }
        else if( strcasecmp("F81", model_string) == 0 ){
            mod = new_F81();
        }
        else if( strcasecmp("HKY", model_string) == 0 ){
            mod = new_HKY();
        }
        else if( strcasecmp("GTR", model_string) == 0 ){
            mod = new_GTR();
        }
        else if( strcasecmp("UREV", model_string) == 0 ){
            mod = new_UnrestrictedNucleotideModel();
        }
        else if( strcasecmp("NONSTAT", model_string) == 0 ){
            mod = new_NONSTATNucleotideModel();
        }
        else if( model_string[0] == '0' ){
            if(strcasecmp("00000", model_string) == 0 && frequencies_user != NULL && !unequal_frequencies){
                free(model_string);
                model_string = String_clone("JC69");
                mod = new_JC69();
            }
            else if( strcasecmp("01001", model_string) == 0 && frequencies_user != NULL && !unequal_frequencies){
                free(model_string);
                model_string = String_clone("K80");
                mod = new_K80();
            }
            else if(strcasecmp("00000", model_string) == 0){
                free(model_string);
                model_string = String_clone("F81");
                mod = new_F81();
            }
            else if( strcasecmp("01001", model_string) == 0){
                free(model_string);
                model_string = String_clone("HKY");
                mod = new_HKY();
            }
            else if( strcasecmp("01234", model_string) == 0){
                free(model_string);
                model_string = String_clone("GTR");
                mod = new_GTR();
            }
            else{
                mod = new_ReversibleNucleotideModel(model_string);
                if( !frequencies_fixed ){
                    ReversibleNucleotideModel_estimate_freqs(mod);
                }
            }
        }
        else{
            assert(0);
        }
        
        fprintf( stdout, "\nModel: %s\n", model_string );
        StringBuffer_append_strings(info, 3, "Model: ", model_string, "\n\n");
        
        /*
         * FREQUENCIES
         */
        if(frequencies_user != NULL){
            memcpy(mod->_freqs, frequencies_user, sizeof(double)*4);
        }
        else if(strcasecmp("JC69", model_string) != 0 && strcasecmp("K80", model_string) != 0 && strcasecmp("NONSTAT", model_string) != 0 ){
            empirical_frequencies(seqs, mod->_freqs);
        }
        
        /*
         * RATES
         */
        if ( rates_user != NULL ) {
            unsigned rateCount = 0;
            double *rates = String_to_double_array( rates_user, ',', &rateCount);
            if( rateCount == 1 ){
                fprintf(stdout, "Kappa (user): %f\n", *rates);
            }
            SubstitutionModel_set_rates(mod,rates);
            free(rates);
        }
        else if( Parameters_count(mod->rates) == 5 ){
            double rates[5];
            Sequence_rel_rates_empirical(seqs, rates);
            SubstitutionModel_set_rates(mod,rates);
        }
        else if( Parameters_count(mod->rates) == 1 ){
            double rate = Sequence_kappa_empirical(seqs, mod->_freqs);
            Parameters_set_value(mod->rates, 0, rate);
        }
    
    }
    else if ( dataType->type == DATA_TYPE_CODON ) {
        
        if( strcasecmp("GY94", model_string) == 0 ){
            mod = new_GY94(genetic_code);
        }
        else if( strcasecmp("MG94", model_string) == 0 ){
            mod = new_MG94(genetic_code);
        }
        
        double freqs_nuc[4];
        double *freqs = dvector(matrixDimension);
        int freq_est = 0;
        
        
        if( frequencies_string_user != NULL ){
            if( strcasecmp("F1x4", frequencies_string_user) == 0 ){
                freq_est = 1;
            }
            else if( strcasecmp("F3x4", frequencies_string_user) == 0 ){
                freq_est = 2;
            }
            else {
                printf("Condon frequency unknwon (F1x3 or F3x4)");
                exit(1);
            }
        }
        
        fprintf( stdout, "\nModel: %s\n", model_string );
        StringBuffer_append_strings(info, 3, "Model: ", model_string, "\n\n");
        
        printf("Nucleotide frequencies ");
        if( freq_est == 1 ){
            printf("F1x4\n");
            empirical_nuc_frequencies( seqs, freqs_nuc);
            printf(" A: %f\n", freqs_nuc[0] );
            printf(" C: %f\n", freqs_nuc[1] );
            printf(" G: %f\n", freqs_nuc[2] );
            printf(" T: %f\n", freqs_nuc[3] );
            printf("\n");
            
            calculate_F1x4_from_nuc_freq(dataType->genetic_code, freqs_nuc, freqs);
            memcpy(mod->_freqs, freqs, sizeof(double)*matrixDimension);
        }
        else if( freq_est == 2 ){
            printf("F3x4\n");
            double nuc_freq[12];
            calculate_nuc_freq_codon( seqs, nuc_freq);
            printf("    1        2        3        Total\n");
            printf(" A: %f %f %f %f\n", nuc_freq[0], nuc_freq[4], nuc_freq[8],  ((nuc_freq[0]+nuc_freq[4]+nuc_freq[8])/3) );
            printf(" C: %f %f %f %f\n", nuc_freq[1], nuc_freq[5], nuc_freq[9],  ((nuc_freq[1]+nuc_freq[5]+nuc_freq[9])/3) );
            printf(" G: %f %f %f %f\n", nuc_freq[2], nuc_freq[6], nuc_freq[10], ((nuc_freq[2]+nuc_freq[6]+nuc_freq[10])/3) );
            printf(" T: %f %f %f %f\n", nuc_freq[3], nuc_freq[7], nuc_freq[11], ((nuc_freq[3]+nuc_freq[7]+nuc_freq[11])/3) );
            printf("\n");
            
            calculate_F3x4_from_nuc_freq(dataType->genetic_code, nuc_freq, freqs);
            
            memcpy(mod->_freqs, freqs, sizeof(double)*matrixDimension);
            
            freqs_nuc[0] = (nuc_freq[0]+nuc_freq[4]+nuc_freq[8])/3;
            freqs_nuc[1] = (nuc_freq[1]+nuc_freq[5]+nuc_freq[9])/3;
            freqs_nuc[2] = (nuc_freq[2]+nuc_freq[6]+nuc_freq[10])/3;
            freqs_nuc[3] = 1 - (freqs_nuc[0]+freqs_nuc[1]+freqs_nuc[2]);
        }
        free(freqs);
        
        /*
         * RATES
         */
        if( rates_user != NULL ){
            unsigned rateCount = 0;
            double *rates = String_to_double_array( rates_user, ',', &rateCount);
            SubstitutionModel_set_rates(mod,rates);
            
            if(rateCount == 2){
                fprintf(stdout, "Kappa (user): %f\n", rates[0]);
                fprintf(stdout, "Omega (user): %f\n", rates[1]);
            }
            else if(rateCount == 3) {
                fprintf(stdout, "Kappa (user): %f\n", rates[0]);
                fprintf(stdout, "Alpha (user): %f\n", rates[1]);
                fprintf(stdout, "Beta  (user): %f\n", rates[2]);
            }
            else{
                error("Too many rate parameters [-r]\n");
            }
            free(rates);
        }
    }
    else if ( dataType->type == DATA_TYPE_AMINO_ACID ) {
        fprintf( stdout, "\nModel: %s\n", model_string );
        StringBuffer_append_strings(info, 3, "Model: ", model_string, "\n\n");
        
        if( strcasecmp("DAYHOFF", model_string) == 0 ){
            mod = new_DAYHOFF();
        }
        else if( strcasecmp("LG", model_string) == 0 ){
            mod = new_LG();
        }
        else if( strcasecmp("WAG", model_string) == 0 ){
            mod = new_WAG();
        }
        else {
            assert(0);
        }
        
        if ( frequencies_string_user == NULL) {
            empirical_frequencies(seqs, mod->_freqs);
        }
    }
    else if ( dataType->type == DATA_TYPE_GENERIC) {
        fprintf( stdout, "\nModel: Generic\n" );
        StringBuffer_append_string(info, "Model: Generic\n\n");
        
        // empirical
        unsigned *rateIndexes = NULL;
        
        if( strcasecmp(model_string,"ER") == 0 ){
            rateIndexes = uivector(matrixDimension*matrixDimension);
        }
        else if( strcasecmp(model_string,"SYM") == 0 ){
            rateIndexes = uivector(matrixDimension*matrixDimension);
            int count = 0;
            for ( int i = 0; i < matrixDimension; i++) {
                rateIndexes[i*matrixDimension+i] = 0;
                for ( int j = i+1; j < matrixDimension; j++) {
                    rateIndexes[i*matrixDimension+j] = rateIndexes[j*matrixDimension+i] = count++;
                }
            }
        }
        else {
            unsigned len = 0;
            rateIndexes = String_to_uint_array(model_string, ',', &len);
            if( matrixDimension*matrixDimension != len ){
                fprintf(stderr, "Number of states (%d) does not match the vector length (%d) \n", matrixDimension, len);
                exit(1);
            }
        }
        
        mod = new_GeneralModel(rateIndexes, matrixDimension);
        
        free(rateIndexes);
        
        if(frequencies_user != NULL){
            memcpy(mod->_freqs, frequencies_user, sizeof(double)*matrixDimension);
        }
        else if(!frequencies_unknown){
            empirical_generic_frequencies(seqs, mod->_freqs);
        }
        if(!normalize_q){
            mod->normalize = false;
            mod->need_update = true;
        }
        if(frequencies_unknown){
            for ( int i = 0; i < matrixDimension; i++ ) {
                mod->_freqs[i] = 1.0;
            }
            mod->need_update = true;
        }
        
        
        if ( rates_user != NULL ) {
            unsigned rateCount = 0;
            double *rates = String_to_double_array(rates_user, ',', &rateCount);
            SubstitutionModel_set_rates(mod, rates);
            free(rates);
        }
    }
    else {
        assert(0);
    }
    
    if(frequencies_user) free(frequencies_user);
    
    free_Sequences(seqs);
    seqs = NULL;
	
	
	/*************************************************************************************************
	 **************************************** Create site model
	 *************************************************************************************************/
	
	SiteModel *sm = NULL;
	

    if(pinv == 0.0 && use_pinv){
        pinv = fmax(0.01, 1.0 - ((double)polymorphisms/sp->nsites));
    }

	bool use_gamma = false;
    
    if ( pinv > 0 && rate_category_count > 1 ) {
        use_gamma = true;
        sm = new_GammaPinvSiteModel( mod, pinv, alpha, rate_category_count );
    }
    else if ( rate_category_count > 1 ) {
        if(sitemodel_string != NULL){
            if ( strcasecmp(sitemodel_string, "gammaquad") == 0 ) {
                use_gamma = true;
                sm = new_GammaLaguerreSiteModel(mod, alpha, rate_category_count);
            }
            else if(strcasecmp(sitemodel_string, "discrete") == 0){
                sm = new_DiscreteSiteModel( mod, rate_category_count );
            }
            else{
                error("Rate heterogeneity not valid (discrete or gammaquad)\n");
            }
        }
        else{
            sm = new_GammaSiteModel(mod, alpha, rate_category_count);
            use_gamma = true;
        }
    }
    else if ( pinv > 0 ) {
        sm = new_PinvSiteModel(mod, pinv);
        use_pinv = true;
    }
    else {
        sm = new_SiteModel(mod);
    }
	
	
    
	fprintf(stdout, "Frequencies:\n");
    if( dataType->type == DATA_TYPE_GENERIC ){
        for ( int i = 0; i < matrixDimension; i++) {
            printf(" %s: %f\n", dataType->state_string(dataType, i), mod->_freqs[i]);
        }
    }
    else {
        print_frequencies(stdout, mod);
    }
	
    if( mod->rates != NULL ){
        fprintf(stdout, "\nRates:\n");
        if( dataType->type == DATA_TYPE_GENERIC ){
            for ( int i = 0; i < Parameters_count(mod->rates); i++ ) {
                fprintf(stdout,"%d %f\n", i, Parameters_value(mod->rates, i));
            }
//            int index = 0;
//            for ( int i = 0; i < matrixDimension; i++) {
//                for ( int j = i+1; j < matrixDimension; j++, index++) {
//                    if( Parameters_count(mod->rates) == index ){
//                        printf(" %s %s: 1\n", sp->datatype->state_string(sp->datatype, i), sp->datatype->state_string(sp->datatype, j));
//                    }
//                    else {
//                        printf(" %s %s: %f\n", sp->datatype->state_string(sp->datatype, i), sp->datatype->state_string(sp->datatype, j), Parameters_value(mod->rates, mod->model[index]));
//                    }
//                }
//            }
        }
        else {
            print_rates(stdout, mod);
        }
    }
	
	if ( sm->pinv != NULL ) {
		fprintf(stdout, "\nProportion of invariant sites: pinv = %f\n", Parameter_value(sm->pinv) );
	}
	
	if ( sm->shape != NULL ) {
		fprintf(stdout, "\nGamma model: shape = %f (%d categories)\n\n", Parameter_value(sm->shape), (sm->pinv == NULL ? sm->cat_count : sm->cat_count-1) );
	}
	
    
    if(approximation > 0) printf("approximation %d\n", approximation);
    // nucleotide: SSE for lower and upper likelihood calculations WITHOUT threads. Threads ONLY for lower calculation WITHOUT SSE.
    // codon:      SSE for lower likelihood calculations WITHOUT threads. Threads ONLY for lower calculation WITHOUT SSE.
    
    // Threads are only used if SSE is not used
    if( use_sse ){
        tlk_threads = 1;
    }
    
    
	/*************************************************************************************************
	 **************************************** Compute rate-free likelihood
	 *************************************************************************************************/

	SingleTreeLikelihood *tlk = new_SingleTreeLikelihood( tree, sm, sp, NULL );
	
	time_t start_time;
	time_t beginning_of_time;
	time_t end_time;
	double diff_time;
	
	fprintf(stdout, "\n\n***************************************************\n" );
	fprintf(stdout, "* Rate-free estimation\n\n" );
	
	time(&start_time);
	beginning_of_time = start_time;
	
	tlk->opt.bl.optimize             = !bl_fixed;
    tlk->opt.bl.tolx = 0.00000001;
    //tlk->opt.bl.tolx = 0.0001;
    //tlk->opt.bl.method = OPT_CG_PR;
	tlk->opt.freqs.optimize          = !frequencies_fixed;
	//tlk->opt.freqs.method  = OPT_CG_PR;
    //if( Parameters_count(mod->rates) > 1 ) tlk->opt.relative_rates.method = OPT_CG_PR; //Brent optmize better in the flub simulated data
    tlk->opt.relative_rates.optimize = !rates_fixed;
	tlk->opt.pinv.optimize           = use_pinv && !pinv_fixed;
	tlk->opt.gamma.optimize          = use_gamma && !alpha_fixed;
    tlk->opt.topology_optimize       = false;
    tlk->opt.verbosity = verbosity;

    
    if ( topology_optimization_algorithm != NULL ) {
        tlk->opt.topology_optimize = true;
        if( strcasecmp(topology_optimization_algorithm, "nni") == 0 ){
            tlk->opt.topology_alogrithm = TREE_SEARCH_NNI;
        }
        else if( strcasecmp(topology_optimization_algorithm, "nnni") == 0 ){
            tlk->opt.topology_alogrithm = TREE_SEARCH_NNNI;
        }
        else if( strcasecmp(topology_optimization_algorithm, "spr") == 0 ){
            tlk->opt.topology_alogrithm = TREE_SEARCH_SPR;
        }
        else {
            tlk->opt.topology_alogrithm = TREE_SEARCH_NNI;
        }
        printf("topology algorithm %d\n", tlk->opt.topology_alogrithm);
    }
	
	if ( mod->modeltype == GY94 || mod->modeltype == WAG || mod->modeltype == DAYHOFF || mod->modeltype == LG ) {
		tlk->opt.freqs.optimize = false;
	}
    if ( mod->modeltype == WAG || mod->modeltype == DAYHOFF || mod->modeltype == LG ) {
        tlk->opt.relative_rates.optimize = false;
    }
    
    SingleTreeLikelihood_set_nthreads(tlk, tlk_threads);
    SingleTreeLikelihood_enable_SSE(tlk, use_sse);
    SingleTreeLikelihood_use_upper(tlk, use_upper);
    if( use_upper ){
        fprintf(stdout, "\nUse upper for likelihood calculation\n");
    }
    
    if( SingleTreeLikelihood_SSE(tlk) ){
        fprintf(stdout, "Use SSE for likelihood calculation\n");
    }
    
    if(ambiguity){
        SingleTreeLikelihood_use_ambiguity(tlk);
    }
    
    fprintf(stdout, "\nNot optimized LnL = %f\n\n", tlk->calculate(tlk) );
	
	fprintf(stdout, "Estimation of branch length and substitution matrix parameters\n\n" );
    
    if( tlk->opt.topology_optimize && use_parsimony ){
        time_t t1;
        time_t t2;
        time(&t1);
        
        TopologyOptimizer *topology = new_TopologyOptimizer( tlk, TREE_SEARCH_PARSIMONY_SPR );
        
        if ( nthreads > 1 ) {
            TopologyOptimizer_set_nthreads( topology, nthreads);
        }

        topology->optimize(topology);
        free_TopologyOptimizer(topology);
        SingleTreeLikelihood_update_all_nodes(tlk);
        
        time(&t2);
        
        diff_time = difftime(t2, t1);
        fprintf(stdout, "\nTotal runtime ");
        print_pretty_time(stdout, diff_time);
        fprintf(stdout, "\n");
    }

//    if( mod->modeltype == GY94 ){
//        printf("Optimizing HKY model (kappa: %f)\n",Parameters_value(mod->rates, 0));
//        SubstitutionModel *hky = new_HKY_with_values(freqs_nuc, Parameters_value(mod->rates, 0));
//        SiteModel *sm2 = new_GammaSiteModel(hky, 0.5,4);
//        
//        seqs->datatype->type = DATA_TYPE_NUCLEOTIDE;
//        SitePattern *sp2 = new_SitePattern( seqs );
//        SingleTreeLikelihood *tlk2 = new_SingleTreeLikelihood(tlk->tree, sm2, sp2, NULL);
//        
//        tlk2->opt.bl.optimize             = true;
//        tlk2->opt.bl.tolx = 0.00000001;
//        tlk2->opt.freqs.optimize          = true;
//        if( Parameters_count(hky->rates) > 1 ) tlk2->opt.relative_rates.method = OPT_CG_PR; //Brent optmize better in the flub simulated data
//        tlk2->opt.relative_rates.optimize = true;
//        tlk2->opt.pinv.optimize           = use_pinv;
//        tlk2->opt.gamma.optimize          = use_gamma;
//        
//        double lk = optimize_singletreelikelihood(tlk2);
//        
//        fprintf(stdout, "\nOptimized LnL = %f\n\n", lk );
//        Parameters_set_value(mod->rates, 0, Parameters_value(hky->rates, 0));
//        
//        free_SingleTreeLikelihood_share2(tlk2, true, false, false, false);
//        printf("Optimizing GY94 model\n");
//    }
    
    
    
//    SubstitutionModel *jc = new_JC69();
//    tlk->sm->m = jc;
//    tlk->opt.max_rounds = 1;
//    optimize_singletreelikelihood(tlk);
//    
//    tlk->sm->m = mod;
//    tlk->opt.max_rounds = 1000;
//    free_SubstitutionModel_share(jc, false, false);
    
    tlk->opt.topology_threads = nthreads;
    
    
    tlk->opt.interruptible = true;
    
    time(&start_time);
    
    double lk = optimize_singletreelikelihood(tlk);
    
    tlk->opt.interruptible = false;
    
	double clock_free_lnl = lk;

	fprintf(stdout, "\nOptimized LnL = %f\n", lk );
    
    double aicc = AICc(lk, SingleTreeLikelihood_df_count(tlk), sp->nsites);
    double aic  = AIC(lk, SingleTreeLikelihood_df_count(tlk));
    double bic  = BIC(lk, SingleTreeLikelihood_df_count(tlk), sp->nsites);
	
	StringBuffer_append_format(info, "Log-likelihood: %f\n", lk);
    
    StringBuffer_append_format(info, "AIC:  %f df: %d\n", aic, SingleTreeLikelihood_df_count(tlk));
    StringBuffer_append_format(info, "AICc: %f df: %d\n", aicc, SingleTreeLikelihood_df_count(tlk));
    StringBuffer_append_format(info, "BIC:  %f df: %d\n\n", bic, SingleTreeLikelihood_df_count(tlk));
	
	if ( mod->dtype == DATA_TYPE_CODON ) {
		fprintf(stdout, "\nRates:\n");
		print_rates(stdout, mod);
		
		StringBuffer_append_string(info, "Relative rates:\n");
		bufferize_rates(info, mod);
	}
	else if( mod->dtype == DATA_TYPE_NUCLEOTIDE){
		if ( mod->modeltype != JC69 ) {
			fprintf(stdout, "Frequencies:\n");
			print_frequencies(stdout, mod);
			
			fprintf(stdout, "\nRates:\n");
			print_rates(stdout, mod);
			
			StringBuffer_append_string(info, "Frequencies:\n");
			bufferize_frequencies(info, mod);
			
			StringBuffer_append_string(info, "Relative rates:\n");
			bufferize_rates(info, mod);
		}
	}
    else if( mod->dtype == DATA_TYPE_GENERIC ){
        fprintf(stdout, "Rates:\n");
        StringBuffer_append_string(info, "\nRates\n");
        
        for ( int i = 0; i < Parameters_count(mod->rates); i++ ) {
            fprintf(stdout,"%d %f\n", i, Parameters_value(mod->rates, i));
            StringBuffer_append_format(info, "%d %f\n", i, Parameters_value(mod->rates, i));
        }
        
//        fprintf(stdout, "\nRate indexes:\n");
//        for ( int i = 0; i < matrixDimension; i++) {
//            StringBuffer_append_format(info, ",%s", sp->datatype->state_string(sp->datatype, i));
//        }
//        StringBuffer_append_string(info, "\n");
//        int index = 0;
//        for ( int i = 0; i < matrixDimension; i++) {
//            StringBuffer_append_format(info, "%s", sp->datatype->state_string(sp->datatype, i));
//            for ( int j = 0; j < matrixDimension; j++,index++) {
//                if( i == j ){
//                    StringBuffer_append_string(info, "*");
//                }
//                StringBuffer_append_format(info, "%s", sp->datatype->state_string(sp->datatype, i));
//                
//                if( Parameters_count(mod->rates) == index ){
//                    printf(" %s %s: 1\n", sp->datatype->state_string(sp->datatype, i), sp->datatype->state_string(sp->datatype, j));
//                }
//                else {
//                    printf(" %s %s: %f\n", sp->datatype->state_string(sp->datatype, i), sp->datatype->state_string(sp->datatype, j), Parameters_value(mod->rates, mod->model[index]));
//                }
//            }
//        }
    }
	
	if ( sm->pinv != NULL ) {
		fprintf(stdout, "\nProportion of invariant sites: pinv = %f\n", Parameter_value(sm->pinv) );
		StringBuffer_append_format(info, "P-inv: %f\n", Parameter_value(sm->pinv));
	}
	
	if ( sm->shape != NULL ) {
		fprintf(stdout, "\nGamma model: shape = %f (%d categories)\n\n", Parameter_value(sm->shape), (sm->pinv == NULL ? sm->cat_count : sm->cat_count-1) );
		StringBuffer_append_format(info, "Gamma model: shape = %f (%d categories)\n", Parameter_value(sm->shape), (sm->pinv == NULL ? sm->cat_count : sm->cat_count-1) );
	}
	
	
	FILE *freerarte_file = fopen(freerate_filename,"w");
	assert(freerarte_file);
	Tree_print_newick(freerarte_file, tree, true);
	fclose(freerarte_file);
	
	
	time(&end_time);
	diff_time = difftime(end_time, start_time);
	fprintf(stdout, "\nTotal runtime ");
	print_pretty_time(stdout, diff_time);
	
    
    if(asr && sm->type != SITEMODEL_NONE ){
        StringBuffer_empty(buffer);
        StringBuffer_append_strings(buffer, 2, output_stem, ".posteriors");
        calculate_posteriors_sites(tlk, buffer->c);
    }
    
    StringBuffer_empty(buffer);
    StringBuffer_append_strings(buffer, 2, output_stem, ".info");
    FILE *info_file = fopen(buffer->c,"w");
	assert(info_file);
	fprintf(info_file, "%s", info->c);
	fclose(info_file);
    
    if ( (bootstrap > 0 || jackknife > 0 ) && !run_greedy_local && !run_custom && !run_strict ) {
        run_bootstrap(tlk, output_stem, nthreads, bootstrap>0, (bootstrap>0? bootstrap: jackknife), bca);
    }
    
    
    if( asr ){
        printf("Inferring ancestral sequences...");
        asr_marginal(tlk);
        Sequences *sss = SitePattern_to_Sequences( tlk->sp );
        StringBuffer_empty(buffer);
        StringBuffer_append_strings(buffer, 2, output_stem, ".asr.fa");
        Sequences_save_fasta(sss, buffer->c);
        
        printf("done\n");
        
        Node **nodes = Tree_nodes(tlk->tree);
        for ( int i = 0; i < Tree_node_count(tlk->tree); i++ ) {
            tlk->mapping[i] = get_sequence_index(tlk->sp, nodes[i]->name);
        }
        tlk->mapping[Tree_root(tree)->id] = -1;
        
        sm->need_update = true;
    }
    

    if( qsearch != NULL  ){
        if( strcasecmp(qsearch, "ga") == 0 ){
            run_searchq(tlk, output_stem, ga_population_size, ga_generations, ga_max_no_improvement, nthreads );
        }
        else if( strcasecmp(qsearch, "lasso") == 0 ){
            qsearch_lasso(tlk);
        }
    }
    
    if( hessian  ){
        int nodeCount = Tree_node_count(tlk->tree);
        double *mat = dvector(nodeCount*nodeCount);
        calculate_hessian_branches(tlk, mat);
        int k = 0;
        Node *root = Tree_root(tlk->tree);
        Node *unwantedNode = Node_right(root);
        
        double *mat2 = dvector((nodeCount-2)*(nodeCount-2));
        for (int i = 0; i < nodeCount; i++) {
            Node *n1 = Tree_node(tlk->tree, i);
            if(n1 == root || n1 == unwantedNode) continue;
            
            for (int j = 0; j < nodeCount; j++) {
                Node *n2 = Tree_node(tlk->tree, j);
                if(n2 == root || n2 == unwantedNode) continue;
                
                mat2[k++] = mat[i*nodeCount+j];
            }
        }
        
        printf("Diagonal hessian\n");
        printf("name,branch\n");
        Node **nodes  = Tree_get_nodes(tlk->tree, POSTORDER);
        for (int i = 0; i < nodeCount; i++) {
            Node *n1 = nodes[i];
            if(n1 == root || n1 == unwantedNode) continue;
            printf("%s,%e\n", n1->name, mat[Node_id(n1)*nodeCount+Node_id(n1)] );
        }
        printf("\n\n");
        
        
        for (int i = 0; i < nodeCount-2; i++) {
            for (int j = 0; j < nodeCount-2; j++) {
                mat2[i*(nodeCount-2)+j] = -mat2[i*(nodeCount-2)+j];
            }
        }
        
        inverse2(mat2, nodeCount-2);
        for (int i = 0; i < nodeCount-2; i++) {
            for (int j = 0; j < nodeCount-2; j++) {
                printf("%e%s",mat2[i*(nodeCount-2)+j], (j == nodeCount-3 ? "\n":","));
            }
        }
    }
    
	/*************************************************************************************************
	 ****************************** Compute derivatives for Taylor expansion *************************
	 *************************************************************************************************/
    
    if( approximation != TREELIKELIHOOD_APPROXIMATION_NONE ){
        
        init_likelihood_approximation(tlk, approximation);
    }
   
    
    BranchModel *bm = NULL;
    SingleTreeLikelihood_use_upper(tlk, false); // does not work well with clocks
    
	/*************************************************************************************************
	 ********************************* Compute strict clock likelihood *******************************
	 *************************************************************************************************/
    double rate_guess = 0;
    
    // always run except when we restart the discrete GA
    if ( run_strict ){
        fprintf(stdout, "\n\n***************************************************\n" );
        fprintf(stdout, "* Strict clock estimation\n\n" );
        
        time(&start_time);
        
        // Codon does not support upper and SSE
        if( use_sse ){
            SingleTreeLikelihood_enable_SSE(tlk, true);
        }
        if( SingleTreeLikelihood_SSE(tlk) ){
            fprintf(stdout, "Use SSE for likelihood calculation\n\n");
        }
        
        
        tlk->opt.bl.optimize             = false;
        tlk->opt.freqs.optimize          = false;
        tlk->opt.relative_rates.optimize = false;
        tlk->opt.pinv.optimize           = false;
        tlk->opt.gamma.optimize          = false;
        tlk->opt.heights.optimize        = true;
        tlk->opt.rates.optimize          = true;
        tlk->opt.topology_optimize       = false;
        tlk->opt.verbosity = 1;
        tlk->opt.max_rounds = 1000;
        

        rate_guess = init_times(tlk, forward, clock_rate_guess, tree_unrooted);
        
        // we could use values spaced on different log scales 0.000001  0.00001 0.0001 0.001 0.01
        if( rate_guess < 0){
            
        }

        
        bm = new_StrictClock(tree);
        SingleTreeLikelihood_set_BranchModel( tlk, bm, false );
        
        bm->set( bm, 0, rate_guess );
        
        lk = tlk->calculate(tlk);
        
        
        fprintf(stdout, "Not optimized\nLnL = %f\n", lk );
        fprintf(stdout, "Root height %f\n\n", Node_height(Tree_root(tree)));
        
        FILE *testFile = fopen(strict_filename,"w");
        assert(testFile);
        
        Tree_print_nexus_header_figtree_Taxa(testFile, tlk->tree);
        Tree_print_nexus_header_figtree_BeginTrees(testFile, tlk->tree);
        
        
    	// This is the tree with times guessed with rate approx and optimized branch lengths
        /*fprintf(testFile, "tree TREE0 [&lk=%f] = [&R] ", clock_strict_lnl);
        
        Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "%e", clock_strict_rate);
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            Node_empty_annotation(nodes[i]);
            if( !Node_isroot(nodes[i]) ){
                Node_set_annotation(nodes[i], "rate", buffer->c);
            }
        }
        Tree_print_nexus_with_annotation(testFile, tlk->tree);
        fprintf(testFile, "\n");
        fflush(testFile);*/
        
        fprintf(stdout, "\nEstimation of divergence times under a strict clock\n\n" );
        
        
        tlk->approx = approximation;

        if(approximation != TREELIKELIHOOD_APPROXIMATION_NONE && approximation != TREELIKELIHOOD_APPROXIMATION_DERIVATIVE_DIAGONAL ){
            tlk->opt.verbosity = 0;
        }

        
        
        if(approximation != TREELIKELIHOOD_APPROXIMATION_NONE){
            
            lk = optimize_singletreelikelihood(tlk);
            
            fprintf(stdout, "\nApproximated LnL = %f\n", lk );
            SingleTreeLikelihood_update_all_nodes(tlk);
            fprintf(stdout, "Recalculated LnL = %f\n", tlk->calculate(tlk) );
            tlk->approx = TREELIKELIHOOD_APPROXIMATION_NONE;
            lk = optimize_singletreelikelihood(tlk);
            tlk->approx = approximation;
        }
        else{
            bm->set(bm, 0,0.0057);
            SingleTreeLikelihood_update_all_nodes(tlk);
            tlk->opt.interruptible = true;
//            tlk->opt.heights.optimize = false;
//            lk = optimize_singletreelikelihood(tlk);
//            tlk->opt.rates.optimize = true;
            lk = optimize_singletreelikelihood(tlk);
            
            tlk->opt.interruptible = false;
        }
        
        fprintf(stdout, "\nLnL = %f\n", lk );
        fprintf(stdout, "Rate = %e\n", Parameters_value( bm->rates, 0) );
        fprintf(stdout, "Root height = %f", Node_height(Tree_root(tree)));
        if ( forward && Tree_dated(tree) ) {
            fprintf(stdout, " (%f)", (get_earliest_date(tree) - Node_height(Tree_root(tree))) );
        }
        fprintf(stdout, "\n");
        
        fprintf(testFile, "tree TREE0 [&lk=%f] = [&R] ", lk);

        Node **nodes = Tree_get_nodes(tlk->tree, POSTORDER);
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "%e", Parameters_value( bm->rates, 0));
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            Node_empty_annotation(nodes[i]);
            if( !Node_isroot(nodes[i]) ){
                Node_set_annotation(nodes[i], "rate", buffer->c);
            }
        }
        Tree_print_nexus_with_annotation(testFile, tlk->tree);
        fprintf(testFile, "\nEnd;");
        
        
        StringBuffer_append_string(info, "\nStrict clock\n");
        StringBuffer_append_format(info, "Log-likelihood: %f\n", lk);
        
        double aicc = AICc(lk, SingleTreeLikelihood_df_count(tlk), sp->nsites);
        double aic  = AIC(lk, SingleTreeLikelihood_df_count(tlk));
        double bic  = BIC(lk, SingleTreeLikelihood_df_count(tlk), sp->nsites);
        
        
        StringBuffer_append_format(info, "AIC:  %f df: %d\n", aic, SingleTreeLikelihood_df_count(tlk));
        StringBuffer_append_format(info, "AICc: %f df: %d\n", aicc, SingleTreeLikelihood_df_count(tlk));
        StringBuffer_append_format(info, "BIC:  %f df: %d\n\n", bic, SingleTreeLikelihood_df_count(tlk));
        
        StringBuffer_append_format(info, "Substitution rate: %e\n", Parameters_value( bm->rates, 0));
        StringBuffer_append_format(info, "Root height: %f", Node_height(Tree_root(tree)) );
        if ( forward && Tree_dated(tree) ) {
            StringBuffer_append_format(info, " (%f)", (get_earliest_date(tree) - Node_height(Tree_root(tree))) );
        }
        
        append_simultron(tlk, model_string, output_stem, info);
        
        StringBuffer_append_string(info, "\n]\n");
        fprintf(testFile, "\n\n%s\n\n", info->c);
        
        fprintf(stdout, "Correlation with branch length %f\n\n", BranchModel_correlation_distance(bm) );
        
        if( Tree_dated(tree) ){
            double p = LRT(lk, clock_free_lnl, 3, Tree_tip_count(tree)); // dfs not right but the difference is ok
            fprintf(stdout, "LRT H0: strict clock against H1: no clock\n  p:%e (df: %d)\n\n", p, (Tree_tip_count(tree)-3) );
        }
        else{
            //TODO
        }
        
        fclose(testFile);
            
        time(&end_time);
        diff_time = difftime(end_time, start_time);
        fprintf(stdout, "\nTotal runtime ");
        print_pretty_time(stdout, diff_time);
        
        
        fflush(stdout);
        
        double ci_alpha = 0.05;
        if( true ){
            double **cis = ConfidenceInterval_rates(tlk, ci_alpha);
            fprintf(stdout, "Rate confidence interval using profile likelihood at %.2f level %e - %e\n", ci_alpha, cis[0][0], cis[0][1] );
            free_dmatrix(cis, 1);
        }
        
        fflush(stdout);
        
        if ( (bootstrap>0 || jackknife > 0) && !run_greedy_local && !run_custom ) {
            run_bootstrap(tlk, output_stem, nthreads, bootstrap>0, (bootstrap>0? bootstrap: jackknife), bca);
        }
    }

    
    /*************************************************************************************************/
    
    SingleTreeLikelihood_set_nthreads(tlk,1);
	
    /*************************************************************************************************
	 ************************************ Compute custom clock likelihood ****************************
	 *************************************************************************************************/
    
    if( run_custom ){
        fprintf(stdout, "\n\n***************************************************\n" );
		fprintf(stdout, "* Custom clock estimation\n\n" );
        
        time(&start_time);
        
        run_custom_likelihood( tlk, output_stem, nthreads, bootstrap, jackknife, bca, forward );
        
        time(&end_time);
        diff_time = difftime(end_time, start_time);
        fprintf(stdout, "\nTotal runtime ");
        print_pretty_time(stdout, diff_time);
        
        fflush(stdout);
    }
    
    
    /*************************************************************************************************/
    
    
    tlk->approx = approximation;
    
    if(approximation != TREELIKELIHOOD_APPROXIMATION_NONE && approximation != TREELIKELIHOOD_APPROXIMATION_DERIVATIVE_DIAGONAL ){
        tlk->opt.verbosity = 0;
    }
    
    
	/*************************************************************************************************
	 ****************************** Compute greedy local clock likelihood ****************************
	 *************************************************************************************************/
	
    if ( run_greedy_local ) {
        fprintf(stdout, "\n\n***************************************************\n" );
        fprintf(stdout, "* Greedy local clock estimation\n\n" );
        
        time(&start_time);
        
        run_local_greedy_likelihood( tlk, output_stem, ic, nthreads, forward, ic_sample_size, (batch ? 1: 2) );
        
        time(&end_time);
        diff_time = difftime(end_time, start_time);
        fprintf(stdout, "\nTotal runtime ");
        print_pretty_time(stdout, diff_time);
        
        if ( bootstrap > 0 ) {
            SingleTreeLikelihood_bootstrap_strict_local_openmp(tlk, bootstrap, output_stem, false, nthreads);
        }
        
        fflush(stdout);
    }
	
	
	/*************************************************************************************************
	 ******************************** Compute GA local clock likelihood ******************************
	 *************************************************************************************************/
    
	if( run_ga_local ){
		fprintf(stdout, "\n\n***************************************************\n" );
		fprintf(stdout, "* Genetic algorithm local clock estimation\n\n" );
		time(&start_time);
        
        run_local_ga_likelihood( tlk, output_stem, ga_population_size, ga_generations, ga_max_no_improvement, nthreads, ic, forward );
        
        time(&end_time);
        diff_time = difftime(end_time, start_time);
        fprintf(stdout, "\nTotal runtime ");
        print_pretty_time(stdout, diff_time);
	}
	

	
	/*************************************************************************************************
	 ******************************** Compute GA discrete clock likelihood ***************************
	 *************************************************************************************************/
	
	
	
	if( run_ga_discrete ){
        fprintf(stdout, "\n\n***************************************************\n" );
        fprintf(stdout, "* Genetic algorithm discrete rate classes estimation\n\n" );
		time(&start_time);
        
        run_discrete_ga_likelihood( tlk , output_stem, ga_population_size, ga_generations, ga_max_no_improvement, nthreads, ic, forward, clock_categories, approximation );
        
        time(&end_time);
        diff_time = difftime(end_time, start_time);
        fprintf(stdout, "\nTotal runtime ");
        print_pretty_time(stdout, diff_time);
	}
	
    time(&end_time);
	diff_time = difftime(end_time, beginning_of_time);
	fprintf(stdout, "\nTotal runtime ");
	print_pretty_time(stdout, diff_time);
	fprintf(stdout, "\n");
	
	
	/*************************************************************************************************
	 **************************************** Free memory ********************************************
	 *************************************************************************************************/

    free(argsparser);
    
    if(seq_file) free(seq_file);
    if(tree_file) free(tree_file);
    if(fix) free(fix);
    if(sitemodel_string) free(sitemodel_string);
    if(clock) free(clock);
    if(markov_states) free(markov_states);
    if(frequencies_string_user) free(frequencies_string_user);
    if(rates_user) free(rates_user);
    if(qsearch) free(qsearch);
    if(topology_optimization_algorithm) free(topology_optimization_algorithm);
    if(clock_algorithm) free(clock_algorithm);
    if(model_string_user)free(model_string_user);
    free(output_stem);
    free(ic);
    
    free(model_string);
	free(freerate_filename);
	free(strict_filename);
	
	free_StringBuffer(buffer);	
	free_StringBuffer(info);
	free_SingleTreeLikelihood(tlk);
    free_DataType(dataType);
}
