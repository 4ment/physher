/*
 *  bootstrap.c
 *  physher
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <strings.h>

#include "phyc/tree.h"
#include "phyc/treeio.h"
#include "phyc/phyboot.h"
#include "phyc/args.h"
#include "phyc/descriptivestats.h"
#include "phyc/matrix.h"
#include "phyc/boot.h"

static void _fixed_topology_strict( Tree *tree, const char *file, const char *output_file, double ci, const char *jackfile );
static void _fixed_topology( Tree *tree, const char *file, const char *output_file, double ci, const char *jackfile );

int main(int argc, char* argv[]){
    
    char *boot_file = NULL;
    char *outfile = NULL;
    char *jackknifefile = NULL; // File containing jackkinfed trees for BCa
    
    Tree *target = NULL;
    bool containBL = false;
    
    double ci = 0.95; // 95% confidence interval

    bool success = false;
    char *targetfile = NULL;
    FILE *pFile = NULL;
    char *tree = NULL;
    
    if ( argc == 1 || args_get_boolean(argc, argv, "-h") || args_get_boolean(argc, argv, "--help") ) {
        printf("\nCalculate confidence intervals using the output tree files from Physher\n\n");
        printf("Command options\n\n");
        printf(" -i bootstrap trees (*.strict.boot.trees) [REQUIRED]\n");
        printf(" -t maximmum likelihood tree (*.strict.tree) [REQUIRED]\n");
        printf(" -o output tree file [REQUIRED]\n");
        printf(" -j jackknife trees for BCa (*.strict.jack.trees)\n");
        printf(" -c confidence interval [default: 0.95]\n");
        printf("\n");
        printf("Example\n\n");
        printf("./bootstrap -i file.strict.boot.trees -t file.strict.tree -o file.strict.boot.tree -c 0.95\n\n");
        exit(0);
    }
	   
    if ( args_contains(argc, argv, "-i") ) {
        boot_file  = args_get_string( argc, argv, "-i");
    }
    else {
        error("Need an input file containing bootstrap trees [-i]\n");
    }
    
	if ( args_contains(argc, argv, "-o") ) {
        outfile  = args_get_string( argc, argv, "-o");
    }
    else {
        error("Need an output file name [-o]\n");
    }
    
    if ( args_contains(argc, argv, "-j") ) {
        jackknifefile  = args_get_string( argc, argv, "-j");
    }
    
	// Needed to get Normal confidence intervals (need MLE)
	if ( args_contains(argc, argv, "-t")) {
		targetfile = args_get_string( argc, argv, "-t");
		if ( file_exists( targetfile) ) {
			tree = readTree(targetfile);
			target = new_Tree(tree, containBL);
			free(tree);
		}
		else {
			fprintf(stderr, "Could not read target file: %s\n", targetfile);
			exit(1);
		}
		free(targetfile);
	}
    else {
        error("Need the file name of the MLE tree [-t]\n");
    }
    
	if ( args_contains(argc, argv, "-c")) {
        success = false;
		ci = args_get_double(argc, argv, "-c", &success);
		if ( !success ) {
			error("Invalid confidence interval [0-1]\n");
		}
	}
    
//    bool local = false;
//    if ( args_contains(argc, argv, "-C")) {
//        char *type = args_get_string( argc, argv, "-C");
//        if ( strcasecmp(type, "local") == 0 ) {
//            local = true;
//        }
//        free(type);
//    }
//    
//    if( local ){
//        _fixed_topology( target, infile, outfile, ci, jackknifefile);
//    }
//    else {
//        _fixed_topology_strict( target, infile, outfile, ci, jackknifefile);
//    }
    
    if( !file_exists(boot_file) ){
        fprintf(stderr,"file does not exist: %s\n",boot_file);
        exit(1);
    }
    
    Phyboot_annotate_fixed_topology(target, boot_file, ci, jackknifefile);
    
    pFile = fopen(outfile,"w");
    assert(pFile);
    Tree_print_nexus_header_figtree_Taxa(pFile, target);
    Tree_print_nexus_header_figtree_BeginTrees(pFile, target);
    fprintf(pFile, "tree TREE0 = [&R] ");
    Tree_print_nexus_with_annotation( pFile, target );
    fprintf(pFile, "\nEnd;");
    fclose(pFile);
	
	free(boot_file);
	free(outfile);
    if( jackknifefile != NULL ) free(jackknifefile);
	if (target != NULL ) free_Tree(target);
	
	return 0;
}


/*void _fixed_topology_strict( Tree *tree, const char *file, const char *output_file, double ci, const char *jackfile ){
    
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    
    //    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
    //        printf("%s %f\n", Node_name(nodes[i]), Node_height(nodes[i]));
    //    }
    int nboot = 0;
    
    printf("\nCalculating confidence intervals using percentile and normal bootstrap methods [strict]\n");
    
    fprintf(stderr, "\nReading trees... ");
    double **matrix = Phyboot_read_strict_rate_and_heights(file, &nboot);
    fprintf(stderr, "done\n\n");
    
    fprintf(stderr, ". Substitution rate\n");
    
    qsort(matrix[Tree_node_count(tree)], nboot, sizeof(double), qsort_asc_dvector);
    
    double rate_median = dmedian_ordered(matrix[Tree_node_count(tree)], nboot);
    double rate_mean   = dmean(matrix[Tree_node_count(tree)], nboot);
    
    
    double qlb = (1.0-ci)/2.0;
    double qub = (1.0-ci)/2.0+ci;
    
    
    // Percentile confidence intervals
    double rate_p_low  = dpercentile_ordered(matrix[Tree_node_count(tree)], nboot, qlb);
    double rate_p_high = dpercentile_ordered(matrix[Tree_node_count(tree)], nboot, qub);
    
    StringBuffer *buffer  = new_StringBuffer(10);
    StringBuffer *buffer2 = new_StringBuffer(10);
    
    StringBuffer_append_format(buffer, "%e", rate_mean);
    char *rate_mean_str = String_clone(buffer->c);
    
    StringBuffer_empty(buffer);
    StringBuffer_append_format(buffer, "%e", rate_median);
    char *rate_median_str = String_clone(buffer->c);
    
    StringBuffer_empty(buffer);
    StringBuffer_append_format(buffer, "{%e,%e}", rate_p_low,rate_p_high);
    char *rate_percentile_CI_str = String_clone(buffer->c);
    
    StringBuffer_empty(buffer);
    StringBuffer_append_format(buffer, "rate_%.0f%%_CI", (ci*100));
    char *rate_percentile_CI_desc_str = String_clone(buffer->c);
    
    
    fprintf(stderr, "    Mean:    %e\n",  rate_mean);
    fprintf(stderr, "    Median:  %e\n", rate_median );
    fprintf(stderr, "    Percentile %.0f%% CI: [%e %e]\n", (ci*100), rate_p_low, rate_p_high);
    
    double clock_strict_rate = Node_get_double_from_info(nodes[0],"rate="); // Get rate MLE
    fprintf(stdout, "    Normal     %.0f%% CI: [%e %e]\n", (ci*100), boot_ci_norm2(matrix[Tree_node_count(tree)], nboot, clock_strict_rate, qlb), boot_ci_norm2(matrix[Tree_node_count(tree)], nboot, clock_strict_rate,qub));
    
    
    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
        qsort(matrix[i], nboot, sizeof(double), qsort_asc_dvector);
        
        
        //height
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "%e", Node_height(nodes[i]));
        Node_set_annotation(nodes[i], "height", buffer->c);
        
        if( !Node_isleaf(nodes[i]) ){
            double median = dmedian_ordered(matrix[i], nboot);
            double mean   = dmean(matrix[i], nboot);
            double p_low  = dpercentile_ordered(matrix[i], nboot, qlb);
            double p_high = dpercentile_ordered(matrix[i], nboot, qub);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", mean);
            Node_set_annotation(nodes[i], "height_mean", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", median);
            Node_set_annotation(nodes[i], "height_median", buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", p_low,p_high);
            StringBuffer_append_format(buffer2, "height_perc_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]), qlb), boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]), qub) );
            StringBuffer_append_format(buffer2, "height_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            if( Node_isroot(nodes[i]) ){
                
                fprintf(stdout, ". Root age\n");
                fprintf(stdout, "    Mean:    %e\n",  mean);
                fprintf(stdout, "    Median:  %e\n", median );
                fprintf(stdout, "    Percentile %.0f%% CI: [%e %e]\n", (ci*100), p_low, p_high);
                fprintf(stdout, "    Normal     %.0f%% CI: [%e %e]\n", (ci*100), boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]), qlb), boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]),qub));
            }
        }
        
        // rate
        if( !Node_isroot(nodes[i]) ){
            Node_set_annotation(nodes[i], "rate_mean", rate_mean_str);
            
            Node_set_annotation(nodes[i], "rate_median", rate_median_str);
            
            Node_set_annotation(nodes[i], rate_percentile_CI_desc_str, rate_percentile_CI_str);
        }
    }
    free(rate_mean_str);
    free(rate_median_str);
    free(rate_percentile_CI_str);
    free(rate_percentile_CI_desc_str);
    
    if( jackfile != NULL && file_exists(jackfile) ){
        printf("\nCalculating confidence intervals using bias-corrected and accelerated bootstrap method (BCa)\n");
        fprintf(stdout, "\nReading trees from %s... ", jackfile);
        int count = 0;// number of sitepatterns or trees in the file
        double ** jack1 = Phyboot_read_strict_rate_and_heights(jackfile, &count);
        double *weights = Phyboot_jackknife_weights( jackfile, &count );
        fprintf(stdout, "done\n\n");
        
        
        double rate_p_low  = bootci_BCa_weighted(jack1[Tree_node_count(tree)], count, weights, matrix[Tree_node_count(tree)], nboot, clock_strict_rate, qlb);
        double rate_p_high = bootci_BCa_weighted(jack1[Tree_node_count(tree)], count, weights, matrix[Tree_node_count(tree)], nboot, clock_strict_rate, qub);
        
        fprintf(stdout, ". Substitution rate\n");
        fprintf(stdout, "    BCa %.0f%% CI: [%e,%e]\n", (ci*100), rate_p_low, rate_p_high);
        
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            
            //height
            if( !Node_isleaf(nodes[i]) ){
                double p_low  = bootci_BCa_weighted(jack1[i], count, weights, matrix[i], nboot, Node_height(nodes[i]), qlb);
                double p_high = bootci_BCa_weighted(jack1[i], count, weights, matrix[i], nboot, Node_height(nodes[i]), qub);
                
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", p_low,p_high);
                StringBuffer_append_format(buffer2, "height_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
                
                if( Node_isroot(nodes[i]) ){
                    fprintf(stdout, ". Root age\n");
                    fprintf(stdout, "    BCa %.0f%% CI: [%e %e]\n", (ci*100), p_low, p_high);
                }
            }
            
            // rate
            if( !Node_isroot(nodes[i]) ){
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", rate_p_low,rate_p_high);
                StringBuffer_append_format(buffer2, "rate_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            }
            
        }
        
        free_dmatrix(jack1, Tree_node_count(tree)+1);
        free(weights);
    }
    
    free_dmatrix(matrix, Tree_node_count(tree)+1);
    
    
    FILE *testFile = fopen(output_file,"w");
    assert(testFile);
    Tree_print_nexus_header_figtree_Taxa(testFile, tree);
    Tree_print_nexus_header_figtree_BeginTrees(testFile, tree);
    fprintf(testFile, "tree TREE0 = [&R] ");
    Tree_print_nexus_with_annotation( testFile, tree );
    fprintf(testFile, "\nEnd;");
    fclose(testFile);
    free_StringBuffer(buffer2);
    
    free_StringBuffer(buffer);
    printf("\n");
    
}


void _fixed_topology( Tree *tree, const char *file, const char *output_file, double ci, const char *jackfile ){
    
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    
    //    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
    //        printf("%s %f\n", Node_name(nodes[i]), Node_height(nodes[i]));
    //    }
    int nboot = 0;
    
    printf("\nCalculating confidence intervals using percentile and normal bootstrap methods [local]\n");
    
    fprintf(stderr, "\nReading trees... ");
    double **matrix = Phyboot_read_rate_and_heights(file, &nboot);
    fprintf(stderr, "done\n\n");
    
    double qlb = (1.0-ci)/2.0;
    double qub = (1.0-ci)/2.0+ci;
    
    
    // Percentile confidence intervals
    double rate_p_low  = dpercentile_ordered(matrix[Tree_node_count(tree)], nboot, qlb);
    double rate_p_high = dpercentile_ordered(matrix[Tree_node_count(tree)], nboot, qub);
    
    StringBuffer *buffer  = new_StringBuffer(10);
    StringBuffer *buffer2 = new_StringBuffer(10);
    
    int nNodes = Tree_node_count(tree);
    
    for ( int i = 0; i < nNodes; i++ ) {
        qsort(matrix[i], nboot, sizeof(double), qsort_asc_dvector);          // heights
        qsort(matrix[nNodes+i], nboot, sizeof(double), qsort_asc_dvector);   // rates
        
        //height
        StringBuffer_empty(buffer);
        StringBuffer_append_format(buffer, "%e", Node_height(nodes[i]));
        Node_set_annotation(nodes[i], "height", buffer->c);
        
        if( !Node_isleaf(nodes[i]) ){
            double median = dmedian_ordered(matrix[i], nboot);
            double mean   = dmean(matrix[i], nboot);
            double p_low  = dpercentile_ordered(matrix[i], nboot, qlb);
            double p_high = dpercentile_ordered(matrix[i], nboot, qub);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", mean);
            Node_set_annotation(nodes[i], "height_mean", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", median);
            Node_set_annotation(nodes[i], "height_median", buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", p_low,p_high);
            StringBuffer_append_format(buffer2, "height_perc_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]), qlb), boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]), qub) );
            StringBuffer_append_format(buffer2, "height_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            if( Node_isroot(nodes[i]) ){
                
                fprintf(stdout, ". Root age\n");
                fprintf(stdout, "    Mean:    %e\n",  mean);
                fprintf(stdout, "    Median:  %e\n", median );
                fprintf(stdout, "    Percentile %.0f%% CI: [%e %e]\n", (ci*100), p_low, p_high);
                fprintf(stdout, "    Normal     %.0f%% CI: [%e %e]\n", (ci*100), boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]), qlb), boot_ci_norm2(matrix[i], nboot, Node_height(nodes[i]),qub));
            }
        }
        
        // rates
        if( !Node_isroot(nodes[i]) ){
            
            double median = dmedian_ordered(matrix[nNodes+i], nboot);
            double mean   = dmean(matrix[nNodes+i], nboot);
            double p_low  = dpercentile_ordered(matrix[nNodes+i], nboot, qlb);
            double p_high = dpercentile_ordered(matrix[nNodes+i], nboot, qub);
            
            StringBuffer_empty(buffer);
            double rate = Node_get_double_from_info(nodes[i], "rate=");
            StringBuffer_append_format(buffer, "%e", rate);
            Node_set_annotation(nodes[i], "rate", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", mean);
            Node_set_annotation(nodes[i], "rate_mean", buffer->c);
            
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%e", median);
            Node_set_annotation(nodes[i], "rate_median", buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", p_low,p_high);
            StringBuffer_append_format(buffer2, "rate_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            StringBuffer_empty(buffer2);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "{%e,%e}", boot_ci_norm2(matrix[nNodes+i], nboot, rate, qlb), boot_ci_norm2(matrix[nNodes+i], nboot, rate,qub) );
            StringBuffer_append_format(buffer2, "rate_perc_%.0f%%_CI", (ci*100));
            Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            
            // count the number of times a local clock was assigned to this node
            int nlocal = 0;
            for ( int j = 0; j < nboot; j++ ) {
                nlocal += matrix[nNodes*2+i][j];
            }
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%d", nlocal);
            Node_set_annotation(nodes[i], "local_count", buffer->c);
            
        }
    }
    
    // not checked
    if( jackfile != NULL && file_exists(jackfile) ){
        printf("\nCalculating confidence intervals using bias-corrected and accelerated bootstrap method (BCa)\n");
        fprintf(stdout, "\nReading trees from %s... ", jackfile);
        int count = 0;// number of sitepatterns or trees in the file
        double ** jack1 = Phyboot_read_rate_and_heights(jackfile, &count);
        double *weights = Phyboot_jackknife_weights( jackfile, &count );
        fprintf(stdout, "done\n\n");
    
        
        for ( int i = 0; i < Tree_node_count(tree); i++ ) {
            
            //height
            if( !Node_isleaf(nodes[i]) ){
                double p_low  = bootci_BCa_weighted(jack1[i], count, weights, matrix[i], nboot, Node_height(nodes[i]), qlb);
                double p_high = bootci_BCa_weighted(jack1[i], count, weights, matrix[i], nboot, Node_height(nodes[i]), qub);
                
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", p_low,p_high);
                StringBuffer_append_format(buffer2, "height_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
                
                if( Node_isroot(nodes[i]) ){
                    fprintf(stdout, ". Root age\n");
                    fprintf(stdout, "    BCa %.0f%% CI: [%e %e]\n", (ci*100), p_low, p_high);
                }
            }
            
            // rate
            if( !Node_isroot(nodes[i]) ){
                StringBuffer_empty(buffer2);
                StringBuffer_empty(buffer);
                StringBuffer_append_format(buffer, "{%e,%e}", rate_p_low,rate_p_high);
                StringBuffer_append_format(buffer2, "rate_%.0f%%_BCa_CI", (ci*100));
                Node_set_annotation(nodes[i], buffer2->c, buffer->c);
            }
            
        }
        
        free_dmatrix(jack1, Tree_node_count(tree)+1);
        free(weights);
    }
    
    free_dmatrix(matrix, Tree_node_count(tree)*3);
    
    
    FILE *testFile = fopen(output_file,"w");
    assert(testFile);
    Tree_print_nexus_header_figtree_Taxa(testFile, tree);
    Tree_print_nexus_header_figtree_BeginTrees(testFile, tree);
    fprintf(testFile, "tree TREE0 = [&R] ");
    Tree_print_nexus_with_annotation( testFile, tree );
    fprintf(testFile, "\nEnd;");
    fclose(testFile);
    free_StringBuffer(buffer2);
    
    free_StringBuffer(buffer);
    printf("\n");
    
}*/



