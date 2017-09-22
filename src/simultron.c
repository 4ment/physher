/*
 *  simultron.c
 *  physher
 *
 *  Created by Mathieu Fourment on 1/11/2013.
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
#include <time.h>
#include <strings.h>
#include <assert.h>

#include "phyc/tree.h"
#include "phyc/treeio.h"
#include "phyc/args.h"
#include "phyc/branchmodel.h"
#include "phyc/sequence.h"
#include "phyc/sequenceio.h"
#include "phyc/sitemodel.h"
#include "phyc/physim.h"
#include "phyc/random.h"
#include "phyc/matrix.h"
#include "phyc/lognormal.h"
#include "phyc/exponential.h"
#include "phyc/descriptivestats.h"
#include "phyc/geneticcode.h"

#include "phyc/substmodel.h"
#include "phyc/jc69.h"
#include "phyc/K80.h"
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

char * Node_get_string_from_info2( const Node *node, const char *str ){
    char *value = NULL;
    if ( node->info != NULL ) {
        char *ptr = node->info;
        if ( String_contains_str(ptr, str) ) {
            StringBuffer *buffer = new_StringBuffer(20);
            ptr += String_index_of_str(ptr, str) + strlen(str);
            bool range = false;
            if(*ptr=='{') range = true;
            
            if( range ){
                while( *ptr != '}' ){
                    StringBuffer_append_char(buffer, *ptr);
                    ptr++;
                }
                StringBuffer_append_char(buffer, *ptr);
            }
            else {
                while( *ptr != ',' && *ptr != ']' ){
                    
                    StringBuffer_append_char(buffer, *ptr);
                    ptr++;
                }
            }
            if (buffer->length !=0 ) {
                value = String_clone(buffer->c);
            }
            
            
            free_StringBuffer(buffer);
        }
    }
    return value;
}

bool array_of_string_contains(const char *str, const char *array[], int count){
    for ( int i = 0; i < count; i++) {
        if( strcasecmp(str, array[i]) == 0 ){
            return true;
        }
    }
    return false;
}

int main( int argc, char *argv[] ){
    if ( argc == 1 || args_get_boolean(argc, argv, "-h") || args_get_boolean(argc, argv, "--help") ) {
        printf("\nSimulate nucleotide data sets given a phylogeny\n\n");
        printf("Command options\n\n");
        printf(" -i input tree [REQUIRED]\n");
        printf(" -m substitution model: GTR, HKY, JC69, or K80 [REQUIRED]\n");
        printf(" -o output sequence file [REQUIRED]\n");
        printf(" -F sequence format: fasta,nexus or phylip]\n");
        printf(" -a rate heterogeneity: shape parameter of gamma distribution [default: 0.3]\n");
        printf(" -c rate heterogeneity: number of categories [default: 1]\n");
        printf(" -i proportion of invariant sites [default: 0]\n");
        printf(" -l length alignment [default: 1000]\n");
        printf(" -r substitution rate parameters(e.g. for GTR: 0.2,0.1,0.2,0.3,0.5)\n");
        printf(" -f Nucleotide frequencies (e.g. 0.25,0.25,0.25,0.25)\n");
        printf(" -s tree scaler. Multiply every branch by this value\n");
        printf(" -R random seed value\n");
        printf(" -d tree , logn\n");
        printf("    tree: rates in the tree file are used to multiply branches\n");
        printf("    logn: rates are drawn from a lognormal distrubiton and branches are multiplied with these rates\n");
        printf(" -M mean of the lognormal distribution\n");
        printf("\n");
        printf("Example\n\n");
        printf("./simultron -i file.tree -o data.fa -F fasta -m HKY -r 3 -l 500\n");
    }
    
    char * output = NULL;
    char *outputtree = NULL;
    Tree *tree = NULL;
    int length = 1000;
    double *rates = NULL;
    int cat = 1;
    double shape = 0.3;
    double scaler = 1;
    double pinv = 0;
    bool ancestral = false;
    int format = 1; // 0 is fasta, 1 is nexus, 2 should be phylip
    double min,max,mean,median;
    int gc = 0;//genetic code default: universal
    double *freq_nuc = NULL;
    double *freq_codon = NULL;
    DataType *dataType = NULL;
    char *model_string = NULL;
    
    const char *nucleotide_models[] = {"JC69", "K80", "HKY", "GTR", "UREV", "NONSTAT"};
    const char *amino_acid_models[] = {"DAYHOFF", "LG", "WAG"};
    const char *codon_models[]      = {"GY94","MG94"};
    
    // Determine the data type from model or number of states
    unsigned matrixDimension = 0;
    
    
    unsigned long seed = time(NULL);
    
    StringBuffer *buffer = new_StringBuffer(100);
    
    StringBuffer_append_string(buffer, "[! Generated by simultron\n\n");
    
    for ( int i = 0; i < argc; i++ ) {
        StringBuffer_append_string(buffer, argv[i]);
        if(i != argc-1 ) StringBuffer_append_string(buffer, " ");
    }
    StringBuffer_append_string(buffer, "\n\n");
    
    if ( args_contains(argc, argv, "-i") ) {
        char * inputree = args_get_string(argc, argv, "-i");
        char *treestring =  readTree( inputree );
        tree = new_Tree( treestring, true );
        
        free(treestring);
        free(inputree);
    }
    else {
        error("Need an input tree [-i]\n");
    }
    
    if ( args_contains(argc, argv, "-o") ) {
        output = args_get_string(argc, argv, "-o");
    }
    else {
        error("Need an output file name [-o]\n");
    }
    
    if ( args_contains(argc, argv, "-t") ) {
        outputtree = args_get_string(argc, argv, "-t");
    }
    
    if ( args_contains(argc, argv, "-F") ) {
        char *format_str = args_get_string(argc, argv, "-F");
        if ( strcasecmp(format_str, "fasta") == 0 ) {
            format = 0;
        }
        if ( strcasecmp(format_str, "nexus") == 0 ) {
            format = 1;
        }
        if ( strcasecmp(format_str, "phylip") == 0 ) {
            format = 2;
        }
        free(format_str);
    }
    
    if ( args_contains(argc, argv, "-A") ) {
        ancestral = true;
    }
    
    if ( args_contains(argc, argv, "-g") ) {
        bool success = false;
        gc = args_get_int(argc, argv, "-g", &success );
        if( !success || gc < 0 || gc > 14 ){
            error("Could not read genetic code [-g]");
        }
    }
    
    // define dataType
    if ( args_contains(argc, argv, "--states") ) {
        char ** states = String_to_string_array( args_get_string(argc, argv, "--states"), ',', &matrixDimension );
        dataType = new_GenericDataType(matrixDimension, (const char**)states);
        for ( int i = 0; i < matrixDimension; i++ ) {
            free(states[i]);
        }
        free(states);
        
        if( args_contains(argc, argv, "-m") ){
            char *temp = args_get_string(argc, argv, "-m");
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
                model_string = String_clone(temp);
            }
        }
        else{
            model_string = String_clone("ER");
        }
    }
    else if( args_contains(argc, argv, "-m") ){
        char *model_type = args_get_string(argc, argv, "-m");
        if( strcmp(model_type, "01234") == 0 ){
            strcpy(model_type, "GTR");
        }
        else if( strcmp(model_type, "00000") == 0 ){
            strcpy(model_type, "JC69");
        }
        else if( strcmp(model_type, "01001") == 0 ){
            strcpy(model_type, "HKY");
        }
        
        model_string = String_clone(model_type);
        free(model_type);
        
        if( array_of_string_contains(model_string, nucleotide_models, sizeof(nucleotide_models) / sizeof(nucleotide_models[0])) ){
            dataType = new_NucleotideDataType();
        }
        else if( array_of_string_contains(model_string, amino_acid_models, sizeof(amino_acid_models) / sizeof(amino_acid_models[0])) ){
            dataType = new_AminoAcidDataType();
            fprintf(stdout, "Genetic code: %s (%d)\n", GENETIC_CODE_NAMES[dataType->genetic_code], dataType->genetic_code );
        }
        else if( array_of_string_contains(model_string, codon_models, sizeof(codon_models) / sizeof(codon_models[0])) ){
            dataType = new_CodonDataType(gc);
        }
        else if(strlen(model_string) == 5 && model_string[0] == '0'){
            dataType = new_NucleotideDataType();
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
    
    
    SubstitutionModel *m = NULL;
    
    if ( dataType->type == DATA_TYPE_NUCLEOTIDE ) {
        
        if( strcasecmp("JC69", model_string) == 0 ){
            m = new_JC69();
        }
        else if( strcasecmp("K80", model_string) == 0 ){
            m = new_K80();
        }
        else if( strcasecmp("HKY", model_string) == 0 ){
            m = new_HKY();
        }
        else if( strcasecmp("GTR", model_string) == 0 ){
            m = new_GTR();
        }
        else if( strcasecmp("UREV", model_string) == 0 ){
            m = new_UnrestrictedNucleotideModel();
        }
        else if( strcasecmp("NONSTAT", model_string) == 0 ){
            m = new_NONSTATNucleotideModel();
        }
        else if( model_string[0] == '0' ){
            m = new_ReversibleNucleotideModel(model_string);
        }
        
        /*
         * FREQUENCIES
         */
        // Should only be specified for GTR and HKY but not NONSTAT, JC69 and K80
        if ( args_contains(argc, argv, "-f") ) {
            int nfreqs = 0;
            double *freqs = args_get_double_array(argc, argv, "-f", ',', &nfreqs);
            memcpy(m->_freqs, freqs, sizeof(double)*4);
            free(freqs);
        }
        
        /*
         * RATES
         */
        if ( args_contains(argc, argv, "-r") ) {
            int rateCount = 0;
            double *rates = args_get_double_array(argc, argv, "-r", ',', &rateCount);
            SubstitutionModel_set_rates(m,rates);
            free(rates);
        }
        
    }
    else if ( dataType->type == DATA_TYPE_CODON ) {
        // FREQUENCIES
        if ( args_contains(argc, argv, "-f") ) {
            int nfreqs = 0;
            double *freqs = args_get_double_array(argc, argv, "-f", ',', &nfreqs);
            memcpy(m->_freqs, freqs, sizeof(double)*nfreqs);
            free(freqs);
        }
        
        // RATES
        if( args_contains(argc, argv, "-r") ){
            int rateCount = 0;
            double *rates = args_get_double_array(argc, argv, "-r", ',', &rateCount);
            SubstitutionModel_set_rates(m,rates);
            free(rates);
        }
    }
    else if ( dataType->type == DATA_TYPE_AMINO_ACID ) {
        if( strcasecmp("DAYHOFF", model_string) == 0 ){
            m = new_DAYHOFF();
        }
        else if( strcasecmp("LG", model_string) == 0 ){
            m = new_LG();
        }
        else if( strcasecmp("WAG", model_string) == 0 ){
            m = new_WAG();
        }
        else {
            assert(0);
        }
        
        if ( args_contains(argc, argv, "-f") ) {
            int freqCount = 0;
            double *freqs = args_get_double_array(argc, argv, "-f", ',', &freqCount);
            memcpy(m->_freqs, freqs, sizeof(double)*20);
            free(freqs);
        }
    }
    else if ( dataType->type == DATA_TYPE_GENERIC) {
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
        
        m = new_GeneralModel(rateIndexes, matrixDimension);
        
        free(rateIndexes);
        
        if ( args_contains(argc, argv, "-f") ) {
            int freqCount = 0;
            double *freqs = args_get_double_array(argc, argv, "-f", ',', &freqCount);
            memcpy(m->_freqs, freqs, sizeof(double)*dataType->stateCount);
            free(freqs);
        }
        
        // special
//        mod->normalize = false;
//        for ( int i = 0; i < matrixDimension; i++ ) {
//            mod->_freqs[i] = 1.0;
//        }
//        mod->need_update = true;
        //end
        
        if ( args_contains(argc, argv, "-r") ) {
            int rateCount = 0;
            double *rates = args_get_double_array(argc, argv, "-r", ',', &rateCount);
            SubstitutionModel_set_rates(m, rates);
            free(rates);
        }
    }
    else {
        assert(0);
    }
    
    SiteModel *sm = NULL;
    
    if ( args_contains(argc, argv, "-c") ) {
        bool success = false;
        cat = args_get_int(argc, argv, "-c", &success );
        if( !success || cat <= 0 ){
            error("Could not read the number of categories [-c]");
        }
    }
    
    if ( args_contains(argc, argv, "-a") ) {
        bool success = false;
        shape = args_get_double(argc, argv, "-a", &success );
        if( !success || shape <= 0 ){
            error("Could not read alpha [-a]");
        }
    }
    if ( args_contains(argc, argv, "-I") ) {
        bool success = false;
        pinv = args_get_double(argc, argv, "-I", &success );
        if( !success || pinv < 0 ){
            error("Could not read alpha [-I]");
        }
    }
    
    if ( pinv == 0 && cat == 1 ){
        sm = new_SiteModel(m);
    }
    else if( pinv > 0 ){
        if( cat > 0 ){
            sm = new_GammaPinvSiteModel( m, pinv, shape, cat );
        }
        else {
            sm = new_PinvSiteModel(m, pinv);
        }
    }
    else {
        sm = new_GammaSiteModel(m, shape, cat);
    }
    
    
    
    
    
    if ( args_contains(argc, argv, "-l") ) {
        bool success = false;
        length = args_get_int(argc, argv, "-l", &success );
        if( !success || length <= 0 ){
            error("Could not read the length of the sequence [-l]");
        }
    }
    if ( args_contains(argc, argv, "-s") && args_contains(argc, argv, "-d") ) {
        error("Use either -s or -d\n");
        
    }
    
    if ( args_contains(argc, argv, "-s") ) {
        bool success = false;
        scaler = args_get_double(argc, argv, "-s", &success );
        if( !success || scaler <= 0 ){
            error("Could not read the tree length scaler [-s]");
        }
        
        //Tree_print_newick(stdout, tree, true);
        Tree_scale_distance(tree, scaler);
        //Tree_print_newick(stdout, tree, true);
    }
    
    if ( args_contains(argc, argv, "-R") ) {
        bool success = false;
        seed = args_get_int(argc, argv, "-R", &success );
    }
    init_genrand(seed);
    
    if ( args_contains(argc, argv, "-d") ) {
        
        char *rate_dist = args_get_string(argc, argv, "-d" );
        
        double *rates = dvector(Tree_node_count(tree)-1);
        Node **nodes = Tree_get_nodes(tree, POSTORDER);
        
        if ( strcasecmp(rate_dist, "tree") == 0 ) {
            for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
                rates[i] = Node_get_double_from_info(nodes[i], "rate=");
            }
            
            StringBuffer_append_format(buffer, "\nRates from tree file\n");
        }
        else {
            if ( strcasecmp(rate_dist, "logn") == 0 || strcasecmp(rate_dist, "exp") == 0 ) {

                char* mean_str = args_get_string(argc, argv, "-M");
                if( mean_str == NULL ) error("Could not read the mean of the lognormal distribution of categories [-M]");
                
                unsigned l = 0;
                double* means = String_to_double_array(mean_str, ',', &l);
                for (int i = 0; i < l; i++) {
                    if( means[i] <= 0 ) error("Mean cannot be negative [-M]");
                }
                free(mean_str);
                
                double* logsigmas = NULL;
                if(strcasecmp(rate_dist, "logn") == 0){
                    char* logsigma_str = args_get_string(argc, argv, "-S");
                    if( logsigma_str == NULL ) error("Could not read the standard deviation of the distribution of categories [-S]");
                    unsigned l2 = 0;
                    logsigmas = String_to_double_array(logsigma_str, ',', &l2);
                    for (int i = 0; i < l2; i++) {
                        if( logsigmas[i] <= 0 ) error("logsigma cannot be negative [-M]");
                    }
                    free(logsigma_str);
                    
                    if(l!=l2){
                        error("-M != -S");
                    }
                }
                
                if(l > 1){
                    int* counts = calloc(l, sizeof(int));
                    for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                        if( !Node_isroot(nodes[i]) ){
                            counts[Node_get_int_from_info(nodes[i], "class=")]++;
                        }
                    }
                    double* temp_rates = dvector(Tree_node_count(tree)-1);
        
                    for (int i = 0; i < l; i++) {
                        if(strcasecmp(rate_dist, "logn") == 0){
                            lognormal_discretize(log(means[i]), logsigmas[i], temp_rates, counts[i]);
                            StringBuffer_append_format(buffer, "\n%d Lognormal distribution mean = %f log(stdev) = %f", i, means[i], logsigmas[i]);
                        }
                        else{
                            exponential_discretize(means[i], temp_rates, counts[i]);
                            StringBuffer_append_format(buffer, "\n%d Exponential distribution mean = %f", i, means[i]);
                        }
                        randomize_dvector( temp_rates, counts[i]);
                        int k = 0;
                        for ( int j = 0; j < Tree_node_count(tree); j++ ) {
                            int class = Node_get_int_from_info(nodes[j], "class=");
                            if( !Node_isroot(nodes[j]) && class == i){
                                rates[j] = temp_rates[k];
                                k++;
                            }
                        }
                    }
                    StringBuffer_append_string(buffer, "\n");
                    free(temp_rates);
                    free(counts);
                }
                else{
                    if(strcasecmp(rate_dist, "logn") == 0){
                        lognormal_discretize(log(means[0]), logsigmas[0], rates, Tree_node_count(tree)-1);
                    }
                    else{
                        exponential_discretize(means[0], rates, Tree_node_count(tree)-1);
                    }
                    randomize_dvector( rates, Tree_node_count(tree)-1);
                }
                free(means);
                if(logsigmas != NULL) {
                    free(logsigmas);
                }
            }
            else if ( strcasecmp(rate_dist, "exp") == 0 ) {
                bool success = false;
                double mean = args_get_double(argc, argv, "-M", &success );
                if( !success || mean <= 0 ) error("Could not read the mean of the exponential distribution of categories [-M]");
                
                exponential_discretize(mean, rates, Tree_node_count(tree)-1);
                randomize_dvector( rates, Tree_node_count(tree)-1 );
                
                StringBuffer_append_format(buffer, "\nExponential distribution mean = %f\n", mean);
            }
            else {
                error("Coud not read the type of distribution -d exp or -d logn\n");
            }
            StringBuffer *info = new_StringBuffer(100);
            for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                Node_empty_annotation(nodes[i]);
                if( !Node_isroot(nodes[i]) ){
                    StringBuffer_empty(info);
                    StringBuffer_append_format(info, "%e", rates[i]);
                    Node_set_annotation(nodes[i], "rate", info->c);
                    
                    StringBuffer_empty(info);
                    int class = Node_get_int_from_info(nodes[i], "class=");
                    StringBuffer_append_format(info, "%d", class);
                    Node_set_annotation(nodes[i], "class", info->c);
                }
            }
            
            // setup heights
            for ( int i = 0; i < Tree_node_count(tree); i++ ) {
                double height = Node_get_double_from_info(nodes[i], "height=");
                Node_set_height(nodes[i], height);
                
                StringBuffer_empty(info);
                StringBuffer_append_format(info, "%f", height);
                Node_set_annotation(nodes[i], "height", info->c);
                
                char *cal = Node_get_string_from_info2(nodes[i], "cal_height=");
                if( cal != NULL ){
                    Node_set_annotation(nodes[i], "cal_height", cal);
                    free(cal);
                }
            }
            free_StringBuffer(info);
        }
        
        
        
        for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
            Node_set_distance(nodes[i], Node_distance(nodes[i])*rates[i]);
        }
        
        
        min = dmin_vector(rates, Tree_node_count(tree)-1);
        max = dmax_vector(rates, Tree_node_count(tree)-1);
        mean = dmean(rates, Tree_node_count(tree)-1);
        median = dmedian(rates, Tree_node_count(tree)-1);
        
        StringBuffer_append_format(buffer, "\nRate min = %f max = %f mean = %f median = %f\n", min,max,mean,median);
        
        free(rates);
        free(rate_dist);
    }
    
    
    StringBuffer_append_string(buffer, "\nFrequencies:\n");
    if( m->dtype == DATA_TYPE_NUCLEOTIDE ){
        bufferize_frequencies(buffer, m);
    }
    else if( m->dtype == DATA_TYPE_AMINO_ACID ){
        bufferize_aa_frequencies(buffer, m);
    }
    else if( m->dtype == DATA_TYPE_CODON ){
        bufferize_codon_frequencies(buffer, m);
    }
    
    StringBuffer_append_string(buffer, "\nRelative rates:\n");
    bufferize_rates(buffer, m);
    
    if( sm->pinv != NULL ) StringBuffer_append_format(buffer, "\nP-inv: %f\n", Parameter_value(sm->pinv));
    if( sm->shape != NULL ) StringBuffer_append_format(buffer, "\nGamma model: shape = %f (%d categories)\n", Parameter_value(sm->shape), sm->cat_count );
    
    Sequences *sim = Sequence_simulate(tree, sm, NULL, dataType, length, ancestral);
    
    if(sm->cat_count > 1 ){
        StringBuffer_append_string(buffer, "\nRate categories:\n");
        init_genrand(seed);
        for ( int i = 0; i < length; i++ ) {
            int rateClass = roulette_wheel(sm->get_proportions(sm), sm->cat_count);
            StringBuffer_append_format(buffer, "%d,", rateClass);
        }
        StringBuffer_chop(buffer);
        StringBuffer_append_char(buffer, '\n');
        
        StringBuffer_append_string(buffer, "Rates:\n");
        for ( int i = 0; i < sm->cat_count; i++ ) {
            StringBuffer_append_format(buffer, "%f,", sm->cat_rates[i] );
        }
        StringBuffer_chop(buffer);
        StringBuffer_append_char(buffer, '\n');
    }
    
    
    if(scaler != 1 ) StringBuffer_append_format(buffer, "\nBranch scaler: %f\n", scaler);
    StringBuffer_append_format(buffer, "\nRandom seed: %lu\n", seed);
    
    
    StringBuffer_append_char(buffer, ']');
    
    switch (format) {
        case 0:
            Sequences_save_fasta(sim, output);
            break;
        case 1:
            Sequences_save_nexus_with_comment(sim, output, buffer->c);
            break;
        case 3:
            Sequences_save_phylip(sim, output);
            break;
        default:
            assert(0);
    }
    
    if( outputtree != NULL ){
        FILE *testFile = fopen(outputtree,"w");
        fprintf(testFile, "#NEXUS\n");
        fprintf(testFile, "%s\n\n", buffer->c);
        Tree_print_nexus_taxa_block(testFile, tree);
        Tree_print_nexus_header_figtree_BeginTrees(testFile, tree);
        
        fprintf(testFile, "tree TREE = [&R] ");
        
        Tree_print_nexus_with_annotation(testFile, tree);
        fprintf(testFile, "\nEnd;");
        fclose(testFile);
    }
    
    int polymorphisms = 0;
    for ( int i = 0; i < sim->length; i++ ) {
        char c = sim->seqs[0]->seq[i];
        int j = 1;
        for ( ; j < sim->size; j++ ) {
            if( c != sim->seqs[j]->seq[i] ) break;
        }
        if( j != sim->size ){
            polymorphisms++;
        }
    }
    Node **nodes = Tree_get_nodes(tree, POSTORDER);
    double tree_length = 0;
    for ( int i = 0; i < Tree_node_count(tree)-1; i++ ) {
        tree_length += Node_distance(nodes[i]);
        //printf("%f %f\n", tree_length,Node_distance(nodes[i]));
    }
    printf("Tree length %f\n", tree_length);
    fprintf(stdout, "Number of polymorphic sites: %d/%d (%f)\n", polymorphisms, length,((double)polymorphisms/length) );
    
    free(freq_nuc);
    if ( freq_codon != NULL ) {
        free(freq_codon);
    }
    free(model_string);
    free_Sequences(sim);
    free(output);
    free_SiteModel(sm);
    free(rates);
    free_StringBuffer(buffer);
    free_Tree(tree);
    
    if( outputtree != NULL ){
        free(outputtree);
    }
    
    return 0;
}

