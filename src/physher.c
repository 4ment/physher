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

#ifndef DISABLED_CONFIG_HEADER
#include "phyc/PhyCConfig.h"
#endif

#include "phyc/parsimony.h"
#include "phyc/treelikelihood.h"
#include "phyc/treeio.h"
#include "phyc/geneticcode.h"
#include "phyc/random.h"
#include "phyc/args.h"
#include "phyc/filereader.h"
#include "phyc/compoundmodel.h"
#include "phyc/distmodel.h"
#include "phyc/mjson.h"
#include "phyc/logger.h"
#include "phyc/vb.h"
#include "phyc/mcmc.h"
#include "phyc/mmcmc.h"
#include "phyc/hessian.h"
#include "phyc/is.h"
#include "phyc/nest.h"
#include "phyc/cpo.h"
#include "phyc/laplace.h"
#include "phyc/bridge.h"
#include "phyc/marginal.h"
#include "phyc/predictive.h"
#include "phyc/mc.h"
#include "phyc/physim.h"


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

void create_json_distance(json_node* root, char* algorithm, char* seq_file, char* model_string, bool nucleotide){

	json_node* jtree = create_json_node_object(root, "model");
	add_json_node(root, jtree);

	add_json_node_string(jtree, "id", "tree");
	add_json_node_string(jtree, "type", "tree");

	json_node* jinittree = create_json_node_object(jtree, "init");
	add_json_node(jtree, jinittree);


	add_json_node_string(jinittree, "algorithm", algorithm);

	if ( model_string != NULL ) {
		if( strcasecmp(model_string, "K2P") == 0 ){
			add_json_node_string(jinittree, "model", "k2p");
		}
		else if( strcasecmp(model_string, "JC69") == 0 ){
			add_json_node_string(jinittree, "model", "jc69");
		}
		else if( strcasecmp(model_string, "K83") == 0 ){
			add_json_node_string(jinittree, "model", "k83");
		}
		else {
			add_json_node_string(jinittree, "model", "uncorrected");
		}
	}
	else {
		add_json_node_string(jinittree, "model", "uncorrected");
	}

	// sitepattern
	json_node* jsitepatterns = create_json_node_object(jinittree, "sitepattern");
	add_json_node(jinittree, jsitepatterns);
	add_json_node_string(jsitepatterns, "id", "patterns");
	add_json_node_string(jsitepatterns, "type", "sitepattern");
	add_json_node_string(jsitepatterns, "datatype", (nucleotide ? "nucleotide": "aa"));

	// alignment
	json_node* jsequences = create_json_node_object(jsitepatterns, "alignment");
	add_json_node(jsitepatterns, jsequences);
	add_json_node_string(jsequences, "id", "sequences");
	add_json_node_string(jsequences, "type", "alignment");
	add_json_node_string(jsequences, "file", seq_file);
	add_json_node_string(jsequences, "datatype", (nucleotide ? "nucleotide": "aa"));

	// physher
	json_node* jphysher = create_json_node_array(root, "physher");
	add_json_node(root, jphysher);

	// logger
	json_node* jlogger = create_json_node_object(jphysher, NULL);
	add_json_node(jphysher, jlogger);
	add_json_node_string(jlogger, "id", "log");
	add_json_node_string(jlogger, "type", "logger");
	add_json_node_string(jlogger, "tree", "&tree");
}

json_node* create_json_file(int argc, char* argv[]){
	json_node* jroot = create_json_node(NULL);
	jroot->node_type = MJSON_OBJECT;
	
	bool overwrite = false;
	char *seq_file  = NULL;
	char *tree_file = NULL;
	char *output_stem_user = NULL;
	
	char* clock = NULL;
	char* clock_algorithm = NULL;
	int nexus_index = -1;
	int genetic_code = 0;
	char* data_type = NULL;
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
	
	
	int approximation = 0;
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
	char* output_stem = NULL;
	char* model_string = NULL;
	StringBuffer *buffer = new_StringBuffer(100);
	StringBuffer_empty(buffer);
	
	struct argsparser_option options[] = {
		{ARGS_OPTION_STRING,  'i', "sequences",    "input.sequences", &seq_file, "Input alignment file"},
		{ARGS_OPTION_STRING,  't', "tree",         "input.tree", &tree_file, "Input tree file"},
		{ARGS_OPTION_STRING,  'o', "stem",         "output.stem", &output_stem_user, "Output stem file"},
		
		{ARGS_OPTION_INTEGER, 0,   "gc",           "sequences.geneticcode", &genetic_code, "Genetic Code"},
		{ARGS_OPTION_INTEGER, 0,   "dt",           "sequences.type", &data_type, "Data type (nucleotide or aa)"},
		
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
		
		{ARGS_OPTION_STRING,  'F', "fix",          "treelikelihood.fix", &fix, "Fix d: branch length, i: invariant, a: alpha, f: frequencies, r: rates"},
		{ARGS_OPTION_BOOLEAN, 0,   "sse",          "treelikelihood.sse", &use_sse, "Use SSE [default true]"},
		{ARGS_OPTION_BOOLEAN, 0,   "upper",        "treelikelihood.upper", &use_upper, "Use upper likelihood"},
		{ARGS_OPTION_INTEGER, 0,   "lk-threads",   "treelikelihood.threads", &tlk_threads, "Number of threads for likelihood calculation"},
		{ARGS_OPTION_INTEGER, 'A', "approx",       "treelikelihood.approximation", &approximation, "Input alignment file"},
		
		{ARGS_OPTION_FLAG,    'U', "unrooted",     "tree.unrooted", &tree_unrooted, "Input tree is rooted"},
		{ARGS_OPTION_DOUBLE,  's', "scaler",       "tree.scaler", &tree_scaler, "Scale input tree"},
		{ARGS_OPTION_STRING,  'O', "treeopt",      "tree.topolology.optimize", &topology_optimization_algorithm, "Optimize topology nni or spr (experimental)"},
		{ARGS_OPTION_BOOLEAN,  0,  "parsimony",    "tree.topolology.optimize.parsimony", &use_parsimony, "Quick optimizaton of tree topology using parsimony before ML optimization"},
		
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
		if (argc == 2){
			if(strcasecmp(argv[1], "-h") == 0 || strcasecmp(argv[1], "--help") == 0){
				help  = true;
			}
			else if(strcasecmp(argv[1], "--hh") == 0){
				help2 = true;
			}
		}
		
		if(!help && !help2){
			if(!file_exists(argv[argc-1])){
				fprintf(stderr, "Cannot read input file: %s\n", argv[argc-1]);
				exit(1);
			}
			argsparser_parse_file(argsparser, argv[argc-1]);
		}
		else{
			argsparser_parse(argsparser, argv, argc);
		}
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
		fprintf("Git:\nBranch: %s\nCommit: %s\n", GIT_BRANCH, GIT_COMMIT_HASH);
		fprintf(stdout, "\n\n");
#endif
		printf("\n");
		printf("Example\n\n");
		printf("%s -i alignment.fa -o stem -m GTR -c 4\n\n", program_name);
		exit(0);
	}
	
	if(distance != NULL){
		bool nucleotide = true;
		if (data_type != NULL && strcasecmp(data_type, "aa") == 0) {
			nucleotide = false;
		}
		create_json_distance(jroot, distance, seq_file, model_string_user, nucleotide);
		goto CLEANUP;
	}

	json_node* jinit = create_json_node_object(jroot, "init");
	add_json_node(jroot, jinit);

	if (fix != NULL) {
		alpha_fixed = String_contains(fix, 'a');
		pinv_fixed = String_contains(fix, 'i');
		bl_fixed = String_contains(fix, 'd');
		frequencies_fixed = String_contains(fix, 'f');
		rates_fixed = String_contains(fix, 'r');
	}
	
	if (output_stem_user != NULL) {
		output_stem = output_stem_user;
	}
	else {
		output_stem = String_clone(seq_file);
	}
	
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
	add_json_node_size_t(jinit, "seed", seed);
	
	json_node* jtlk = create_json_node_object(jroot, "treelikelihood");
	add_json_node(jroot, jtlk);
	
	json_node* jsitepatterns = create_json_node_object(jtlk, "sitepattern");
	json_node* jsitemodel = create_json_node_object(jtlk, "sitemodel");
	add_json_node(jtlk, jsitepatterns);
	add_json_node(jtlk, jsitemodel);
	
	if(frequencies_unknown){
		add_json_node_size_t(jtlk, "root_frequencies", 1);
	}
	
	add_json_node_string(jsitepatterns, "id", "patterns");
	add_json_node_string(jsitepatterns, "type", "sitepattern");
	//	add_json_node_string(jsitepatterns, "datatype", "nucleotide");
	
	json_node* jsequences = create_json_node_object(jsitepatterns, "sequences");
	add_json_node(jsitepatterns, jsequences);
	add_json_node_string(jsequences, "id", "sequences");
	add_json_node_string(jsequences, "type", "alignment");
	//	add_json_node_string(jsequences, "datatype", "nucleotide");
	add_json_node_string(jsequences, "file", seq_file);
	if( nexus_index >= 0 ){
		add_json_node_size_t(jsequences, "index", nexus_index);
	}
	
	const char* nucleotide_models[] = {"JC69", "K80", "F81", "HKY", "GTR", "UREV", "NONSTAT"};
	const char* amino_acid_models[] = {"DAYHOFF", "LG", "WAG"};
	const char* codon_models[]      = {"GY94", "MG94"};
	
	/*************************************************************************************************
	 ***************************************** Create DataType ***************************************
	 *************************************************************************************************/
	
	datatype dataType = DATA_TYPE_NUCLEOTIDE;
	
	// Determine the data type from model or number of states
	unsigned matrixDimension = 0;
	
	
	if( markov_states != NULL ){
		char ** states = String_to_string_array( markov_states, ',', &matrixDimension );
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
		dataType = DATA_TYPE_GENERIC;
	}
	else if( model_string_user != NULL ){
		model_string = String_clone(model_string_user);
		
		if(strlen(model_string) == 5 && model_string[0] == '0'){
			add_json_node_string(jsitepatterns, "datatype", "nucleotide");
			add_json_node_string(jsequences, "datatype", "nucleotide");
			matrixDimension = 4;
		}
		else if( array_of_string_contains(model_string, nucleotide_models, sizeof(nucleotide_models) / sizeof(nucleotide_models[0])) ){
			add_json_node_string(jsitepatterns, "datatype", "nucleotide");
			add_json_node_string(jsequences, "datatype", "nucleotide");
			matrixDimension = 4;
		}
		else if( array_of_string_contains(model_string, amino_acid_models, sizeof(amino_acid_models) / sizeof(amino_acid_models[0])) ){
			add_json_node_string(jsitepatterns, "datatype", "aa");
			add_json_node_string(jsequences, "datatype", "aa");
			dataType = DATA_TYPE_AMINO_ACID;
			matrixDimension = 20;
		}
		else if( array_of_string_contains(model_string, codon_models, sizeof(codon_models) / sizeof(codon_models[0])) ){
			add_json_node_string(jsitepatterns, "datatype", "codon");
			add_json_node_string(jsitepatterns, "code", GENETIC_CODE_NAMES[genetic_code]);
			add_json_node_string(jsequences, "datatype", "codon");
			matrixDimension = NUMBER_OF_CODONS[genetic_code];
			dataType = DATA_TYPE_CODON;
		}
		else {
			fprintf(stderr,"Does not recognize the model type (%s)\n", model_string);
			exit(1);
		}
	}
	else {
		error("No substitution model\n");
	}
	
	/*************************************************************************************************
	 ****************************************** Create tree ******************************************
	 *************************************************************************************************/
	
	json_node* jtree = create_json_node_object(jtlk, "tree");
	add_json_node(jtlk, jtree);
	
	add_json_node_string(jtree, "id", "tree");
	add_json_node_string(jtree, "type", "tree");
	add_json_node_string(jtree, "parameters", "tree.distances");
	
	if( tree_file != NULL ){
		char *treestring =  readTree( tree_file );
		add_json_node_string(jtree, "newick", treestring);
		free(treestring);
	}
	else {
		json_node* jinittree = create_json_node_object(jtree, "init");
		add_json_node(jtree, jinittree);
		add_json_node_string(jinittree, "sitepattern", "&patterns");
		add_json_node_string(jinittree, "parameters", "tree.distances");
		
		if( dataType == DATA_TYPE_AMINO_ACID ){
			add_json_node_string(jinittree, "model", "kimura");
		}
		else{
			add_json_node_string(jinittree, "model", "uncorrected");
		}
	}
	
	if ( tree_scaler > 0 ) {
		add_json_node_double(jtree, "scaler", tree_scaler);
	}
	
	//	check_aln_tree(tree, sp);
	
	/*************************************************************************************************
	 *********************************** Create substitution model ***********************************
	 *************************************************************************************************/
	
	json_node* jsubstmodel = create_json_node_object(jsitemodel, "substitutionmodel");
	add_json_node(jsitemodel, jsubstmodel);
	
	add_json_node_string(jsubstmodel, "id", "substitutionmodel");
	add_json_node_string(jsubstmodel, "type", "substitutionmodel");
	json_node* pp = get_json_node(jsitepatterns, "datatype");
	// Could be a ref or a native type such as nucleotide
	if(pp->node_type == MJSON_STRING){
		add_json_node_string(jsubstmodel, "datatype", pp->value);
	}
	
	double *frequencies = NULL;
	bool equal_frequencies = (frequencies_string_user != NULL && strlen(frequencies_string_user) == 1 && tolower(frequencies_string_user[0]) == 'e');
	
	if(frequencies_string_user != NULL){
		unsigned nfreqs = 0;
		frequencies = String_to_double_array( frequencies_string_user, ',', &nfreqs);
		assert(nfreqs==matrixDimension);
	}
	else{
		frequencies = dvector(matrixDimension);
		for (int i = 0; i < matrixDimension; i ++) {
			frequencies[i] = 1.0/matrixDimension;
		}
	}
	size_t rateCount = 0;
	
	if ( dataType == DATA_TYPE_NUCLEOTIDE ) {
		// custom model
		if( model_string[0] == '0' ){
			if(strcasecmp("00000", model_string) == 0 && equal_frequencies){
				free(model_string);
				model_string = String_clone("JC69");
			}
			else if( strcasecmp("01001", model_string) == 0 && equal_frequencies){
				free(model_string);
				model_string = String_clone("K80");
			}
			else if(strcasecmp("00000", model_string) == 0){
				free(model_string);
				model_string = String_clone("F81");
			}
			else if( strcasecmp("01001", model_string) == 0){
				free(model_string);
				model_string = String_clone("HKY");
			}
			else if( strcasecmp("01234", model_string) == 0){
				free(model_string);
				model_string = String_clone("GTR");
			}
		}
		
		if( strcasecmp("GTR", model_string) == 0){
			rateCount = 5;
		}
		else if( strcasecmp("HKY", model_string) == 0 || strcasecmp("K80", model_string) == 0){
			rateCount = 1;
		}
		
		add_json_node_string(jsubstmodel, "model", model_string);
		
		if(strcasecmp("NONSTAT", model_string) != 0){
			json_node* jfreqs = create_json_node_object(jsubstmodel, "frequencies");
			add_json_node(jsubstmodel, jfreqs);
			
			add_json_node_string(jfreqs, "id", "freqs");
			add_json_node_string(jfreqs, "type", "simplex");
			add_json_node_array_double(jfreqs, "values", frequencies, 4);
			
			if(!equal_frequencies){
				//TODO: init with this empirical_frequencies(seqs, freqs);
			}
		}
		else{
			rateCount = 11;
		}
		
		/*
		 * RATES
		 */
		if(rates_user != NULL || rateCount > 0){
			json_node* jrates = create_json_node_object(jsubstmodel, "rates");
			add_json_node(jsubstmodel, jrates);
			
			double *rates = NULL;
			char** names = NULL;
			
			if ( rates_user != NULL ) {
				unsigned rcount = 0;
				rates = String_to_double_array( rates_user, ',', &rcount);
				if( rateCount == 1 ){
					fprintf(stdout, "Kappa (user): %f\n", *rates);
				}
				names = (char**)malloc(sizeof(char*)*rateCount);
				StringBuffer* str = new_StringBuffer(10);
				for (int i = 0; i < rateCount; i++) {
					StringBuffer_empty(str);
					StringBuffer_append_format(str, "r%d", i);
					names[i] = StringBuffer_tochar(str);
				}
				free_StringBuffer(str);
			}
			else if( rateCount == 5 ){
				rateCount = 5;
				rates = dvector(rateCount);
				names = (char**)malloc(sizeof(char*)*rateCount);
				char* params[5] = {"ac", "ag", "at", "cg", "ct"};
				for (int i = 0; i < rateCount; i++) {
					names[i] = String_clone(params[i]);
					rates[i] = 1;
				}
			}
			else if( rateCount == 1 ){
				rateCount = 1;
				rates = dvector(rateCount);
				names = (char**)malloc(sizeof(char*));
				names[0] = String_clone("kappa");
				rates[0] = 1;
			}
			
			for (int i = 0; i < rateCount; i++) {
				json_node* jrate = create_json_node_parameter(jrates, names[i], rates[i], 0, INFINITY);
				add_json_node(jrates, jrate);
				free(names[i]);
			}
			free(names);
			free(rates);
		}
	}
	else if ( dataType == DATA_TYPE_CODON ) {
		add_json_node_string(jsubstmodel, "model", model_string);
		add_json_node_string(jsitepatterns, "code", GENETIC_CODE_NAMES[genetic_code]);
		matrixDimension = NUMBER_OF_CODONS[genetic_code];
		double* frequencies = dvector(matrixDimension);
		
		json_node* jfreqs = create_json_node_object(jsubstmodel, "frequencies");
		add_json_node(jsubstmodel, jfreqs);
		
		add_json_node_string(jfreqs, "id", "freqs");
		add_json_node_string(jfreqs, "type", "simplex");
		add_json_node_array_double(jfreqs, "values", frequencies, matrixDimension);
		free(frequencies);
	}
	else if ( dataType == DATA_TYPE_AMINO_ACID ) {
		add_json_node_string(jsubstmodel, "model", model_string);
		
		if ( frequencies_string_user != NULL) {
			unsigned nfreqs = 0;
			frequencies = String_to_double_array( frequencies_string_user, ',', &nfreqs);
			
			json_node* jfreqs = create_json_node_object(jsubstmodel, "frequencies");
			add_json_node(jsubstmodel, jfreqs);
			
			add_json_node_string(jfreqs, "id", "freqs");
			add_json_node_string(jfreqs, "type", "simplex");
			add_json_node_array_double(jfreqs, "values", frequencies, matrixDimension);
		}
		rates_fixed = true;
		frequencies_fixed = true;
	}
	else if ( dataType == DATA_TYPE_GENERIC) {
		// empirical
		unsigned *rateIndexes = NULL;
		
		if(frequencies_string_user != NULL){
			unsigned nfreqs = 0;
			frequencies = String_to_double_array( frequencies_string_user, ',', &nfreqs);
			assert(nfreqs==matrixDimension);
		}
		else{
			frequencies = dvector(matrixDimension);
			for (int i = 0; i < matrixDimension; i ++) {
				frequencies[i] = 1.0/matrixDimension;
			}
		}
		
		json_node* jfreqs = create_json_node_object(jsubstmodel, "frequencies");
		add_json_node(jsubstmodel, jfreqs);
		
		add_json_node_string(jfreqs, "id", "freqs");
		add_json_node_string(jfreqs, "type", "simplex");
		add_json_node_array_double(jfreqs, "values", frequencies, matrixDimension);
		
		rateCount = 1;
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
			rateCount = count;
		}
		else {
			unsigned len = 0;
			rateIndexes = String_to_uint_array(model_string, ',', &len);
			if( matrixDimension*matrixDimension != len ){
				fprintf(stderr, "Number of states (%d) does not match the vector length (%d) \n", matrixDimension, len);
				exit(1);
			}
			rateCount = umax_vector(rateIndexes, len) + 1;
		}
		
		json_node* discrete_node = create_json_node_object(jsubstmodel, "model");
		add_json_node(jsubstmodel, discrete_node);
		add_json_node_array_unsigned(discrete_node, "values", rateIndexes, matrixDimension*matrixDimension);
		if(!normalize_q){
			add_json_node_bool(jsubstmodel, "normalize", normalize_q);
		}
		free(rateIndexes);
		
		json_node* jrates = create_json_node_object(jsubstmodel, "rates");
		add_json_node(jsubstmodel, jrates);
		
		if ( rates_user != NULL ) {
			unsigned len = 0;
			double *rates = String_to_double_array(rates_user, ',', &len);
			for (int i = 0; i < len; i++) {
				json_node* jrate = create_json_node_parameter(jrates, "a", rates[i], 0, INFINITY);
				add_json_node(jrates, jrate);
			}
			free(rates);
		}
		else{
			for (int i = 0; i < rateCount; i++) {
				json_node* jrate = create_json_node_parameter(jrates, "a", 1.0, 0, INFINITY);
				add_json_node(jrates, jrate);
			}
		}
	}
	else {
		assert(0);
	}
	
	if(frequencies != NULL) free(frequencies);
	
	/*************************************************************************************************
	 **************************************** Create site model
	 *************************************************************************************************/
	
	add_json_node_string(jsitemodel, "id", "sitemodel");
	add_json_node_string(jsitemodel, "type", "sitemodel");
	
	if ( pinv > 0 || rate_category_count > 1 ) {
		json_node* jrates = create_json_node_object(jsitemodel, "rates");
		add_json_node(jsitemodel, jrates);
		if (pinv > 0) {
			json_node* jrate = create_json_node_object(jrates, "invariant");
			add_json_node(jrates, jrate);
			add_json_node_string(jrate, "id", "pinv");
			add_json_node_string(jrate, "type", "parameter");
			add_json_node_double(jrate, "value", pinv);
			add_json_node_double(jrate, "lower", 0);
			add_json_node_string(jrate, "upper", "infinity");
		}
		if(rate_category_count > 1){
			//			add_json_node_string(jsitemodel, "distribution", "gamma");
			add_json_node_size_t(jsitemodel, "categories", rate_category_count);
			if ( sitemodel_string != NULL && strcasecmp(sitemodel_string, "gammaquad") == 0 ) {
				add_json_node_string(jsitemodel, "discretization", "laguerre");
			}
			json_node* jrate = create_json_node_object(jrates, "alpha");
			add_json_node(jrates, jrate);
			add_json_node_string(jrate, "id", "alpha");
			add_json_node_string(jrate, "type", "parameter");
			add_json_node_double(jrate, "value", alpha);
			add_json_node_double(jrate, "lower", 0);
			add_json_node_string(jrate, "upper", "infinity");
		}
		//TODO: init
		//		if(pinv == 0.0 && use_pinv){
		//			pinv = fmax(0.01, 1.0 - ((double)polymorphisms/sp->nsites));
		//		}
	}
	
	/*************************************************************************************************
	 **************************************** Compute rate-free likelihood
	 *************************************************************************************************/
	
	json_node* jphysher = create_json_node_array(jroot, "physher");
	add_json_node(jroot, jphysher);
	
	if (use_parsimony) {
		json_node* jopt_topo_parsimony = create_json_node_object(jphysher, NULL);
		add_json_node(jphysher, jopt_topo_parsimony);
		add_json_node_string(jopt_topo_parsimony, "id", "optopopars");
		add_json_node_string(jopt_topo_parsimony, "type", "optimizer");
		add_json_node_string(jopt_topo_parsimony, "algorithm", "topology");
		add_json_node_string(jopt_topo_parsimony, "model", "&treelikelihood");
		add_json_node_string(jopt_topo_parsimony, "move", "spr");
	}
	
	json_node* joptmeta = create_json_node_object(jphysher, NULL);
	add_json_node(jphysher, joptmeta);
	
	add_json_node_string(joptmeta, "id", "metaopt");
	add_json_node_string(joptmeta, "type", "optimizer");
	add_json_node_string(joptmeta, "algorithm", "meta");
	add_json_node_string(joptmeta, "model", "&treelikelihood");
	add_json_node_double(joptmeta, "min", 1);
	add_json_node_double(joptmeta, "max", 10000);
	add_json_node_double(joptmeta, "precision", 0.001);
	
	json_node* joptmetalist = create_json_node_array(joptmeta, "list");
	add_json_node(joptmeta, joptmetalist);
	
	if(!bl_fixed){
		json_node* jopt_bl = create_json_node_object(joptmetalist, NULL);
		add_json_node(joptmetalist, jopt_bl);
		add_json_node_string(jopt_bl, "id", "optbl");
		add_json_node_string(jopt_bl, "type", "optimizer");
		add_json_node_string(jopt_bl, "algorithm", "serial");
		add_json_node_string(jopt_bl, "model", "&treelikelihood");
		add_json_node_string(jopt_bl, "treelikelihood", "&treelikelihood");
	}
	
	if(!frequencies_fixed || !rates_fixed){
		json_node* jopt_mat = create_json_node_object(joptmetalist, NULL);
		add_json_node(joptmetalist, jopt_mat);
		add_json_node_string(jopt_mat, "id", "optmat");
		add_json_node_string(jopt_mat, "type", "optimizer");
		add_json_node_string(jopt_mat, "algorithm", "serial");
		add_json_node_string(jopt_mat, "model", "&treelikelihood");
		add_json_node_string(jopt_mat, "parameters", "[\"$frequencies\"]");
	}
	
	if(!pinv_fixed){
		json_node* jopt_pinv = create_json_node_object(joptmetalist, NULL);
		add_json_node(joptmetalist, jopt_pinv);
		add_json_node_string(jopt_pinv, "id", "optpinv");
		add_json_node_string(jopt_pinv, "type", "optimizer");
		add_json_node_string(jopt_pinv, "algorithm", "brent");
		add_json_node_string(jopt_pinv, "model", "&treelikelihood");
		add_json_node_string(jopt_pinv, "parameters", "&invariant");
	}
	
	if(!alpha_fixed){
		json_node* jopt_alpha = create_json_node_object(joptmetalist, NULL);
		add_json_node(joptmetalist, jopt_alpha);
		add_json_node_string(jopt_alpha, "id", "optalpha");
		add_json_node_string(jopt_alpha, "type", "optimizer");
		add_json_node_string(jopt_alpha, "algorithm", "brent");
		add_json_node_string(jopt_alpha, "model", "&treelikelihood");
		add_json_node_string(jopt_alpha, "parameters", "&alpha");
	}
	
	if ( topology_optimization_algorithm != NULL ) {
		json_node* jopt_topo = create_json_node_object(joptmetalist, NULL);
		add_json_node(joptmetalist, jopt_topo);
		add_json_node_string(jopt_topo, "id", "optopo");
		add_json_node_string(jopt_topo, "type", "optimizer");
		add_json_node_string(jopt_topo, "algorithm", "topology");
		add_json_node_string(jopt_topo, "model", "&treelikelihood");
		add_json_node_string(jopt_topo, "treelikelihood", "&treelikelihood");
		add_json_node_string(jopt_topo, "move", topology_optimization_algorithm);
		add_json_node_size_t(jopt_topo, "threads", nthreads);
	}
	
CLEANUP:

	json_tree_print(jroot);
	
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
	if(output_stem)free(output_stem);
	if(ic)free(ic);
	
	if(model_string)free(model_string);
	
	if(buffer)free_StringBuffer(buffer);
	return jroot;
}

int main(int argc, char* argv[]){
	time_t start_time;
	time_t beginning_of_time;
	time_t end_time;
	double diff_time;
	
	time(&start_time);
	beginning_of_time = start_time;
	
	json_node* json = NULL;

    if (argc == 1) {
        fprintf(stdout, "\n%s:\n\n", get_program_name(argv));
        printf("Fourment M and Holmes EC. Novel non-parametric models to estimate evolutionary rates and divergence times from heterochronous sequence data.\n");
        printf("BMC Evolutionary Biology 14:163, 2014\n\n");
		
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
        fprintf(stdout, "Git:\nBranch: %s\nCommit: %s\n", GIT_BRANCH, GIT_COMMIT_HASH);
        fprintf(stdout, "\n\n");
#endif
		return 0;
	}
	else if(argc > 2){
		json = create_json_file(argc, argv);
	}
	else{
		char* content = load_file(argv[1]);
		printf("Reading file %s\n", argv[1]);
		printf("done\n\n");

		json = create_json_tree(content);
		free(content);
	}

	long seeed = time(NULL);
	
	
	Hashtable* hash2 = new_Hashtable_string(100);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	json_node* run_node = get_json_node(json, "physher");
	json_node* init_node = get_json_node(json, "init");
	
	if (init_node != NULL) {
		seeed = get_json_node_value_size_t(init_node, "seed", seeed);
	}

	printf("seed: %ld\n", seeed);
	init_genrand(seeed);
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(r, seeed);
	char* rand_key = "RANDOM_GENERATOR!@";
	Hashtable_add(hash2, rand_key, r);
	
	size_t model_count = json->child_count;
	if (run_node != NULL) model_count--; // maybe physher
	if (init_node != NULL) model_count--; // maybe init
	
	Model** models = NULL;
	if(model_count > 0){
		models = calloc(model_count, sizeof(Model*));
		size_t index = 0;
		for (int i = 0; i < json->child_count; i++) {
			json_node* child = json->children[i];
			if (strcasecmp(child->key, "physher") == 0 || strcasecmp(child->key, "init") == 0) continue;
			
			json_node* type_node = get_json_node(child, "type");
			char* id = get_json_node_value_string(child, "id");
			
			if (strcasecmp((char*)type_node->value, "compound") == 0) {
				models[index] = new_CompoundModel_from_json(child, hash2);
			}
			else if(strcasecmp((char*)type_node->value, "variational") == 0){
				models[index] = new_Variational_from_json(child, hash2);
			}
			else if(strcasecmp((char*)type_node->value, "distribution") == 0){
				models[index] = new_DistributionModel_from_json(child, hash2);
			}
			else if(strcasecmp((char*)type_node->value, "treelikelihood") == 0){
				models[index] = new_TreeLikelihoodModel_from_json(child, hash2);
			}
			else if(strcasecmp((char*)type_node->value, "parsimony") == 0){
				models[index] = new_ParsimonyModel_from_json(child, hash2);
			}
			else if(strcasecmp((char*)type_node->value, "tree") == 0){
				models[index] = new_TreeModel_from_json(child, hash2);
			}
			
			Hashtable_add(hash2, id, models[index]);
			index++;
		}
	}
	
	if(run_node != NULL)
	for (int i = 0; i < run_node->child_count; i++) {
        json_node* child = run_node->children[i];
        char* type = get_json_node_value_string(child, "type");
        bool ignore = get_json_node_value_bool(child, "ignore", false);
        
        if(ignore) continue;
        
		if (strcasecmp(type, "optimizer") == 0) {
			Optimizer* opt = new_Optimizer_from_json(child, hash2);
			double logP;
			opt_optimize(opt, NULL, &logP);
			free_Optimizer(opt);
		}
		else if (strcasecmp(type, "logger") == 0) {
			struct Logger* logger = new_logger_from_json(child, hash2);
			logger->log(logger);
			free_Logger(logger);
		}
		else if(strcasecmp(type, "mcmc") == 0){
			MCMC* mcmc = new_MCMC_from_json(child, hash2);
			mcmc->run(mcmc);
			mcmc->free(mcmc);
		}
		else if(strcasecmp(type, "mmcmc") == 0){
			MMCMC* mmcmc = new_MMCMC_from_json(child, hash2);
			mmcmc->run(mmcmc);
			mmcmc->free(mmcmc);
		}
		else if(strcasecmp(type, "hessian") == 0){
			Hessian* hessian = new_Hessian_from_json(child, hash2);
			hessian->calculate(hessian);
			print_hessian(hessian);
			hessian->free(hessian);
		}
		else if(strcasecmp(type, "cpo") == 0){
			CPO* cpo = new_CPO_from_json(child, hash2);
			cpo->calculate(cpo);
			cpo->free(cpo);
		}
		else if(strcasecmp(type, "vbis") == 0 || strcasecmp(type, "is") == 0){
			ImportanceSampler* mvb = new_ImportanceSampler_from_json(child, hash2);
			printf("Marginal likelihood using IS: %f\n", mvb->calculate(mvb));
			mvb->free(mvb);
		}
		else if(strcasecmp(type, "nest") == 0){
			NEST* nest = new_NEST_from_json(child, hash2);
			nest->run(nest);
			free_NEST(nest);
		}
		else if(strcasecmp(type, "laplace") == 0){
			Laplace* laplace = new_Laplace_from_json(child, hash2);
			laplace->calculate(laplace);
			laplace->free(laplace);
		}
		else if(strcasecmp(type, "bridgesampling") == 0){
			BridgeSampling* bridge = new_BridgeSampling_from_json(child, hash2);
			bridge->run(bridge);
			bridge->free(bridge);
		}
		else if(strcasecmp(type, "marginallikelihood") == 0){
			MarginaLikelihood* margl = new_MarginaLikelihood_from_json(child, hash2);
			margl->run(margl);
			margl->free(margl);
		}
		else if(strcasecmp(type, "mc") == 0){
			MC* mc = new_MonteCarlo_from_json(child, hash2);
			mc->calculate(mc);
			mc->free(mc);
		}
		else if(strcasecmp(type, JSON_PREDICTIVE) == 0){
			Predictive* predictive = new_Predictive_from_json(child, hash2);
			predictive->calculate(predictive);
			predictive->free(predictive);
		}
		else if(strcasecmp(type, JSON_SIMULTRON) == 0){
			printf("Simulating sequences...\n");
			SimulateSequences_from_json(child, hash2);
		}
	}
	if(models != NULL){
		for (int i = 0; i < model_count; i++) {
			models[i]->free(models[i]);
		}
		free(models);
	}
	free_Hashtable(hash2);
	
	//json_tree_to_string(json);
	json_free_tree(json);
	gsl_rng_free(r);
	
	time(&end_time);
	diff_time = difftime(end_time, start_time);
	fprintf(stdout, "\nTotal runtime ");
	print_pretty_time(stdout, diff_time);
	return 0;
}
