//
//  physhercmd.c
//  physher
//
//  Created by Mathieu Fourment on 5/12/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "physhercmd.h"

#include <time.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <ctype.h>

#include "hashtable.h"
#include "filereader.h"
#include "mjson.h"
#include "args.h"
#include "random.h"
#include "matrix.h"
#include "utils.h"

#include "datatype.h"
#include "geneticcode.h"
#include "treeio.h"


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

json_node* create_json_distance(Hashtable* options){
	
	char* algorithm = Hashtable_get(options, "distance");
	char* seq_file = Hashtable_get(options, "sequences");
	char* model_string = Hashtable_get(options, "model");
	char* datatype = Hashtable_get(options, "datatype");
	bool nucleotide = true;
	if (strlen(datatype) != 0 && strcasecmp(datatype, "aa") == 0) {
		nucleotide = false;
	}
	
	json_node* jroot = create_json_node(NULL);
	jroot->node_type = MJSON_OBJECT;
	
	json_node* jtree = create_json_node_object(jroot, "model");
	add_json_node(jroot, jtree);
	
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
	json_node* jphysher = create_json_node_array(jroot, "physher");
	add_json_node(jroot, jphysher);
	
	// logger
	json_node* jlogger = create_json_node_object(jphysher, NULL);
	add_json_node(jphysher, jlogger);
	add_json_node_string(jlogger, "id", "log");
	add_json_node_string(jlogger, "type", "logger");
	add_json_node_string(jlogger, "tree", "&tree");
	
	return jroot;
}

void create_json_datatype(Hashtable* options, Hashtable* nodes){
	
	const char* nucleotide_models[] = {"JC69", "K80", "F81", "HKY", "SYM", "GTR", "UREV", "NONSTAT"};
	const char* amino_acid_models[] = {"DAYHOFF", "LG", "WAG"};
	const char* codon_models[]      = {"GY94", "MG94"};
	
	char* markov_states = Hashtable_get(options, "states");
	char* model_string = Hashtable_get(options, "model");
	
	json_node* jdatatype = NULL;
	
	if( strlen(markov_states) != 0 ){
		unsigned matrixDimension;
		char ** states = String_to_string_array( markov_states, ',', &matrixDimension );
		jdatatype = create_json_node_object(NULL, "datatype");
		add_json_node_string(jdatatype, "id", "datatype");
		add_json_node_string(jdatatype, "type", "datatype");
		json_node* jstates = create_json_node(jdatatype);
		add_json_node(jdatatype, jstates);
		jstates->node_type = MJSON_ARRAY;
		jstates->key = String_clone("states");
		
		for ( int i = 0; i < matrixDimension; i++ ) {
			add_json_node_string(jstates, NULL, states[i]);
			free(states[i]);
		}
		free(states);
	}
	else if( strlen(model_string) != 0 ){
		if( (strlen(model_string) == 5 && model_string[0] == '0') ||
		   array_of_string_contains(model_string, nucleotide_models, sizeof(nucleotide_models) / sizeof(nucleotide_models[0]), false) ){
			jdatatype = create_json_node(NULL);
			jdatatype->key = String_clone("datatype");
			jdatatype->value = String_clone("nucleotide");
			jdatatype->node_type = MJSON_STRING;
		}
		else if( array_of_string_contains(model_string, amino_acid_models, sizeof(amino_acid_models) / sizeof(amino_acid_models[0]), false) ){
			jdatatype = create_json_node(NULL);
			jdatatype->key = String_clone("datatype");
			jdatatype->value = String_clone("aa");
			jdatatype->node_type = MJSON_STRING;
		}
		else if( array_of_string_contains(model_string, codon_models, sizeof(codon_models) / sizeof(codon_models[0]), false) ){
			jdatatype = create_json_node(NULL);
			jdatatype->key = String_clone("datatype");
			jdatatype->value = String_clone("codon");
			jdatatype->node_type = MJSON_STRING;
		}
		else {
			fprintf(stderr,"Does not recognize the model (%s)\n", model_string);
			exit(1);
		}
	}
	else {
		error("No substitution model\n");
	}
	
	Hashtable_add(nodes, "jdatatype", jdatatype);
}

void create_json_phylo_tree(Hashtable* options, Hashtable* nodes){
	
	json_node* jtlk = Hashtable_get(nodes, "jtlk");
	json_node* jpatterns = Hashtable_get(nodes, "jsitepatterns");
	json_node* joptmetalist = Hashtable_get(nodes, "joptmetalist");
	
	json_node* jtree = create_json_node_object(jtlk, "tree");
	add_json_node(jtlk, jtree);
	
	add_json_node_string(jtree, "id", "tree");
	add_json_node_string(jtree, "type", "tree");
	add_json_node_string(jtree, "parameters", "tree.distances");
	
	char* tree_file = Hashtable_get(options, "tree");
	if( strlen(tree_file) != 0 ){
		char *treestring =  readTree( tree_file );
		add_json_node_string(jtree, "newick", treestring);
		free(treestring);
	}
	else {
		json_node* jinittree = create_json_node_object(jtree, "init");
		add_json_node(jtree, jinittree);
		add_json_node_string(jinittree, "sitepattern", "&patterns");
		add_json_node_string(jinittree, "parameters", "tree.distances");
		add_json_node_string(jinittree, "algorithm", "NJ");
		
		json_node* jdatatype = get_json_node(jpatterns, "datatype");
		
		if( jdatatype->node_type == MJSON_STRING && strcasecmp(jdatatype->value, "aa") == 0){
			add_json_node_string(jinittree, "model", "kimura");
		}
		else{
			add_json_node_string(jinittree, "model", "uncorrected");
		}
	}
	
	int tree_scaler = atoi(Hashtable_get(options, "scaler"));
	if ( tree_scaler > 0 ) {
		add_json_node_double(jtree, "scaler", tree_scaler);
	}
	
	bool bl_fixed = false;
	if (strlen(Hashtable_get(options, "fix")) != 0) {
		char* fix = Hashtable_get(options, "fix");
		bl_fixed = String_contains(fix, 'd');
	}
	
	if(!bl_fixed){
		json_node* jopt_bl = create_json_node_object(joptmetalist, NULL);
		add_json_node(joptmetalist, jopt_bl);
		add_json_node_string(jopt_bl, "id", "optbl");
		add_json_node_string(jopt_bl, "type", "optimizer");
		add_json_node_string(jopt_bl, "algorithm", "serial");
		add_json_node_string(jopt_bl, "model", "&treelikelihood");
		add_json_node_string(jopt_bl, "treelikelihood", "&treelikelihood");
	}
	
	if (strlen(Hashtable_get(options, "treeopt")) > 0) {
		json_node* jopt_topo = create_json_node_object(joptmetalist, NULL);
		add_json_node(joptmetalist, jopt_topo);
		add_json_node_string(jopt_topo, "id", "opttopo");
		add_json_node_string(jopt_topo, "type", "optimizer");
		add_json_node_string(jopt_topo, "algorithm", "topology");
		add_json_node_string(jopt_topo, "model", "&treelikelihood");
		add_json_node_string(jopt_topo, "move", Hashtable_get(options, "treeopt"));
	}
	
	Hashtable_add(nodes, "jtree", jtree);
}

void create_json_substitution_model(Hashtable* options, Hashtable* nodes){
	json_node* jsitemodel = Hashtable_get(nodes, "jsitemodel");
	json_node* jsitepatterns = Hashtable_get(nodes, "jsitepatterns");
	json_node* joptmetalist = Hashtable_get(nodes, "joptmetalist");
	
	json_node* jdatatype = get_json_node(jsitepatterns, "datatype");
	
	
	json_node* jopt_mat = create_json_node_object(joptmetalist, NULL);
	add_json_node_string(jopt_mat, "id", "optmat");
	add_json_node_string(jopt_mat, "type", "optimizer");
	add_json_node_string(jopt_mat, "algorithm", "serial");
	add_json_node_string(jopt_mat, "model", "&treelikelihood");
	
	json_node* jparams = create_json_node(jopt_mat);
	add_json_node(jopt_mat, jparams);
	jparams->node_type = MJSON_ARRAY;
	jparams->key = String_clone("parameters");
	
	char* fix = Hashtable_get(options, "fix");
	bool frequencies_fixed = String_contains(fix, 'f');
	bool rates_fixed = String_contains(fix, 'r');
	
	char* datatype = NULL;
	int genetic_code = atoi(Hashtable_get(options, "genetic-code")); // 0 is default
	
	json_node* jsubstmodel = create_json_node_object(jsitemodel, "substitutionmodel");
	add_json_node(jsitemodel, jsubstmodel);
	
	add_json_node_string(jsubstmodel, "id", "substitutionmodel");
	add_json_node_string(jsubstmodel, "type", "substitutionmodel");
	
	int matrixDimension = 4;
	enum datatype dt;
	
	if(jdatatype->node_type == MJSON_STRING){
		datatype = (char*)jdatatype->value;
		add_json_node_string(jsubstmodel, "datatype", datatype);
		if (strcasecmp(datatype, "aa") == 0) {
			matrixDimension = 20;
			dt = DATA_TYPE_AMINO_ACID;
		}
		else if (strcasecmp(datatype, "codon") == 0) {
			matrixDimension = NUMBER_OF_CODONS[genetic_code];
			dt = DATA_TYPE_CODON;
		}
		else if (strcasecmp(datatype, "nucleotide") == 0) {
			matrixDimension = 4;
			dt = DATA_TYPE_NUCLEOTIDE;
		}
		else{
			error(" create_json_substitution_model not done yet\n");
		}
	}
	else{
		dt = DATA_TYPE_GENERIC;
		json_node* jstates = get_json_node(jdatatype, "states");
		matrixDimension = jstates->child_count;
	}
	
	char* frequencies_string_user = Hashtable_get(options, "frequencies");
	char* rates_user = Hashtable_get(options, "rates");
	
	bool equal_frequencies = (strlen(frequencies_string_user) != 0 && strlen(frequencies_string_user) == 1 && tolower(frequencies_string_user[0]) == 'e');

	double *frequencies = NULL;
	if(strlen(frequencies_string_user) != 0){
		unsigned nfreqs = 0;
		frequencies = String_to_double_array( frequencies_string_user, ',', &nfreqs);
	}
	else{
		frequencies = dvector(matrixDimension);
		for (int i = 0; i < matrixDimension; i ++) {
			frequencies[i] = 1.0/matrixDimension;
		}
	}
	size_t rateCount = 0;
	char* model_string = String_clone(Hashtable_get(options, "model"));
	
//	StringBuffer* buffer = new_StringBuffer(10);
	if (dt == DATA_TYPE_NUCLEOTIDE) {
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
		
		if( strcasecmp("GTR", model_string) == 0 || strcasecmp("SYM", model_string) == 0){
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
			
			if (!frequencies_fixed && strcasecmp("K80", model_string) != 0 && strcasecmp("SYM", model_string) != 0 && strcasecmp("JC69", model_string) != 0) {
				add_json_node_string(jparams, NULL, "$freqs");
			}
			
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
		if(strlen(rates_user) != 0 || rateCount > 0){
			json_node* jrates = create_json_node_object(jsubstmodel, "rates");
			add_json_node(jsubstmodel, jrates);
			
			double *rates = NULL;
			char** names = NULL;
			
			if ( strlen(rates_user) != 0 ) {
				unsigned rcount = 0;
				rates = String_to_double_array( rates_user, ',', &rcount);
				rateCount = rcount;
				
				if( rateCount == 1 ){
					names = (char**)malloc(sizeof(char*));
					names[0] = String_clone("kappa");
				}
				else if( rateCount == 5 ){
					names = (char**)malloc(sizeof(char*)*rateCount);
					char* params[5] = {"ac", "ag", "at", "cg", "ct"};
					for (int i = 0; i < rateCount; i++) {
						names[i] = String_clone(params[i]);
					}
				}
				else{
					names = (char**)malloc(sizeof(char*)*rateCount);
					StringBuffer* str = new_StringBuffer(10);
					for (int i = 0; i < rateCount; i++) {
						StringBuffer_empty(str);
						StringBuffer_append_format(str, "r%d", i);
						names[i] = StringBuffer_tochar(str);
					}
					free_StringBuffer(str);
				}
			}
			else if( rateCount == 5 ){
				rates = dvector(rateCount);
				names = (char**)malloc(sizeof(char*)*rateCount);
				char* params[5] = {"ac", "ag", "at", "cg", "ct"};
				for (int i = 0; i < rateCount; i++) {
					names[i] = String_clone(params[i]);
					rates[i] = 1;
				}
			}
			else if( rateCount == 1 ){
				rates = dvector(rateCount);
				names = (char**)malloc(sizeof(char*));
				names[0] = String_clone("kappa");
				rates[0] = 1;
			}
			
			StringBuffer* buffer = new_StringBuffer(10);
			for (int i = 0; i < rateCount; i++) {
				create_json_node_parameter(jrates, names[i], rates[i], 0, INFINITY);
				if (!rates_fixed) {
					StringBuffer_empty(buffer);
					StringBuffer_append_format(buffer, "&%s", names[i]);
					add_json_node_string(jparams, NULL, buffer->c);
				}
				free(names[i]);
			}
			free_StringBuffer(buffer);
			free(names);
			free(rates);
		}
	}
	else if ( dt == DATA_TYPE_CODON) {
		int genetic_code = atoi(Hashtable_get(options, "genetic-code"));
		add_json_node_string(jsubstmodel, "model", model_string);
		add_json_node_string(jsitepatterns, "code", GENETIC_CODE_NAMES[genetic_code]);
		matrixDimension = NUMBER_OF_CODONS[genetic_code];
		
		json_node* jfreqs = create_json_node_object(jsubstmodel, "frequencies");
		add_json_node(jsubstmodel, jfreqs);
		
		add_json_node_string(jfreqs, "id", "freqs");
		add_json_node_string(jfreqs, "type", "simplex");
		add_json_node_array_double(jfreqs, "values", frequencies, matrixDimension);
	}
	else if ( dt == DATA_TYPE_AMINO_ACID ) {
		add_json_node_string(jsubstmodel, "model", model_string);
		
		if ( strlen(frequencies_string_user) != 0) {
			json_node* jfreqs = create_json_node_object(jsubstmodel, "frequencies");
			add_json_node(jsubstmodel, jfreqs);
			
			add_json_node_string(jfreqs, "id", "freqs");
			add_json_node_string(jfreqs, "type", "simplex");
			add_json_node_array_double(jfreqs, "values", frequencies, matrixDimension);
		}
//		rates_fixed = true;
//		frequencies_fixed = true;
	}
	else{
		
		add_json_node_string(jsubstmodel, "model", "generic");
		
		// empirical
		unsigned *rateIndexes = NULL;
		
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
		
		json_node* discrete_node = create_json_node_object(jsubstmodel, "structure");
		add_json_node(jsubstmodel, discrete_node);
		add_json_node_string(discrete_node, "id", "discrete");
		add_json_node_string(discrete_node, "type", "discrete");
		add_json_node_array_unsigned(discrete_node, "values", rateIndexes, matrixDimension*matrixDimension);
		bool normalize_q = atoi(Hashtable_get(options, "q-normalize"));
		if(!normalize_q){
			add_json_node_bool(jsubstmodel, "normalize", false);
		}
		free(rateIndexes);
		
		json_node* jrates = create_json_node_object(jsubstmodel, "rates");
		add_json_node(jsubstmodel, jrates);
		
		StringBuffer* buffer = new_StringBuffer(10);
		if ( strlen(rates_user) != 0 ) {
			unsigned len = 0;
			double *rates = String_to_double_array(rates_user, ',', &len);
			for (int i = 0; i < len; i++) {
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "r%d", i);
				create_json_node_parameter(jrates, buffer->c, rates[i], 0, INFINITY);
			}
			free(rates);
		}
		else{
			for (int i = 0; i < rateCount; i++) {
				StringBuffer_empty(buffer);
				StringBuffer_append_format(buffer, "r%d", i);
				create_json_node_parameter(jrates, buffer->c, 1.0, 0, INFINITY);
				if (!rates_fixed) {
					StringBuffer_empty(buffer);
					StringBuffer_append_format(buffer, "&r%d", i);
					add_json_node_string(jparams, NULL, buffer->c);
				}
			}
		}
		free_StringBuffer(buffer);
	}
	
	if(frequencies != NULL) free(frequencies);
	free(model_string);
	
	if(jparams->child_count > 0){
		add_json_node(joptmetalist, jopt_mat);
	}
	else{
		json_free_tree(jopt_mat);
	}
	
	Hashtable_add(nodes, "jsubstmodel", jsubstmodel);
}

void create_json_site_model(Hashtable* options, Hashtable* nodes){
	json_node* joptmetalist = Hashtable_get(nodes, "joptmetalist");
	json_node* jtlk = Hashtable_get(nodes, "jtlk");
	
	json_node* jsitemodel = create_json_node_object(jtlk, "sitemodel");
	add_json_node(jtlk, jsitemodel);
	
	add_json_node_string(jsitemodel, "id", "sitemodel");
	add_json_node_string(jsitemodel, "type", "sitemodel");

	json_node* dt = Hashtable_get(nodes, "jdatatype");
	if (dt->node_type == MJSON_STRING) {
		add_json_node_string(jsitemodel, "datatype", dt->value);
	}
	else{
		StringBuffer* buffer = new_StringBuffer(10);
		StringBuffer_append_format(buffer, "&%", get_json_node_value_string(dt, "id"));
		add_json_node_string(jsitemodel, "datatype", buffer->c);
		free_StringBuffer(buffer);
	}
	
	double pinv = atof(Hashtable_get(options, "invariant"));
	int rate_category_count = atoi(Hashtable_get(options, "cat"));
	int alpha = atoi(Hashtable_get(options, "alpha"));
	bool alpha_fixed = false;
	bool pinv_fixed = false;
	if (strlen(Hashtable_get(options, "fix")) != 0) {
		char* fix = Hashtable_get(options, "fix");
		alpha_fixed = String_contains(fix, 'a');
		pinv_fixed = String_contains(fix, 'i');
	}
	
	if ( pinv > 0 || rate_category_count > 1 ) {
		json_node* jrates = create_json_node_object(jsitemodel, "rates");
		add_json_node(jsitemodel, jrates);
		if (pinv > 0) {
			create_json_node_parameter(jrates, "invariant", pinv, 0, INFINITY);
			if (!pinv_fixed) {
				json_node* jopt_pinv = create_json_node_object(joptmetalist, NULL);
				add_json_node(joptmetalist, jopt_pinv);
				add_json_node_string(jopt_pinv, "id", "optpinv");
				add_json_node_string(jopt_pinv, "type", "optimizer");
				add_json_node_string(jopt_pinv, "algorithm", "brent");
				add_json_node_string(jopt_pinv, "model", "&treelikelihood");
				add_json_node_string(jopt_pinv, "parameters", "&invariant");
			}
		}
		if(rate_category_count > 1){
			//			add_json_node_string(jsitemodel, "distribution", "gamma");
			add_json_node_size_t(jsitemodel, "categories", rate_category_count);
			char* sitemodel_string = Hashtable_get(options, "het");
			if ( strlen(sitemodel_string) != 0 && strcasecmp(sitemodel_string, "gammaquad") == 0 ) {
				add_json_node_string(jsitemodel, "discretization", "laguerre");
			}
			create_json_node_parameter(jrates, "alpha", alpha, 0, INFINITY);
			if(!alpha_fixed){
				json_node* jopt_alpha = create_json_node_object(joptmetalist, NULL);
				add_json_node(joptmetalist, jopt_alpha);
				add_json_node_string(jopt_alpha, "id", "optalpha");
				add_json_node_string(jopt_alpha, "type", "optimizer");
				add_json_node_string(jopt_alpha, "algorithm", "brent");
				add_json_node_string(jopt_alpha, "model", "&treelikelihood");
				add_json_node_string(jopt_alpha, "parameters", "&alpha");
			}
		}
		//TODO: init
		//		if(pinv == 0.0 && use_pinv){
		//			pinv = fmax(0.01, 1.0 - ((double)polymorphisms/sp->nsites));
		//		}
	}
	
	Hashtable_add(nodes, "jsitemodel", jsitemodel);
}

void create_json_site_pattern(Hashtable* options, Hashtable* nodes){
	json_node* jtlk = Hashtable_get(nodes, "jtlk");
	
	json_node* jsitepatterns = create_json_node_object(jtlk, "sitepattern");
	add_json_node(jtlk, jsitepatterns);
	
	add_json_node_string(jsitepatterns, "id", "patterns");
	add_json_node_string(jsitepatterns, "type", "sitepattern");
	
	json_node* jsequences = create_json_node_object(jsitepatterns, "alignment");
	add_json_node(jsitepatterns, jsequences);
	Hashtable_add(nodes, "jsequences", jsequences);
	
	add_json_node_string(jsequences, "id", "sequences");
	add_json_node_string(jsequences, "type", "alignment");
	
	char* seq_file = Hashtable_get(options, "sequences");
	int nexus_index = atoi(Hashtable_get(options, "nexus-index"));
	add_json_node_string(jsequences, "file", seq_file);
	if( nexus_index >= 0 ){
		add_json_node_size_t(jsequences, "index", nexus_index);
	}
	
	json_node* dt = Hashtable_get(nodes, "jdatatype");
	if (dt->node_type == MJSON_STRING) {
		add_json_node_string(jsequences, "datatype", dt->value);
		add_json_node_string(jsitepatterns, "datatype", dt->value);
	}
	else{
		StringBuffer* buffer = new_StringBuffer(10);
		StringBuffer_append_format(buffer, "&%s", get_json_node_value_string(dt, "id"));
		add_json_node(jsitepatterns, dt);
		dt->parent = jsitepatterns;
		add_json_node_string(jsequences, "datatype", buffer->c);
		free_StringBuffer(buffer);
	}
	
	Hashtable_add(nodes, "jsitepatterns", jsitepatterns);
	Hashtable_add(nodes, "jsequences", jsequences);
}

void create_json_node_physher(Hashtable* nodes){
	json_node* jroot = Hashtable_get(nodes, "jroot");
	json_node* jphysher = create_json_node_object(jroot, "physher");
	add_json_node(jroot, jphysher);
	jphysher->node_type = MJSON_ARRAY;
	Hashtable_add(nodes, "jphysher", jphysher);
}

void create_json_node_meta_opt(Hashtable* nodes){
	json_node* jphysher = Hashtable_get(nodes, "jphysher");
	
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
	
	Hashtable_add(nodes, "joptmeta", joptmeta);
	Hashtable_add(nodes, "joptmetalist", joptmetalist);
}

Hashtable * extract_arguments(int argc, char* argv[]){
	
	struct argsparser2_option options[] = {
		{ARGS_OPTION_STRING,  'i', "sequences",    NULL, "Input alignment file"},
		{ARGS_OPTION_STRING,  't', "tree",         NULL, "Input tree file"},
		{ARGS_OPTION_STRING,  'o', "stem",         NULL, "Output stem file"},
		
		{ARGS_OPTION_INTEGER, 'g', "genetic-code", String_clone("0"), "Genetic Code"},
		{ARGS_OPTION_INTEGER, 'd', "datatype",     NULL, "Data type (nucleotide or aa)"},
		
		{ARGS_OPTION_STRING,  'm', "model",        NULL, "Susbtitution model"},
		{ARGS_OPTION_STRING,  0,   "states",       NULL, "State space of Markov process. For nucleotide --states A,C,G,T"},
		{ARGS_OPTION_STRING,  'f', "frequencies",  NULL, "Frequencies as an array or 'e' for equal frequencies. For nucleotide -f 0.2,0.3,0.4,0.1"},
		{ARGS_OPTION_STRING,  'r', "rates",        NULL, "Relative rates of the susbtitution matrix"},
		{ARGS_OPTION_BOOLEAN, 0,   "q-normalize",  String_clone("1"), "Normalize rate matrix"},
		{ARGS_OPTION_FLAG,    0,   "f-unknown",    String_clone("0"), "Set frequencies to 1.0"},
		{ARGS_OPTION_STRING,  0,   "q-search",     NULL, "Find best rate matrix (ga)"},
		
		{ARGS_OPTION_INTEGER, 'c', "cat",          String_clone("1"), "Number of rate categories for gamma distribution"},
		{ARGS_OPTION_STRING,  'H', "het",          NULL, "discrete or gammaquad"},
		{ARGS_OPTION_DOUBLE,  'a', "alpha",        String_clone("0.5"), "Value of the alpha parameter of the gamma distribution"},
		{ARGS_OPTION_DOUBLE,  'I', "invariant",    String_clone("0"), "Value of the proportion of invariant sites"},
		{ARGS_OPTION_FLAG,    0,   "ps",           String_clone("0"), "Caclulate posterior estimates of rates at each site"},
		
		{ARGS_OPTION_STRING,  'F', "fix",          NULL, "Fix d: branch length, i: invariant, a: alpha, f: frequencies, r: rates"},
//		{ARGS_OPTION_BOOLEAN, 0,   "sse",          "treelikelihood.sse", &use_sse, "Use SSE [default true]"},
//		{ARGS_OPTION_BOOLEAN, 0,   "upper",        "treelikelihood.upper", &use_upper, "Use upper likelihood"},
//		{ARGS_OPTION_INTEGER, 0,   "lk-threads",   "treelikelihood.threads", &tlk_threads, "Number of threads for likelihood calculation"},
//		{ARGS_OPTION_INTEGER, 'A', "approx",       "treelikelihood.approximation", &approximation, "Input alignment file"},
//		
//		{ARGS_OPTION_FLAG,    'U', "unrooted",     "tree.unrooted", &tree_unrooted, "Input tree is rooted"},
		{ARGS_OPTION_DOUBLE,  's', "scaler",       String_clone("-1"), "Scale input tree"},
		{ARGS_OPTION_STRING,  'O', "treeopt",      NULL, "Optimize topology nni or spr (experimental)"},
//		{ARGS_OPTION_BOOLEAN,  0,  "parsimony",    "tree.topolology.optimize.parsimony", &use_parsimony, "Quick optimizaton of tree topology using parsimony before ML optimization"},
//
//		{ARGS_OPTION_STRING,  'C', "clock",        "clock", &clock, "Clock type: strict, local, discrete"},
//		{ARGS_OPTION_FLAG,    0,   "forward",      "clock.forward", &forward, "Time is forward"},
//		{ARGS_OPTION_DOUBLE,  0,   "clock-rate",   "clock.rate", &clock_rate_guess, "A rate guess"},
//		{ARGS_OPTION_STRING,  'S', "clock-search", "clock.algorithm", &clock_algorithm, "Algorithm for local and discrete clock: ga or greedy or exhaustive"},
//		{ARGS_OPTION_INTEGER, 0,   "clock-cat",    "clock.discrete.cat", &clock_categories, "Number of discrete rate categories along phylogeny"},
//		
//		// GA
//		{ARGS_OPTION_INTEGER, 0,   "ga-pop",       "ga.popsize", &ga_population_size, "Genetic algorithm population size"},
//		{ARGS_OPTION_INTEGER, 0,   "ga-gen",       "ga.ngen", &ga_generations, "Genetic algorithm number of generations"},
//		{ARGS_OPTION_INTEGER, 0,   "ga-no-improv", "ga.maxnoimprovement", &ga_max_no_improvement, "Genetic algorithm number of generation without improvment before stopping"},
//		
//		{ARGS_OPTION_STRING,  0,   "ic",           "ic", &ic, "Information criterion (AIC, AICc, BIC)"},
//		{ARGS_OPTION_INTEGER, 0,   "ic-ss",        "ic.samplesize", &ic_sample_size, "Sample size for information criterion"},
//		
//		{ARGS_OPTION_INTEGER, 'b', "bootstrap",    "resampling.bootstrap", &bootstrap, "Number of bootstrap replicates"},
//		{ARGS_OPTION_FLAG,    0,   "bca",          "resampling.bootstrap.bca", &bca, "Use BCA bootstrap"},
//		{ARGS_OPTION_INTEGER, 'j', "jackknife",    "resampling.jackknife", &jackknife, "Jackknife"},
//		
//		// for GA, greedy
//		{ARGS_OPTION_INTEGER, 'T', "nthreads",     "nthreads", &nthreads, "Number of threads for GA and bootstrap algorithms"},
//		
//		{ARGS_OPTION_INTEGER, 'V', "verbose",      "verbosity", &verbosity, "Verbosity"},
		{ARGS_OPTION_LONG,    'R', "seed",         String_clone("-1"), "Random seed"},
//		{ARGS_OPTION_FLAG,    0,   "asr",          "asr", &asr, "Ancestral sequence reconstruction"},
//		{ARGS_OPTION_FLAG,    0,   "double-matrix","distance.matrix.double.precision", &use_double_distance, "Use double precision for distance matrix"},
//		{ARGS_OPTION_FLAG,    0,   "batch",        "batch", &batch, "Use batch mode"},
//		{ARGS_OPTION_FLAG,    0,   "hessian",      "derivative.hessian", &hessian, "Calculate Hessian"},
		
		{ARGS_OPTION_STRING,  'D', "distance",      NULL, "NJ or UPGMA"},
		
		
		{ARGS_OPTION_INTEGER, 0,   "nexus-index",   String_clone("-1"), "Index of tree in a multi-tree nexus file"},
		
		{ARGS_OPTION_FLAG,  '!',   "overwrite",     NULL, "Overwrite output files"},
		{ARGS_OPTION_FLAG,    0,   "dry",           NULL, "Only print the json file"},
		{ARGS_OPTION_FLAG,  'h',   "help",          NULL, "Print help"}
	};
	
	args_parser2* argsparser = argsparser2_init(options, sizeof(options)/sizeof(struct argsparser2_option));
	Hashtable* hash = argsparser2_parse(argsparser, argv, argc);
	
	if (strlen(Hashtable_get(hash, "help")) != 0) {
		argsparser2_help(argsparser);
		exit(0);
	}
	
	argsparser2_free(argsparser);
	
	return hash;
}

json_node* create_json_file(int argc, char* argv[]){
	
	Hashtable* options = extract_arguments(argc, argv);
	
	if (strlen(Hashtable_get(options, "distance")) != 0) {
		json_node* jroot = create_json_distance(options);
		free_Hashtable(options);
		return jroot;
	}
	
	Hashtable* nodes = new_Hashtable_string(10);
	hashtable_set_key_ownership(nodes, false);
	hashtable_set_value_ownership(nodes, false);
	
	json_node* jroot = create_json_node(NULL);
	jroot->node_type = MJSON_OBJECT;
	Hashtable_add(nodes, "jroot", jroot);
	
	json_node* jinit = create_json_node_object(jroot, "init");
	add_json_node(jroot, jinit);
	
	long seed = atol(Hashtable_get(options, "seed"));
	if(seed < 0) seed = time(NULL);
	init_genrand(seed);
	
	add_json_node_size_t(jinit, "seed", seed);
	
	json_node* jtlk = create_json_node_object(jroot, "model");
	add_json_node_string(jtlk, "id", "treelikelihood");
	add_json_node_string(jtlk, "type", "treelikelihood");
	add_json_node(jroot, jtlk);
	Hashtable_add(nodes, "jtlk", jtlk);
	
	bool frequencies_unknown = atoi(Hashtable_get(options, "f-unknown"));
	if(frequencies_unknown){
		add_json_node_size_t(jtlk, "root_frequencies", 1);
	}
	
	
	create_json_node_physher(nodes);
	create_json_node_meta_opt(nodes);
	
	create_json_datatype(options, nodes);
	
	create_json_site_pattern(options, nodes);
	
	create_json_phylo_tree(options, nodes);
	
	create_json_site_model(options, nodes);
	
	create_json_substitution_model(options, nodes);
	
	json_node* jdatatype = Hashtable_get(nodes, "jdatatype");
	// not included the jroot tree
	if (jdatatype->node_type == MJSON_STRING) {
		json_free_tree(jdatatype);
	}
	free_Hashtable(nodes);
	free_Hashtable(options);

	return jroot;
}
