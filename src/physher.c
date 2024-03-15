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

#include "phyc/PhyCConfig.h"

#include "phyc/parsimony.h"
#include "phyc/treelikelihood.h"
#include "phyc/random.h"
#include "phyc/args.h"
#include "phyc/filereader.h"
#include "phyc/compoundmodel.h"
#include "phyc/distmodelfactory.h"
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
#include "phyc/demographicmodels.h"
#include "phyc/asr.h"
#include "phyc/ppsites.h"
#include "phyc/cat.h"
#include "phyc/sbn.h"
#include "phyc/checkpoint.h"
#include "phyc/optimizer.h"

#include "phyc/physhercmd.h"


int main(int argc, const char* argv[]){
	time_t start_time;
	time_t beginning_of_time;
	time_t end_time;
	double diff_time;
	
	time(&start_time);
	beginning_of_time = start_time;
	
	json_node* json = NULL;
	
    if (argc == 1 || args_contains(argc, (char**)argv, "--help") || args_contains(argc, (char**)argv, "-h")) {
		fprintf(stdout, "\nphysher version %d.%d.%s", PHYC_VERSION_MAJOR, PHYC_VERSION_MINOR, PHYC_PATCH_LEVEL);
		if(strcmp(GIT_HEAD_TAG, "") == 0 && strcmp(GIT_COMMIT_HASH, "") != 0){
			fprintf(stdout, "-%s built on", GIT_COMMIT_HASH);
		}
		// Build from a directory without git
		else if(strcmp(GIT_COMMIT_HASH, "") == 0){
			fprintf(stdout, " compiled on");
        }
        else{
			fprintf(stdout, " released on\n");
        }
		fprintf(stdout, " %s\n", __DATE__);
		fprintf(stdout, "Compiled with:");
		bool flag = false;
		if(PHYC_SSE_ENABLED){
			fprintf(stdout, " %s", PHYC_SSE_LEVEL);
			flag = true;
		}
		if(PHYC_OPENMP_ENABLED){
			fprintf(stdout, "%s OpenMP", (flag ? ",": ""));
			flag = true;
		}
		if(PHYC_PTHREAD_ENABLED){
			fprintf(stdout, "%s PThread", (flag ? ",": ""));
		}
		fprintf(stdout, "\n");
		fprintf(stdout, "Developped by: Mathieu Fourment\n");
		fprintf(stdout, "Project hosted on: https://github.com/4ment/physher\n");

        // printf("Fourment M and Holmes EC. Novel non-parametric models to estimate evolutionary rates and divergence times from heterochronous sequence data.\n");
        // printf("BMC Evolutionary Biology 14:163, 2014\n\n");

		fprintf(stdout, "\nusage: physher [--help] [--seed SEED] input-file-name\n");
		fprintf(stdout, "\npositional arguments:\n");
  		fprintf(stdout, "  input-file-name     JSON configuration file\n\n");
		fprintf(stdout, "options:\n");
  		fprintf(stdout, "  --help              show this help message and exit\n");
  		fprintf(stdout, "  --seed              SEED  initialize seed\n");
		return 0;
	}
	else if(array_of_string_contains("--dry", argv, argc, true)){
		json = create_json_file(argc, argv);
		json_tree_print(json);
		printf("\n");
		json_free_tree(json);
		exit(0);
	}
	else{
		char* content = load_file(argv[argc-1]);
		printf("Reading file %s\n", argv[argc-1]);
		printf("done\n\n");

		json = create_json_tree(content);
		free(content);
		
		// remove json_nodes containing "ignore": true
		json_prune_ignored(json);
        json_prune_underscored(json);
	}

	long seeed = args_get_long(argc, (char**)argv, "--seed", time(NULL));
	
	Hashtable* hash2 = new_Hashtable_string(100);
	hashtable_set_key_ownership( hash2, false );
	hashtable_set_value_ownership( hash2, false );
	
	json_node* run_node = get_json_node(json, "physher");
	json_node* init_node = get_json_node(json, "init");
	
	// seed provided through command line overrides seed specified in json file
	if (init_node != NULL && !args_contains(argc, (char**)argv, "--seed")) {
		seeed = get_json_node_value_size_t(init_node, "seed", seeed);
	}
	else{
		printf("seed: %ld\n", seeed);
	}
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
			
			model_t model_type = check_model(type_node->value);
			if (model_type < 0) {
				fprintf(stderr, "model type not recognized: %s\n\n", (char*)type_node->value);
				fprintf(stderr, "Possible models:\n");
				for (int j = 0; j < sizeof(model_type_strings)/sizeof(model_type_strings[0]); j++) {
					fprintf(stderr, "%s\n", model_type_strings[j]);
				}
				exit(2);
			}
			
			if (model_type == MODEL_COMPOUND) {
				models[index] = new_CompoundModel_from_json(child, hash2);
			}
			else if(model_type == MODEL_VARIATIONAL){
				models[index] = new_Variational_from_json(child, hash2);
			}
			else if(model_type == MODEL_DISTRIBUTION){
				models[index] = new_DistributionModel_from_json(child, hash2);
			}
			else if(model_type == MODEL_TREELIKELIHOOD){
				models[index] = new_TreeLikelihoodModel_from_json(child, hash2);
			}
			else if(model_type == MODEL_PARSIMONY){
				models[index] = new_ParsimonyModel_from_json(child, hash2);
			}
			else if(model_type == MODEL_TREE){
				models[index] = new_TreeModel_from_json(child, hash2);
			}
			else if(model_type == MODEL_COALESCENT){
				models[index] = new_CoalescentModel_from_json(child, hash2);
			}
			
			Hashtable_add(hash2, id, models[index]);
			index++;
		}
	}
	
	if(run_node != NULL)
	for (int i = 0; i < run_node->child_count; i++) {
        json_node* child = run_node->children[i];
        char* type = get_json_node_value_string(child, "type");
        
		if (strcasecmp(type, "optimizer") == 0) {
			Optimizer* opt = new_Optimizer_from_json(child, hash2);
			char* checkpoint = args_get_string(argc, (char**)argv, "-c");
			if(checkpoint != NULL){
				checkpoint_apply(checkpoint, opt_parameters(opt));
			}
			double logP;
			opt_optimize(opt, NULL, &logP);
			free_Optimizer(opt);
		}
		else if (strcasecmp(type, "logger") == 0) {
			struct Logger* logger = new_logger_from_json(child, hash2);
			logger->log(logger);
			free_Logger(logger);
		}
        else if (strcasecmp(type, "dumper") == 0) {
            struct Dumper* dumper = new_Dumper_from_json(child, hash2);
            dumper->dump(dumper);
            dumper->free(dumper);
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
		else if(strcasecmp(type, "sbn") == 0){
			SBN* sbn = new_SBN_from_json(child, hash2);
		}
		else if(strcasecmp(type, "asr") == 0){
			asr_calculator_from_json(child, hash2);
		}
		else if(strcasecmp(type, "ppsite") == 0){
			posteriors_sites_calculator_from_json(child, hash2);
		}
		else if(strcasecmp(type, "cat") == 0){
			cat_estimator_from_json(child, hash2);
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
