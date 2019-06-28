//
//  calculator.c
//  physher
//
//  Created by Mathieu Fourment on 28/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "calculator.h"

#include <strings.h>

#include "asr.h"
#include "ppsites.h"
#include "cat.h"

void calculatorModel_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {"whattodo"};
	json_check_required(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	char* what_node = get_json_node_value_string(node, "whattodo");
	
	if(strcasecmp(what_node, "asr") == 0){
		asr_marginal_calculator_from_json(node, hash);
	}
	else if(strcasecmp(what_node, "ppsite") == 0){
		posteriors_sites_calculator_from_json(node, hash);
	}
	else if(strcasecmp(what_node, "cat") == 0){
		cat_estimator_from_json(node, hash);
	}
}
