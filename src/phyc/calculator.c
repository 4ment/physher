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

void calculatorModel_from_json(json_node* node, Hashtable* hash){
	char* what_node = get_json_node_value_string(node, "whattodo");
	
	if(strcasecmp(what_node, "asr") == 0){
		asr_marginal_calculator_from_json(node, hash);
	}
	else if(strcasecmp(what_node, "ppsite") == 0){
		posteriors_sites_calculator_from_json(node, hash);
	}
}
