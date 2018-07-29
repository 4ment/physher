//
//  calculator.h
//  physher
//
//  Created by Mathieu Fourment on 28/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef calculator_h
#define calculator_h

#include "parameters.h"
#include "hashtable.h"

void calculatorModel_from_json(json_node* node, Hashtable* hash);

#endif /* calculator_h */
