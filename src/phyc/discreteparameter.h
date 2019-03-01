//
//  discreteparameter.h
//  physher
//
//  Created by Mathieu Fourment on 1/03/2019.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#ifndef discreteparameter_h
#define discreteparameter_h

#include "parameters.h"

struct _DiscreteParameter;
typedef struct _DiscreteParameter DiscreteParameter;

struct _DiscreteParameter{
	unsigned* values;
	unsigned* stored_values;
	size_t length;
	
	DiscreteParameter* (*clone)(DiscreteParameter*);
	void (*free)(DiscreteParameter*);
	
	void (*set_value)( DiscreteParameter*, int, unsigned );
	void (*set_values)( DiscreteParameter*, const unsigned* );
	
	ListenerList* listeners;
};

DiscreteParameter * new_DiscreteParameter( size_t dim );

DiscreteParameter * new_DiscreteParameter_with_postfix( const char *postfix, size_t dim );

DiscreteParameter * new_DiscreteParameter_with_values( const unsigned* values, size_t dim );

DiscreteParameter * new_DiscreteParameter_with_postfix_values( const char *postfix, const unsigned* values, size_t dim );

Model * new_DiscreteParameterModel( const char* name, DiscreteParameter *dp );

Model* new_DiscreteParameterModel_from_json(json_node* node, Hashtable* hash);

#endif /* discreteparameter_h */
