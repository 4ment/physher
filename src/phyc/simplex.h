//
//  simplex.h
//  physher
//
//  Created by Mathieu Fourment on 7/11/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef simplex_h
#define simplex_h

#include <stdio.h>
#include "parameters.h"
#include "mjson.h"


struct _Simplex;
typedef struct _Simplex Simplex;


struct _Simplex{
	size_t K;
	Parameter* parameter; // K-1 parameters
	double* values; // K double
	double* stored_values; // K double
	const double* (*get_values)(Simplex*);
	double (*get_value)(Simplex*, int);
	void (*set_values)(Simplex*, const double*);
	
	void (*set_parameter_value)(Simplex*, int, double);
	void (*gradient)(Simplex*, size_t, double*);
	bool need_update;
};

Simplex* new_Simplex_with_values(const char* name, const double *x, size_t K);

Simplex* new_Simplex_with_parameter(const char* name, Parameter* parameter);

Simplex* new_Simplex(const char* name, size_t K);

void free_Simplex(Simplex* simplex);

Simplex* clone_Simplex(const Simplex* simplex);

Model * new_SimplexModel( const char* name, Simplex *simplex );

Model* new_SimplexModel_from_json(json_node*node, Hashtable*hash);

Parameter* new_SimplexParameter_from_json(json_node*node, Hashtable*hash);

void Simplex_use_stan_transform(Simplex* simplex, bool use_stan);

#endif /* simplex_h */
