/*
 *  parameters.h
 *  PhyC
 *
 *  Created by Mathieu Fourment on 11/15/10.
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

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "utils.h"

#include <float.h>

#include "mstring.h"
#include "parser.h"
#include "model.h"

#include "mjson.h"
#include "hashtable.h"

#define PARAMETER_TINY 1.0e-25
#define PARAMETER_ZERO 0.
#define PARAMETER_ZERO_PLUS TINY
#define PARAMETER_ONE 1.
#define PARAMETER_ONE_MINUS (1-TINY)
#define PARAMETER_POSITIVE_INFINTY INFINITY
#define PARAMETER_NEGATIVE_INFINTY (-INFINITY)


struct _Constraint;
typedef struct _Constraint Constraint;

struct _Parameter;
typedef struct _Parameter Parameter;

struct _Parameters;
typedef struct _Parameters Parameters;

struct _DiscreteParameter;
typedef struct _DiscreteParameter DiscreteParameter;


struct _ListenerList;
typedef struct _ListenerList ListenerList;

struct _Listeners;
typedef struct _Listener Listener;

struct _Model;
typedef struct _Model Model;


struct _Parameter{
	char *name;
	int id;
	double value;
	Constraint *cnstr;
	bool estimate;
	ListenerList *listeners;
	int refCount;
};

struct _DiscreteParameter{
	char *name;
	int id;
	unsigned* values;
	unsigned length;
	
	DiscreteParameter* (*clone)(DiscreteParameter*);
	void (*free)(DiscreteParameter*);
	
	void (*set_value)( DiscreteParameter*, int, unsigned );
	void (*set_values)( DiscreteParameter*, const unsigned* );
	ListenerList *listeners;
	int refCount;
};

#pragma mark -
#pragma mark Constraint

Constraint * new_Constraint( const double lower, const double upper );

void free_Constraint( Constraint *c );

Constraint * clone_Constraint( Constraint *cnstr );

StringBuffer * Constraint_bufferize( StringBuffer *buffer, Constraint *c );

void * Constraint_SML_to_object( SMLNode node );

bool Constraint_lower_fixed( const Constraint *c );

bool Constraint_upper_fixed( const Constraint *c );


void Constraint_set_lower_fixed( Constraint *c, const bool fixed );

void Constraint_set_upper_fixed( Constraint *c, const bool fixed );

double Constraint_upper( const Constraint *c );

double Constraint_lower( const Constraint *c );

void Constraint_set_upper( Constraint *c, const double upper );

void Constraint_set_lower( Constraint *c, const double lower );

void Constraint_set_bounds( Constraint *constr, const double lower, const double upper );

double Constraint_fupper( const Constraint *c );

double Constraint_flower( const Constraint *c );

void Constraint_set_fupper( Constraint *c, const double fupper );

void Constraint_set_flower( Constraint *c, const double flower );

void compare_constraint( const Constraint *c1, const Constraint *c2 );

#pragma mark -
#pragma mark Parameter

Parameter * new_Parameter( const char *name, const double value, Constraint *constr );

Parameter * new_Parameter_with_postfix( const char *name, const char *postfix, const double value, Constraint *constr );

Parameter* new_Parameter_from_json(json_node* node, Hashtable* hash);

void free_Parameter( Parameter *p );

Parameter * clone_Parameter( Parameter *p );

char * Parameter_stringify( Parameter *p );

StringBuffer * Parameter_bufferize( StringBuffer *buffer, Parameter *p );

StringBuffer * Parameter_bufferize_parameter_set( StringBuffer *buffer, Parameter *p );

void * Parameter_SML_to_object( SMLNode node );

void * Parameter_SML_to_object_with_postfix( SMLNode node, const char *postfix );

char * Parameter_name( const Parameter *p );

void Parameter_set_name( Parameter *p, const char *name );

void Parameter_set_value( Parameter *p, const double value );

double Parameter_value( const Parameter *p );

double check_value( Constraint *cnstr, double value );


bool Parameters_is_at_boundry( const Parameters *p, const int index, double precision );

Constraint * Parameter_constraint( const Parameter *p );

bool Parameter_estimate( const Parameter *p );

void Parameter_set_estimate( Parameter *p, const bool estimate );

double Parameter_upper( const Parameter *p );

double Parameter_lower( const Parameter *p );

void Parameter_set_upper( Parameter *p, const double value );

void Parameter_set_lower( Parameter *p, const double value );

void Parameter_set_bounds( Parameter *p, const double lower, const double upper );


void Parameter_warn( const Parameter *p, const double value );

void Parameter_die( const Parameter *p, const double value );

void Parameter_print( const Parameter *p );

void compare_parameter( const Parameter *p1, const Parameter *p2 );

#pragma mark -
#pragma mark Parameters

Parameters * new_Parameters( const size_t capacity );

Parameters * new_Parameters_from_json(json_node* node, Hashtable* hash);

void free_Parameters( Parameters *ps );

Parameters * clone_Parameters( Parameters *p );

char * Parameters_stringify( Parameters *ps );

StringBuffer * Parameters_SML_bufferize( StringBuffer *buffer, Parameters *ps );

Parameters * Parameters_SML_to_object( SMLNode node );


Parameter * Parameters_at( const Parameters *p, const size_t index );

void Parameters_add(Parameters *ps, Parameter *p);

void Parameters_move( Parameters *ps, Parameter *p);

void Parameters_add_free_parameters(Parameters *dst, const Parameters *src);

void Parameters_add_parameters(Parameters *dst, const Parameters *src);

void Parameters_set_name( Parameters *p, const size_t index, const char *name );

char * Parameters_name( const Parameters *p, const size_t index );

size_t Parameters_count( const Parameters *p );

size_t Parameters_capacity( const Parameters *p );

void Parameters_set_value( Parameters *p, const int index, const double value );

void Parameters_set_all_value( Parameters *p, const double value );

double Parameters_value( const Parameters *p, const int index );

bool Parameters_estimate( const Parameters *p, const int index );

void Parameters_set_estimate( Parameters *p, const bool estimate, const int index );

Constraint * Parameters_constraint( const Parameters *p, const int index );

double Parameters_upper( const Parameters *p, const int index );

double Parameters_lower( const Parameters *p, const int index );

void Parameters_set_upper( Parameters *p, const int index, const double value );

void Parameters_set_lower( Parameters *p, const int index, const double value );

void Parameters_set_bounds( Parameters *p, const int index, const double lower, const double upper );

void Parameters_remove( Parameters *params, size_t index );

void Parameters_pop( Parameters *params );


Parameters * get_sub_parameters( Parameters *p, const int start, const int end );

Parameters * Parameters_optimizable( Parameters *p, int **map, const char postfix[] );

void Parameters_store_value( const Parameters *p, double *store );

void Parameters_restore_value( Parameters *p, const double *store );

int Parameters_count_optimizable( const Parameters *p, const char postfix[] );

Parameters * pack_parameters(Parameters *ps, Parameter *p);


void Parameters_print( Parameters *ps );

void compare_parameters( const Parameters *p1, const Parameters *p2 );

void Parameters_swap( Parameter **a, Parameter **b );

void Parameters_swap_index( Parameters *ps, unsigned a, unsigned b );

void Parameters_sort_from_ivector( Parameters *p, int *s );

#pragma mark -
#pragma mark DiscreteParameter

DiscreteParameter * new_DiscreteParameter( const char *name, size_t dim );

DiscreteParameter * new_DiscreteParameter_with_postfix( const char *name, const char *postfix, size_t dim );

DiscreteParameter * new_DiscreteParameter_with_values( const char *name, const unsigned* values, size_t dim );

DiscreteParameter * new_DiscreteParameter_with_postfix_values( const char *name, const char *postfix, const unsigned* values, size_t dim );

DiscreteParameter* new_DiscreteParameter_from_json(json_node* node, Hashtable* hash);

struct _ListenerList {
	Model** models;
	int count;
	int capacity;
	void (*free)( ListenerList*);
	void (*fire)( ListenerList*, Model*, int );
	void (*add)( ListenerList*, Model* );
	void (*remove)( ListenerList*, Model* );
	void (*removeAll)( ListenerList*);
};

struct _Model {
	void *obj; // pointer to model
	char *name;
	char* type;
	void* data;
	double (*logP)( Model * );
	double (*dlogP)( Model *, const Parameter* );
	void (*gradient)( Model *, double* );
	Model* (*clone)( Model *, Hashtable* );
	void (*free)( Model * );
	void (*update)( Model *, Model *, int );
	void (*get_free_parameters)(Model*, Parameters*);
	
	ListenerList *listeners;
	int ref_count;
};


#pragma mark -

Model * new_Model( const char *type, const char *name, void *obj );

void free_Model( Model *model );

double Model_first_derivative( Model *model, Parameter* parameter, double eps );

#pragma mark -

ListenerList * new_ListenerList( const unsigned capacity );

#pragma mark -

void get_parameters_references(json_node* node, Hashtable* hash, Parameters* parameters);



#endif
