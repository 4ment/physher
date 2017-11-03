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
	bool ownconstraint;
#ifdef LISTENERS
	ListenerList *listeners;
	int refCount;
#endif
};

struct _DiscreteParameter{
	char *name;
	int id;
	unsigned* values;
	unsigned length;
	
	DiscreteParameter* (*clone)(DiscreteParameter*);
//	void (*free)(DiscreteParameter*);
	
	void (*set_value)( DiscreteParameter*, int, unsigned );
	void (*set_values)( DiscreteParameter*, const unsigned* );
#ifdef LISTENERS
	ListenerList *listeners;
	int refCount;
#endif
};

#pragma mark -
#pragma mark Constraint

Constraint * new_Constraint( const double lower, const double upper );

void free_Constraint( Constraint *c );

Constraint * clone_Constraint( Constraint *cnstr );

StringBuffer * Constraint_bufferize( StringBuffer *buffer, Constraint *c );

void * Constraint_SML_to_object( SMLNode node );

void remove_contraint( Parameter *p );

void add_contraint( Parameter *p, Constraint *constr, bool ownconstraint );

bool Constraint_fixed( const Constraint *c );

bool Constraint_lower_fixed( const Constraint *c );

bool Constraint_upper_fixed( const Constraint *c );

void Constraint_set_fixed( Constraint *c, const bool fixed );

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

Parameter * new_Parameter_with_postfix_and_ownership( const char *name, const char *postfix, const double value, Constraint *constr, bool owner );

void free_Parameter( Parameter *p );

Parameter * clone_Parameter( Parameter *p, bool duplicateconstraint );

char * Parameter_stringify( Parameter *p );

StringBuffer * Parameter_bufferize( StringBuffer *buffer, Parameter *p );

StringBuffer * Parameter_bufferize_parameter_set( StringBuffer *buffer, Parameter *p );

void * Parameter_SML_to_object( SMLNode node );

void * Parameter_SML_to_object_with_postfix( SMLNode node, const char *postfix );

char * Parameter_name( const Parameter *p );

void Parameter_set_value( Parameter *p, const double value );

double Parameter_value( const Parameter *p );

double check_value( Constraint *cnstr, double value );


bool Parameters_is_at_boundry( const Parameters *p, const int index, double precision );

Constraint * Parameter_constraint( const Parameter *p );

bool Parameter_fixed( const Parameter *p );

void Parameter_set_fixed( Parameter *p, const bool fixed );

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

Parameters * new_Parameters( const int n );

Parameters * new_Parameters_with_contraint( const int capacity, Constraint *cnstr, const char *name );

void free_Parameters( Parameters *ps );

void free_Parameters_soft( Parameters *ps );

Parameters * clone_Parameters( Parameters *p, bool ownconstraint );

char * Parameters_stringify( Parameters *ps );

StringBuffer * Parameters_SML_bufferize( StringBuffer *buffer, Parameters *ps );

Parameters * Parameters_SML_to_object( SMLNode node );


void Parameters_set( Parameters *ps, const int index, Parameter *p );

Parameter * Parameters_at( const Parameters *p, const int index );

void Parameters_add(Parameters *ps, Parameter *p);

void Parameters_move( Parameters *ps, Parameter *p);

void Parameters_add_parameters(Parameters *dst, const Parameters *src);

char * Parameters_name( const Parameters *p, const int index );

int Parameters_count( const Parameters *p );

int Parameters_capacity( const Parameters *p );

void Parameters_set_value( Parameters *p, const int index, const double value );

void Parameters_set_all_value( Parameters *p, const double value );

double Parameters_value( const Parameters *p, const int index );



bool Parameters_fixed( const Parameters *p, const int index );

void Parameters_set_fixed( Parameters *p, const bool fixed, const int index );

Constraint * Parameters_constraint( const Parameters *p, const int index );

double Parameters_upper( const Parameters *p, const int index );

double Parameters_lower( const Parameters *p, const int index );

void Parameters_set_upper( Parameters *p, const int index, const double value );

void Parameters_set_lower( Parameters *p, const int index, const double value );

void Parameters_set_bounds( Parameters *p, const int index, const double lower, const double upper );

void Parameters_remove( Parameters *params, int index );

void Parameters_pop( Parameters *params );

void Parameters_pop_soft( Parameters *params );

bool update_matching_parameters( Parameters *params, const Parameters *new, const double precision );

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

DiscreteParameter * new_DiscreteParameter( const char *name, int dim );

DiscreteParameter * new_DiscreteParameter_with_postfix( const char *name, const char *postfix, int dim );

DiscreteParameter * new_DiscreteParameter_with_values( const char *name, const unsigned* values, int dim );

DiscreteParameter * new_DiscreteParameter_with_postfix_values( const char *name, const char *postfix, const unsigned* values, int dim );



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
	void* data;
	double (*logP)( Model * );
	double (*dlogP)( Model *, const Parameter* );
	Model* (*clone)( Model *, Hashtable* );
	void (*free)( Model * );
	void (*update)( Model *, Model *, int );
	
	Parameters *parameters; // pointers to the parameters of the model
	ListenerList *listeners;
	int ref_count;
};


#pragma mark -

Model * new_Model( const char *name, void *obj );

void free_Model( Model *model );

#pragma mark -

ListenerList * new_ListenerList( const unsigned capacity );

void free_DiscreteParameter( DiscreteParameter *p );

#endif
