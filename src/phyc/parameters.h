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

typedef enum model_t{
	MODEL_BRANCHMODEL=0,
	MODEL_COALESCENT,
	MODEL_COMPOUND,
	MODEL_DISCRETE_PARAMETER,
	MODEL_DISTRIBUTION,
	MODEL_LAPLACE,
	MODEL_PARSIMONY,
	MODEL_SIMPLEX,
	MODEL_SITEMODEL,
	MODEL_SUBSTITUTION,
	MODEL_TREE,
	MODEL_TREE_TRANSFORM,
	MODEL_TREELIKELIHOOD,
    MODEL_VARIATIONAL,
    MODEL_VARIATIONAL_TREELIKELIHOOD
}model_t;

static const char* model_type_strings[] = {
	"branchmodel",
	"coalescent",
	"compound",
	"discreteparameter",
	"distribution",
	"laplace",
	"parsimony",
	"simplex",
	"sitemodel",
	"substitutionmodel",
	"tree",
	"treetransform",
	"treelikelihood",
    "variational",
    "variationaltreelikelihood"
};

struct _Constraint;
typedef struct _Constraint Constraint;

struct _Parameter;
typedef struct _Parameter Parameter;

struct _Parameters;
typedef struct _Parameters Parameters;

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
	double stored_value;
	Constraint *cnstr;
	bool estimate;
	ListenerList *listeners;
	int refCount;
	model_t model; // model it belongs to
};


#pragma mark -
#pragma mark Constraint

Constraint * new_Constraint( const double lower, const double upper );

void free_Constraint( Constraint *c );

Constraint * clone_Constraint( Constraint *cnstr );

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

#pragma mark -
#pragma mark Parameter

Parameter * new_Parameter( const char *name, const double value, Constraint *constr );

Parameter * new_Parameter_with_postfix( const char *name, const char *postfix, const double value, Constraint *constr );

Parameter* new_Parameter_from_json(json_node* node, Hashtable* hash);

Parameters* new_MultiParameter_from_json(json_node* node, Hashtable* hash);

void free_Parameter( Parameter *p );

Parameter * clone_Parameter( Parameter *p );

json_node* Parameter_to_json(Parameter* parameter, json_node* parent);

json_node* Parameters_to_json(Parameters* parameters, json_node* parent);

char * Parameter_name( const Parameter *p );

void Parameter_set_name( Parameter *p, const char *name );

void Parameter_set_value( Parameter *p, const double value );

void Parameter_set_value_quietly( Parameter *p, const double value );

void Parameter_fire(Parameter *p);

double Parameter_value( const Parameter *p );

void Parameter_store(Parameter *p);

void Parameter_restore(Parameter *p);

void Parameter_restore_quietly(Parameter *p);

bool Parameter_changed(Parameter *p);

double check_value( Constraint *cnstr, double value );


bool Parameters_is_at_boundry( const Parameters *p, const size_t index, double precision );

Constraint * Parameter_constraint( const Parameter *p );

bool Parameter_estimate( const Parameter *p );

void Parameter_set_estimate( Parameter *p, const bool estimate );

double Parameters_fupper( const Parameters *p, const size_t index );

double Parameters_flower( const Parameters *p, const size_t index );

double Parameter_upper( const Parameter *p );

double Parameter_lower( const Parameter *p );

void Parameter_set_upper( Parameter *p, const double value );

void Parameter_set_lower( Parameter *p, const double value );

void Parameter_set_bounds( Parameter *p, const double lower, const double upper );


void Parameter_warn( const Parameter *p, const double value );

void Parameter_die( const Parameter *p, const double value );

void Parameter_print( const Parameter *p );

#pragma mark -
#pragma mark Parameters

Parameters * new_Parameters( const size_t capacity );

Parameters * new_Parameters_with_name( const char* name, const size_t capacity );

Parameters * new_Parameters_from_json(json_node* node, Hashtable* hash);

void free_Parameters( Parameters *ps );

Parameters * clone_Parameters( Parameters *p );

void Parameters_set_name2(Parameters* ps, const char* name);

const char* Parameters_name2(const Parameters* ps);

char* Parameters_mutable_name2(Parameters* ps);

Parameter * Parameters_at( const Parameters *p, const size_t index );

void Parameters_add(Parameters *ps, Parameter *p);

void Parameters_move( Parameters *ps, Parameter *p);

void Parameters_add_free_parameters(Parameters *dst, const Parameters *src);

void Parameters_add_parameters(Parameters *dst, const Parameters *src);

void Parameters_set_name( Parameters *p, const size_t index, const char *name );

char * Parameters_name( const Parameters *p, const size_t index );

size_t Parameters_count( const Parameters *p );

size_t Parameters_capacity( const Parameters *p );

void Parameters_set_value( Parameters *p, const size_t index, const double value );

void Parameters_set_value_quietly( Parameters *p, const size_t index, const double value );

void Parameters_set_values( Parameters *p, const double* values );

void Parameters_set_all_value( Parameters *p, const double value );

void Parameters_fire( Parameters *p );

double Parameters_value( const Parameters *p, const size_t index );

void Parameters_store(Parameters* ps);

bool Parameters_estimate( const Parameters *p, const size_t index );

void Parameters_set_estimate( Parameters *p, const bool estimate, const size_t index );

Constraint * Parameters_constraint( const Parameters *p, const size_t index );

double Parameters_upper( const Parameters *p, const size_t index );

double Parameters_lower( const Parameters *p, const size_t index );

void Parameters_set_upper( Parameters *p, const size_t index, const double value );

void Parameters_set_lower( Parameters *p, const size_t index, const double value );

void Parameters_set_bounds( Parameters *p, const size_t index, const double lower, const double upper );

void Parameters_remove( Parameters *params, size_t index );

void Parameters_removeAll( Parameters *params);

void Parameters_pop( Parameters *params );


Parameters * get_sub_parameters( Parameters *p, const int start, const int end );

Parameters * Parameters_optimizable( Parameters *p, int **map, const char postfix[] );

void Parameters_store_value( const Parameters *p, double *store );

void Parameters_restore_value( Parameters *p, const double *store );

int Parameters_count_optimizable( const Parameters *p, const char postfix[] );

Parameters * pack_parameters(Parameters *ps, Parameter *p);


void Parameters_print( Parameters *ps );

void Parameters_swap( Parameter **a, Parameter **b );

void Parameters_swap_index( Parameters *ps, unsigned a, unsigned b );

void Parameters_sort_from_ivector( Parameters *p, int *s );

void check_constraint(Parameter* rate, double lower, double upper, double flower, double fupper);

void check_constraints(Parameters* rates, double lower, double upper, double flower, double fupper);

#pragma mark -

struct _ListenerList {
	Model** models;
	int count;
	int capacity;
	void (*free)( ListenerList*);
	void (*fire)( ListenerList*, Model*, int );
	void (*fire_restore)( ListenerList*, Model*, int );
	void (*add)( ListenerList*, Model* );
	void (*remove)( ListenerList*, Model* );
	void (*removeAll)( ListenerList*);
};

model_t check_model(const char* type);

struct _Model {
	void *obj; // pointer to model
	char *name;
	model_t type;
	void* data;
	double (*logP)( Model * );
	double (*full_logP)( Model * );
	double (*dlogP)( Model *, const Parameter* );
	double (*d2logP)( Model *, const Parameter* );
	double (*ddlogP)(Model*, const Parameter*, const Parameter*);
	void (*gradient)( Model *, const Parameters*, double* );
	void (*prepare_gradient)( Model *, const Parameters* );
	Model* (*clone)( Model *, Hashtable* );
	void (*free)( Model * );
	void (*update)( Model *, Model *, int );
	void (*handle_restore)( Model *, Model *, int );
	void (*reset)(Model*);
	void (*sample)(Model*, double*, double* logP);
	double (*sample_evaluate)(Model*);
	
	ListenerList *listeners;
	bool need_update;
	int ref_count;
	
	void(*store)(Model*);
	void(*restore)(Model*);
	double lp;
	double storedLogP;
	bool samplable; // model is a distribution that can sampled directly
	void (*print)(Model*, FILE*);
    void (*jsonize)(Model*, json_node*);
};


#pragma mark -

Model * new_Model( model_t type, const char *name, void *obj );

void free_Model( Model *model );

double Model_first_derivative( Model *model, Parameter* parameter, double eps );

double Model_second_derivative( Model *model, Parameter* parameter, double* first, double eps );

double Model_mixed_derivative( Model *model, Parameter* p1, Parameter* p2 );

#pragma mark -

ListenerList * new_ListenerList( const unsigned capacity );

#pragma mark -

void get_parameters_from_node(json_node* node, Hashtable* hash, Parameters* parameters);

bool get_parameter_list_from_node(json_node* node, Parameters* parameters);

void get_parameters_references(json_node* node, Hashtable* hash, Parameters* parameters);

void get_parameters_references2(json_node* node, Hashtable* hash, Parameters* parameters, const char* tag);

void get_parameters_slice(const char* ref, Parameters* parameters, Hashtable* hash);

#endif
