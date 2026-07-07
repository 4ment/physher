// Copyright (C) 2010-2026 Mathieu Fourment
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "utils.h"

#include <float.h>

#include "mstring.h"

#include "mjson.h"
#include "hashtable.h"
#include "model.h"
// #include "transforms.h"

#define PARAMETER_TINY 1.0e-25
#define PARAMETER_ZERO 0.
#define PARAMETER_ZERO_PLUS TINY
#define PARAMETER_ONE 1.
#define PARAMETER_ONE_MINUS (1-TINY)
#define PARAMETER_POSITIVE_INFINTY INFINITY
#define PARAMETER_NEGATIVE_INFINTY (-INFINITY)

struct _Transform;
typedef struct _Transform Transform;

struct _Transform {
    size_t dim;            // dimension of x
    Parameter *parameter;  // y
    double lower;
    double upper;
    // y = f(x)
    void (*transform)(const double *, double *, size_t, double, double);
    // x = f^{-1}(y)
    void (*inverse_transform)(double *, const double *, size_t, double, double);
    // dL/dy = dL/dx * dx/dy
    void (*backward_inverse_transform)(double *, const double *, const double *, size_t,
                                       double, double);
    // dx/dy = df^{-1}(y)/dy
    double (*inverse_transform_log_det_jacobian)(double *, const double *, size_t,
                                                 double, double);
    double (*inverse_transform_gradient_log_det_jacobian)(double *, const double *,
                                                          size_t, double, double);
    void (*inverse_transform_jacobian)(double *, const double *, size_t, double,
                                       double);
    void (*get)(Transform *, double *);
    void (*set)(Transform *, const double *);
    void (*backward)(Transform *, const double *);  //  dx/dy
    void (*jacobian)(Transform *, double *);        //  dx/dy
    double (*log_det_jacobian)(Transform *);
    double (*gradient_log_det_jacobian)(Transform *);
};

Transform *new_Transform_with_parameter(const char *type, double lower, double upper,
                                        Parameter *parameter);

Transform *new_SimplexTransform_with_parameter(const char *type, Parameter *parameter);

void free_Transform(Transform* transform);

struct _Parameter{
	char *name;
	int id;
	double *value;
	double *stored_value;
	bool stored;
	size_t dim;
	bool simplex;
	Constraint *cnstr;
	bool estimate;
	ListenerList *listeners;
	int refCount;
	model_t model; // model it belongs to
	double* grad;
	// Model* model_obj;
	Transform *transform;
	void (*update)(Parameter *, Model *, Parameter *, int);
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

Parameter * new_Parameter2( const char *name, const double *value, size_t dim, Constraint *constr );

Parameter * new_Parameter_with_postfix( const char *name, const char *postfix, const double value, Constraint *constr );
Parameter * new_Parameter_with_postfix2( const char *name, const char *postfix, const double* value, size_t dim, Constraint *constr );
Parameter * new_Parameter_full( const char *name, double value, size_t dim, Constraint *constr);

Parameter* new_Parameter_from_json(json_node* node, Hashtable* hash);

Parameters* new_MultiParameter_from_json(json_node* node, Hashtable* hash);

void free_Parameter( Parameter *p );

Parameter * clone_Parameter( Parameter *p );

json_node* Parameter_to_json(Parameter* parameter, json_node* parent);

json_node* Parameters_to_json(Parameters* parameters, json_node* parent);

char * Parameter_name( const Parameter *p );

void Parameter_set_name( Parameter *p, const char *name );

size_t Parameter_size(const Parameter* p);

void Parameter_set_value( Parameter *p, const double value );
void Parameter_set_values( Parameter *p, const double *values );
void Parameter_set_value_at( Parameter *p, const double value, size_t index);


void Parameter_grad_mul_inverse_transform(Parameter* p);
void Parameter_grad_add_log_det_inverse_transform(Parameter* p);
void Parameter_zero_grad(Parameter*p);

void Parameter_set_value_quietly( Parameter *p, const double value );
void Parameter_set_values_quietly( Parameter *p, const double* value );
void Parameter_set_value_at_quietly( Parameter *p, const double value, size_t index);

void Parameter_fire(Parameter *p, int index);

double Parameter_value(const Parameter *p);
double Parameter_value_at(const Parameter *p, size_t index);
const double *Parameter_values(const Parameter *p);

void Parameter_set_model( Parameter *p, model_t model );

void Parameter_save_to( const Parameter *p, double *dst );

void Parameter_store(Parameter *p);

void Parameter_restore(Parameter *p);

void Parameter_restore_quietly(Parameter *p);

void Parameter_accept(Parameter* p);

bool Parameter_changed(Parameter *p);

double check_value( Constraint *cnstr, double value );


bool Parameters_is_at_boundry( const Parameters *p, const size_t index, double precision );

Constraint * Parameter_constraint( const Parameter *p );

Parameter *Parameters_depends(const Parameters *parameters, const Parameter *x);

bool Parameter_estimate( const Parameter *p );

void Parameter_set_estimate( Parameter *p, const bool estimate );

double Parameters_fupper( const Parameters *p, const size_t index );

double Parameters_flower( const Parameters *p, const size_t index );

double Parameter_upper( const Parameter *p );

double Parameter_lower( const Parameter *p );

void Parameter_set_upper( Parameter *p, const double value );

void Parameter_set_lower( Parameter *p, const double value );

void Parameter_set_bounds( Parameter *p, const double lower, const double upper );

void Parameter_allocate_grad(Parameter* p);

#pragma mark -
#pragma mark Parameters

Parameters * new_Parameters( const size_t capacity );

Parameters * new_Parameters_with_name( const char* name, const size_t capacity );

Parameters * new_Parameters_from_json(json_node* node, Hashtable* hash);

void free_Parameters( Parameters *ps );

void free_Parameters_weak( Parameters *ps );

Parameters * clone_Parameters( Parameters *p );

void Parameters_set_name2(Parameters* ps, const char* name);

const char* Parameters_name2(const Parameters* ps);

char* Parameters_mutable_name2(Parameters* ps);

Parameter * Parameters_at( const Parameters *p, const size_t index );

void Parameters_add(Parameters *ps, Parameter *p);

void Parameters_add_recursively(Parameters *ps, Parameter *p);

void Parameters_move( Parameters *ps, Parameter *p);

void Parameters_add_free_parameters(Parameters *dst, const Parameters *src);

void Parameters_add_parameters(Parameters *dst, const Parameters *src);

void Parameters_add_parameters_recursively(Parameters *dst, const Parameters *src);

void Parameters_add_listener(Parameters *parameters, Model *model);

void Parameters_set_name( Parameters *p, const size_t index, const char *name );

char * Parameters_name( const Parameters *p, const size_t index );

size_t Parameters_count( const Parameters *p );

size_t Parameters_capacity( const Parameters *p );

void Parameters_set_value( Parameters *p, const size_t index, const double value );

void Parameters_set_value_quietly( Parameters *p, const size_t index, const double value );

void Parameters_set_values( Parameters *p, const double* values );

void Parameters_set_values_quietly( Parameters *p, const double* values );

void Parameters_set_all_value( Parameters *p, const double value );

double Parameters_value( const Parameters *p, const size_t index );

void Parameters_store(Parameters* ps);

void Parameters_restore(Parameters *p);

void Parameters_accept(Parameters *p);

bool Parameters_estimate( const Parameters *p, const size_t index );

void Parameters_set_estimate( Parameters *p, const bool estimate, const size_t index );

Constraint * Parameters_constraint( const Parameters *p, const size_t index );

double Parameters_upper( const Parameters *p, const size_t index );

double Parameters_lower( const Parameters *p, const size_t index );

void Parameters_set_upper( Parameters *p, const size_t index, const double value );

void Parameters_set_lower( Parameters *p, const size_t index, const double value );

void Parameters_set_bounds( Parameters *p, const size_t index, const double lower, const double upper );

double Parameter_fupper(const Parameter *p);

double Parameter_flower(const Parameter *p);

void Parameters_remove( Parameters *params, size_t index );

void Parameters_removeAll( Parameters *params);

void Parameters_pop( Parameters *params );


Parameters * get_sub_parameters( Parameters *p, const int start, const int end );

Parameters * Parameters_optimizable( Parameters *p, int **map, const char postfix[] );

void Parameters_store_value( const Parameters *p, double *store );

void Parameters_restore_value( Parameters *p, const double *store );

int Parameters_count_optimizable( const Parameters *p, const char postfix[] );

Parameters * pack_parameters(Parameters *ps, Parameter *p);

void Parameters_swap( Parameter **a, Parameter **b );

void Parameters_swap_index( Parameters *ps, unsigned a, unsigned b );

void Parameters_sort_from_ivector( Parameters *p, int *s );

void check_constraint(Parameter* rate, double lower, double upper, double flower, double fupper);

void check_constraints(Parameters* rates, double lower, double upper, double flower, double fupper);

bool Parameters_contains(const Parameters *parameters, const Parameter *p);

void Parameters_zero_grad(Parameters *parameters);

size_t Parameters_size(const Parameters *ps);

#pragma mark -

void *safe_get_reference_parameter(const char *ref, Hashtable *hash,
                                   const char *parent);

void *safe_get_reference_model(const char *ref, Hashtable *hash,
                               const char *parent);

bool safe_is_reference(const char *ref, const char *parent);

#pragma mark -

Parameter * new_ParameterModel( const char *name, const double* value, size_t dim, Constraint *constr, Model* model);

#pragma mark -

void get_parameters_from_node(json_node* node, Hashtable* hash, Parameters* parameters);

bool get_parameter_list_from_node(json_node* node, Parameters* parameters);

void get_parameters_references(json_node* node, Hashtable* hash, Parameters* parameters);

void get_parameters_references2(json_node* node, Hashtable* hash, Parameters* parameters, const char* tag);

void get_parameters_slice(const char* ref, Parameters* parameters, Hashtable* hash);

void get_parameter_reference(const char* ref, Hashtable* hash, Parameters* parameters);

void grab_parameters(json_node *node, Hashtable *hash, Parameters *parameters);

#endif
