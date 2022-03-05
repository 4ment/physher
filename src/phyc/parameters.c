/*
 *  parameters.c
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

#include "parameters.h"

#include <assert.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <math.h>

#include "mstring.h"
#include "utils.h"
#include "matrix.h"

#include "simplex.h"


struct _Constraint{
	double lower;
	double upper;
	bool lower_fixed;
	bool upper_fixed;
	double flower;
	double fupper;
};

struct _Parameters{
	char* name;
	Parameter **list;
	size_t count;
	size_t capacity;
};


#pragma mark -
#pragma mark Constraint

Constraint * new_Constraint( const double lower, const double upper ){
	Constraint *c = (Constraint *)malloc( sizeof(Constraint) );
	assert(c);
	c->lower = lower;
	c->upper = upper;
	c->lower_fixed = false;
	c->upper_fixed = false;
	
	c->flower = lower;
	c->fupper = upper;
	//( ( lower == DBL_MIN && upper == DBL_MAX ) ? true : false);
	return c;
}

void free_Constraint( Constraint *c ){
	free(c);
	c = NULL;
}

Constraint * clone_Constraint( Constraint *cnstr ){
	Constraint *new = new_Constraint(cnstr->lower, cnstr->upper);
	new->lower_fixed = cnstr->lower_fixed;
	new->upper_fixed = cnstr->upper_fixed;
	return new;
}

void Constraint_set_bounds( Constraint *constr, const double lower, const double upper ){
	if( constr != NULL ){
		constr->lower = lower;
		constr->upper = upper;
		constr->flower = lower;
		constr->fupper = upper;
	}
}

bool Constraint_lower_fixed( const Constraint *c ){
	return c->lower_fixed;
}

bool Constraint_upper_fixed( const Constraint *c ){
	return c->upper_fixed;
}

void Constraint_set_lower_fixed( Constraint *c, const bool fixed ){
	c->lower_fixed = fixed;
}


void Constraint_set_upper_fixed( Constraint *c, const bool fixed ){
	c->upper_fixed = fixed;
}


double Constraint_upper( const Constraint *c ){
	return c->upper;
}

double Constraint_lower( const Constraint *c ){
	return c->lower;
}

double Constraint_fupper( const Constraint *c ){
	return c->fupper;
}

double Constraint_flower( const Constraint *c ){
	return c->flower;
}

void Constraint_set_upper( Constraint *c, const double upper ){
	c->upper = upper;
}

void Constraint_set_lower( Constraint *c, const double lower ){
	c->lower = lower;
}

void Constraint_set_fupper( Constraint *c, const double fupper ){
	c->fupper = fupper;
}

void Constraint_set_flower( Constraint *c, const double flower ){
	c->flower = flower;
}

#pragma mark -
#pragma mark Parameter

// Owns the constraint by default
Parameter * new_Parameter( const char *name, const double value, Constraint *constr ){
	return new_Parameter_with_postfix( name, "", value, constr );
}

Parameter * new_Parameter_with_postfix( const char *name, const char *postfix, const double value, Constraint *constr ){
	Parameter *p = (Parameter *)malloc( sizeof(Parameter) );
	assert(p);
	size_t name_len    = strlen(name);
	size_t postfix_len = strlen(postfix);
	
	p->name = (char*)malloc( sizeof(char)*(name_len+postfix_len+2) );
	assert(p->name);
	strcpy(p->name, name);
	if ( postfix_len > 0 ) {
		p->name[name_len] = '.';
		strcpy(p->name+name_len+1, postfix);
	}
	p->value = value;
	p->stored_value = value;
	p->cnstr = constr;
	p->estimate = true;
	p->id = 0;
	p->listeners = new_ListenerList(1);
	p->refCount = 1;
	p->model = -1;
	return p;
}

Parameter* new_Parameter_from_json(json_node* node, Hashtable* hash){
	char* allowed[] = {
		"flower",
		"fupper",
		"lower",
		"upper",
		"value"
	};
	json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
	
	if (node->node_type == MJSON_STRING) {
		char* ref = (char*)node->value;
		Parameter*p = Hashtable_get(hash, ref+1);
		p->refCount++;
		return p;
	}
	double value = 0;
	json_node* value_node = get_json_node(node, "value");
	if(value_node != NULL){
		value = atof((char*)value_node->value);
	}
	double lower = -INFINITY;
	double upper = INFINITY;
	json_node* lower_node = get_json_node(node, "lower");
	if(lower_node != NULL && strcasecmp((char*)lower_node->value, "-infinity") != 0){
		lower = atof((char*)lower_node->value);
		
	}
	json_node* upper_node = get_json_node(node, "upper");
	if(upper_node != NULL && strcasecmp((char*)upper_node->value, "infinity") != 0){
		upper = atof((char*)upper_node->value);
	}
	Constraint* cnstr = new_Constraint(lower, upper);
	lower_node = get_json_node(node, "flower");
	if(lower_node != NULL && strcasecmp((char*)lower_node->value, "-infinity") != 0){
		cnstr->flower = atof((char*)lower_node->value);
		
	}
	upper_node = get_json_node(node, "fupper");
	if(upper_node != NULL && strcasecmp((char*)upper_node->value, "infinity") != 0){
		cnstr->fupper = atof((char*)upper_node->value);
	}
	json_node* id_node = get_json_node(node, "id");
	return new_Parameter((char*)id_node->value, value, cnstr);
}

Parameters* new_MultiParameter_from_json(json_node* node, Hashtable* hash){
    char* allowed[] = {
        "dimension",
        "flower",
        "fupper",
        "lower",
        "upper",
        "values"
    };
    json_check_allowed(node, allowed, sizeof(allowed)/sizeof(allowed[0]));
    
    if (node->node_type == MJSON_STRING) {
        char* ref = (char*)node->value;
        Parameters*p = Hashtable_get(hash, ref+1);
        // refCount of parameters inside p are NOT incremented
        return p;
    }
    
    size_t dim = get_json_node_value_size_t(node, "dimension", 0);
    Parameters* parameters = new_Parameters(dim);
    
    json_node* lower_node = get_json_node(node, "lower");
    json_node* upper_node = get_json_node(node, "upper");
    char* id = get_json_node_value_string(node, "id");
    Parameters_set_name2(parameters, id);
    double lower = -INFINITY;
    double upper = INFINITY;
    if(lower_node != NULL && strcasecmp((char*)lower_node->value, "-infinity") != 0){
        lower = atof((char*)lower_node->value);
        
    }
    if(upper_node != NULL && strcasecmp((char*)upper_node->value, "infinity") != 0){
        upper = atof((char*)upper_node->value);
        
    }
    json_node* values = get_json_node(node, "values");
    size_t K = values->child_count;
    if (dim == 0) {
        dim = K;
    }
    if(dim < K){
        fprintf(stderr, "%s - dimension attribute (%zu) cannot be smaller than the number of values (%zu)\n", id, dim, K);
        exit(2);
    }
    
    StringBuffer* buffer = new_StringBuffer(10);
    int i = 0;
    while (i != dim) {
        for(int j = 0; j < values->child_count; j++){
            if(i == dim) break;
            Constraint* cnstr = new_Constraint(lower, upper);
            json_node* child = values->children[j];
            double value = atof((char*)child->value);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%s.%d", id, i+1);
            Parameters_move(parameters, new_Parameter(buffer->c, value, cnstr));
            i++;
        }
    }
    free_StringBuffer(buffer);
    return parameters;
}

Parameters * new_Parameters_from_json(json_node* node, Hashtable* hash){
	Parameters* parameters = new_Parameters(node->child_count);
	bool found = get_parameter_list_from_node(node, parameters);
	if (!found) {
		get_parameters_from_node(node, hash, parameters);
	}
	return parameters;
}

Parameter * clone_Parameter( Parameter *p ){
	Parameter *pnew = NULL;
	Constraint *cnstr = NULL;

	if( p->cnstr != NULL ){
		cnstr = clone_Constraint(p->cnstr);
	}
	
	pnew = new_Parameter( p->name, p->value,  cnstr );//pnew ref_count is reset to 1
	pnew->estimate = p->estimate;
	pnew->id = p->id;
	pnew->model = p->model;
    pnew->stored_value = p->stored_value;
	return pnew;
}

void free_Parameter( Parameter *p ){
	if(p == NULL) return;
	if(p->refCount == 1){
		free(p->name);
		p->listeners->free(p->listeners);
		if( p->cnstr != NULL ) free_Constraint(p->cnstr);
		free(p);
	}
	else{
		p->refCount--;
	}
}

json_node* Parameter_to_json(Parameter* parameter, json_node* parent){
    json_node* jnode = create_json_node_object(parent, parameter->name);
    add_json_node(parent, jnode);
    add_json_node_string(jnode, "id", parameter->name);
    add_json_node_string(jnode, "type", "parameter");
    add_json_node_double(jnode, "value", parameter->value);
    if(!isinf(parameter->cnstr->lower))
        add_json_node_double(jnode, "lower", parameter->cnstr->lower);
    if(!isinf(parameter->cnstr->upper))
        add_json_node_double(jnode, "upper", parameter->cnstr->upper);
//    add_json_node_double(node, "flower", parameter->cnstr->flower);
//    add_json_node_double(node, "fupper", parameter->cnstr->fupper);
    return jnode;
}

json_node* Parameters_to_json(Parameters* parameters, json_node* parent){
    json_node* jnode = create_json_node_object(parent, Parameters_name2(parameters));
    add_json_node(parent, jnode);
    add_json_node_string(jnode, "id", Parameters_name2(parameters));
    add_json_node_string(jnode, "type", "parameter");
    add_json_node_size_t(jnode, "dimension", Parameters_count(parameters));
    if(!isinf(Parameters_lower(parameters, 0)))
        add_json_node_double(jnode, "lower", Parameters_lower(parameters, 0));
    if(!isinf(Parameters_upper(parameters, 0)))
        add_json_node_double(jnode, "upper", Parameters_upper(parameters, 0));
    double* values = dvector(Parameters_count(parameters));
    for (int i = 0; i < Parameters_count(parameters); i++) {
        values[i] = Parameters_value(parameters, i);
    }
    add_json_node_array_double(jnode, "values", values, Parameters_count(parameters));
    free(values);
    return jnode;
}

void Parameter_set_name( Parameter *p, const char *name ){
	if ( p->name != NULL ) {
		free(p->name);
	}
	p->name = String_clone(name);
}

char * Parameter_name( const Parameter *p ){
	return p->name;
}

void Parameter_set_value( Parameter *p, const double value ){
	p->value = value;
	p->listeners->fire(p->listeners, NULL, p->id);
}


void Parameter_set_value_quietly( Parameter *p, const double value ){
	p->value = value;
}

void Parameter_fire(Parameter *p){
	p->listeners->fire(p->listeners, NULL, p->id);
}

double Parameter_value( const Parameter *p ){
	return p->value;
}

void Parameter_store(Parameter *p){
	p->stored_value = p->value;
}

void Parameter_restore(Parameter *p){
	if (p->stored_value != p->value) {
		p->value = p->stored_value;
		p->listeners->fire_restore(p->listeners, NULL, p->id);
	}
}

void Parameter_restore_quietly(Parameter *p){
	p->value = p->stored_value;
}

bool Parameter_changed(Parameter *p){
	return p->value != p->stored_value;
}

void assign_value( Parameter *p, const double value ){
	if ( p->cnstr == NULL ) 
		p->value = value;
	else{
		if ( value >= p->cnstr->lower &&  value <= p->cnstr->upper ) p->value = value;
		else if( value < p->cnstr->lower ) p->value = p->cnstr->lower;
		else if( value > p->cnstr->upper ) p->value = p->cnstr->upper;
		else {
			fprintf(stderr, "\n\n");
			Parameter_print(p);
			error("parameters.c:assign_value:should not be here\n");
		}
	}
#ifdef ASSIGN_DEBUG
	Parameter_die(p, value);
#endif
}

Constraint * Parameter_constraint( const Parameter *p ){
	return p->cnstr;
}

bool Parameter_estimate( const Parameter *p ){
	return p->estimate;
}

void Parameter_set_estimate( Parameter *p, const bool estimate ){
	p->estimate = estimate;
}

double Parameter_upper( const Parameter *p ){
	return Constraint_upper( p->cnstr );
}

double Parameter_lower( const Parameter *p ){
	return Constraint_lower( p->cnstr );	
}

void Parameter_set_upper( Parameter *p, const double value ){
	Constraint_set_upper( p->cnstr, value );
}

void Parameter_set_lower( Parameter *p, const double value ){
	Constraint_set_lower( p->cnstr, value );	
}

void Parameter_set_bounds( Parameter *p, const double lower, const double upper ){
	if ( p->cnstr == NULL ) {
		p->cnstr = new_Constraint(lower, upper);
	}
	else {
		Constraint_set_bounds( p->cnstr, lower, upper );
	}
}


void Parameter_warn( const Parameter *p, const double value ){
	if ( p->cnstr != NULL && (value < p->cnstr->lower || value > p->cnstr->upper ) ){
		fprintf(stderr, "Parameter out of bound\n%s\n",__FILE__);
		Parameter_print(p);
		fprintf(stderr, "\n");
	}
}

void Parameter_die( const Parameter *p, const double value ){
	if ( p->cnstr != NULL && (value < p->cnstr->lower || value > p->cnstr->upper ) ){
		fprintf(stderr, "Parameter out of bound\n%s\n",__FILE__);
		Parameter_print(p);
		fprintf(stderr, "\n");
		exit(1);
	}
}

double check_value( Constraint *cnstr, double value ){
	if ( cnstr == NULL ) {
		return value;
	}
	else{
		if ( value >= cnstr->lower &&  value <= cnstr->upper ) return value;
#ifdef ASSIGN_DEBUG
		fprintf(stderr, "Parameter out of bound\n%s\n",__FILE__);
		fprintf(stderr, "%f [%f - %f]\n", value, cnstr->lower, cnstr->upper);
		exit(0);
#endif
		if( value < cnstr->lower ) return cnstr->lower;
		else if( value > cnstr->upper ) return cnstr->upper;
		fprintf(stderr, "lower %f value = %f upper = %f\n", cnstr->lower ,value, cnstr->upper);
		error("parameters.c:check_value:should not be here\n");
	}
	return -1;
}



void Parameter_print( const Parameter *p ){
	if(p->cnstr != NULL)fprintf(stderr, "%s %e [%e - %e] estimate = %d lower_fixed = %d upper_fixed = %d\n",p->name, p->value, p->cnstr->lower, p->cnstr->upper, p->estimate, p->cnstr->lower_fixed, p->cnstr->upper_fixed);
	else fprintf(stderr, "%s %e estimate = %d\n",p->name, p->value, p->estimate);
}

#pragma mark -
#pragma mark Parameters

Parameters * new_Parameters( const size_t capacity ){
    Parameters *ps = (Parameters *)malloc( sizeof(Parameters) );
    assert(ps);
    ps->name = NULL;
    ps->list = (Parameter **)malloc( capacity * sizeof(Parameter *) );
    assert(ps->list);
    for (int i = 0; i < capacity; i++) {
        ps->list[i] = NULL;
    }
    ps->count = 0;
    ps->capacity = capacity;
    return ps;
}


Parameters * new_Parameters_with_name( const char* name, const size_t capacity ){
    Parameters* ps = new_Parameters(capacity);
    ps->name = String_clone(name);
    return ps;
}

void free_Parameters( Parameters *ps ){
	if(ps == NULL) return;
	
	for (int i = 0; i < ps->count; i++) {
		free_Parameter(ps->list[i]);
	}
	if(ps->name != NULL) free(ps->name);
	free(ps->list);
	free(ps);
	ps = NULL;
}

// The cloned parameters can have different capacities but same count of course
Parameters * clone_Parameters( Parameters *p ){
	Parameters * clone = new_Parameters( p->capacity );
	for (int i = 0; i < p->count; i++) {
		Parameters_move(clone, clone_Parameter( p->list[i] ) );
	}
	if(p->name != NULL){
		clone->name = String_clone(p->name);
	}
	return clone;
}

void Parameters_set_name2(Parameters* ps, const char* name){
	if(ps->name != NULL){
		free(ps->name);
	}
	ps->name = String_clone(name);
}

const char* Parameters_name2(const Parameters* ps){
	return ps->name;
}

char* Parameters_mutable_name2(Parameters* ps){
	return ps->name;
}

Parameter * Parameters_at( const Parameters *p, const size_t index ){
	return p->list[index];
}

void Parameters_add( Parameters *ps, Parameter *p){
	if( ps->count == ps->capacity ){
		ps->list = realloc(ps->list, (ps->capacity+1) * sizeof(Parameter*) );
		assert(ps->list);
		ps->capacity++;
	}
	ps->list[ps->count++] = p;
	p->refCount++;
}

// add and do not increment refCount (move ownership)
void Parameters_move( Parameters *ps, Parameter *p){
	if( ps->count == ps->capacity ){
		ps->list = realloc(ps->list, (ps->capacity+1) * sizeof(Parameter*) );
		assert(ps->list);
		ps->capacity++;
	}
	ps->list[ps->count++] = p;
}

void Parameters_add_parameters(Parameters *dst, const Parameters *src){
	if(Parameters_count(src) == 0) return;
	if( dst->count+src->count > dst->capacity ){
		dst->list = realloc(dst->list, (dst->count+src->count) * sizeof(Parameter*) );
		assert(dst->list);
		dst->capacity = dst->count+src->count;
	}
	for (int i = 0; i < src->count; i++) {
		Parameters_add( dst, src->list[i] );
	}
}

void Parameters_add_free_parameters(Parameters *dst, const Parameters *src){
	if(Parameters_count(src) == 0) return;
	if( dst->count+src->count > dst->capacity ){
		dst->list = realloc(dst->list, (dst->count+src->count) * sizeof(Parameter*) );
		assert(dst->list);
		dst->capacity = dst->count+src->count;
	}
	for (int i = 0; i < src->count; i++) {
		if(src->list[i]->estimate){
			Parameters_add( dst, src->list[i] );
		}
	}
}

Parameters * pack_parameters(Parameters *ps, Parameter *p){
	ps->list = realloc(ps->list, ps->count * sizeof(Parameter*) );
	assert(ps->list);
	ps->capacity = ps->count;
	return ps;
}



void Parameters_set_name( Parameters *p, const size_t index, const char *name ){
	assert( index < p->count);
	Parameter_set_name(p->list[index], name);
}

char * Parameters_name( const Parameters *p, const size_t index ){
	return Parameter_name(p->list[index]);
}

size_t Parameters_count( const Parameters *p ){
    if( p == NULL ) return 0;
	return p->count;
}


size_t Parameters_capacity( const Parameters *p ){
	return p->capacity;
}

void Parameters_set_value( Parameters *p, const size_t index, const double value ){
	assert( index < p->count);
	Parameter_set_value(p->list[index], value);
}

void Parameters_set_value_quietly( Parameters *p, const size_t index, const double value ){
	assert( index < p->count);
	Parameter_set_value_quietly(p->list[index], value);
}

void Parameters_set_values( Parameters *p, const double* values ){
	for (int i = 0; i < Parameters_count(p); ++i) {
		p->list[i]->value = values[i];
	}
	Parameter_fire(p->list[0]);
}

void Parameters_fire( Parameters *p ){
	for (int i = 0; i < Parameters_count(p); ++i) {
		Parameter_fire(p->list[i]);
	}
}

void Parameters_set_all_value( Parameters *p, const double value ){
	for (int i = 0; i < Parameters_count(p); ++i) {
		Parameter_set_value(p->list[i], value);
	}
}

double Parameters_value( const Parameters *p, const size_t index ){
	return Parameter_value(p->list[index]);
}

void Parameters_store(Parameters* ps){
	for (int i = 0; i < Parameters_count(ps); i++) {
		Parameter_store(Parameters_at(ps, i));
	}
}

bool Parameters_estimate( const Parameters *p, const size_t index ){
	return Parameter_estimate(p->list[index] );
}

void Parameters_set_estimate( Parameters *p, const bool estimate, const size_t index ){
	Parameter_set_estimate(p->list[index], estimate );
}

Constraint * Parameters_constraint( const Parameters *p, const size_t index ){
	return Parameter_constraint( Parameters_at(p, index) );
}

double Parameters_fupper( const Parameters *p, const size_t index ){
	return Constraint_fupper( Parameters_at(p, index)->cnstr );
}

double Parameters_flower( const Parameters *p, const size_t index ){
	return Constraint_flower( Parameters_at(p, index)->cnstr );
}

double Parameters_upper( const Parameters *p, const size_t index ){
	return Parameter_upper( Parameters_at(p, index) );
}

double Parameters_lower( const Parameters *p, const size_t index ){
	return Parameter_lower( Parameters_at(p, index) );
}

void Parameters_set_upper( Parameters *p, const size_t index, const double value ){
	Parameter_set_upper( Parameters_at(p, index), value );
}

void Parameters_set_lower( Parameters *p, const size_t index, const double value ){
	Parameter_set_lower( Parameters_at(p, index), value );	
}

void Parameters_set_bounds( Parameters *p, const size_t index, const double lower, const double upper ){
	Parameter_set_bounds( Parameters_at(p, index), lower, upper );
}


bool Parameters_is_at_boundry( const Parameters *p, const size_t index, double precision ){
    return ( Parameter_value(p->list[index]) - Parameters_lower(p, index) < precision
            || Parameters_upper(p, index) - Parameter_value(p->list[index]) < precision );
    
}

// not cloned
Parameters * get_sub_parameters( Parameters *p, const int start, const int end ){
	Parameters * sub = new_Parameters( end - start + 1 );
	for (int i = start; i <= end; i++) {
		Parameters_add(sub, p->list[i]);
	}
	return sub;
}

void Parameters_print( Parameters *ps ){
	fprintf(stderr, "==================================\n");
	for (int i = 0; i < ps->count; i++) {
		Parameter_print(ps->list[i]);
	}
	fprintf(stderr, "----------------------------------\n");
}

// The parameters are not cloned
// assume that .height are consecutive, followed by .bl...
Parameters * Parameters_optimizable( Parameters *p, int **map, const char postfix[] ){
	if( p == NULL ) return NULL;
	
	int count = Parameters_count_optimizable(p, postfix);
	
	if( count == 0 ) return NULL;

	*map = malloc( sizeof(int) * count );
	assert(map);
	
	Parameters * sub = new_Parameters( count );
	count = 0;

	int index = 0; 

	for ( int i = 0; i < p->count; i++) {
		if ( postfix != NULL ) {
			if( index == 0 && strstr(p->list[i]->name, postfix) != NULL )
				index = i;
		}
			
		
		if( p->list[i]->estimate ){
			if ( postfix != NULL ){
				if(strstr(p->list[i]->name, postfix) == NULL ) continue;
			}
			
			Parameters_add(sub, p->list[i]);
			(*map)[count++] = i - index;
		}
	}
	return sub;
}

void Parameters_store_value( const Parameters *p, double *store ){
	for (int i = 0; i < Parameters_count(p); i++) {
		store[i] = Parameters_value(p, i);
	}
}

void Parameters_restore_value( Parameters *p, const double *store ){
	for (int i = 0; i < Parameters_count(p); i++) {
		Parameters_set_value(p, i, store[i]);
	}
}


int Parameters_count_optimizable( const Parameters *p, const char postfix[] ){
	if( p == NULL ) return 0;
	int count = 0;
	for ( int i = 0; i < p->count; i++) {
		if( p->list[i]->estimate ){
			if ( postfix != NULL ){
				if(strstr(p->list[i]->name, postfix) != NULL ) count++;
			}
			else count++;
		}
	}
	return count;
}

void Parameters_remove( Parameters *params, size_t index ){
	assert(index < params->count);
	if( params->count !=0 ){
		free_Parameter( params->list[index] );
        
        for( size_t i = index+1; i < params->count; i++ ){
            params->list[i-1] = params->list[i];
        }
		params->list[params->count-1] = NULL;
		params->count--;
	}
}

void Parameters_removeAll( Parameters *params){
	for( size_t i = 0; i < params->count; i++ ){
		free_Parameter( params->list[i] );
		params->list[i] = NULL;
	}
	params->count = 0;
}

void Parameters_pop( Parameters *params ){
	if( params->count !=0 ){
		free_Parameter( params->list[params->count-1] );
		params->list[params->count-1] = NULL;
		params->count--;
	}
}

void Parameters_swap( Parameter **a, Parameter **b ){
	Parameter *p = *a;
	*a = *b;
	*b = p;
}

void Parameters_swap_index( Parameters *ps, unsigned a, unsigned b ){
	Parameter *p = ps->list[a];
	ps->list[a] = ps->list[b];
	ps->list[b] = p;
}

void Parameters_sort_from_ivector( Parameters *p, int *s ){
	bool done = false;
	size_t size = p->count;
	while ( !done ) {
		done = true;
		for ( int i = 0 ; i < size-1 ; i++ ) {
			if ( s[i] > s[i+1] ) {
				done = false;
				swap_int( &s[i], &s[i+1] );
				Parameters_swap_index(p, i, i+1);
			}
		}
		size--;
	}
}

void check_constraint(Parameter* rate, double lower, double upper, double flower, double fupper){
	if (isinf(Constraint_lower(rate->cnstr))) {
		Constraint_set_lower(rate->cnstr, lower);
		Constraint_set_flower(rate->cnstr, flower);
	}
	if (isinf(Constraint_upper(rate->cnstr))) {
		Constraint_set_upper(rate->cnstr, upper);
		Constraint_set_fupper(rate->cnstr, fupper);
	}
}

void check_constraints(Parameters* rates, double lower, double upper, double flower, double fupper){
	for (int i = 0; i <Parameters_count(rates); i++) {
		check_constraint(Parameters_at(rates, i), lower, upper, flower, fupper);
	}
}

#pragma mark -

static void dummyUpdate( Model *self, Model *model, int index ){}
static double _logP(Model *model){return 0;}
static double _fulllogP(Model *model){return model->logP(model);}
static double _dlogP(Model *model, const Parameter* p){return 0;}
static double _d2logP(Model *model, const Parameter* p){return 0;}
static double _ddlogP(Model *self, const Parameter* p1, const Parameter* p2){return 0;}
static void _dummy_prepare_gradient(Model *model, const Parameters* ps){ }
static void _dummy_reset(Model* m){}
static void _dummy_restore(Model* m){}
static void _dummy_store(Model* m){}
static void _dummy_sample(Model* m, double* samples, double* logP){
	fprintf(stderr, "Cannot sample from model %s\n", m->name);
	exit(2);
}
static double _dummy_sample_evaluate(Model* m){
	fprintf(stderr, "Cannot sample and evaluate from model %s\n", m->name);
	exit(2);
}

static void _dummy_jsonize(Model* m, json_node* node){
    fprintf(stderr, "jsonize function not implemented for model %s\n", m->name);
    exit(2);
}

#pragma mark -


Model * new_Model( model_t type, const char *name, void *obj ){
	Model *model = (Model*)malloc(sizeof(Model));
	assert(model);
	model->name = String_clone(name);
	model->type = type;
	model->obj = obj;
	model->logP = _logP;
	model->full_logP = _fulllogP;
	model->dlogP = _dlogP;
	model->d2logP = _d2logP;
	model->ddlogP = _ddlogP;
	model->update  = dummyUpdate;
	model->handle_restore  = dummyUpdate;
	model->free = free_Model;
	model->data = NULL;
	model->listeners = new_ListenerList(1);
	model->need_update = true;
	model->clone = NULL;
	model->ref_count = 1;
	model->reset = _dummy_reset;
	model->restore = _dummy_restore;
	model->store = _dummy_store;
	model->logP = 0;
	model->lp = 0;
	model->sample = _dummy_sample;
	model->sample_evaluate = _dummy_sample_evaluate;
	model->samplable = false;
	model->print = NULL;
	model->prepare_gradient = _dummy_prepare_gradient;
    model->jsonize = _dummy_jsonize;
	return model;
}

void free_Model( Model *model ){
	assert(model->ref_count >= 1);
	if(model->ref_count == 1){
		free(model->name);
		model->listeners->free(model->listeners);
		free(model);
	}
	else{
		model->ref_count--;
	}
	
}

model_t check_model(const char* type){
	int len = sizeof(model_type_strings)/sizeof(model_type_strings[0]);
	for (int i = 0; i < len; i++) {
		if (strcasecmp(model_type_strings[i], type) == 0) {
			return i;
		}
	}
	return -1;
}

#pragma mark -

static void _free_ListenerList( ListenerList *listeners ){
	free(listeners->models);
	free(listeners);
}

static void _ListenerList_fire( ListenerList *listeners, Model* model, int index){
	for ( int i = 0; i < listeners->count; i++ ) {
		listeners->models[i]->update( listeners->models[i], model, index );
	}
}

static void _ListenerList_fire_restore( ListenerList *listeners, Model* model, int index){
	for ( int i = 0; i < listeners->count; i++ ) {
		listeners->models[i]->handle_restore( listeners->models[i], model, index );
	}
}

static void _ListenerList_remove( ListenerList *listeners, Model* model ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		if ( listeners->models[i] == model ) {
			break;
		}
	}
	if ( i == listeners->count) {
		return;
	}
	i++;
	for ( ; i < listeners->count; i++ ) {
		listeners->models[i-1] = listeners->models[i];
	}
	listeners->models[listeners->count-1] = NULL;
	listeners->count--;
}

static void _ListenerList_remove_all( ListenerList *listeners ){
	int i = 0;
	for ( ; i < listeners->count; i++ ) {
		listeners->models[i] = NULL;
	}
	listeners->count = 0;
}

static void _ListenerList_add( ListenerList *listeners, Model *model ){
	if ( listeners->count == listeners->capacity) {
		listeners->capacity++;
		listeners->models = realloc(listeners->models, listeners->capacity*sizeof(Model*));
	}
	listeners->models[listeners->count] = model;
	listeners->count++;
}
ListenerList * new_ListenerList( const unsigned capacity ){
	ListenerList *listeners = (ListenerList*)malloc( sizeof(ListenerList));
	assert(listeners);
	listeners->capacity = (capacity == 0 ? 1 : capacity);
	listeners->count = 0;
	listeners->models = (Model**)malloc( listeners->capacity * sizeof(Model*));
	assert(listeners->models);
	listeners->free = _free_ListenerList;
	listeners->add = _ListenerList_add;
	listeners->remove = _ListenerList_remove;
	listeners->removeAll = _ListenerList_remove_all;
	listeners->fire = _ListenerList_fire;
	listeners->fire_restore = _ListenerList_fire_restore;
	return listeners;
}

double Model_mixed_derivative( Model *model, Parameter* p1, Parameter* p2 ) {
	double eps = 0.001;
	
	double v1 = Parameter_value(p1);
	double e1 = eps * v1;
	double v2 = Parameter_value(p2);
	double e2 = eps * v2;
	
	// + +
	Parameter_set_value(p1, v1 + e1);
	Parameter_set_value(p2, v2 + e2);
	double pp = model->logP(model);
	
	// - -
	Parameter_set_value(p1, v1 - e1);
	Parameter_set_value(p2, v2 - e2);
	double mm = model->logP(model);
	
	// + -
	Parameter_set_value(p1, v1 + e1);
	Parameter_set_value(p2, v2 - e2);
	double pm = model->logP(model);
	
	// - +
	Parameter_set_value(p1, v1 - e1);
	Parameter_set_value(p2, v2 + e2);
	double mp = model->logP(model);
	
	Parameter_set_value(p1, v1);
	Parameter_set_value(p2, v2);
	
	return (pp + mm - pm - mp) / (4.0*e1*e2);
}

double Model_second_derivative( Model *model, Parameter* parameter, double* first, double eps ) {
	//double eps = 0.00001;
	
	double lnl = model->logP(model);
	
	double v = Parameter_value(parameter);
	double e = eps * v;
	
	// +
	Parameter_set_value(parameter, v + e);
	double p = model->logP(model);
	
	// -
	Parameter_set_value(parameter, v - e);
	double m = model->logP(model);
	
	Parameter_set_value(parameter, v);
	
	if (first != NULL) {
		*first = (p - m)/(2.0*e);
	}
	
	return (p + m -2*lnl)/(e*e);
}

// First derivative using central differences
// eps = 1.0e-8;
double Model_first_derivative( Model *model, Parameter* parameter, double eps ) {	
	double v = Parameter_value(parameter);
	double e = eps * v;
	
	Parameter_set_value(parameter, v + e);
	double pp = model->logP(model);
	
	Parameter_set_value(parameter, v - e);
	double mm = model->logP(model);
	
	Parameter_set_value(parameter, v);
	
	return (pp - mm)/(2.0*e);
}

#pragma mark -

void get_parameters_slice(const char* ref, Parameters* parameters, Hashtable* hash){
	char* copy = String_clone(ref);
	char* start = copy;
	while (*start != '['){
		start++;
	}
	*start = '\0';
	Parameters* ps = Hashtable_get(hash, copy);
	start++;
	start[strlen(start)-1] = '\0';
	int begin = 0;
	int end = Parameters_count(ps);
	int inc = 1;
	int c = 0;
	char** chars = String_split_char(start, ':', &c);
	if(c > 1){
		if(strlen(chars[0]) > 0) begin = atoi(chars[0]);
		if(strlen(chars[1]) > 0){
			int temp = atoi(chars[1]);
			if(temp < 0) end = end + temp;
		}
		if(c == 3){
			inc = atoi(chars[2]);
		}
		// [:] [1:] [:3] [1:2] [1::1] [:4:2] [::2]
		if(inc > 0){
			for (int i = begin; i < end; i+=inc) {
				Parameters_add(parameters, Parameters_at(ps, i));
			}
		}
		// [:4:-1] [::-1] [::-2]
		else{
			for (int i = end-1; i >= begin; i+=inc) {
				Parameters_add(parameters, Parameters_at(ps, i));
			}
		}
	}
	else{
		// simple indexing
		int index = atoi(start);
		if (index < 0) {
			index = Parameters_count(ps) + index;
		}
		Parameters_add(parameters, Parameters_at(ps, index));
	}
	free_cmatrix(chars, c);
	free(copy);
}

void get_multi_parameter_from_node(json_node* node, Parameters* parameters){
    size_t dim = get_json_node_value_size_t(node, "dimension", 0);
    json_node* lower_node = get_json_node(node, "lower");
    json_node* upper_node = get_json_node(node, "upper");
    char* id = get_json_node_value_string(node, "id");
    double lower = -INFINITY;
    double upper = INFINITY;
    if(lower_node != NULL && strcasecmp((char*)lower_node->value, "-infinity") != 0){
        lower = atof((char*)lower_node->value);
        
    }
    if(upper_node != NULL && strcasecmp((char*)upper_node->value, "infinity") != 0){
        upper = atof((char*)upper_node->value);
        
    }
    json_node* values = get_json_node(node, "values");
    size_t K = values->child_count;
    if (dim == 0) {
        dim = K;
    }
    if(dim < K){
        fprintf(stderr, "%s - dimension attribute (%zu) cannot be smaller than the number of values (%zu)\n", id, dim, K);
        exit(2);
    }
    
    StringBuffer* buffer = new_StringBuffer(10);
    int i = 0;
    while (i != dim) {
        for(int j = 0; j < values->child_count; j++){
            if(i == dim) break;
            Constraint* cnstr = new_Constraint(lower, upper);
            json_node* child = values->children[j];
            double value = atof((char*)child->value);
            StringBuffer_empty(buffer);
            StringBuffer_append_format(buffer, "%s.%d", id, i+1);
            Parameters_move(parameters, new_Parameter(buffer->c, value, cnstr));
            i++;
        }
    }
    free_StringBuffer(buffer);
}

void get_parameters_from_node(json_node* node, Hashtable* hash, Parameters* parameters){
	if(node->node_type == MJSON_ARRAY){
		for (int i = 0; i < node->child_count; i++) {
			json_node* child = node->children[i];
			char* ref = (char*)child->value;
			// it's a ref
			if (child->node_type == MJSON_STRING) {
				if (ref[0] == '&') {
					Parameter* p = Hashtable_get(hash, ref+1);
					Parameters_add(parameters, p);
					
				}
				// tree
				else if (ref[0] == '%') {
					// slicing
					if (ref[strlen(ref)-1] == ']') {
						get_parameters_slice(ref+1, parameters, hash);					}
					else{
						Parameters* ps = Hashtable_get(hash, ref+1);
						Parameters_add_parameters(parameters, ps);
					}
					Parameters_set_name2(parameters, ref+1);
				}
				// simplex
				else if (ref[0] == '$') {
					Model* msimplex = Hashtable_get(hash, ref+1);
					Simplex* simplex = msimplex->obj;
					Parameters_add_parameters(parameters, simplex->parameters);
				}
			}
			// it's a value
			else if(child->node_type == MJSON_PRIMITIVE){
				double v = atof((char*)child->value);
				Parameters_move(parameters, new_Parameter("anonymous", v, NULL));
			}
			else{
				exit(1);
			}
		}
	}
	// it's a ref
	else if(node->node_type == MJSON_STRING){
		char* ref = (char*)node->value;
		if (ref[0] == '&') {
			Parameter* p = Hashtable_get(hash, ref+1);
			Parameters_add(parameters, p);
		}
		else if (ref[0] == '%') {
			// slicing
			if (ref[strlen(ref)-1] == ']') {
				get_parameters_slice(ref+1, parameters, hash);
			}
			else{
				Parameters* ps = Hashtable_get(hash, ref+1);
				Parameters_add_parameters(parameters, ps);
			}
			Parameters_set_name2(parameters, ref+1);
		}
		// simplex
		else if (ref[0] == '$') {
			Model* msimplex = Hashtable_get(hash, ref+1);
			Simplex* simplex = msimplex->obj;
			Parameters_add_parameters(parameters, simplex->parameters);
		}
	}
	else if(node->node_type == MJSON_OBJECT){
        if(get_json_node(node, "type") != NULL && get_json_node(node, "id") != NULL){
            json_node* p_node_dimension = get_json_node(node, "dimension");
            if(p_node_dimension != NULL){
                get_multi_parameter_from_node(node, parameters);
            }
            else{
                Parameter* p = new_Parameter_from_json(node, hash);
                Parameters_move(parameters, p);
            }
            return;
        }
		for(int i = 0; i < node->child_count; i++){
			json_node* p_node = node->children[i];
            json_node* p_node_dimension = get_json_node(p_node, "dimension");
            if(p_node_dimension != NULL){
                get_multi_parameter_from_node(p_node, parameters);
            }
            else{
                Parameter* p = new_Parameter_from_json(p_node, hash);
                Parameters_move(parameters, p);
            }
		}
	}
	else{
        fprintf(stderr, "Do not recognize node type of %s", node->key);
		exit(1);
	}
}

bool get_parameter_list_from_node(json_node* node, Parameters* parameters){
	// node can be multidimensional: "parameters":{"dimension": 4, "values": [2], "lower": 0, "upper": 100}
	for (int i = 0; i < node->child_count; i++) {
		json_node* child = node->children[i];
		if(strcasecmp(child->key, "dimension") == 0 && child->node_type == MJSON_PRIMITIVE){
            get_multi_parameter_from_node(node, parameters);
			return true;
		}
	}
	return false;
}

void get_parameters_references(json_node* node, Hashtable* hash, Parameters* parameters){
	bool found = get_parameter_list_from_node(node, parameters);
	if (found) return;
	
	get_parameters_references2(node, hash, parameters, "parameters");
}

void get_parameters_references2(json_node* node, Hashtable* hash, Parameters* parameters, const char* tag){
    json_node* x_node = get_json_node(node, tag);
	get_parameters_from_node(x_node, hash, parameters);
}


